#' Train and Evaluate Many Regression Models for Predicting Viral Load or CD4 Counts
#'
#' This function builds, trains, and evaluates a set of statistical learning models for predicting viral load or CD4 counts. 
#' It implements multiple pre-processing options (simple, normalized, full quadratic) and model types (MARS, neural network, KNN). 
#' The best model is selected based on RMSE.
#'
#' @param target A character string specifying the column name of the target variable to predict.
#' @param pliegues An integer specifying the number of folds for cross-validation.
#' @param repeticiones An integer specifying the number of times the cross-validation should be repeated.
#' @param rejilla An integer specifying the number of grid search iterations for tuning hyperparameters.
#' @param semilla An integer specifying the seed for random number generation to ensure reproducibility.
#' @param data A data frame containing the predictors and the target variable.
#'
#' @return A list containing two elements: \code{predictions} (a vector of predicted values for the target variable) 
#' and \code{RMSE} (the root mean square error of the best model).
#' @export
#'
#' @examples
#' \donttest{
#' library(tidyverse)
#' library(baguette)
#' library(kernlab)
#' library(kknn)
#' library(ranger)
#' library(rules)
#' library(glmnet)
#' # Define the function to impute values in the undetectable range
#' set.seed(123)
#' impute_undetectable <- function(column) {
#' ifelse(column <= 40,
#'       rexp(sum(column <= 40), rate = 1/13) + 1,
#'             column)
#'             }
#' # Apply the function to all vl columns using purrr's map_dfc
#' library(viraldomain)
#' data("viral", package = "viraldomain")
#' viral_imputed <- viral |>
#' mutate(across(starts_with("vl"), ~impute_undetectable(.x)))
#' traindata <- viral_imputed
#' target <- "cd_2022"
#' viralvars <- c("vl_2019", "vl_2021", "vl_2022")
#' logbase <- 10
#' pliegues <- 5
#' repeticiones <- 2
#' rejilla <- 2
#' semilla <- 123
#' viralpreds(target, pliegues, repeticiones, rejilla, semilla, traindata)
#' }
viralpreds <- function(target, pliegues, repeticiones, rejilla, semilla, data) {
  set.seed(semilla)
  workflow_set <- workflowsets::workflow_set(
    preproc = list(
      simple = workflows::workflow_variables(outcomes = tidyselect::all_of(target), predictors = tidyselect::everything()), 
      normalized = recipes::step_normalize(recipes::recipe(stats::as.formula(paste(target, "~ .")), data = data), recipes::all_predictors()), 
      full_quad = recipes::step_interact(recipes::step_poly(recipes::step_normalize(recipes::recipe(stats::as.formula(paste(target, "~ .")), data = data), recipes::all_predictors()), recipes::all_predictors()), ~all_predictors():all_predictors())
    ), 
    models = list(
      MARS = parsnip::set_mode(parsnip::set_engine(parsnip::mars(prod_degree = parsnip::tune(), num_terms = parsnip::tune(), prune_method = parsnip::tune()), "earth"), "regression"),
      neural_network = parsnip::set_mode(parsnip::set_engine(parsnip::mlp(hidden_units = parsnip::tune(), penalty = parsnip::tune(), epochs = parsnip::tune()), "nnet", MaxNWts = 2600), "regression"), 
      KNN = parsnip::set_mode(parsnip::set_engine(parsnip::nearest_neighbor(neighbors = parsnip::tune(), dist_power = parsnip::tune(), weight_func = parsnip::tune()), "kknn"), "regression")
    )
  )
  
  workflow_results <- workflowsets::workflow_map(
    workflow_set,
    seed = semilla,
    resamples = rsample::vfold_cv(rsample::training(rsample::initial_split(data)), v = pliegues, repeats = repeticiones),
    grid = rejilla,
    control = tune::control_grid(save_pred = TRUE, parallel_over = "everything", save_workflow = TRUE)
  )
  
  best_workflow <- workflow_results |>
    dplyr::mutate(rmse = purrr::map_dbl(result, ~ .x |>
                                          tune::show_best(metric = "rmse") |>
                                          dplyr::slice(1) |>
                                          dplyr::pull(mean))) |>
    dplyr::slice_min(rmse, with_ties = FALSE) |>
    dplyr::pull(wflow_id)
  
  best_metrics <- workflow_results |>
    dplyr::mutate(metrics = purrr::map(result, ~ .x |>
                                         tune::show_best(metric = c("rmse")))) |>
    dplyr::filter(wflow_id == best_workflow) |>
    dplyr::pull(metrics)
  
  best_model <- tune::select_best(
    workflowsets::extract_workflow_set_result(workflow_results, id = best_workflow),
    metric = "rmse"
  )
  
  final_workflow <- tune::finalize_workflow(
    workflowsets::extract_workflow(workflow_set, id = best_workflow),
    best_model
  )
  
  fitted_model <- parsnip::fit(final_workflow, data = data)
  predictions <- stats::predict(fitted_model, new_data = data)
  
  return(list(predictions = predictions$.pred, RMSE = best_metrics[[1]]$mean[1]))
}

utils::globalVariables(c("result", "rmse", "wflow_id", "metrics"))