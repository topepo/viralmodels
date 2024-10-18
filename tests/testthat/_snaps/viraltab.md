# viraltab() works

    Code
      print(viraltab(traindata, semilla, target, viralvars, logbase, pliegues,
        repeticiones, rejilla))
    Message
      i Creating pre-processing data to finalize unknown parameter: mtry
    Output
                          wflow_id              .config .metric   mean std_err n
      1                  simple_rf Preprocessor1_Model1    rmse 113.35    0.91 2
      2                  simple_rf Preprocessor1_Model1     rsq   0.79    0.05 2
      3         simple_CART_bagged Preprocessor1_Model1    rmse 117.26    2.16 2
      4         simple_CART_bagged Preprocessor1_Model1     rsq   0.77    0.05 2
      5             normalized_KNN Preprocessor1_Model1    rmse 127.32    6.92 2
      6             normalized_KNN Preprocessor1_Model1     rsq   0.72    0.00 2
      7              full_quad_KNN Preprocessor1_Model1    rmse 149.64    5.18 2
      8              full_quad_KNN Preprocessor1_Model1     rsq   0.67    0.03 2
      9  normalized_neural_network Preprocessor1_Model1    rmse 157.17   17.18 2
      10 normalized_neural_network Preprocessor1_Model1     rsq   0.65    0.07 2
      11       normalized_SVM_poly Preprocessor1_Model1    rmse 166.73   15.43 2
      12       normalized_SVM_poly Preprocessor1_Model1     rsq   0.69    0.01 2
      13             simple_Cubist Preprocessor1_Model1    rmse 179.49   26.23 2
      14             simple_Cubist Preprocessor1_Model1     rsq   0.68    0.01 2
      15     normalized_SVM_radial Preprocessor1_Model1    rmse 233.47    2.09 2
      16     normalized_SVM_radial Preprocessor1_Model1     rsq   0.65    0.00 2
      17      full_quad_linear_reg Preprocessor1_Model1    rmse 984.45  822.60 2
      18      full_quad_linear_reg Preprocessor1_Model1     rsq   0.27    0.25 2
               preprocessor            model rank
      1  workflow_variables      rand_forest    1
      2  workflow_variables      rand_forest    1
      3  workflow_variables         bag_tree    2
      4  workflow_variables         bag_tree    2
      5              recipe nearest_neighbor    3
      6              recipe nearest_neighbor    3
      7              recipe nearest_neighbor    4
      8              recipe nearest_neighbor    4
      9              recipe              mlp    5
      10             recipe              mlp    5
      11             recipe         svm_poly    6
      12             recipe         svm_poly    6
      13 workflow_variables     cubist_rules    7
      14 workflow_variables     cubist_rules    7
      15             recipe          svm_rbf    8
      16             recipe          svm_rbf    8
      17             recipe       linear_reg    9
      18             recipe       linear_reg    9

