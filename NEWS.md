# viralmodels 1.3.2

* Added a `NEWS.md` file to track changes to the package.

## Major changes

* Added `obs-by-obs` functionality to the `viralpreds` function, allowing row-by-row predictions in addition to normal full-dataset predictions. This is controlled by the new `prediction_type` parameter.

* Added `batch` functionality to the `viralpreds` function, allowing smaller size batches predictions in addition to normal full-dataset predictions. This is controlled by the new `prediction_type` parameter.
