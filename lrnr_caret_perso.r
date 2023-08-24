library(R6)

Lrnr_caret_correct <- R6Class(
  classname = "Lrnr_caret", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(algorithm,
                          metric = NULL,
                          trControl = caret::trainControl(method = "cv"),
                          ...) {
      params <- list(
        method = algorithm,
        metric = metric,
        ...
      )
      # provide two ways for users to specify trControl
      ## 1. Pass the method to `method`
      ## 2. Pass a list of trainControl arguments to `trControl`
      if (typeof(trControl) == "list") {
        params$trControl <- sl3:::call_with_args(caret::trainControl, trControl)
      } else {
        stop("Specified trControl type is unsupported in Lrnr_caret.")
      }
      super$initialize(params = params)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "categorical", "wrapper"),
    .train = function(task){
      # set type
      outcome_type <- 'binomial'
      train_type <- "classification"

      # load args
      args <- self$params
      
      # data
      args$x <- as.matrix(task$X)
      args$y <- task$Y
      
      # metric
      if (is.null(args$metric)) {
        args$metric <- ifelse(train_type == "regression", "RMSE", "Accuracy")
      }
      
      # fit
      fit_object <- sl3:::call_with_args(caret::train, args, keep_all = TRUE)
      return(fit_object)
    },
    .predict = function(task) {
      outcome_type <- self$training_outcome_type
      if (outcome_type$type == "continuous") {
        predict_type <- "regression"
      } else if (outcome_type$type %in% c("binomial", "categorical")) {
        predict_type <- "classification"
      } else {
        stop("Specified outcome type is unsupported in Lrnr_caret.")
      }
      
      if (predict_type == "regression") {
        predictions <- stats::predict(
          private$.fit_object,
          newdata = task$X, type = "raw"
        )
      } else {
        predictions <- stats::predict(
          private$.fit_object,
          newdata = task$X, type = "prob"
        )[, 2]
      }
      predictions <- as.numeric(predictions)
      return(predictions)
    },
    .required_packages = c("caret")
  )
)

