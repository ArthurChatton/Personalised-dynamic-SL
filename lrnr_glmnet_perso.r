library(R6)
Lrnr_glmnet_perso <- R6Class(
  classname = "Lrnr_glmnet",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lambda = NULL, type.measure = "deviance",
                          nfolds = 10, alpha = 1, nlambda = 100,
                          use_min = TRUE, stratify_cv = FALSE, ...) {
      super$initialize(params = args_to_list(), ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical",
      "weights", "ids"
    ),
    .train = function(task) {
      args <- self$params
      
      outcome_type <- self$get_outcome_type(task)
      
      if (is.null(args$family)) {
        args$family <- outcome_type$glm_family()
      }
      
      # specify data
      args$x <- as.matrix(task$X)
      args$y <- outcome_type$format(task$Y)
      
      if (task$has_node("weights")) {
        args$weights <- task$weights
      }
      
      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      
      if (task$has_node("id")) {
        args$foldid <- origami::folds2foldvec(task$folds)
      }
      
      if (args$stratify_cv) {
        if (outcome_type$type == "binomial" & is.null(args$foldid)) {
          folds <- origami::make_folds(
            n = length(args$y), strata_ids = args$y, fold_fun = folds_vfold,
            V = as.integer(args$nfolds)
          )
          args$foldid <- origami::folds2foldvec(folds)
        } else {
          warning(
            "stratify_cv is TRUE; but inner cross-validation folds cannot ",
            "be stratified. Either the outcome is not binomial, or foldid ",
            "has already been established (user specified foldid upon ",
            "initializing the learner, or it was set according to task id's)."
          )
        }
      }
      
      fit_object <- sl3:::call_with_args(
        glmnet::glmnet, args,
        other_valid = names(formals(glmnet::glmnet)),
        ignore = c("use_min", "stratify_cv")
      )
      fit_object$glmnet.fit$call <- NULL
      return(fit_object)
    },
    .predict = function(task) {
      args <- list(
        object = private$.fit_object, newx = as.matrix(task$X), type = "response"
      )
      
      # set choice regularization penalty
      if (self$params$use_min) {
        args$s <- "lambda.min"
      } else {
        args$s <- "lambda.1se"
      }
      
      if (task$has_node("offset")) {
        if (private$.fit_object$glmnet.fit$offset) {
          args$newoffset <- task$offset
        } else {
          warning(
            "Prediction task has offset, but an offset was not included in ",
            "the task for training the glmnet learner. The prediction task's ",
            "offset will not be considered for prediction."
          )
        }
      }
      
      # get predictions via S3 method
      predictions <- do.call(stats::predict, args)
      
      # reformat predictions based on outcome type
      if (private$.training_outcome_type$type == "categorical") {
        cat_names <- dimnames(predictions)[[2]]
        # predictions is a 3-dim matrix, convert to 2-dim matrix
        dim(predictions) <- dim(predictions)[1:2]
        colnames(predictions) <- cat_names
        # pack predictions in a single column
        predictions <- pack_predictions(predictions)
      }
      return(predictions)
    },
    .required_packages = c("glmnet", "origami")
  )
)

