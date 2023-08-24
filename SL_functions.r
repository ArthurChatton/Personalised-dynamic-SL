####################################################""""
### Implementation of POSL (Malenica et al., 2022, arXiv)
###################################################

#1) user-supplied input (data & Lrnr)

library(sl3) # learners
library(origami) #cross-validation


#2) Functions call in cSL()

`%nin%` <- Negate(`%in%`)


# time-weighted meta-learner (NNLS), needed because Lrnr_nnls does not be weighted (yet?) through the wrapper Lrnr_ts_weights. It's therefore hardcoded as in the code of I. Malenica and R. Phillips.
get_weights_nnls <- function(pred, observed, weights, convex){
  fit.nnls <- nnls::nnls(sqrt(weights) * as.matrix(pred), 
                         sqrt(weights) * as.numeric(observed))
  initCoef <- coef(fit.nnls)
  initCoef[is.na(initCoef)] <- 0.0
  
  if(convex) {
    # normalize so sum(coef) = 1 if possible
    if (sum(initCoef) > 0) {
      coef <- initCoef / sum(initCoef)
    } else {
      coef <- initCoef
    }
  } else {
    coef <- initCoef
  }
  return(as.numeric(coef))
}

# to compute time-dpt weights, weights_control is a list of the time-dependant parameters, and times is the predicted session
process_weights <- function(weights_control, times){
  # intialize weights of 1 for all losses
  weights <- rep(1, length(times))
  
  # update weights based on weights_control list
  if (!is.null(weights_control$window)) {
    window <- max(times) - weights_control$window
    weights <- weights * ifelse(times <= window, 0, 1)
  }
  
  if (!is.null(weights_control$rate_decay)) {
    lags <- max(times) - times
    if (!is.null(weights_control$delay_decay)) {
      lags_delayed <- lags - weights_control$delay_decay
      lags <- ifelse(lags_delayed < 0, 0, lags_delayed)
    }
    weights <- weights * (1 - weights_control$rate_decay)^lags
  }
  return(weights)
}

# to bound the NLL
bound <- function(preds, bounds = 0.001) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  preds_bounded <- pmin(pmax(preds, lower), upper)
  return(preds_bounded)
}

#function to be CVed
cv_forecasts <- function(fold, data, stack_ind, stack_hist, outcome, covar, session, historical=TRUE, type='continuous', lab){
  
  print(fold$v)
  
  # Get training and validation data
  train_data <- training(data)  
  valid_data <- validation(data) 
  valid_size <- length(valid_data) # required by cross_validate()
  
  #remove predictor with only one value (i.e., time-fixed predictors) for convergence
  covar_ind <- names(train_data[,covar])[apply(train_data[,covar], 2, FUN= function(x) length(unique(x))!=1)]

  
  #individual learners' predictions
  #Define task for training and validation
  train_task_ind <- make_sl3_Task(
    data = train_data,
    outcome = outcome,
    covariates = covar_ind,
    outcome_type=type
  )

  valid_task_ind <- make_sl3_Task( 
    data = valid_data,
    outcome = outcome,
    covariates = covar_ind,
    outcome_type=type
  )
  
  
  ## specific part for my application, not useful in other context
  if(length(covar_ind)<5){ # 5 being the number of var selected by rf
    lrn_train <- stack_ind$params$learners[[1]]$train(train_task_ind)
    pred <- c(lrn_train$predict(valid_task_ind), lrn_train$predict(valid_task_ind), stack_ind$params$learners[[3]]$train(train_task_ind)$predict(valid_task_ind))
  }else{
  ## end of the specific part  
  
  
  
   # learners training
    lrn_train <- stack_ind$train(train_task_ind)
    
    # prediction on validation set
    pred <- lrn_train$predict(valid_task_ind)
  } #remove } if not in the specific case of my application 
  
  if(historical){
    
    valid_task_hist <- make_sl3_Task( #needed for historical learners predictions, since the covariates arg differs
      data = valid_data,
      outcome = outcome,
      covariates = covar,
      outcome_type=type
    )
    
    # predictions from historical learners
    pred_hist <- unlist(lapply(1:length(stack_hist), function(x) stack_hist[[x]]$predict(valid_task_hist))) 
    
    pred <- c(pred, pred_hist)
    
  }
  
  #loss for the dSL selector (square root for type='continuous', negative loglikelihood for type='binomial')
  if(type=='continuous'){
    sqe <- sapply(pred, function(x){
      sum((valid_data[,outcome] - x)^2)
    })
  }else{
    sqe <- sapply(pred, function(x){
      -1 * ifelse(valid_data[,outcome] == 1, log(bound(x)), log(1 - bound(x)))
    })
  }
  

  
  output <- list(
    meta_db = data.frame( #meta-level dataset
      fold = fold$v, 
      outcome = valid_data[,outcome],
      data.frame(as.list(pred))
    ) |> setNames(nm=c('fold', 'outcome', lab)), 
    sqe = data.frame( #estimator-specific loss for dSL selector
      fold = fold$v,
      data.frame(as.list(sqe))
    ) |> setNames(nm=c('fold', lab)),
    pred = data.frame( #estimator-specific predictions
      fold = fold$v,
      data.frame(as.list(pred))
    ) |> setNames(nm=c('fold', lab))
  )
  
  if(any(is.na(names(lrn_train$learner_fits)))){ #remove folds where lrnrs don't fit
    next
  }else{
    return(output)
  }  
  
}

lab_stack <- function(hist){ #return the names of the different candidate learners
  
  
  dig <- function(stk){
    while(!is.null(stk$name) && (stk$name=="Stack" | length(grep("Pipeline", stk$name))>0)){
      stk <- stk$params$learners
    }
    return(stk)
  }
  
  label <- NULL
  for(i in 1:length(hist)){
    
    for(j in 1:length(tmp <- hist[[i]]$params$learners)){
      
      lab_tmp <- NULL
      stk <- tmp[[j]] |> dig() 
      if(length(stk)==2){
        if(stk[[2]]$name=="Stack"){
          lab_tmp <- stk[[1]]$name
          stk <- stk[[2]] |> dig()
        }
        stk <- stk |> lapply(function(x) x$name) |> unlist() 
      }else{
        if(is.list(stk)){
          stk <- unlist(sapply(stk, function(x) x$name))
        }else{
          stk <- stk$name
        }
      }
      
      if(!is.null(lab_tmp)) stk <- paste(stk, lab_tmp, sep="_")
      if(!is.null(names(hist)[[i]])) stk <- paste(names(hist)[[i]], stk, sep="_")
      
      label <- c(label, stk)
      
    }
    
  }
  
  return(label)                      
}



# 3) clustered SL function

cSL <- function(data, newdata, outcome, type="continuous", covar, session, stack_ind, stack_hist, tw_arg=list(window=NULL,delay_decay=NULL,rate_decay=NULL), cv_rolling="O", fw=10, b=1, cv_save=NA){ #fw is the minimal size of the training set (in term of obs) and b is the delay between each fold
  
  
  ########## Format ############
  
  if(!cv_rolling %in% c("W", "B", "O")) stop("Wrong value for the argument 'cv_rolling', must be B, O or W.")
  
  if(!is.list(stack_hist)){ #one can supply several historical stacks (e.g., with different individuals / temporal windows), this line deals when only one stack is supplied
    stack_hist <- list(stack_hist)
  }
  
  label <- c(
              paste0("ind_", lab_stack(list(stack_ind))),
              paste0("hist_", lab_stack(stack_hist))
  )
  
  ########## Cross-validation ###########
  
  
  if(cv_rolling=="O"){
    folds <- folds_rolling_origin(
      dim(data)[1],
      first_window = fw, validation_size = 1, gap = 0, batch = b
    )
  }else{
    if(cv_rolling=="W"){
      folds <- folds_rolling_window(
        dim(data)[1],
        window_size = fw, validation_size = 1, gap = 0, batch = b
      )
    }else{ #both CV schemes
      folds <- list(O=folds_rolling_origin(
                        dim(data)[1],
                        first_window = fw, batch = b,
                        validation_size = 1, gap = 0
                    ),
                    W=folds_rolling_window(
                      dim(data)[1],
                      window_size = fw, batch = b,
                      validation_size = 1, gap = 0
                    )
                )
    }
  }
  
  if(any(!is.na(cv_save))){
    new_folds <- setdiff(unlist(lapply(folds[["W"]], "[[", 'v')), cv_save$meta_db$fold)
  }else{
    new_folds <-  unlist(lapply(folds[["W"]], "[[", 'v')) 
  }


  
  if(cv_rolling=="B"){
    
    CV1 <- cross_validate(
              cv_fun = cv_forecasts, fold = folds[["O"]][new_folds], data = data, stack_ind=stack_ind, stack_hist=stack_hist, outcome=outcome, type=type, covar=covar, session=session, historical=FALSE, lab=paste0(lab_stack(list(stack_ind)), '_indO'),
              use_future = FALSE
            ) 
          
    CV2 <- cross_validate(
              cv_fun = cv_forecasts, fold = folds[["W"]][new_folds], data = data, stack_ind=stack_ind, stack_hist=stack_hist, outcome=outcome, type=type, covar=covar, session=session, lab=c(paste0(lab_stack(list(stack_ind)), '_indW'), paste0(lab_stack(stack_hist), '_hist')),
              use_future = FALSE
            ) 
    
    if(length(CV1$errors$index)>0)  warning(paste0("At least one individual learner did not converge in the following origin-rolling cross-validation folds: ", paste(CV1$errors$index, collapse=", "), ". These folds were not considered."), immediate. = TRUE)
    
    if(length(CV2$errors$index)>0)  warning(paste0("At least one individual learner did not converge in the following window-rolling cross-validation folds: ", paste(CV2$errors$index, collapse=", "), ". These folds were not considered."), immediate. = TRUE)
    
    CV1$errors <- CV2$errors <- NULL
    
    CV <- Map(merge, CV1, CV2, all=F, by="fold"); rm(CV1); rm(CV2) #doesnot work if CV1 and CV2 are of different lengths due to non-convergent learner
    
    CV$meta_db <- CV$meta_db[,names(CV$meta_db) != "outcome.y"]
  
    label <- c(paste0(grep("ind_", label, value=TRUE), "_OriginCV"), paste0(grep("ind_", label, value=TRUE), "_WindowCV"), grep("hist_", label, value=TRUE))
    
  }else{
    
    CV <- cross_validate(
      cv_fun = cv_forecasts, fold = folds, data = data, stack_ind=stack_ind, stack_hist=stack_hist, outcome=outcome, type=type, covar=covar, session=session,
      use_future = FALSE
    ) 
    
    if(!is.null(CV$errors))  warning(paste0("At least one learner did not converge in the following cross-validation folds: ", paste(CV$errors$index, collapse=", "), ". These folds were not considered."), immediate. = TRUE)
    
  }
  
  if(any(!is.na(cv_save))) CV <- Map(rbind, cv_save, CV)
  
  
  ########## Meta-learner ###########
  
  tw <- process_weights(weights_control=tw_arg, times=CV$meta_db$fold)
  
  
  # lrnr specific weights
  lrn_weights_convex <- get_weights_nnls(pred=CV$meta_db[,-(1:2)], observed=CV$meta_db$outcome, weights=tw, convex=TRUE)
  lrn_weights <- get_weights_nnls(pred=CV$meta_db[,-(1:2)], observed=CV$meta_db$outcome, weights=tw, convex=FALSE)
  
  
  
  ########## Prediction on newdata ###########
  
  covar_ind <- names(data[,covar])[apply(data[,covar], 2, FUN= function(x) length(unique(x))!=1)]
  
  # task for predictions for the new session
  task_ind <- make_sl3_Task( #for individual lrnrs' training
    data = data, 
    outcome = outcome,
    covariates = covar_ind,
    outcome_type=type
  )
  
  newtask_ind <- make_sl3_Task( #for individual lrnrs' prediction
    data = newdata, 
    outcome = outcome,
    covariates = covar_ind,
    outcome_type=type
  )
  
  newtask_hist <- make_sl3_Task( #for historical lrnrs' prediction
    data = newdata, 
    outcome = outcome,
    covariates = covar,
    outcome_type = type
  )
  
  if(length(covar_ind)<5){ # 5 being the number of var selected by rf
    pred_lrnr <- c(stack_ind$params$learners[[1]]$train(task_ind)$predict(newtask_ind) |> rep(2) |> unlist(), stack_ind$params$learners[[3]]$train(task_ind)$predict(newtask_ind))
    pred_lrnr <- unlist(c(rep(pred_lrnr, ifelse(cv_rolling=="B",2,1)),
                          unlist(lapply(1:length(stack_hist), function(x) stack_hist[[x]]$predict(newtask_hist)))
    ))
  }else{
  
  pred_lrnr <- unlist(c(
    rep(stack_ind$train(task_ind)$predict(newtask_ind), ifelse(cv_rolling=="B",2,1)),
    unlist(lapply(1:length(stack_hist), function(x) stack_hist[[x]]$predict(newtask_hist)))
    ))
  }
  
  
  pred_esl <- t(as.matrix(pred_lrnr)) %*% lrn_weights 
  pred_esl_convex <- t(as.matrix(pred_lrnr)) %*% lrn_weights_convex
  
  
  ########## Output ###########
  
  
  output <- list(
    eSL_convex=pred_esl_convex, #eSL final prediction
    eSL = pred_esl,
    dSL=pred_lrnr[which.min(colSums(tw*CV$sqe[,-1]))] |> setNames(nm=label[which.min(colSums(tw*CV$sqe[,-1]))]), #dSL final prediction labelled with the corresponding lrnr
    weights_convex=lrn_weights_convex |> setNames(nm=label), #weights in eSL
    weights=lrn_weights |> setNames(nm=label),
    lrnr=pred_lrnr |> setNames(nm=label), #lrnr specific predictions
    loss=colSums(tw*CV$sqe)[-1] |> setNames(nm=label), # CV error for each fold and learner
    dsl_lrnr = label[which.min(colSums(tw*CV$sqe[,-1]))], # lrnr selected for the dSL
    CV = CV #meta level data set to pass to the next loop/prediction
    )
  
  return(output)
}



# 4) Forward validation


run_csl <- function(ind, data, outcome, type, covar, stack_hist, stack_ind, fw=10){ 
  
  ## Pre-processing (data splitting according to \{ind} and historical training)
  ind_data <- data[data$id==ind,]
  hist_data <- data[data$id!=ind,]
  
  lrn_hist <- stack_hist$train(task=
                                 make_sl3_Task(
                                   data = hist_data, 
                                   outcome = outcome,
                                   outcome_type=type,
                                   covariates = covar
                                 )
  )
  
  ## POSL looped over each session
  
  N <- sort(ind_data$obs)[11:length(ind_data$obs)]
  
  prediction <- omega <- omega_convex <- dsl_lrnr <- NULL
  cv_save <- NA
  
  for(s in N[-length(N)]){  
    
    if(nrow(ind_data[ind_data$obs<=s,]) <= fw) next #pass if missing sessions at beginning because CV not possible
    
    if(exists('csl_res')) cv_save <- csl_res$CV
    
    
    s_valid <- ind_data$obs[1+which(sort(ind_data$obs)==s)]
    
    csl_res <- cSL(data=ind_data[ind_data$obs<=s,], 
                   newdata=ind_data[ind_data$obs==s_valid,],
                   outcome=outcome, covar=covar, session="obs",
                   stack_ind=stack_ind, 
                   stack_hist=lrn_hist, 
                   cv_rolling = "B", 
                   tw_arg=list(window=NULL, delay_decay=5, rate_decay=0.1),
                   fw=fw,
                   type=type,
                   cv_save=cv_save
    ) 
    
    prediction <- rbind(prediction, c(ind_data[ind_data$obs==s_valid,outcome], csl_res$eSL_convex, csl_res$eSL, csl_res$dSL, csl_res$lrnr))
    omega <- rbind(omega, csl_res$weights)
    omega_convex <- rbind(omega_convex, csl_res$weights_convex)
    dsl_lrnr <- c(dsl_lrnr, csl_res$dsl_lrnr)
    
  }
  
  ## Initial prediction at t1 with only the historical lrnrs (no individual trajectory)
  
  ipred <- lrn_hist$predict(task=make_sl3_Task( 
    data = ind_data[which.min(ind_data$obs),], #in case of unordered data
    outcome = outcome,
    outcome_type=type,
    covariates = covar
  )
  )
  
  cnames <- c('True', 'eSL_convex', 'eSL', 'dSL', colnames(prediction)[5:ncol(prediction)])
  
  prediction <- rbind(c(ind_data[which.min(ind_data$obs),outcome], rep(NA, ncol(prediction)-length(ipred)-1), unlist(ipred)), prediction) 
  
  
  ## Output
  
  output <- list(
    pred = prediction,
    weights = omega,
    weights_convex= omega_convex,
    dsl_lrnr = dsl_lrnr
  )
  
  dimnames(output$pred)[[2]] <- cnames
  
  return(output)
  
} 
