source("SL_functions.r") #POSL functions
source("lrnr_glmnet_perso.r") #glmnet learner adapted because the native sl3 function use cv.glmnet() rather than glmnet()

#import db here

#####
## Ind - Hist - Tun dbs
#####

set.seed(291122)
draw <- sample(unique(db$id), 173, replace=FALSE)
db_tune <- db[db$id %in% draw,]
db_train <- db[db$id %in% setdiff(unique(db$id), draw),] 


####
## Tuning parameters definition
####

library(caret)
library(sl3)
library(ranger)

task_tune <- make_sl3_Task(
  data = db_tune,
  outcome = outcome,
  covariates = covar
)


tun_ranger <- Lrnr_caret$new(algorithm = "ranger", name = "rf_autotune")$train(task_tune)
tun_ranger$fit_object$bestTune
# mtry  splitrule min.node.size
#   20 extratrees             5

tun_lasso <- Lrnr_caret$new(algorithm="glmnet", name="lasso_autotune", tuneGrid=expand.grid(
  .alpha=1,
  .lambda=seq(0, 1, by = 0.001)))$train(task_tune)
tun_lasso$fit_object$bestTune
# lambda = 0.004 

tun_ridge <- Lrnr_caret$new(algorithm="glmnet", name="ridge_autotune", tuneGrid=expand.grid(
  .alpha=0,
  .lambda=seq(0, 1, by = 0.001)))$train(task_tune)
tun_ridge$fit_object$bestTune
# lambda = 0.285

# tun_mars <- Lrnr_caret$new(algorithm = "earth", name = "mars_autotune")$train(task_tune)
# tun_mars$fit_object$bestTune #does not work...
train(x=db_tune[,covar], y=db_tune[,outcome], method="earth", tuneGrid=expand.grid(degree = 1:2, nprune = (1:10) * 2))
# nprune = 6 and degree = 2

tun_xgb <- Lrnr_caret$new(algorithm = "xgbTree", name = "xgb_autotune")$train(task_tune)
tun_xgb$fit_object$bestTune
# nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
# 150         3     0.4     0              0.8                1         1


####
## Learners definition
####

# initialisation of the learners
lrn_glm <- Lrnr_glm$new(name="lm")
lrn_mean <- Lrnr_mean$new(name="mean")
lrn_lasso <- Lrnr_glmnet_perso$new(alpha = 1, lambda = 0.004, name="lasso", family="gaussian")
lrn_ridge <- Lrnr_glmnet_perso$new(alpha = 0, lambda = 0.285, name="ridge", family="gaussian")
lrn_mars <- Lrnr_earth$new(nprune = 20, degree = 2, name="mars")
lrn_xgb <- Lrnr_xgboost$new(nrounds=150, max_depth=3, eta=0.4, gamma=0, colsample_bytree=0.8, min_child_weight=1, subsample=1, name="xgb")


lrn_rf_importance <- Lrnr_ranger$new(importance = "impurity_corrected", splitrule="extratrees", min.node.size=5) #top 5 predictors only, mtry unspecified due to convergence issue
RFscreen_top5 <- Lrnr_screener_importance$new(
  learner = lrn_rf_importance, num_screen = 5, name="rf"
)

# stack historical learners together and we will apply the functions directly on the stack (similar to the apply family)
stack <- Stack$new(lrn_glm, lrn_mars, lrn_ridge, lrn_lasso, lrn_xgb)
RFscreen_top5_stack <- Pipeline$new(RFscreen_top10, stack)


stack_hist <- stack_ind <- Stack$new(stack, RFscreen_top5_stack, lrn_mean)

####
## Sequential cSL
####

seeds <- read.table("seeds.txt",header=F)$V1


runned <- list.files(getwd(), pattern="prediction") |> stringr::str_extract("\\d+") |> as.numeric() # recover id of individuals already predicted (in case of a potential bug)

torun <-  unique(db_train$id)[unique(db_train$id) %nin% runned] # remaining individuals 


library(doParallel)  
n_cores <- detectCores() - 2  
cl <- makeCluster(n_cores, type="PSOCK", outfile="")  # for linux, use type="FORK" instead (reduce computing time)
registerDoParallel(cl)  
result <- foreach(i=torun, .packages=c("sl3", "origami"), .errorhandling = "pass") %dopar%{
  set.seed(seeds[i]);
  res <- run_csl(ind=i, data=db_train, outcome=outcome, type="continuous", covar=covar, stack_hist=stack_hist, stack_ind=stack_ind);
  write.csv2(res$pred, paste0("prediction_", i, ".csv"), row.names = FALSE);
  write.csv2(res$weights, paste0("weights_", i, ".csv"), row.names = FALSE);
  write.csv2(res$dsl_lrnr, paste0("dsl_", i, ".csv"), row.names = FALSE);
  write.csv2(res$weights_convex, paste0("wconvex_", i, ".csv"), row.names = FALSE);
  print(i)
} 
stopCluster(cl)

result # id 51 and 343 doesn't run due to too few sessions (cannot train individual learners)


