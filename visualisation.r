library(ggplot2)
library(ggpubr)

colind <- ggokabeito::palette_okabe_ito(1)
colhist <- ggokabeito::palette_okabe_ito(8)
colsl <- ggokabeito::palette_okabe_ito(7)

#####
# Results recovering
#####

#setwd("F:/cSL/Binary")
setwd("F:/cSL/Final/all sessions")

files <- list.files(getwd(),, pattern="prediction") 
files2 <- list.files(getwd(), pattern="weights") 
files3 <- list.files(getwd(),, pattern="wconvex") 
files4 <- list.files(getwd(), pattern="dsl_") 

extract <- function(files){
  pred <- c()
  for(f in files){
    
    db <- read.csv2(f, h=T)
    db$id <- f |> stringr::str_extract("\\d+")
    db$t <- 1:nrow(db)
    
    pred <- rbind(pred, db)
    
  }
  return(pred)
}

pred <- extract(files) # predictions
poid <- extract(files2) # non-convex eSl weights
pconvex <- extract(files3) # convex eSL weights
dsl <- extract(files4) # candidate learner selected in the dSL


#### 
# Data manadgement for ggplot2 
####

mef_poid <- function(poid, binary=FALSE){
  
  n_lrnr <- ncol(poid)-2
  
  poid2 <- sapply(colnames(poid)[1:n_lrnr],
                  function(x){
                    nnsum <- 100 * sum(poid[,x]>0) / length(poid[,x])
                    rwgt <- sum(poid[,x])
                    c(nnsum, rwgt) 
                  } ,
                  simplify=T
  ) |>  t() |> as.data.frame()
  

  poid2$Type <- c(rep("Individual ROCV", n_lrnr/3), rep("Individual RWCV", n_lrnr/3), rep("Historical", n_lrnr/3))
  
  if(binary){
    lrnrs <- c("logistic model", "MARS", "ridge", "XGBoost")
  }else{
    lrnrs <- c("logistic model", "MARS", "ridge", "lasso", "XGBoost")
  }
  
  poid2$Learner2 <- rep(c(lrnrs, paste("RF + ", lrnrs), "mean"),3)
  
  
  poid2$Learner2 <- paste(substr(poid2$Type, 12, 15), poid2$Learner2)
  
  poid2$Type <- substr(poid2$Type, 1, 10)
  
  return(poid2)
}

mef_dsl <- function(dsl, binary=FALSE){
  
  poid2 <- table(dsl$x) |> prop.table() |> as.data.frame()
  
  poid2$Type <- ifelse(grepl("hist", poid2$Var1), "Historical",
                       ifelse(grepl("Origin", poid2$Var1), "Individual ROCV",
                              "Individual RWCV"))
  
  

  poid2$Learner2 <- ifelse(grepl("lasso", poid2$Var1), 'lasso',
                           ifelse(grepl("lm", poid2$Var1), ifelse(binary, 'logistic model','linear model'),
                                  ifelse(grepl("mars", poid2$Var1), 'MARS',
                                         ifelse(grepl("mean", poid2$Var1), 'mean',
                                                ifelse(grepl("ridge", poid2$Var1), 'ridge',
                                                        'XGBoost'
                                                       )))))
  
  poid2$Learner2 <- ifelse(grepl("rf", poid2$Var1), paste('RF +', poid2$Learner2), poid2$Learner2) 
  
  poid2$Learner2 <- paste(substr(poid2$Type, 12, 15), poid2$Learner2)
  
  poid2$Type <- substr(poid2$Type, 1, 10)
  
  return(poid2)
}

mef_pred <- function(pred, metric='mse', cluster='id', outlier='NA', binary=FALSE){ #ape = mediane de l'erreur de prediction absolue
  
  n_lrnr <- ncol(pred)-2
  
  if(outlier=="bounds"){
    for(l in names(pred)[2:n_lrnr]){
      if(binary){
        pred[,l] <- ifelse(pred[,l]<0, 0, ifelse(pred[,l]>1, 1, pred[,l]))
      }else{
        pred[,l] <- ifelse(pred[,l]<0, 0, ifelse(pred[,l]>50, 50, pred[,l]))
      }
    }
  }else if(outlier=="exclude"){
    for(l in names(pred)[2:n_lrnr]){
      pred[,l] <- pred[pred[,l]>=0 & pred[,l]<=50,l]
    }
  }
  
  pred2 <- pred[pred$t>1,] 
pred2 <- sapply(unique(pred2[,cluster]),
                  function(i) sapply(names(pred2)[2:n_lrnr],
                                     function(l){
                                       if(metric=="mse"){
                                         mean((pred2$True[pred2[,cluster]==i] - pred2[pred2[,cluster]==i,l])^2)
                                       }else if(metric=='mcal'){
                                         if(binary){
                                           if(length(unique(pred2$True[pred2[,cluster]==i]))==1){
                                             NA
                                           }else{
                                             glm(pred2$True[pred2[,cluster]==i]~offset(pred2[pred2[,cluster]==i,l]), family=binomial)$coef
                                           }
                                           
                                         }else{
                                           lm(pred2$True[pred2[,cluster]==i]~offset(pred2[pred2[,cluster]==i,l]))$coef
                                         }
                                         
                                       }else if(metric=="wcal"){
                                         if(binary){
                                           if(length(unique(pred2$True[pred2[,cluster]==i]))==1){
                                             NA
                                           }else{
                                            glm(pred2$True[pred2[,cluster]==i]~pred2[pred2[,cluster]==i,l], family=binomial)$coef[2]
                                           }
                                         }else{
                                           lm(pred2$True[pred2[,cluster]==i]~pred2[pred2[,cluster]==i,l])$coef[2]
                                         }
                                       }else if(metric=="ape"){
                                         median(abs(pred2$True[pred2[,cluster]==i] - pred2[pred2[,cluster]==i,l]))
                                       }
                                       
                                       },
                                     simplify=T
                  ),
                  simplify=T) |> t() |> as.data.frame()
  
  
  pred2 <- reshape(pred2, direction="long",
                   varying = names(pred2),
                   v.names="pred",
                   times = names(pred2),
                   timevar="Learner")
  
  pred2$CV <- ifelse(pred2$Learner %in% unique(grep("Origin", pred2$Learner, value=TRUE)) , "Rolling-Origin" , ifelse(pred2$Learner %in% unique(grep("Window", pred2$Learner, value=TRUE)) , "Rolling-Window" , "None"))
  pred2$Type <- ifelse(pred2$Learner %in% unique(grep("ind", pred2$Learner, value=TRUE)) , "Individual" , ifelse(pred2$Learner %in% unique(grep("hist", pred2$Learner, value=TRUE)) , "Historical", "POSL"))
  pred2$Type[pred2$Type=="Individual"] <- ifelse(pred2$CV[pred2$Type=="Individual"]=='Rolling-Origin', "Individual", "Individual (RWCV)")
  
  
  pred2$Learner2 <- ifelse(pred2$Learner %in% unique(grep("lm_rf", pred2$Learner, value=TRUE)), 'RF + linear model', 
                           ifelse(pred2$Learner %in% unique(grep("lm", pred2$Learner, value=TRUE)), 'linear model', ifelse(pred2$Learner %in% unique(grep("mars_rf", pred2$Learner, value=TRUE)), 'RF + MARS', ifelse(pred2$Learner %in% unique(grep("mars", pred2$Learner, value=TRUE)), 'MARS', ifelse(pred2$Learner %in% unique(grep("ridge_rf", pred2$Learner, value=TRUE)), 'RF + ridge', ifelse(pred2$Learner %in% unique(grep("ridge", pred2$Learner, value=TRUE)), 'ridge', ifelse(pred2$Learner %in% unique(grep("lasso_rf", pred2$Learner, value=TRUE)), 'RF + lasso', ifelse(pred2$Learner %in% unique(grep("lasso", pred2$Learner, value=TRUE)), 'lasso', ifelse(pred2$Learner %in% unique(grep("xgb_rf", pred2$Learner, value=TRUE)), 'RF + XGBoost', ifelse(pred2$Learner %in% unique(grep("xgb", pred2$Learner, value=TRUE)), 'XGBoost', ifelse(pred2$Learner %in% unique(grep("mean", pred2$Learner, value=TRUE)), 'mean', ifelse(pred2$Learner %in% unique(grep("dSL", pred2$Learner, value=TRUE)), 'dSL', ifelse(pred2$Learner %in% unique(grep("convex", pred2$Learner, value=TRUE)), 'convex eSL', 'non-convex eSL')        ))))))))))))
  

  return(pred2)
}

mef_pred_sd <- function(pred, metric='mse', cluster='id', outlier='NA', binary=FALSE){
  
  n_lrnr <- ncol(pred)-2
  
  if(outlier=="bounds"){
    for(l in names(pred)[2:n_lrnr]){
      if(binary){
        pred[,l] <- ifelse(pred[,l]<0, 0, ifelse(pred[,l]>1, 1, pred[,l]))
      }else{
        pred[,l] <- ifelse(pred[,l]<0, 0, ifelse(pred[,l]>50, 50, pred[,l]))
      }
    }
  }else if(outlier=="exclude"){
    for(l in names(pred)[2:n_lrnr]){
      pred[,l] <- pred[pred[,l]>=0 & pred[,l]<=50,l]
    }
  }
  
  pred2 <- pred[pred$t>1,]
  pred2 <- sapply(unique(pred2[,cluster]),
                  function(i){
                    print(i)
                    if(nrow(pred2[pred2[,cluster]==i,])==1) rep(NA, length(names(pred2)[2:n_lrnr])) else
                    sapply(names(pred2)[2:n_lrnr],
                                     function(l){
                                       print(l)
                                       if(metric=="mse"){
                                         sd((pred2$True[pred2[,cluster]==i] - pred2[pred2[,cluster]==i,l])^2)
                                       }else if(metric=='mcal'){
                                          if(binary){
                                           if(length(unique(pred2$True[pred2[,cluster]==i]))==1){
                                             NA
                                           }else{
                                             summary(glm(pred2$True[pred2[,cluster]==i]~offset(pred2[pred2[,cluster]==i,l]), family=binomial))$coef[2]
                                           }
                                         }else{
                                           summary(lm(pred2$True[pred2[,cluster]==i]~offset(pred2[pred2[,cluster]==i,l])))$coef[2]
                                         }
                                       }else if(metric=="wcal"){
                                         if(l=='hist_mean') NA else{
                                           if(binary){
                                           if(length(unique(pred2$True[pred2[,cluster]==i]))==1){
                                             NA
                                           }else{
                                             summary(glm(pred2$True[pred2[,cluster]==i]~offset(pred2[pred2[,cluster]==i,l]), family=binomial))$coef[2]
                                           }
                                         }else{
                                           summary(lm(pred2$True[pred2[,cluster]==i]~pred2[pred2[,cluster]==i,l]))$coef[2,2]
                                         }
                                         }
                                       }else if(metric=="ape"){
                                         sd(abs(pred2$True[pred2[,cluster]==i] - pred2[pred2[,cluster]==i,l]))
                                       }
                                      
                                     },
                                     simplify=T
                  )},
                  simplify=T)  |> t() |> as.data.frame()
  
  
  pred2 <- reshape(pred2, direction="long",
                   varying = names(pred2),
                   v.names="pred",
                   times = names(pred2),
                   timevar="Learner")
  
  pred2$CV <- ifelse(pred2$Learner %in% unique(grep("Origin", pred2$Learner, value=TRUE)) , "Rolling-Origin" , ifelse(pred2$Learner %in% unique(grep("Window", pred2$Learner, value=TRUE)) , "Rolling-Window" , "None"))
  pred2$Type <- ifelse(pred2$Learner %in% unique(grep("ind", pred2$Learner, value=TRUE)) , "Individual" , ifelse(pred2$Learner %in% unique(grep("hist", pred2$Learner, value=TRUE)) , "Historical", "POSL"))
  pred2$Type[pred2$Type=="Individual"] <- ifelse(pred2$CV[pred2$Type=="Individual"]=='Rolling-Origin', "Individual", "Individual (RWCV)")
  
  
  pred2$Learner2 <- ifelse(pred2$Learner %in% unique(grep("lm_rf", pred2$Learner, value=TRUE)), ifelse(binary, 'RF + logistic model','RF + linear model'), 
                           ifelse(pred2$Learner %in% unique(grep("lm", pred2$Learner, value=TRUE)), ifelse(binary, 'logistic model','linear model'), ifelse(pred2$Learner %in% unique(grep("mars_rf", pred2$Learner, value=TRUE)), 'RF + MARS', ifelse(pred2$Learner %in% unique(grep("mars", pred2$Learner, value=TRUE)), 'MARS', ifelse(pred2$Learner %in% unique(grep("ridge_rf", pred2$Learner, value=TRUE)), 'RF + ridge', ifelse(pred2$Learner %in% unique(grep("ridge", pred2$Learner, value=TRUE)), 'ridge', ifelse(pred2$Learner %in% unique(grep("lasso_rf", pred2$Learner, value=TRUE)), 'RF + lasso', ifelse(pred2$Learner %in% unique(grep("lasso", pred2$Learner, value=TRUE)), 'lasso', ifelse(pred2$Learner %in% unique(grep("xgb_rf", pred2$Learner, value=TRUE)), 'RF + XGBoost', ifelse(pred2$Learner %in% unique(grep("xgb", pred2$Learner, value=TRUE)), 'XGBoost', ifelse(pred2$Learner %in% unique(grep("mean", pred2$Learner, value=TRUE)), 'mean', ifelse(pred2$Learner %in% unique(grep("dSL", pred2$Learner, value=TRUE)), 'dSL', ifelse(pred2$Learner %in% unique(grep("convex", pred2$Learner, value=TRUE)), 'convex eSL', 'non-convex eSL')        ))))))))))))
  
  
  return(pred2)
}


mef_ccurv <- function(pred, cluster='id', outlier='NA', binary=FALSE){
  
  n_lrnr <- ncol(pred)-2
  
  if(outlier=="bounds"){
    for(l in names(pred)[2:n_lrnr]){
      pred[,l] <- ifelse(pred[,l]<8, 8, ifelse(pred[,l]>40, 40, pred[,l]))
    }
  }else if(outlier=="exclude"){
    for(l in names(pred)[2:n_lrnr]){
      pred[,l] <- pred[pred[,l]>=0 & pred[,l]<=50,l]
    }
  }
  
  pred2 <- pred[pred$t>1,]
  
  
  pred2 <- reshape(pred2, direction="long",
                   varying = names(pred2)[2:n_lrnr], #37 instead of 38:39
                   v.names="pred",
                   times = names(pred2)[2:n_lrnr],
                   timevar="Learner", 
                   idvar=c("id", "t"))
  
  
  pred2$CV <- ifelse(pred2$Learner %in% unique(grep("Origin", pred2$Learner, value=TRUE)) , "Rolling-Origin" , ifelse(pred2$Learner %in% unique(grep("Window", pred2$Learner, value=TRUE)) , "Rolling-Window" , "None"))
  pred2$Type <- ifelse(pred2$Learner %in% unique(grep("ind", pred2$Learner, value=TRUE)) , "Individual" , ifelse(pred2$Learner %in% unique(grep("hist", pred2$Learner, value=TRUE)) , "Historical", "POSL"))
  pred2$Type[pred2$Type=="Individual"] <- ifelse(pred2$CV[pred2$Type=="Individual"]=='Rolling-Origin', "Individual", "Individual (RWCV)")
  
  
  pred2$Learner2 <- ifelse(pred2$Learner %in% unique(grep("lm_rf", pred2$Learner, value=TRUE)), ifelse(binary, 'RF + logistic model','RF + linear model'), 
                           ifelse(pred2$Learner %in% unique(grep("lm", pred2$Learner, value=TRUE)), ifelse(binary, 'logistic model','linear model'), ifelse(pred2$Learner %in% unique(grep("mars_rf", pred2$Learner, value=TRUE)), 'RF + MARS', ifelse(pred2$Learner %in% unique(grep("mars", pred2$Learner, value=TRUE)), 'MARS', ifelse(pred2$Learner %in% unique(grep("ridge_rf", pred2$Learner, value=TRUE)), 'RF + ridge', ifelse(pred2$Learner %in% unique(grep("ridge", pred2$Learner, value=TRUE)), 'ridge', ifelse(pred2$Learner %in% unique(grep("lasso_rf", pred2$Learner, value=TRUE)), 'RF + lasso', ifelse(pred2$Learner %in% unique(grep("lasso", pred2$Learner, value=TRUE)), 'lasso', ifelse(pred2$Learner %in% unique(grep("xgb_rf", pred2$Learner, value=TRUE)), 'RF + XGBoost', ifelse(pred2$Learner %in% unique(grep("xgb", pred2$Learner, value=TRUE)), 'XGBoost', ifelse(pred2$Learner %in% unique(grep("mean", pred2$Learner, value=TRUE)), 'mean', ifelse(pred2$Learner %in% unique(grep("convex", pred2$Learner, value=TRUE)), 'convex eSL', ifelse(pred2$Learner %in% unique(grep("eSL", pred2$Learner, value=TRUE)), 'non-convex eSL','dSL')))        )))))))))) 
  
  return(pred2)
}


####
# Plotting function
####



plot_poid <- function(poid2, dsl2, ABC="AUTO"){

  pperc <- ggpubr::ggdotchart(data=poid2, x="Learner2", y="V1",
                              color = "Type", fill="Type",                               # Color by groups
                              palette = c(colhist, colind), # Custom color palette
                              sorting = "descending",                       # Sort value in descending order
                              rotate = TRUE,                                # Rotate vertically
                              dot.size = 2,                                 # Large dot size
                              ggtheme = theme_pubr(legend=c(0.8, 0.2)),                       # ggplot2 theme
                              ylab="Presence in eSL (%)", xlab="", legend.title=''
  ) + font("xy.title", face = "bold") + geom_segment(aes(y=0, yend=sort(poid2$V1), xend=poid2$Learner2[order(poid2$V1)], color=poid2$Type[order(poid2$V1)])) 
  
  
  pcum <- ggpubr::ggdotchart(data=poid2, x="Learner2", y="V2",
                             color = "Type", fill="Type",                               # Color by groups
                             #add="segment",
                             palette = c(colhist, colind), # Custom color palette
                             sorting = "descending",                       # Sort value in descending order
                             rotate = TRUE,                                # Rotate vertically
                             dot.size = 2,                                 # Large dot size
                             ggtheme = theme_pubr(legend=c(0.8, 0.2)),                       # ggplot2 theme
                             ylab="Cumulative weights in eSL", xlab="", legend.title=''
  ) + font("xy.title", face = "bold") + geom_segment(aes(y=0, yend=sort(poid2$V2), xend=poid2$Learner2[order(poid2$V2)], color=poid2$Type[order(poid2$V2)])) 
  
  pdsl <- ggpubr::ggdotchart(data=dsl2, x="Learner2", y="Freq",
                             color = "Type", fill="Type",                               # Color by groups
                             #add="segment",
                             palette = c(colhist, colind), # Custom color palette
                             sorting = "descending",                       # Sort value in descending order
                             rotate = TRUE,                                # Rotate vertically
                             dot.size = 2,                                 # Large dot size
                             ggtheme = theme_pubr(legend=c(0.8, 0.2)),                       # ggplot2 theme
                             ylab="Presence in dSL (%)", xlab="", legend.title=''
  ) + font("xy.title", face = "bold") + geom_segment(aes(y=0, yend=sort(dsl2$Freq), xend=dsl2$Learner2[order(dsl2$Freq)], color=dsl2$Type[order(dsl2$Freq)])) 
  
  
  ggarrange(pperc, pcum, pdsl, ncol=3, common.legend = T, labels=ABC)

}

plot_pred <- function(pred2, metric='mse', ylim=c(0, 100)){
  
  yint <- ifelse(metric=='mse', median(pred2[pred2$Learner2=='convex eSL', metric]), 
                 ifelse(metric=='mcal', 0, ifelse(metric=='ape', median(pred2[pred2$Learner2=='non-convex eSL', metric], 1), 1))
                 )
  
  ylab <- ifelse(metric=='mse', 'MSE', 
                 ifelse(metric=='mcal', 'Calibration intercept', ifelse(metric=='ape', 'MdAE', 'Calibration slope'))
                 )
  
  ggboxplot(pred2, x = 'Learner2', y = metric, fill = 'Type', color = 'Type', alpha=0.5,orientation="horiz",
            ggtheme = theme_pubr(legend=c(0.9,0.9), x.text.angle = 45), ylab=ylab, xlab="", legend.title='',  outlier.shape = NA, 
            palette = c(colhist, colind, colsl)
  ) + font("xy.title", face = "bold") + geom_hline(yintercept=yint, linetype="dashed") + coord_cartesian(ylim = ylim)
  
}

plot_pred_evol <- function(dbt, metric, col=c("#56B4E9", "#CC79A7", "#0072B2"), smooth=FALSE, colblind='n', se=TRUE){
  
  ylab <- ifelse(metric=='mse', 'MSE', 
                 ifelse(metric=='mcal', 'Calibration intercept', 
                        ifelse(metric=='ape', 'MdAE', 'Calibration slope'))
  )
  
  col <- switch(colblind,
                'deu' = colorspace::deutan(col),
                'pro' = colorspace::protan(col),
                'tri' = colorspace::tritan(col), 
                col
                )
  
  if(metric=='ape'){
    g <- ggplot(data=dbt, aes(y=ape, x=id))
  }else if(metric=='mse'){
    g <- ggplot(data=dbt, aes(y=mse, x=id))
  }else if(metric=='mcal'){
    g <- ggplot(data=dbt, aes(y=mcal, x=id))  +
      geom_hline(yintercept=0, linetype="dashed")
  }else{
    g <- ggplot(data=dbt, aes(y=wcal, x=id)) +
      geom_hline(yintercept=1, linetype="dashed")
  } 
  
  if(smooth){
    g <- g + geom_smooth(aes(color=Learner2, fill=Learner2), size=1, se=se, method='loess', alpha=0.25) +
      scale_fill_manual(values=col)
  }else{
    g <- g +
      geom_line(aes(color=Learner2), linewidth=.75)
  } 
  
  g <- g +
    theme_pubr() +
    scale_color_manual(values=col) +
    xlab('Time') + ylab(ylab) + font("xy.title", face = "bold")
  
   ggpubr::ggpar(g, legend.title="")
}


plot_poid_evol <- function(poid, ABC="AUTO"){
  
  #cumulative weights
  poid2ind <- sapply(unique(poid$t),
         function(u) apply(poid[poid$t==u,grep("ind_", names(poid), value=T)], 1, function(x) sum(x)) |> unname()
  ) |> lapply(mean) |> unlist()
  
  poid2hist <- sapply(unique(poid$t),
                    function(u) apply(poid[poid$t==u,grep("hist_", names(poid), value=T)], 1, function(x) sum(x)) |> unname()
  ) |> lapply(mean) |> unlist()
  
  #number
  nnind <- sapply(unique(poid$t),
         function(u) apply(poid[poid$t==u,grep("ind_", names(poid), value=T)], 1, function(x) sum(x>0)) |> unname()
  ) |> lapply(mean) |> unlist()
  
  nnhist <- sapply(unique(poid$t),
                  function(u) apply(poid[poid$t==u,grep("hist_", names(poid), value=T)], 1, function(x) sum(x>0)) |> unname()
  ) |> lapply(mean) |> unlist()
  
  poid2 <- data.frame(
    t=c(1:length(poid2ind),1:length(poid2hist)),
    cum=c(poid2ind,poid2hist),
    numb=c(nnind, nnhist),
    Learner=c(rep("Individual", length(poid2ind)), rep("Historical", length(poid2hist)))
  )
  
  
  gcum <- ggplot(data=poid2, aes(y=cum, x=t)) +
    geom_line(aes(color=Learner), linewidth=.75) +
    theme_pubr() +
    scale_color_manual(values=c(colhist,colind)) +
    xlab('Time') + ylab("Learners' cumulative weights") + font("xy.title", face = "bold") +
    theme(legend.position = "none")
  
  gnumb <- ggplot(data=poid2, aes(y=numb, x=t)) +
    geom_line(aes(color=Learner), linewidth=.75) +
    theme_pubr() +
    scale_color_manual(values=c(colhist,colind)) +
    xlab('Time') + ylab("Number of learners") + font("xy.title", face = "bold") + 
    theme(legend.position = "none")
  
  gnumb <- ggpubr::ggpar(gnumb, legend.title="")
  gcum <- ggpubr::ggpar(gcum, legend.title="")
  
  
  ggpubr::ggarrange(gnumb, gcum, ncol=2, labels=ABC)
  
  
  
}

####
## Results analysis
####

### Figure 4

nonconvex <- mef_poid(poid)
convex <- mef_poid(pconvex)

dsl_gg <- mef_dsl(dsl)
dsl_gg$Freq <- dsl_gg$Freq*100

plot_poid(convex, dsl_gg, ABC="AUTO") #"auto" for lowercase letters
p1 <- plot_poid(nonconvex, dsl_gg) # almost identical between nonconvex and convex (inversion ROCv rf + linear model & ROCv rf + mars)

p2 <- plot_poid_evol(poid[poid$t<=500,], ABC=c("D", "E"))

ggarrange(p1,p2, ncol=1, common.legend = T, heights = c(2, 1))


### Figure 5

db <- cbind(mef_pred(pred, 'mse', outlier='bounds'), mef_pred(pred, 'mcal', outlier='bounds')[,2], mef_pred(pred, 'wcal', outlier='bounds')[,2], mef_pred(pred, 'ape', outlier='bounds')[,2],mef_pred_sd(pred, 'mse', outlier='bounds')[,2], mef_pred_sd(pred, 'mcal', outlier='bounds')[,2], mef_pred_sd(pred, 'wcal', outlier='bounds')[,2], mef_pred_sd(pred, 'ape', outlier='bounds')[,2]) 
names(db)[c(2, 7:13)] <- c('mse','mcal', 'wcal', 'ape', 'se.se','mcal.se', 'wcal.se', 'ape.se')


pmse <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='mse', ylim=c(0,75))
pape <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='ape', ylim=c(0,5))
pmcal <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='mcal', ylim=c(-2,2))
pwcal <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='wcal', ylim=c(-0.5,2.5))

ggarrange(pape, pmse, pmcal, pwcal, ncol=2, nrow=2, common.legend = TRUE, labels="AUTO")


### Figure 6

dbt <- cbind(mef_pred(pred, 'mse', cluster='t', outlier='bounds'), mef_pred(pred, 'mcal', cluster='t', outlier='bounds')[,2], mef_pred(pred, 'wcal', cluster='t', outlier='bounds')[,2], mef_pred(pred, 'ape', cluster='t', outlier='bounds')[,2],mef_pred_sd(pred, 'mse', cluster='t', outlier='bounds')[,2], mef_pred_sd(pred, 'mcal', cluster='t', outlier='bounds')[,2], mef_pred_sd(pred, 'wcal', cluster='t', outlier='bounds')[,2], mef_pred_sd(pred, 'ape', cluster='t', outlier='bounds')[,2]) 
names(dbt)[c(2, 7:13)] <- c('mse','mcal', 'wcal', 'ape', 'se.se','mcal.se', 'wcal.se', 'ape.se')


pt1 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_lm_OriginCV')) & dbt$id<=500,], metric='ape', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt1b <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_lm_OriginCV')) & dbt$id<=500,], metric='mse', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt2 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_lm_OriginCV')) & dbt$id<=500,], metric='mcal', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt3 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_lm_OriginCV')) & dbt$id<=500,], metric='wcal', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

ptres <- ggarrange(pt1, pt1b, pt2, pt3, ncol=2, nrow=2, common.legend = TRUE, labels = "AUTO")



### Supplementary Figure 1


pred2 <- mef_ccurv(pred, outlier='bounds')
pred2$Learner2 <- factor(pred2$Learner2, levels=c("convex eSL", "non-convex eSL", "dSL", "linear model", "MARS", "ridge", "lasso", "XGBoost", "RF + linear model", "RF + MARS", "RF + ridge", "RF + lasso", "RF + XGBoost", "mean"))


g <- ggplot(data=pred2[pred2$CV!="Rolling-Window" & pred2$Learner2!="mean",])  + geom_segment(aes(x=8,y=8,xend=40,yend=40), col="black", size=.5, linetype=2) + geom_smooth(aes(x = pred, y = True, color = Type, fill=Type), method = "gam", size = 1, se=T) +
 facet_wrap(Learner2 ~ ., ncol=5) + scale_color_manual(values=c(colhist, colind, colsl)) + ggpubr::font("xy.title", face = "bold") + ggpubr::theme_pubr(legend=c(0.75,0.1)) + scale_fill_manual(values=c(colhist, colind, colsl))
ggpubr::ggpar(g, legend.title="") + coord_cartesian(xlim = c(8, 40), ylim=c(8,40)) + xlab(expression(bold("Predictions"))) + ylab(expression(bold("Observations")))


pred2 <- rbind(
            mef_ccurv(pred[pred$t==2,], outlier='bounds'),
            mef_ccurv(pred[pred$t==102,], outlier='bounds'),
            mef_ccurv(pred[pred$t==202,], outlier='bounds'),
            mef_ccurv(pred[pred$t==352,], outlier='bounds'),
            mef_ccurv(pred[pred$t==502,], outlier='bounds')
          )

### Figure 5, panel E

pred2$t <- as.factor(pred2$t)

levels(pred2$t) <- c("Time 1", "Time 100", "Time 200", "Time 350", "Time 500")

g <- ggplot(data=pred2[pred2$Type=="POSL",])  + geom_segment(aes(x=8,y=8,xend=40,yend=40), col="black", size=.5, linetype=2) + geom_smooth(aes(x = pred, y = True, color = Learner2, fill=Learner2), method = "gam", size = 1, se=T) +
  xlab("Predictions") + ylab("Observations") + facet_grid(cols=vars(t)) + scale_color_manual(values=ggokabeito::palette_okabe_ito(c(2:3,5))) + ggpubr::font("xy.title", face = "bold") + ggpubr::theme_pubr(legend="none") + scale_fill_manual(values=ggokabeito::palette_okabe_ito(c(2:3,5)))
tflex <- ggpubr::ggpar(g, legend.title="") + coord_cartesian(xlim = c(8, 40), ylim=c(8,40)) + xlab(expression(bold("Predictions"))) + ylab(expression(bold("Observations")))

ggarrange(ptres,tflex, labels=c("","E"), heights = c(2, 1), ncol = 1)


### Supplementary Figure 2 

db_auc <- pred[pred$t>1,]

db_auc$btrue <- 1*(db_auc$True>=24)

aucs <- sapply(names(db_auc)[c(2:14, 26:36)],
               function(x) pROC::ci.auc(pROC::roc(db_auc$btrue, db_auc[,x], percent=F, ci=T)),
               simplify=T) |> round(digits = 3)


pred2 <- reshape(data.frame(aucs), direction="long",
                 varying = colnames(aucs),
                 v.names="pred",
                 times = colnames(aucs),
                 timevar="Learner")

pred2$Type <- ifelse(pred2$Learner %in% unique(grep("ind", pred2$Learner, value=TRUE)) , "Individual" , ifelse(pred2$Learner %in% unique(grep("hist", pred2$Learner, value=TRUE)) , "Historical", "POSL"))
pred2$CV <- ifelse(pred2$Learner %in% unique(grep("Origin", pred2$Learner, value=TRUE)) , "Rolling-Origin" , ifelse(pred2$Learner %in% unique(grep("Window", pred2$Learner, value=TRUE)) , "Rolling-Window" , "None"))
pred2$Type[pred2$Type=="Individual"] <- ifelse(pred2$CV[pred2$Type=="Individual"]=='Rolling-Origin', "Individual", "Individual (RWCV)")


pred2$Learner2 <- ifelse(pred2$Learner %in% unique(grep("lm_rf", pred2$Learner, value=TRUE)), 'RF + linear model', 
                         ifelse(pred2$Learner %in% unique(grep("lm", pred2$Learner, value=TRUE)), 'linear model', ifelse(pred2$Learner %in% unique(grep("mars_rf", pred2$Learner, value=TRUE)), 'RF + MARS', ifelse(pred2$Learner %in% unique(grep("mars", pred2$Learner, value=TRUE)), 'MARS', ifelse(pred2$Learner %in% unique(grep("ridge_rf", pred2$Learner, value=TRUE)), 'RF + ridge', ifelse(pred2$Learner %in% unique(grep("ridge", pred2$Learner, value=TRUE)), 'ridge', ifelse(pred2$Learner %in% unique(grep("lasso_rf", pred2$Learner, value=TRUE)), 'RF + lasso', ifelse(pred2$Learner %in% unique(grep("lasso", pred2$Learner, value=TRUE)), 'lasso', ifelse(pred2$Learner %in% unique(grep("xgb_rf", pred2$Learner, value=TRUE)), 'RF + XGBoost', ifelse(pred2$Learner %in% unique(grep("xgb", pred2$Learner, value=TRUE)), 'XGBoost', ifelse(pred2$Learner %in% unique(grep("mean", pred2$Learner, value=TRUE)), 'mean', ifelse(pred2$Learner %in% unique(grep("dSL", pred2$Learner, value=TRUE)), 'dSL', ifelse(pred2$Learner %in% unique(grep("convex", pred2$Learner, value=TRUE)), 'convex eSL', "non-convex eSL")))        ))))))))))  

pred3 <- pred2[pred2$id==1,]   
pred4 <- pred2[pred2$id==3,]   
pred2 <- pred2[pred2$id==2,] 

ggpubr::ggdotchart(data=pred2[pred2$CV!='Rolling-Window',], x = 'Learner2', y = 'pred', color = 'Type', 
                           palette = c(colhist, colind, colsl),
                           sorting = "none",  order=rev(unique(pred2$Learner2)), 
                           rotate = TRUE, shape=18,                         
                           dot.size = 5, ylab="AUROC", xlab="", 
                           ggtheme = theme_pubr(),
                           legend.title = "") + font("xy.title", face = "bold") +  geom_segment(aes(y=pred3$pred[1], yend=pred4$pred[1], x=13, xend=13), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[2], yend=pred4$pred[2], x=12, xend=12), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[3], yend=pred4$pred[3], x=11, xend=11), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[4], yend=pred4$pred[4], x=10, xend=10), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[5], yend=pred4$pred[5], x=9, xend=9), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[6], yend=pred4$pred[6], x=8, xend=8), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[7], yend=pred4$pred[7], x=7, xend=7), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[8], yend=pred4$pred[8], x=6, xend=6), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[9], yend=pred4$pred[9], x=5, xend=5), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[10], yend=pred4$pred[10], x=4, xend=4), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[11], yend=pred4$pred[11], x=3, xend=3), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[12], yend=pred4$pred[12], x=2, xend=2), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[13], yend=pred4$pred[13], x=1, xend=1), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[15], yend=pred4$pred[15], x=10, xend=10), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[16], yend=pred4$pred[16], x=9, xend=9), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[17], yend=pred4$pred[17], x=8, xend=8), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[18], yend=pred4$pred[18], x=7, xend=7), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[19], yend=pred4$pred[19], x=6, xend=6), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[20], yend=pred4$pred[20], x=5, xend=5), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[21], yend=pred4$pred[21], x=4, xend=4), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[22], yend=pred4$pred[22], x=3, xend=3), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[23], yend=pred4$pred[23], x=2, xend=2), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[24], yend=pred4$pred[24], x=1, xend=1), col=colhist, size=1, linetype=1)


### Supplementary Figure 3

pred2 <- pred[pred$t>1,]

names(pred2)[c(2:5,31)] <- c('convex eSL', 'non-convex eSL', 'dSL', 'linear model', 'XGBoost')

benef <- function(x, true, cut, db, all=F){
  tp <- sum(db[,true]>=cut & db[,x]>=cut)
  fp <- sum(db[,true]<cut & db[,x]>=cut)
  p <- prop.table(table(db[,true]>=cut))[2]
  
  net_b <- tp/nrow(db) - fp/nrow(db) * p/(1-p)
  
  if(all) net_b <- sum(db[,true]>=cut)/nrow(db) - sum(db[,true]<cut)/nrow(db) * p/(1-p)
  
  return(net_b)
}

nb <- sapply(names(pred2)[c(2:5,31)], function(x) sapply(10:40, function(y)  benef(x=x, true='True', cut=y, db=pred2), simplify=T)
       , simplify=T) |> as.data.frame()


nb$cut <- 10:40


dca <- reshape(nb, direction="long",
                 varying = names(nb)[1:5],
                 v.names="benef",
                 times = names(nb)[1:5],
                 timevar="Learner")

d <- ggplot(data=dca, aes(x = cut, y = benef, color = Learner)) +
  geom_line(size=1) + 
  xlab("Predicted convection volume (L)") + ylab("Net benefit") + 
  scale_color_manual(values=ggokabeito::palette_okabe_ito(c(2:3,1,5,8))) + ggpubr::font("xy.title", face = "bold") + ggpubr::theme_pubr() + scale_fill_manual(values=ggokabeito::palette_okabe_ito(c(2:3,5))) + scale_x_continuous(breaks=seq(10,40,2))
ggpubr::ggpar(d, legend.title="") + ggpubr::font("xy.title", face = "bold")

  

####
# Summary table
####

tbl_res <- function(pred2, metric){
  res <- sapply(unique(pred2$Learner2),
                    function(l){
                      if(length(grep("SL", l)) > 0){
                        s <- summary(pred2[pred2$Learner2==l & pred2$CV!="Rolling-Window", metric])
                        se <- var(pred2[pred2$Learner2==l & pred2$CV!="Rolling-Window", metric], na.rm=T)
                        output <- c(s[c(3,4)], se, s[c(1,6)]) 
                      }else{
                        s <- summary(pred2[pred2$Learner2==l & pred2$Type=="Historical", metric])
                        se <- var(pred2[pred2$Learner2==l & pred2$Type=="Historical", metric], na.rm=T)
                        s2 <- summary(pred2[pred2$Learner2==l & pred2$Type=="Individual", metric])
                        se2 <- var(pred2[pred2$Learner2==l & pred2$Type=="Individual", metric], na.rm=T)
                        output <- rbind(c(s[c(3,4)], se, s[c(1,6)]), c(s2[c(3,4)], se2, s2[c(1,6)])) 
                      }
                      return(output ) 
                    },
                    simplify=T ) 
  
  rnames <- names(res)[-(1:3)] 
  rnames <- c('convex eSL', 'non-convex eSL', 'dSL', paste(c("Historical", "Individual"), sapply(rnames, function(n) rep(n,2),simplify=F) |> unlist()))
  
  res <- res |> do.call(what=rbind) |> round(digits=2) |> data.frame(row.names = rnames) |> setNames(nm=c('Median', 'Mean', 'Variance', 'Minimum', 'Maximum')) 
  
  return(res)
}

res.ape <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'ape')
res.mse <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'mse')
res.mcal <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'mcal')
res.wcal <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'wcal')


res <- merge(res.ape, res.mse, by = 'row.names', sort=F) |> merge(y=res.mcal, by.y = 'row.names', by.x='Row.names', sort=F) |> merge(y=res.wcal, by.y = 'row.names', by.x='Row.names', sort=F)

data.table::setorder(res, 'Median.x')

resf <- rbind(res[,1:11], res[,c(1,12:21)])

xtable::xtable(resf)


res.apet <- tbl_res(dbt[dbt$Learner2!='mean',], 'ape')
res.mset <- tbl_res(dbt[dbt$Learner2!='mean',], 'mse')
res.mcalt <- tbl_res(dbt[dbt$Learner2!='mean',], 'mcal')
res.wcalt <- tbl_res(dbt[dbt$Learner2!='mean',], 'wcal')


rest <- merge(res.apet, res.mset, by = 'row.names', sort=F) |> merge(y=res.mcalt, by.y = 'row.names', by.x='Row.names', sort=F) |> merge(y=res.wcalt, by.y = 'row.names', by.x='Row.names', sort=F)

data.table::setorder(rest, 'Median.x')

restf <- rbind(rest[,1:11], rest[,c(1,12:21)])

xtable::xtable(restf)






####
## Results Binary outcome (sensitivity analysis)
####

### Supplementary Figure 4

nonconvex <- mef_poid(poid, binary=TRUE)
convex <- mef_poid(pconvex, binary=TRUE)

dsl_gg <- mef_dsl(dsl, binary=TRUE)
dsl_gg$Freq <- dsl_gg$Freq*100

p1 <- plot_poid(convex, dsl_gg, ABC="AUTO") #"auto" for lowercase letters
plot_poid(nonconvex, dsl_gg) 

p2 <- plot_poid_evol(pconvex[pconvex$t<=500,], ABC=c("D", "E"))

ggarrange(p1,p2, ncol=1, common.legend = T, heights = c(2, 1))


###  Supplementary Figure 5

db <- cbind(mef_pred(pred, 'mse', outlier='bounds', binary=TRUE), mef_pred(pred, 'mcal', outlier='bounds', binary=TRUE)[,2], mef_pred(pred, 'wcal', outlier='bounds', binary=TRUE)[,2], mef_pred(pred, 'ape', outlier='bounds', binary=TRUE)[,2],mef_pred_sd(pred, 'mse', outlier='bounds', binary=TRUE)[,2], mef_pred_sd(pred, 'mcal', outlier='bounds', binary=TRUE)[,2], mef_pred_sd(pred, 'wcal', outlier='bounds', binary=TRUE)[,2], mef_pred_sd(pred, 'ape', outlier='bounds', binary=TRUE)[,2]) 
names(db)[c(2, 7:13)] <- c('mse','mcal', 'wcal', 'ape', 'se.se','mcal.se', 'wcal.se', 'ape.se')

db$Learner2[db$Learner2=="linear model"] <- "logistic model"
db$Learner2[db$Learner2=="RF + linear model"] <- "RF + logistic model"


pmse <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='mse', ylim=c(0,0.5))
pape <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='ape', ylim=c(0,0.75))
pmcal <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='mcal', ylim=c(-3,3))
pwcal <- plot_pred(db[db$CV!="Rolling-Window" & db$Learner2!="mean",], metric='wcal', ylim=c(-1,20))

ggarrange(pape, pmse, pmcal, pwcal, ncol=2, nrow=2, common.legend = TRUE, labels="AUTO")


### Supplementary Figure 6

dbt <- cbind(mef_pred(pred, 'mse', cluster='t', outlier='bounds', binary = TRUE), mef_pred(pred, 'mcal', cluster='t', outlier='bounds', binary = TRUE)[,2], mef_pred(pred, 'wcal', cluster='t', outlier='bounds', binary = TRUE)[,2], mef_pred(pred, 'ape', cluster='t', outlier='bounds', binary = TRUE)[,2],mef_pred_sd(pred, 'mse', cluster='t', outlier='bounds', binary = TRUE)[,2], mef_pred_sd(pred, 'mcal', cluster='t', outlier='bounds', binary = TRUE)[,2], mef_pred_sd(pred, 'wcal', cluster='t', outlier='bounds', binary = TRUE)[,2], mef_pred_sd(pred, 'ape', cluster='t', outlier='bounds', binary = TRUE)[,2]) 
names(dbt)[c(2, 7:13)] <- c('mse','mcal', 'wcal', 'ape', 'se.se','mcal.se', 'wcal.se', 'ape.se')

dbt$Learner2[dbt$Learner2=="linear model"] <- "logistic model"
dbt$Learner2[dbt$Learner2=="RF + linear model"] <- "RF + logistic model"



pt1 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_glm_OriginCV')) & dbt$id<=500,], metric='ape', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt1b <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_glm_OriginCV')) & dbt$id<=500,], metric='mse', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt2 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_glm_OriginCV')) & dbt$id<=500,], metric='mcal', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

pt3 <- plot_pred_evol(dbt[(dbt$Type=="POSL" | dbt$Learner%in%c('hist_xgb', 'ind_glm_OriginCV')) & dbt$id<=500,], metric='wcal', smooth=T, col=ggokabeito::palette_okabe_ito(c(2:3,1,5,8)), se=F)

ptres <- ggarrange(pt1, pt1b, pt2, pt3, ncol=2, nrow=2, common.legend = TRUE, labels = "AUTO")



###  Supplementary Figure 7

pred2 <- mef_ccurv(pred, outlier='NA')
pred2$Learner2[pred2$Learner2=="linear model"] <- "logistic model"
pred2$Learner2[pred2$Learner2=="RF + linear model"] <- "RF + logistic model"

pred2$Learner2 <- factor(pred2$Learner2, levels=c("convex eSL", "non-convex eSL", "dSL", "logistic model", "MARS", "ridge", "lasso", "XGBoost", "RF + logistic model", "RF + MARS", "RF + ridge", "RF + lasso", "RF + XGBoost", "mean"))

pred2$pred[pred2$pred>1] <- 1

g <- ggplot(data=pred2[pred2$CV!="Rolling-Window" & pred2$Learner2!="mean",])  + geom_segment(aes(x=0,y=0,xend=1,yend=1), col="black", size=.5, linetype=2) + geom_smooth(aes(x = pred, y = True, color = Type, fill=Type), method = "gam", size = 1, se=T) +
  xlab("Predictions") + ylab("Observations") + facet_wrap(Learner2 ~ ., ncol=4) + scale_color_manual(values=c(colhist, colind, colsl)) + ggpubr::font("xy.title", face = "bold") + ggpubr::theme_pubr(legend=c(0.85,0.15)) + scale_fill_manual(values=c(colhist, colind, colsl))
ggpubr::ggpar(g, legend.title="") + coord_cartesian(xlim = c(0, 1), ylim=c(0,1)) + xlab(expression(bold("Predictions"))) + ylab(expression(bold("Observations")))
#13749/105036=13.1% pour non-convex eSL (seul avec pred >1)

###  Supplementary Figure 6, panel E

pred2 <- rbind(
  mef_ccurv(pred[pred$t==2,], outlier='NA'),
  mef_ccurv(pred[pred$t==102,], outlier='NA'),
  mef_ccurv(pred[pred$t==202,], outlier='NA'),
  mef_ccurv(pred[pred$t==352,], outlier='NA'),
  mef_ccurv(pred[pred$t==502,], outlier='NA')
)
pred2$Learner2[pred2$Learner2=="linear model"] <- "logistic model"
pred2$Learner2[pred2$Learner2=="RF + linear model"] <- "RF + logistic model"


pred2$pred[pred2$pred>1] <- 1


pred2$t <- as.factor(pred2$t)

levels(pred2$t) <- c("Time 1", "Time 100", "Time 200", "Time 350", "Time 500")

g <- ggplot(data=pred2[pred2$Type=="POSL",])  + geom_segment(aes(x=0,y=0,xend=1,yend=1), col="black", size=.5, linetype=2) + geom_smooth(aes(x = pred, y = True, color = Learner2, fill=Learner2), method = "gam", size = 1, se=F) +
  xlab("Predictions") + ylab("Observations") + facet_grid(cols=vars(t)) + scale_color_manual(values=ggokabeito::palette_okabe_ito(c(2:3,5))) + ggpubr::font("xy.title", face = "bold") + ggpubr::theme_pubr(legend="none") + scale_fill_manual(values=ggokabeito::palette_okabe_ito(c(2:3,5)))
tflex <- ggpubr::ggpar(g, legend.title="") + coord_cartesian(xlim = c(0, 1), ylim=c(0,1)) + xlab(expression(bold("Predictions"))) + ylab(expression(bold("Observations")))
ggarrange(ptres,tflex, labels=c("","E"), heights = c(2, 1), ncol = 1)


###  Supplementary Figure 8

db_auc <- pred[pred$t>1,]
db_auc$eSL[db_auc$eSL>1] <- 1

aucs <- sapply(names(db_auc)[c(2:12, 23:30)],
               function(x) pROC::ci.auc(pROC::roc(db_auc$True, db_auc[,x], percent=F, ci=T)),
               simplify=T) |> round(digits = 3)


pred2 <- reshape(data.frame(aucs), direction="long",
                 varying = colnames(aucs),
                 v.names="pred",
                 times = colnames(aucs),
                 timevar="Learner")

pred2$Type <- ifelse(pred2$Learner %in% unique(grep("ind", pred2$Learner, value=TRUE)) , "Individual" , ifelse(pred2$Learner %in% unique(grep("hist", pred2$Learner, value=TRUE)) , "Historical", "POSL"))
pred2$CV <- ifelse(pred2$Learner %in% unique(grep("Origin", pred2$Learner, value=TRUE)) , "Rolling-Origin" , ifelse(pred2$Learner %in% unique(grep("Window", pred2$Learner, value=TRUE)) , "Rolling-Window" , "None"))
pred2$Type[pred2$Type=="Individual"] <- ifelse(pred2$CV[pred2$Type=="Individual"]=='Rolling-Origin', "Individual", "Individual (RWCV)")


pred2$Learner2 <- ifelse(pred2$Learner %in% unique(grep("glm_rf", pred2$Learner, value=TRUE)), 'RF + logistic model', 
                         ifelse(pred2$Learner %in% unique(grep("glm", pred2$Learner, value=TRUE)), 'logistic model', ifelse(pred2$Learner %in% unique(grep("mars_rf", pred2$Learner, value=TRUE)), 'RF + MARS', ifelse(pred2$Learner %in% unique(grep("mars", pred2$Learner, value=TRUE)), 'MARS', ifelse(pred2$Learner %in% unique(grep("ridge_rf", pred2$Learner, value=TRUE)), 'RF + ridge', ifelse(pred2$Learner %in% unique(grep("ridge", pred2$Learner, value=TRUE)), 'ridge', ifelse(pred2$Learner %in% unique(grep("lasso_rf", pred2$Learner, value=TRUE)), 'RF + lasso', ifelse(pred2$Learner %in% unique(grep("lasso", pred2$Learner, value=TRUE)), 'lasso', ifelse(pred2$Learner %in% unique(grep("xgb_rf", pred2$Learner, value=TRUE)), 'RF + XGBoost', ifelse(pred2$Learner %in% unique(grep("xgb", pred2$Learner, value=TRUE)), 'XGBoost', ifelse(pred2$Learner %in% unique(grep("mean", pred2$Learner, value=TRUE)), 'mean', ifelse(pred2$Learner %in% unique(grep("dSL", pred2$Learner, value=TRUE)), 'dSL', ifelse(pred2$Learner %in% unique(grep("convex", pred2$Learner, value=TRUE)), 'convex eSL', "non-convex eSL")))        ))))))))))  

pred3 <- pred2[pred2$id==1,]   
pred4 <- pred2[pred2$id==3,]   
pred2 <- pred2[pred2$id==2,] 

ggpubr::ggdotchart(data=pred2[pred2$CV!='Rolling-Window',], x = 'Learner2', y = 'pred', color = 'Type', 
                   palette = c(colhist, colind, colsl),
                   sorting = "none",  order=rev(unique(pred2$Learner2)), 
                   rotate = TRUE, shape=18,                         
                   dot.size = 5, ylab="AUROC", xlab="", 
                   ggtheme = theme_pubr(),
                   legend.title = "") + font("xy.title", face = "bold") +  geom_segment(aes(y=pred3$pred[1], yend=pred4$pred[1], x=11, xend=11), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[2], yend=pred4$pred[2], x=10, xend=10), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[3], yend=pred4$pred[3], x=9, xend=9), col=colsl, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[4], yend=pred4$pred[4], x=8, xend=8), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[5], yend=pred4$pred[5], x=7, xend=7), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[6], yend=pred4$pred[6], x=6, xend=6), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[7], yend=pred4$pred[7], x=5, xend=5), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[8], yend=pred4$pred[8], x=4, xend=4), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[9], yend=pred4$pred[9], x=3, xend=3), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[10], yend=pred4$pred[10], x=2, xend=2), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[11], yend=pred4$pred[11], x=1, xend=1), col=colind, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[12], yend=pred4$pred[12], x=8, xend=8), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[13], yend=pred4$pred[13], x=7, xend=7), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[14], yend=pred4$pred[14], x=6, xend=6), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[15], yend=pred4$pred[15], x=5, xend=5), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[16], yend=pred4$pred[16], x=4, xend=4), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[17], yend=pred4$pred[17], x=3, xend=3), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[18], yend=pred4$pred[18], x=2, xend=2), col=colhist, size=1, linetype=1) + 
  geom_segment(aes(y=pred3$pred[19], yend=pred4$pred[19], x=1, xend=1), col=colhist, size=1, linetype=1) 


###  Supplementary Figure 9

names(db_auc)[2:3] <- c("convex_eSL", "eSL")


dca_res <- dcurves::dca(formula=reformulate(names(db_auc[,2:4]), 'True'), data=data.frame(db_auc), label=list(convex_eSL="convex eSL", eSL="non-convex eSL", dSL="dSL"))

as_tibble(dca_res) %>%
  dplyr::filter(!is.na(net_benefit)) %>%
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  geom_line(size=1) +
  coord_cartesian(ylim = c(-0.0826128928456505, 0.826128928456505
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values = c('black', 'red', ggokabeito::palette_okabe_ito(c(2,5,3)))) + theme_pubr() +
  ggpubr::font("xy.title", face = "bold")



####
# Summary Table for sensitivity analysis
####

tbl_res <- function(pred2, metric){
  res <- sapply(unique(pred2$Learner2),
                function(l){
                  if(length(grep("SL", l)) > 0){
                    s <- summary(pred2[pred2$Learner2==l & pred2$CV!="Rolling-Window", metric])
                    se <- var(pred2[pred2$Learner2==l & pred2$CV!="Rolling-Window", metric], na.rm=T)
                    output <- c(s[c(3,4)], se, s[c(1,6)]) 
                  }else{
                    s <- summary(pred2[pred2$Learner2==l & pred2$Type=="Historical", metric])
                    se <- var(pred2[pred2$Learner2==l & pred2$Type=="Historical", metric], na.rm=T)
                    s2 <- summary(pred2[pred2$Learner2==l & pred2$Type=="Individual", metric])
                    se2 <- var(pred2[pred2$Learner2==l & pred2$Type=="Individual", metric], na.rm=T)
                    output <- rbind(c(s[c(3,4)], se, s[c(1,6)]), c(s2[c(3,4)], se2, s2[c(1,6)])) 
                  }
                  return(output ) 
                },
                simplify=T ) 
  
  rnames <- names(res)[-(1:3)] 
  rnames <- c('convex eSL', 'non-convex eSL', 'dSL', paste(c("Historical", "Individual"), sapply(rnames, function(n) rep(n,2),simplify=F) |> unlist()))
  
  res <- res |> do.call(what=rbind) |> round(digits=2) |> data.frame(row.names = rnames) |> setNames(nm=c('Median', 'Mean', 'Variance', 'Minimum', 'Maximum')) 
  
  return(res)
}

res.ape <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'ape')
res.mse <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'mse')
res.mcal <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'mcal')
res.wcal <- tbl_res(db[db$Learner2!='mean' & db$Type!="Individual (RWCV)",], 'wcal')


res <- merge(res.ape, res.mse, by = 'row.names', sort=F) |> merge(y=res.mcal, by.y = 'row.names', by.x='Row.names', sort=F) |> merge(y=res.wcal, by.y = 'row.names', by.x='Row.names', sort=F)

data.table::setorder(res, 'Median.x')

resf <- rbind(res[,1:11], res[,c(1,12:21)])

xtable::xtable(resf)

