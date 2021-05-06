
rm(list=ls())

library(dplyr)
library(caret)
library(foreach)
library(SuperLearner)
library(readxl)
library(doParallel)#Load parallel library

no_cores <- 5 # Set number of cores for parallel computing
main_data_path <- "./derivation/main_analysis_data.csv"
validation_file <- "./validation-data-4-5-21.csv"
derivation_file <- "./derivation/results/SL_bootstrap.RData"
retrain=F

my_round <- function(x, num){
  xr <- round(x, num)
  format(xr, nsmall=num)
}

group_testind <- function(groups, all_idxs, k=10){
  trfolds <- caret::groupKFold(groups, k=k)
  lapply(trfolds, function(tr) setdiff(all_idxs, tr))
}


raw <- read.csv(main_data_path,
                header=T)
# raw = raw[-1,]
colnames(raw)[1] <- "Record.ID"


raw <- raw %>% mutate(BI=ifelse(is.na(True.Bacterial.Infection), 
                                0, 
                                True.Bacterial.Infection))

#==========================================================================================#
#============================ step 1: data processing =====================================#
#==========================================================================================#
##Clean data set
library(dplyr)
features.raw <- raw %>%
  mutate(sex= factor(Sex, labels = c("male", "female")),
         insurance = factor(Insurance.Type, labels = c("private", "public", "none")),
         cm_condition = factor(Presence.of.Chronic.Medical.Condition..CMC., labels = c("no", "yes")),
         age = as.numeric(Age),
         gestationage = as.numeric(Gestational.Age),
         appearance = factor(Ill.Appearing, labels = c("well", "ill","unknown")),
         maxtemp = as.numeric(Maximum.Temperature),
         numdaysill = as.numeric(Number.of.Days.of.Illness),
         cough = factor(Cough.Present.During.Illness., labels = c("yes", "no","unknown" )),
         puti = factor(Positive.Urinary.Tract.Inflammation, labels = c("no", "yes", "unknown")),
         BI = factor(BI, labels = c("no", "yes"))
  )

features_nosocial = features.raw %>% select(Record.ID, sex, insurance, cm_condition, age, gestationage, 
                                            appearance, maxtemp, numdaysill, cough, puti,
                                            BI)


dim(features_nosocial)

##Explore data set
missingFeatures_nonsocial <- as.data.frame(apply(features_nosocial, 2, function(x) sum(is.na(x))))
# no missingness

data.features_nosocial <- data.frame(model.matrix(BI ~., data = features_nosocial)[,-1],
                                     BI=features_nosocial$BI)
record_ids <- data.features_nosocial[,1]

data.use_nosocial <- data.features_nosocial[-1]

data.use_nosocial$BI<-as.numeric(data.use_nosocial$BI)-1



if(retrain){
  nb<-500
  bs.dat <- vector("list", length = nb)
  
  set.seed(608991)
  for(i in 1:nb) {
    idxs <- sample(nrow(data.use_nosocial),replace=TRUE)
    bs.dat[[i]]<-idxs
    
  }
  
  Y_nosocial <-as.numeric(data.use_nosocial$BI)
  X_nosocial<-data.use_nosocial[,!names(data.use_nosocial)%in%c("BI")]
  
  # learner=create.Learner("SL.randomForest", tune = list(mtry = seq(6,29,1),ntrees=c(500,1000)))	
  
  rflrn=create.Learner("SL.randomForest", 
                       tune = list(mtry = seq(6,13,2),
                                   ntrees=c(1000)))
  
  SL.lib <- c(
    "SL.earth",
    "SL.gam",
    "SL.glm",
    "SL.glmnet",
    rflrn$names
  )
  
  sl_final = SuperLearner(
    Y = Y_nosocial,
    X = X_nosocial,
    family = binomial(),
    method = "method.AUC",
    SL.library = SL.lib
  )
  
  glm_final = glm(
    Y_nosocial ~ ., data=X_nosocial,
    family = binomial()
  )
  
  save.image("./validation.RData")
} else {
  load("./validation.RData")
}

######################################
#### Validation Set ####
######################################
val_raw <- read.csv(validation_file, header=T)
validation_features <- val_raw %>% 
  rename(
    Sex=Gender,
    Presence.of.Chronic.Medical.Condition..CMC. = CMC,
    Ill.Appearing = Ill.Appearance,
    Maximum.Temperature = Temperature,
    Number.of.Days.of.Illness = Number.of.Days.with.Illness,
    Cough.Present.During.Illness. = Cough.Present,
    Positive.Urinary.Tract.Inflammation = Urinary.tract.inflammation,
    BI = Bacterial.Infection.Present,
    II = Invasive.Bacterial.Infection
  ) %>% 
  filter(!is.na(Insurance.Type) & Insurance.Type != 3) %>% 
  na.omit %>%  # TODO: Check this!
  mutate(sex= factor(Sex, labels = c("male", "female")),
         insurance = factor(Insurance.Type, labels = c("private", "public", "none")),
         cm_condition = factor(Presence.of.Chronic.Medical.Condition..CMC., labels = c("no", "yes")),
         age = as.numeric(Age),
         gestationage = as.numeric(Gestational.Age),
         appearance = factor(Ill.Appearing, levels=c(0,1,2), 
                             labels = c("well", "ill","unknown")),
         maxtemp = as.numeric(Maximum.Temperature),
         numdaysill = as.numeric(Number.of.Days.of.Illness),
         cough = factor(Cough.Present.During.Illness., labels = c("yes", "no","unknown" )),
         puti = factor(Positive.Urinary.Tract.Inflammation, labels = c("no", "yes", "unknown")),
         BI = factor(BI, labels = c("no", "yes")),
         II = factor(II, labels = c("no", "yes"))
  ) %>% 
  select(
    Record.ID, sex, insurance, cm_condition, age, gestationage, 
    appearance, maxtemp, numdaysill, cough, puti,
    BI, II
  )
val_ids <- pull(validation_features, Record.ID)
validation_features <- validation_features %>% select(!one_of("Record.ID"))
validation_df <- data.frame(
  model.matrix(BI ~ . - II, data = validation_features)[,-1],
  BI = validation_features$BI,
  II = validation_features$II
  )


###########################################
#### Check both Invasive Infections/BI ####
###########################################

conf_stat_names <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
cut.pts <- c(.001, .005, .01, .03, .05)

for(i in 1:2){
  invasiveB <- i==2
  # set output files
  tmp_str <- ifelse(invasiveB, "-inv", "")
  
  result_dir <- sprintf("./results%s/", tmp_str)
  
  report_table_file <- paste0(result_dir, "validation_results.csv")
  sens1_table_file <- paste0(result_dir, "8to60_validation_results.csv")
  sens2_table_file <- paste0(result_dir, "30plus_validation_results.csv")
  sens1_roc_file <- paste0(result_dir, "8to60_roc.txt")
  sens2_roc_file <- paste0(result_dir, "30plus_roc.txt")
  
  sens5_table_file <- paste0(result_dir, "wellappearing_validation_results.csv")
  sens5_roc_file <- paste0(result_dir, "wellappearing_roc.txt")
  
  rec_table_file <- paste0(result_dir, "rec_validation_results.csv")
  recommended_roc_file <- paste0(result_dir, "rec_roc.txt")
  
  missed_records_file <- paste0(result_dir, "false_negatives.csv")
  rec_missed_records_file <- paste0(result_dir, "rec_false_negatives.csv")
  low_risk_file <- paste0(result_dir, "low_risk.csv")
  
  auc_plot_file <- paste0(result_dir, "validation_roc.png")
  both_auc_plot_file <- paste0(result_dir, "derivation_validation_roc.png")
  
  
  if(invasiveB) {
    Y_val <- as.numeric(validation_df$II) - 1
  } else {
    Y_val <- as.numeric(validation_df$BI) - 1
  }
  X_val <- validation_df %>% select(-all_of(c("BI","II")))
  # data.use_nosocial[,!names(validation_features)%in%c("BI","Record.ID")]
  
  GLM.hat <- glm_final %>% predict(X_val, type="response")
  SL.hat <- predict(sl_final, newdata=X_val, X=X_nosocial, 
                    Y=Y_nosocial, onlySL=T)$pred %>% 
    as.numeric
  
  val_df <- data.frame(BI=Y_val, GLM=GLM.hat, SL=SL.hat, ID=val_ids)
  
  sens_df <- validation_features %>% 
    mutate(preterm=gestationage <= 36,
           greater30 = (31 <= age & age <= 90 ),
           sens1=(!preterm) & ( 8 <= age & age <= 60 ),
           sens2=(!preterm) & greater30,
           wellappearing=appearance == "well"
           ) %>% 
    dplyr::select(sens1, sens2, greater30, wellappearing) %>% 
    cbind(val_df)
  
  
  tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    
    sl_conf <- val_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- val_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  cl <- parallel::makeCluster(10, "FORK")
  registerDoParallel(cl)
  ci_tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:2000L, .combine=cbind) %dopar% {
      
      set.seed(34598+b)
      boot <- sample(nrow(val_df), replace=T)
      boot_df <- val_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(0.025,0.975)) %>% t %>% 
      apply(1, function(x) paste0("(", paste(my_round(x, 3), collapse=", "), ")"))
  }
  
  stopifnot(all.equal(dim(tbl), dim(ci_tbl)))
  final_tbl <- foreach(a=1:nrow(ci_tbl), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl), .combine=cbind) %do% 
    {
      paste0(my_round(tbl[a,b], 3), " ", ci_tbl[a,b])
    }
  colnames(final_tbl) <- cut.pts
  rownames(final_tbl) <- rownames(tbl)
  write.csv(final_tbl, file=report_table_file)
  
  
  
  
  ### ROC Curves
  library(pROC)
  GLM.roc <- pROC::roc(Y_val ~ GLM.hat)
  SL.roc <- pROC::roc(Y_val ~ SL.hat)
  print(GLM.roc)
  print(pROC::ci.auc(GLM.roc))
  
  print(SL.roc)
  print(pROC::ci.auc(SL.roc))
  
  load(derivation_file)
  train_sl_roc = pROC::roc(Y~SL, data=full_data_pred_df)
  train_lr_roc = pROC::roc(Y~GLM, data=full_data_pred_df)
  
  auc_fun <- function(i) {
    if(i == 1){
      roc = SL.roc
      name_str = "Validation SL"
    } else if(i == 2){
      roc = GLM.roc
      name_str = "Validation LR"
    } else if(i==3){
      roc = train_sl_roc
      name_str = "Derivation SL"
    } else if(i==4){
      roc = train_lr_roc
      name_str = "Derivation LR"
    }
    ci = ci.auc(roc)
    # lims = my_round(ci[-2],3)
    lims = ci[-2]
    auc = ci[2]
    # auc = my_round(ci[2],3)
    
    sprintf(paste0(name_str, " AUC: %0.2f CI: (%0.2f, %0.2f)"), auc, lims[1], lims[2])
  }
  sl.color <- "grey"
  glm.color <- "black"
  auc_legend = sapply(1:2, auc_fun)
  png(auc_plot_file, width=600, height=600, pointsize = 16)
  plot(SL.roc, col=sl.color, lty=1,
       xlab="Specificity (1 - False Positive %)",
       ylab="Sensitivity (True Positive %)")
  plot(GLM.roc, col=glm.color, lty=1,
       add=T)
  legend(x=0.65, y=0.15,
         legend=auc_legend,
         col=c(sl.color, glm.color),
         lty=c(1,1),
         cex=0.9)
  dev.off()
  
  
  auc_legend = sapply(1:4, auc_fun)
  png(both_auc_plot_file, width=600, height=600, pointsize = 16)
  plot(SL.roc, col=sl.color, lty=1,
       xlab="Specificity (1 - False Positive %)",
       ylab="Sensitivity (True Positive %)")
  plot(GLM.roc, col=glm.color, lty=1, add=T)
  plot(train_sl_roc, col=sl.color, lty=3, add=T)
  plot(train_lr_roc, col=glm.color, lty=3, add=T)
  legend(x=0.75, y=0.2,
         legend=auc_legend,
         col=c(sl.color, glm.color, sl.color, glm.color),
         lty=c(1,1,3,3),
         cex=0.99)
  dev.off()
  
  
  # Sensitivity Analysis 1
  tbl1 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    sens1_df <- sens_df %>% filter(sens1)
    sl_conf <- sens1_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- sens1_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  ci_tbl1 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:2000L, .combine=cbind) %dopar% {
      sens1_df <- sens_df %>% filter(sens1)
      set.seed(34598+b)
      boot <- sample(nrow(sens1_df), replace=T)
      boot_df <- sens1_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(0.025,0.975)) %>% t %>% 
      apply(1, function(x) paste0("(", paste(my_round(x, 3), collapse=", "), ")"))
  }
  
  stopifnot(all.equal(dim(tbl1), dim(ci_tbl1)))
  sens_tbl1 <- foreach(a=1:nrow(ci_tbl1), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl1), .combine=cbind) %do% 
    {
      paste0(my_round(tbl1[a,b], 3), " ", ci_tbl1[a,b])
    }
  colnames(sens_tbl1) <- cut.pts
  rownames(sens_tbl1) <- rownames(tbl1)
  write.csv(sens_tbl1, file=sens1_table_file)
  
  
  
  
  ### Sensitivity 1 ROC Curves
  library(pROC)
  capture.output(
    {
      GLM.roc <- pROC::roc(BI ~ GLM, data=sens_df %>% filter(sens1))
      SL.roc <- pROC::roc(BI ~ SL, data=sens_df %>% filter(sens1))
      print(GLM.roc)
      print(pROC::ci.auc(GLM.roc))
      
      print(SL.roc)
      print(pROC::ci.auc(SL.roc))
    },
    file=sens1_roc_file
  )
  
  
  
  
  
  
  # Sensitivity Analysis 2
  tbl2 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    sens2_df <- sens_df %>% filter(sens2)
    sl_conf <- sens2_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- sens2_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  ci_tbl2 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:2000L, .combine=cbind) %dopar% {
      sens2_df <- sens_df %>% filter(sens2)
      set.seed(34598+b)
      boot <- sample(nrow(sens2_df), replace=T)
      boot_df <- sens2_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(0.025,0.975)) %>% t %>% 
      apply(1, function(x) paste0("(", paste(my_round(x, 3), collapse=", "), ")"))
  }
  
  stopifnot(all.equal(dim(tbl2), dim(ci_tbl2)))
  sens_tbl2 <- foreach(a=1:nrow(ci_tbl2), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl2), .combine=cbind) %do% 
    {
      paste0(my_round(tbl2[a,b], 3), " ", ci_tbl2[a,b])
    }
  colnames(sens_tbl2) <- cut.pts
  rownames(sens_tbl2) <- rownames(tbl2)
  write.csv(sens_tbl2, file=sens2_table_file)
  
  
  
  
  ### Sensitivity 2 ROC Curves
  library(pROC)
  capture.output(
    {
      GLM.roc <- pROC::roc(BI ~ GLM, data=sens_df %>% filter(sens2))
      SL.roc <- pROC::roc(BI ~ SL, data=sens_df %>% filter(sens2))
      print(GLM.roc)
      print(pROC::ci.auc(GLM.roc))
      
      print(SL.roc)
      print(pROC::ci.auc(SL.roc))
    },
    file=sens2_roc_file
  )
  
  
  ## recommended approach: <31 => +, >=31 => model
  # ROC
  rec_df <- sens_df %>% mutate(GLM=if_else(greater30, GLM, 1.0),
                               SL =if_else(greater30, SL, 1.0))
  capture.output(
    {
      GLM.roc <- pROC::roc(BI ~ GLM, data=rec_df)
      SL.roc <- pROC::roc(BI ~ SL, data=rec_df)
      print(GLM.roc)
      print(pROC::ci.auc(GLM.roc))
      
      print(SL.roc)
      print(pROC::ci.auc(SL.roc))
    },
    file=recommended_roc_file
  )
  
  # Cut point
  tbl_rec <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    sl_conf <- rec_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- rec_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  ci_tbl_rec <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:2000L, .combine=cbind) %dopar% {
      set.seed(34598+b)
      boot <- sample(nrow(rec_df), replace=T)
      boot_df <- rec_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(0.025,0.975)) %>% t %>% 
      apply(1, function(x) paste0("(", paste(my_round(x, 3), collapse=", "), ")"))
  }
  
  stopifnot(all.equal(dim(tbl_rec), dim(ci_tbl_rec)))
  sens_tbl_rec <- foreach(a=1:nrow(ci_tbl_rec), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl_rec), .combine=cbind) %do% 
    {
      paste0(my_round(tbl_rec[a,b], 3), " ", ci_tbl_rec[a,b])
    }
  colnames(sens_tbl_rec) <- cut.pts
  rownames(sens_tbl_rec) <- rownames(tbl_rec)
  write.csv(sens_tbl_rec, file=rec_table_file)
  
  
  
  
  # Performance on Well-Appearing Infants
  tbl5 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    sens5_df <- sens_df %>% filter(wellappearing)
    sl_conf <- sens5_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- sens5_df %>%
      select(GLM, BI) %>%
      mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    glm_conf <- glm_conf$byClass[conf_stat_names]
    
    ret <- c(sl_conf, glm_conf)
    names(ret) <- lapply(c("SL.", "GLM."),
                         paste0, conf_stat_names) %>% do.call(c, .)
    
    sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
    glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
    sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
    glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
    lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
    names(lr) <- c("SL.LR+", "GLM.LR+",
                   "SL.LR-", "GLM.LR-")
    
    ret <- c(ret, lr)
    
    return(ret)
  }
  
  ci_tbl5 <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:2000L, .combine=cbind) %dopar% {
      sens5_df <- sens_df %>% filter(wellappearing)
      set.seed(34598+b)
      boot <- sample(nrow(sens5_df), replace=T)
      boot_df <- sens5_df[boot,]
      sl_conf <- boot_df %>%
        select(SL, BI) %>%
        mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      sl_conf <- sl_conf$byClass[conf_stat_names]
      
      glm_conf <- boot_df %>%
        select(GLM, BI) %>%
        mutate(GLM=factor(1*(GLM > cut.point), levels=0:1)) %>%
        table %>%
        caret::confusionMatrix(positive="1")
      glm_conf <- glm_conf$byClass[conf_stat_names]
      
      ret <- c(sl_conf, glm_conf)
      names(ret) <- lapply(c("SL.", "GLM."),
                           paste0, conf_stat_names) %>% do.call(c, .)
      
      sl.lrp <- ret["SL.Sensitivity"]/(1-ret["SL.Specificity"])
      glm.lrp <- ret["GLM.Sensitivity"]/(1-ret["GLM.Specificity"])
      sl.lrn <- (1-ret["SL.Sensitivity"])/ret["SL.Specificity"]
      glm.lrn <- (1-ret["GLM.Sensitivity"])/ret["GLM.Specificity"]
      lr <- c(sl.lrp, glm.lrp, sl.lrn, glm.lrn)
      names(lr) <- c("SL.LR+", "GLM.LR+",
                     "SL.LR-", "GLM.LR-")
      
      ret <- c(ret, lr)
      ret[is.nan(ret)] <- 1e3
      
      return(ret)
    } %>% apply(1, quantile, probs=c(0.025,0.975)) %>% t %>% 
      apply(1, function(x) paste0("(", paste(my_round(x, 3), collapse=", "), ")"))
  }
  
  stopifnot(all.equal(dim(tbl5), dim(ci_tbl5)))
  sens_tbl5 <- foreach(a=1:nrow(ci_tbl5), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl5), .combine=cbind) %do% 
    {
      paste0(my_round(tbl5[a,b], 3), " ", ci_tbl5[a,b])
    }
  colnames(sens_tbl5) <- cut.pts
  rownames(sens_tbl5) <- rownames(tbl5)
  write.csv(sens_tbl5, file=sens5_table_file)
  
  ### Well-Appearing ROC Curves
  library(pROC)
  capture.output(
    {
      GLM.roc <- pROC::roc(BI ~ GLM, data=sens_df %>% filter(wellappearing))
      SL.roc <- pROC::roc(BI ~ SL, data=sens_df %>% filter(wellappearing))
      print(GLM.roc)
      print(pROC::ci.auc(GLM.roc))
      
      print(SL.roc)
      print(pROC::ci.auc(SL.roc))
    },
    file=sens5_roc_file
  )
 
  
  
  
  # print which infants were missed at which cut point
  sens_df %>% 
    rename(Outcome=BI) %>% 
    # filter((SL <= .05 | GLM <= .05) & BI == 1) %>% 
    arrange(SL, GLM) %>% 
    mutate(across(one_of(c("SL", "GLM")), ~cut(.x, breaks=c(0, cut.pts, 1),labels=F))) %>% 
    mutate(across(one_of(c("SL", "GLM")), 
                  ~factor(.x, levels=1:(length(cut.pts)+1),
                          labels=c(
                            paste0("Missed ", cut.pts, "+"), 
                            paste0("Not missed below ", cut.pts[length(cut.pts)]))
                  )
    )) %>% 
    rename(`8to60`=sens1,
           `30plus`=sens2) %>% 
    select(ID, Outcome, SL, GLM, `30plus`, `8to60`, wellappearing) %>% 
    write.csv(file=missed_records_file, row.names = F)
  
  
  # print which infants were missed at which cut point
  rec_df %>% 
    rename(Outcome=BI) %>% 
    # filter((SL <= .05 | GLM <= .05) & BI == 1) %>% 
    arrange(SL, GLM) %>% 
    mutate(across(one_of(c("SL", "GLM")), ~cut(.x, breaks=c(0, cut.pts, 1),labels=F))) %>% 
    mutate(across(one_of(c("SL", "GLM")), 
                  ~factor(.x, levels=1:(length(cut.pts)+1),
                          labels=c(
                            paste0("Predicted negative at ", cut.pts, "+"), 
                            paste0("Not predicted negative below ", cut.pts[length(cut.pts)]))
                  )
    )) %>% 
    select(ID, Outcome, SL, GLM) %>% 
    write.csv(file=rec_missed_records_file, row.names = F)
  
  
  # rec_df %>% select(ID, SL, GLM) %>% 
  #   rename(SL_rec=SL, GLM_rec=GLM) %>%
  #   inner_join(val_df, by="ID") %>% 
  #   filter(SL <= .05 | GLM <= .05 | SL_rec <= .05 | GLM_rec <= .05) %>% 
  #   mutate(across(one_of(c("SL", "GLM", "SL_rec", "GLM_rec")), ~as.numeric(sprintf("%0.2f", .x)))) %>% 
  #   select(ID, SL, GLM, SL_rec, GLM_rec) %>% 
  #   write.csv(file=low_risk_file, row.names=F)
  
  stopImplicitCluster()
  
}

  
