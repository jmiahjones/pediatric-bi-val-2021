
rm(list=ls())

library(dplyr)
library(caret)
library(foreach)
library(SuperLearner)
library(readxl)
library(doParallel)#Load parallel library

no_cores <- 5 # Set number of cores for parallel computing
validation_file <- "./data/validation-12-6-21.csv"
# derivation_file <- "./derivation/results/revised_trained.RData"

trained_uti_file <- "./derivation/results/revised_trained_uti.RData"
trained_le_file  <- "./derivation/results/revised_trained_le.RData"

load(trained_uti_file)
glm_uti <- glm_nosocial
sl_uti <- sl_nosocial


load(trained_le_file)
glm_le <- glm_nosocial
sl_le <- sl_nosocial

# load(derivation_file)

######################################
#### Validation Set ####
######################################
val_raw <- read.csv(validation_file, header=T)
validation_features <- val_raw %>% 
  rename(
    Record.ID=recordid,
    Sex=gender,
    Presence.of.Chronic.Medical.Condition..CMC. = cmc,
    Ill.Appearing = illappearing,
    Maximum.Temperature = temp,
    Gestational.Age=gestationalage,
    Number.of.Days.of.Illness = numberdaysillness,
    Cough.Present.During.Illness. = coughpresent,
    Positive.Urinary.Tract.Inflammation = urinarytractinflamm,
    BI = bi,
    II = ibi
  ) %>% 
  na.omit %>%  # TODO: Check this!
  mutate(sex= factor(Sex, labels = c("male", "female")),
         cm_condition = factor(Presence.of.Chronic.Medical.Condition..CMC., labels = c("no", "yes")),
         Age = as.numeric(age),
         gestationage = as.numeric(Gestational.Age),
         appearance = factor(Ill.Appearing, levels=c(0,1,2), 
                             labels = c("well", "ill","unknown")),
         maxtemp = as.numeric(Maximum.Temperature),
         numdaysill = as.numeric(Number.of.Days.of.Illness),
         cough = factor(Cough.Present.During.Illness., labels = c("yes", "no","unknown" )),
         puti = factor(Positive.Urinary.Tract.Inflammation, labels = c("no", "yes", "unknown")),
         BI = factor(BI, labels = c("no", "yes")),
         II = factor(II, labels = c("no", "yes")),
         le = factor(Leukocyte.Esterase.Present, labels = c("no", "yes", "unknown"))
  ) %>% 
  select(
    Record.ID, sex, insurance, cm_condition, age, gestationage, 
    appearance, maxtemp, numdaysill, cough, puti, le,
    BI, II
  )
val_ids <- pull(validation_features, Record.ID)
validation_features <- validation_features %>% select(!one_of("Record.ID"))
validation_df <- data.frame(
  model.matrix(BI ~ . - II, data = validation_features)[,-1],
  BI = validation_features$BI,
  II = validation_features$II
)

glm_final <- glm_nosocial
sl_final <- sl_nosocial

###########################################
########### Functions #####################
###########################################


create_report_table <- function(
  the_df, cut.pts, table_file,
  clNum=20,
  B=2000L, # number of bootstraps
  alpha=0.05, # 1-confidence level
  digit_format="%.3f", # number of digits
  conf_stat_names=c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value')
){
  
  sprintf_ci_format <- paste0("(", digit_format, ", ", digit_format, ")")
  sprintf_est_ci_format <- paste0(digit_format, " %s")
  
  
  if(clNum > 1){
    cl <- parallel::makeCluster(clNum, "FORK")
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }
  
  tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    
    sl_conf <- the_df %>%
      select(SL, BI) %>%
      mutate(SL=factor(1*(SL > cut.point), levels=0:1)) %>%
      table %>%
      caret::confusionMatrix(positive="1")
    sl_conf <- sl_conf$byClass[conf_stat_names]
    
    glm_conf <- the_df %>%
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
  
  ci_tbl <- foreach(cut.point=cut.pts, .combine=cbind) %do% {
    foreach(b=1:B, .combine=cbind) %dopar% {
      
      set.seed(34598+b)
      boot <- sample(nrow(the_df), replace=T)
      boot_df <- the_df[boot,]
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
    } %>% apply(1, quantile, probs=c(alpha/2, 1-(alpha/2))) %>% t %>% 
      apply(1, function(x) sprintf(sprintf_ci_format, x[1], x[2]))
  }
  
  stopifnot(all.equal(dim(tbl), dim(ci_tbl)))
  final_tbl <- foreach(a=1:nrow(ci_tbl), .combine=rbind) %:% 
    foreach(b=1:ncol(ci_tbl), .combine=cbind) %do% 
    {
      sprintf(sprintf_est_ci_format, tbl[a,b], ci_tbl[a,b])
    }
  colnames(final_tbl) <- cut.pts
  rownames(final_tbl) <- rownames(tbl)
  write.csv(final_tbl, file=table_file)
  if(clNum > 1)
    stopCluster(cl)
  
  Sys.sleep(2)
  
}

print_roc_results <- function(
  the_df, ci_format="%.2f, (%.2f, %.2f)"
) {
  GLM.roc <- pROC::roc(BI ~ GLM, data=the_df)
  SL.roc <- pROC::roc(BI ~ SL, data=the_df)
  GLM.ci <- (pROC::ci.auc(GLM.roc))
  SL.ci <- (pROC::ci.auc(SL.roc))
  
  paste0(
    "GLM: ",
    sprintf(ci_format,
            GLM.ci[2], GLM.ci[1], GLM.ci[3]),
    "\nSL: ",
    sprintf(ci_format,
            SL.ci[2], SL.ci[1], SL.ci[3])
  )
  
}


###########################################
#### Check both Invasive Infections/BI ####
###########################################

conf_stat_names <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value")
cut.pts <- c(.001, .005, .01, .03, .05)

for(j in 1:2){
  uti_le_B <- j==2
  uti_le_tmp_str <- ifelse(uti_le_B, "uti", "leukesterase")
  
  if(uti_le_B){
    glm_final <- glm_uti
    sl_final <- sl_uti
  } else {
    glm_final <- glm_le
    sl_final <- sl_le
  }
  
  
  for(i in 1:2){
    invasiveB <- i==2
    # set output files
    tmp_str <- ifelse(invasiveB, "ibi", "bi")
    
    result_dir <- sprintf("./results/%s/%s/", uti_le_tmp_str, tmp_str)
    
    report_table_file <- paste0(result_dir, "validation_results.csv")
    sens1_table_file <- paste0(result_dir, "8to60_validation_results.csv")
    sens2_table_file <- paste0(result_dir, "30plus_validation_results.csv")
    sens5_table_file <- paste0(result_dir, "wellappearing_validation_results.csv")
    rec_table_file <- paste0(result_dir, "rec_validation_results.csv")
    
    roc_file <- paste0(result_dir, "other_auc_results.txt")
    
    
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
    
    GLM.hat <- predict(glm_final, newdata=X_val, X=X_nosocial, 
                       Y=Y_nosocial, onlySL=T)$pred %>% 
      as.numeric
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
    
    ## recommended approach: <31 => +, >=31 => model
    rec_df <- sens_df %>% mutate(GLM=if_else(greater30, GLM, 1.0),
                                 SL =if_else(greater30, SL, 1.0))
    
    # Overall validation report table
    create_report_table(
      val_df, cut.pts=cut.pts, table_file=report_table_file
    )
    
    # sensitivity 1
    sens_df %>% filter(sens1) %>% 
      create_report_table(
        cut.pts=cut.pts, table_file=sens1_table_file
      )
    
    # sensitivity 2
    sens_df %>% filter(sens2) %>% 
      create_report_table(
        cut.pts=cut.pts, table_file=sens2_table_file
      )
    
    # recommended approach
    rec_df %>%
      create_report_table(
        cut.pts=cut.pts, table_file=rec_table_file
      )
    
    # Performance on Well-Appearing Infants
    sens_df %>% filter(wellappearing) %>% 
      create_report_table(
        cut.pts=cut.pts, table_file=sens5_table_file
      )
    
    
    ### ROC Curves
    library(pROC)
    GLM.roc <- pROC::roc(Y_val ~ GLM.hat)
    SL.roc <- pROC::roc(Y_val ~ SL.hat)
    print(GLM.roc)
    print(pROC::ci.auc(GLM.roc))
    
    print(SL.roc)
    print(pROC::ci.auc(SL.roc))
    
    #load(derivation_file)
    #train_sl_roc = pROC::roc(Y~SL, data=full_data_pred_df)
    #train_lr_roc = pROC::roc(Y~GLM, data=full_data_pred_df)
    
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
      
      lims = ci[-2]
      auc = ci[2]
      
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
    
    
    # auc_legend = sapply(1:4, auc_fun)
    # png(both_auc_plot_file, width=600, height=600, pointsize = 16)
    # plot(SL.roc, col=sl.color, lty=1,
    #      xlab="Specificity (1 - False Positive %)",
    #      ylab="Sensitivity (True Positive %)")
    # plot(GLM.roc, col=glm.color, lty=1, add=T)
    # plot(train_sl_roc, col=sl.color, lty=3, add=T)
    # plot(train_lr_roc, col=glm.color, lty=3, add=T)
    # legend(x=0.75, y=0.2,
    #        legend=auc_legend,
    #        col=c(sl.color, glm.color, sl.color, glm.color),
    #        lty=c(1,1,3,3),
    #        cex=0.99)
    # dev.off()
    
    
    
    ### Sensitivity Analysis ROC Curves
    library(pROC)
    cat(
      {
        paste0("8 to 60 AUC Results:\n",
               print_roc_results(sens_df %>% filter(sens1)),
               "\n",
               "------------------\n",
               "30+ AUC Results:\n",
               print_roc_results(sens_df %>% filter(sens2)),
               "\n",
               "------------------\n",
               "Pragmatic AUC Results:\n",
               print_roc_results(rec_df),
               "\n",
               "------------------\n",
               "Well-appearing Infants AUC Results:\n",
               print_roc_results(sens_df %>% filter(wellappearing)),
               "\n"
        )
      }, file=roc_file)
    
    
    
    
    
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
    
    
    stopImplicitCluster()
    
  }
}

