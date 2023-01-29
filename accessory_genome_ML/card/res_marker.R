library(glmnet)
# library(rms)
library(devtools)
library("ggpubr")
# install.packages("tidyverse")
# library(tidyverse)

# library(ROCR)
library(pROC)
library(glmtoolbox)
library(dcurves)
library(gtsummary)
library(dplyr)
library(tidyr)
library(reshape2)
library(randomForest)

# clade info ----
clade_info_2 <- read.csv('E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv', 
                         header = FALSE,
                         row.names = 1
)

# readin ----
res_df <- read.csv("res_mech_filtered.csv",
                   #sep = ",",
                   encoding = "UTF-8",
                   header = FALSE,
                   #row.names = FALSE
) %>% t() %>% as.data.frame()
res_df <- res_df[-c(1), ]
# colnames(res_df) <- res_df[1, ]
# res_df <- res_df[-c(1), ]

# deprecate columns
needed_cols <- c()
already_added <- c()

for (x in colnames(res_df)) {
  if (res_df[1, x] %in% already_added) {
    next
  } else {
    needed_cols <- c(needed_cols, x)
    already_added <- c(already_added, res_df[1, x])
  }
}

length(needed_cols)

needed_df <- res_df[, needed_cols]

# colname form
name_of_col <- c()

for (x in needed_df[1, ]) {
  print(x)
  name_of_col <- c(name_of_col,
                   gsub(' ', '_', x))
}

colnames(needed_df) <- name_of_col
needed_df <- needed_df[-c(1), ]

rownames(needed_df) <- needed_df$gene
needed_df <- needed_df[, -c(1)]

# add clade info ----
needed_df$clade <- rep('1', 1127)

for (x in rownames(clade_info_2)) {
  needed_df[x, "clade"] <- clade_info_2[x, "V2"]
}

# marker talbe ----
marker_table <- matrix(0, 6, 260) %>% data.frame()
rownames(marker_table) <- c('A', 'B1', 'B2', 'D', 'E', 'N9')
colnames(marker_table) <- colnames(needed_df)[1: 260]

# function of training ----
get_marker <- function(in_table, clade_) {
  binary_table <- in_table
  
  # binarize ----
  for (x in seq(1, length(rownames(binary_table)))) {
    if (binary_table[x, 261] == clade_) {
      binary_table[x, 261] <- 1
    } else {
      binary_table[x, 261] <- 0
    }
  }
  
  # to numeric ----
  binary_table <- mutate_all(binary_table,
                             function(x) as.numeric(x))
  # str(binary_table)
  
  # regularization by /10000 ----
  reg_tab <- binary_table
  
  for (x in colnames(reg_tab)[1: (length(colnames(reg_tab)) -1)]) {
    reg_tab[x] <- reg_tab[x]/10000#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # max(reg_tab[x])
  }
  
  zero_rows <- reg_tab[reg_tab["clade"] == 0, ]
  one_rows <- reg_tab[reg_tab["clade"] == 1, ]
  
  train_rows <- c(rownames(one_rows) %>% sample(130),
                  rownames(zero_rows) %>% sample(720)
  )
  
  test_rows <- c()
  
  for (x in rownames(reg_tab)) {
    if (x %in% train_rows) {
      next
    } else {
      test_rows <- c(test_rows, x)
    }
  }
  
  trian_set <- reg_tab[train_rows, ]# contain clade
  test_set <- reg_tab[test_rows, ]
  
  togo_matrix <- as.matrix(trian_set[, 1: (length(colnames(trian_set)) -1)])# matrix for training without answer
  test_matrix <- as.matrix(test_set[, 1: (length(colnames(test_set)) -1)])
  
  # multi-fator n lasso ----
  y <- as.factor(trian_set[, 261])
  
  lasso_test_all <- cv.glmnet(x = togo_matrix,
                              y = y,
                              alpha = 1,
                              family = "binomial",
                              nlambda = 1000,
                              nfolds = 3,
                              type.measure = "deviance"#deviance auc
  )
  
  lasso_best <- glmnet(x = togo_matrix,
                       y = y,
                       alpha = 1,
                       family = "binomial",
                       #nlambda = 1000
                       lambda = lasso_test_all$lambda.min
  )
  
  # ROC ----
  # input and answers
  Test <- test_set[, 1: 260] %>% as.matrix()
  Train <- trian_set[, 1: 260] %>% as.matrix()
  
  # auroc
  
  train_roc <- roc(response = trian_set$clade,
                   predictor = predict(lasso_best,
                                       type = "response",
                                       newx = togo_matrix
                   ) %>% as.matrix()
  )
  
  test_roc <- roc(response = test_set$clade,
                  predictor = predict(lasso_best,
                                      type = "response",
                                      newx = test_matrix
                  ) %>% as.matrix()
  )
  
  
  
  a <- coef(lasso_best)
  
  return(list(train_roc$auc, test_roc$auc, a))
}

# marker make ----
for (x in c('A', 'B1', 'B2', 'D', 'E', 'N9')) {
  training_times <- 50
  auroc_train <- c()
  auroc_test <- c()
  for (y in seq(1, training_times)) {# get auroc of 50 training
    res <- get_marker(in_table = needed_df,
                      clade_ = x
    )
    auroc_train <- c(auroc_train, res[[1]])
    auroc_test <- c(auroc_test, res[[2]])
    
    for (z in seq(1, length(res[[3]]))) {
      if (res[[3]][z] != 0 & res[[3]]@Dimnames[[1]][z] != '(Intercept)') {
        marker_table[x, res[[3]]@Dimnames[[1]][z]] <- marker_table[x,
                                                                   res[[3]]@Dimnames[[1]][z]
        ] + res[[3]][z]
      } else {
        next
      }
    }
  }
  
  # print(x)
  # print(sum(auroc_train)/training_times)
  # print(sum(auroc_test)/training_times)
  
  fileConn<-file("E:/adcadmic/Prof_guo/all_ec2/card/card_marker_log.txt", open = 'a')
  write(paste(x,
              sum(auroc_train)/training_times,
              sum(auroc_test)/training_times, sep = '\t'), 
        fileConn, 
        append = TRUE)
  close(fileConn)
  
}



# marker_table <- apply(marker_table, 2, function(x){x/training_times})
marker_table <- marker_table/training_times

# getwd()
write.csv(marker_table,
          "E:/adcadmic/Prof_guo/all_ec2/card/card_marker.csv", 
          row.names = TRUE)

