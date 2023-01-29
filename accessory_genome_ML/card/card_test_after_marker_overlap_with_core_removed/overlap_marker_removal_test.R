library(glmnet)
library(devtools)
library("ggpubr")
library(pROC)
library(glmtoolbox)
library(dcurves)
library(gtsummary)
library(dplyr)
library(tidyr)
library(reshape2)
library(randomForest)
library(stringr)

# clade info ----
clade_info_2 <- read.csv('E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv', 
                         header = FALSE,
                         row.names = 1
)

# read in ----
ori_table <- read.csv("E:/adcadmic/Prof_guo/all_ec2/card/res_mech_filtered.csv",
                      header = FALSE,
                      # row.names = NA,
)

colnames(ori_table) <- ori_table[1, ]
ori_table <- ori_table[-c(1), -c(1)]

to_keep <- c()
rep_checker <- c()

for (x in rownames(ori_table)) {
  if (ori_table[x, "gene"] %in% rep_checker) {
    next
  } else {
    to_keep <- c(to_keep, x)
    rep_checker <- c(rep_checker, ori_table[x, "gene"])
  }
}


ori_table <- ori_table[to_keep, ]

# transpose
ori_table <- ori_table %>% t() %>% as.data.frame()
colnames(ori_table) <- ori_table[1, ]
ori_table <- ori_table[-c(1), ]

# remove marker overlap with core !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# colnames(ori_table) %>% unique() %>% length()
# 'cpxA' %in% colnames(ori_table)

# [10] "Escherichia coli soxR with mutation conferring antibiotic resistance*"                     
# [11] "Escherichia coli soxS with mutation conferring antibiotic resistance*" 
# [4] "CRP"  

ori_table <- ori_table[, 
                       -which(names(ori_table) %in% c("Escherichia coli soxR with mutation conferring antibiotic resistance*",
                                                      "Escherichia coli soxS with mutation conferring antibiotic resistance*",
                                                      "CRP" , "cpxA",
                                                      'TEM-206'
                                                      )
                              )
                       ]

# sort according to name
# ori_table <- ori_table[rownames(ori_table) %>% sort() %>% as.vector(), ]
str(ori_table)

# the table needed ----
combined_tbl <- ori_table
# "Escherichia coli soxR with mutation conferring antibiotic resistance*" %in% colnames(combined_tbl) # F
# "MCR-1.1" %in% colnames(combined_tbl) # T

combined_tbl$clade <- rep('1', length(rownames(combined_tbl)
)
)
for (x in rownames(combined_tbl)) {
  combined_tbl[x, "clade"] <- clade_info_2[x, "V2"]
  
}


# marker table ----
marker_table <- matrix(0, 6, length(colnames(ori_table)
)
) %>% data.frame()
rownames(marker_table) <- combined_tbl$clade %>% unique()
colnames(marker_table) <- colnames(combined_tbl)[1: length(colnames(ori_table)
)
]

# function of training ----
get_marker <- function(in_table, clade_) {
  binary_table <- in_table
  
  # binarize ----
  for (x in seq(1, length(rownames(binary_table)
  )
  )
  ) {
    if (binary_table[x, length(colnames(in_table))] == clade_) {
      binary_table[x, length(colnames(in_table))] <- 1
    } else {
      binary_table[x, length(colnames(in_table))] <- 0
    }
  }
  
  # to numeric ----
  binary_table <- mutate_all(binary_table,
                             function(x) as.numeric(x))
  
  # regularization by /10000 ----
  reg_tab <- binary_table
  
  for (x in colnames(reg_tab)[1: (length(colnames(reg_tab)) -1)]) {
    reg_tab[x] <- reg_tab[x]/10#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # max(reg_tab[x])
  }
  
  zero_rows <- reg_tab[reg_tab["clade"] == 0, ]
  one_rows <- reg_tab[reg_tab["clade"] == 1, ]
  
  train_rows <- c(rownames(one_rows) %>% sample(150),
                  rownames(zero_rows) %>% sample(680)
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
  y <- as.factor(trian_set[, length(colnames(trian_set)
  )
  ]
  )
  
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
for (x in rownames(marker_table)) {
  training_times <- 1
  auroc_train <- c()
  auroc_test <- c()
  
  for (y in seq(1, training_times)) {# get auroc of 50 training
    res <- get_marker(in_table = combined_tbl,
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
  
  fileConn<-file("E:/adcadmic/Prof_guo/all_ec2/card/card_test_after_marker_overlap_with_core_removed/res_marker_log.txt", 
                 open = 'a')
  write(paste(x,
              sum(auroc_train)/training_times,
              sum(auroc_test)/training_times, 
              sep = '\t'), 
        fileConn, 
        append = TRUE)
  close(fileConn)
  
}

marker_table <- marker_table/training_times

write.csv(marker_table,
          "E:/adcadmic/Prof_guo/all_ec2/card/card_test_after_marker_overlap_with_core_removed/res_marker.csv", 
          row.names = TRUE
          )