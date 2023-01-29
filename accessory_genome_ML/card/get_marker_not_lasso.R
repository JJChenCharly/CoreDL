library(devtools)
library(pROC)
library(dplyr)
library('caTools')

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
ori_table <- ori_table[-c(1) ,]

# sort according to name 
ori_table <- ori_table[rownames(ori_table) %>% sort() %>% as.vector(), ]
str(ori_table)

# the table needed ----
combined_tbl <- ori_table
combined_tbl$clade <- rep('1', length(rownames(combined_tbl)
)
)
for (x in rownames(combined_tbl)) {
  combined_tbl[x, "clade"] <- clade_info_2[x, "V2"]
  
}

# pre treat table ----
col_ <- colnames(combined_tbl)
new_col <- c()

for (x in seq(1, length(col_))) {
  new_col <- c(new_col, paste0("v", as.character(x)
  )
  )
}


for_train <- combined_tbl

colnames(for_train) <- new_col


formula_vec <- c()
colnames(for_train)

for (x in colnames(for_train)) {
  # print(x)
  if (x == "v261") {
    next
  } else {
    formula_vec <- c(formula_vec, x)
  }
  
}

f <- paste("v261~", paste(sprintf('%s', formula_vec), collapse="+"))# no tick

formula_togo <- formula(f)
# formula_togo

# colname hashmap ----
col_hash <- for_train[1, ]
col_hash[1, ] <- colnames(combined_tbl)
rownames(col_hash) <- c('v')
# str(col_hash)

# marker table ----
marker_table <- matrix(0, 6, length(colnames(combined_tbl)
                                    )
                       ) %>% data.frame()
rownames(marker_table) <- combined_tbl$clade %>% unique()
colnames(marker_table) <- colnames(combined_tbl)[1: length(colnames(combined_tbl)
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
                             function(x) as.numeric(x)
                             )
  
  # regularization by /10000 ----
  reg_tab <- binary_table
  
  for (x in colnames(reg_tab)[1: (length(colnames(reg_tab)) -1)]) {
    reg_tab[x] <- reg_tab[x]#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # max(reg_tab[x])
  }
  
  zero_rows <- reg_tab[reg_tab[, "v261"] == 0, ]
  one_rows <- reg_tab[reg_tab[, "v261"] == 1, ]
  
  choose_zero <- sample.split(rownames(zero_rows), SplitRatio = 0.75)
  choose_one <- sample.split(rownames(one_rows), SplitRatio = 0.75)
  
  train_rows <- rbind(zero_rows[choose_zero, ], one_rows[choose_one, ])
  test_rows <- rbind(zero_rows[!choose_zero, ], one_rows[!choose_one, ])
  
  logistic_model <- glm(formula_togo,
                        data = train_rows,
                        family = "binomial"
  )
  
  roc_train <- roc(response = train_rows$v261,
                   predictor = predict(logistic_model,
                                       type = "response",
                                       newdata = train_rows
                                       ) %>% as.matrix()
                   )
  
  roc_test <- roc(response = test_rows$v261,
                   predictor = predict(logistic_model,
                                       type = "response",
                                       newdata = test_rows
                                       ) %>% as.matrix()
                  )
  
  a <- logistic_model$coefficients
  
  return(list(roc_train$auc, roc_test$auc, a)
         )
}


# marker make ----
for (x in rownames(marker_table)) {
  training_times <- 50
  auroc_train <- c()
  auroc_test <- c()
  
  for (y in seq(1, training_times)) {# get auroc of 50 training
    res <- get_marker(in_table = for_train,
                      clade_ = x
    )
    
    auroc_train <- c(auroc_train, res[[1]])
    auroc_test <- c(auroc_test, res[[2]])
    
    for (z in names(res[[3]])) {
      if (!(is.na(res[[3]][z])) & z != '(Intercept)') {
        marker_table[x, col_hash['v', z]] <- marker_table[x,
                                                          col_hash['v', z]
                                                          ] + res[[3]][z]
      } else {
        next
      }
    }
  }
  
  fileConn<-file("E:/adcadmic/Prof_guo/all_ec2/card/res_marker_log.txt", 
                 open = 'a')
  write(paste(x,
              sum(auroc_train)/training_times,
              sum(auroc_test)/training_times, sep = '\t'), 
        fileConn, 
        append = TRUE)
  close(fileConn)
  
}

marker_table <- marker_table/training_times

write.csv(marker_table,
          "E:/adcadmic/Prof_guo/all_ec2/card/res_marker.csv", 
          row.names = TRUE)
