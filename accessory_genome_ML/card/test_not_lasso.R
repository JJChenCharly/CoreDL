library(devtools)
library(pROC)
library(dplyr)

# install.packages('Rtools')
# install.packages('caTools')
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

# glm test ----
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
formula_togo

# binarize ----
for (x in rownames(for_train)) {
  if (for_train[x, "v261"] == 'B2') {
    for_train[x, "v261"] <- 1
  } else {
    for_train[x, "v261"] <- 0
  }
}

for_train <- mutate_all(for_train,
                        function(x) as.numeric(x)
                        )


zero_rows <- for_train[for_train[, "v261"] == 0, ]
one_rows <- for_train[for_train[, "v261"] == 1, ]

choose_zero <- sample.split(rownames(zero_rows), SplitRatio = 0.75)
choose_one <- sample.split(rownames(one_rows), SplitRatio = 0.75)

train_rows <- rbind(zero_rows[choose_zero, ], one_rows[choose_one, ])
test_rows <- rbind(zero_rows[!choose_zero, ], one_rows[!choose_one, ])


logistic_model <- glm(formula_togo,
                      data = train_rows,
                      family = "binomial"
                      )

# predict_res <- predict(logistic_model, 
#                        train_rows, 
#                        type = "response")

# predict_res <- ifelse(predict_res > 0.5, 1, 0)

# table(train_rows$v261, predict_res)



roc_ <- roc(response = test_rows$v261,
            predictor = predict(logistic_model,
                                type = "response",
                                newdata = test_rows
                                ) %>% as.matrix()
            )

roc_$auc

logistic_model$coefficients %>% unname() %>% unlist() %>% max(na.rm = TRUE)

v <- colnames(combined_tbl)
which(v == 'TEM-206')
which(v == 'AcrE')

logistic_model$coefficients['v162']
names_ <- logistic_model$coefficients %>% names()
names_[1]
