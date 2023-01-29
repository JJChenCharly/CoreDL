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

# binary classification ----
B2_no_B2_table <- needed_df

# binarize ----
for (x in seq(1, length(rownames(B2_no_B2_table)))) {
  if (B2_no_B2_table[x, 261] == "D" | B2_no_B2_table[x, 261] == "N9") {
  # if (B2_no_B2_table[x, 261] == "A") {
    B2_no_B2_table[x, 261] <- 1
  } else {
    B2_no_B2_table[x, 261] <- 0
  }
}

# to numeric ----
B2_no_B2_table <- mutate_all(B2_no_B2_table,
                             function(x) as.numeric(x))
# str(B2_no_B2_table)

# regularization by /10000 ----
reg_tab <- B2_no_B2_table

# for (x in colnames(reg_tab)) {
#   if ( reg_tab[, x] %>% sum() == 0) {
#     print(x)
#   }
# }

for (x in colnames(reg_tab)[1: (length(colnames(reg_tab)) -1)]) {
  reg_tab[x] <- reg_tab[x]/1000#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # max(reg_tab[x])
}

non_B2_rows <- reg_tab[reg_tab["clade"] == 0, ]
B2_rows <- reg_tab[reg_tab["clade"] == 1, ]

length(rownames(non_B2_rows))
length(rownames(B2_rows))

train_rows <- c(rownames(B2_rows) %>% sample(280),
                rownames(non_B2_rows) %>% sample(600))
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
y

test_all_lambda <- glmnet(x =togo_matrix,
                          y = y,
                          alpha = 1,
                          family = "binomial",
                          nlambda = 300,
                          # lambda = lasso_test_all$lambda.min
)

# str(test_all)
# print(test_all_lambda)

# test_all_lambda

plot(test_all_lambda, xvar = "lambda")

# cv.glmnet
lasso_test_all <- cv.glmnet(x = togo_matrix,
                            y = y,
                            alpha = 1,
                            family = "binomial",
                            nlambda = 1000,
                            nfolds = 3,
                            type.measure = "deviance"#deviance auc
)

plot(lasso_test_all)
lasso_test_all$lambda.min
# lasso_test_all$glmnet.fit


lasso_best <- glmnet(x = togo_matrix,
                     y = y,
                     alpha = 1,
                     family = "binomial",
                     #nlambda = 1000
                     lambda = lasso_test_all$lambda.min
)
coef(lasso_best)

# ROC ----
# input and answers
Test <- test_set[, 1: 260] %>% as.matrix()
Train <- trian_set[, 1: 260] %>% as.matrix()

# togo_matrix

train_roc <- roc(response = trian_set$clade,
                 predictor = predict(lasso_best,
                                     type = "response",
                                     newx = togo_matrix
                 ) %>% as.matrix()
)
train_roc$auc


test_roc <- roc(response = test_set$clade,
                predictor = predict(lasso_best,
                                    type = "response",
                                    newx = test_matrix
                ) %>% as.matrix()
)

# predictor = predict(lasso_best,
#                     type = "response",
#                     newx = test_matrix
#                     )
# 
# predictor %>% class()

test_roc$auc

ggroc(train_roc)
ggroc(test_roc)

# predict value res distribution ----
lasso.prob <- predict(lasso_best,
                      type="response",
                      newx = Test,
                      # s = 'lambda.min'
)



density_p <- cbind(test_set$clade, lasso.prob) %>% as.data.frame()
# str(density_p)

for (x in rownames(density_p)) {
  if (density_p[x, "V1"] == 0) {
    density_p[x, "V1"] <- "non-B2"
  } else {
    density_p[x, "V1"] <- "B2"
  }
}
str(density_p)

density_p$V1 <- as.factor(density_p$V1)

str(density_p)
levels(density_p$V1)
ggdensity(density_p, 
          x = "s0",
          color = "V1")

# hosmer lemeshow ----
# write our own ----


binned <- reg_tab %>%
  mutate(predprob = predict(lasso_best, type = "response", newx = as.matrix(reg_tab[, 1: 260])), 
         linpred = predict(lasso_best, type = "link", newx = as.matrix(reg_tab[, 1: 260])),
         bin = cut(linpred, breaks = unique(quantile(linpred, (1:25)/26, type = 1)
         )
         )
  )
# binned$predprob <- binned$predprob/max(binned$predprob)

binned <- group_by(binned, bin)

binned <- binned %>% summarise(y       = sum(ifelse(clade == 1, 1, 0)),
                               avgpred = mean(predprob),
                               count   = n()
)
binned <- mutate(binned, se_fit = sqrt(avgpred * (1 - avgpred) / count))

# sum(binned$count)

binned %>%
  ggplot(mapping = aes(x = avgpred, y = y / count)) + 
  geom_point() +
  geom_linerange(mapping = aes(ymin = y / count - 2 * se_fit,
                               ymax = y / count + 2 * se_fit), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Predicted probability", y = "Observed proportion")

(hlstat <- with(binned, sum((y - count * avgpred)^2 / (count * avgpred * (1 - avgpred)
)
)
)
)
nrow(binned)
binned[order(binned$avgpred), ]

pchisq(hlstat, nrow(binned) - 1, lower.tail = FALSE)# good


# which genes are of high correlation ----

length(coef(lasso_best))

a <- coef(lasso_best)
str(a)
a
a[1]
a[2] == 0
length(a)

a@Dimnames[[1]][1]
length(a@Dimnames[[1]])

for (x in seq(1, length(a))) {
  if (a[x] != 0) {
    print(a@Dimnames[[1]][x])
    print(a[x])
  }
}

