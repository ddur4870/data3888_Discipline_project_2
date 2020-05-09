library(GEOquery)  ## go to https://github.com/seandavi/GEOquery for installation details
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)


F1score <- function(mat) {
      TN <- mat[2,2]
      FP <- mat[1,2]
      TP <- mat[1,1]
      FN <- mat[2,1]
      return(2*TP/(2*TP+FP+FN))

}

false_pos <- function(mat) {
      TN <- mat['No','No']
      FP <- mat['Yes','No']
      return(FP/(FP+TN))

}

gse <- read.csv("GSE120396_expression_matrix.txt",header = TRUE,row.names = 'X')
rejection_status <- read.delim("reactions_120396.txt")
rejection_status1 = rejection_status$Status

gse_pca <- prcomp(t(gse))
df_toplot <- data.frame(rejection_status1, 
                        pc1 = gse_pca$x[,1], pc2 = gse_pca$x[,2]  )


get_knn = function(num_folds, k_val, num_repeats){


largevar = apply(gse, 1, var)
ind = which(largevar > quantile(largevar, 0.9))

X = as.matrix(t(gse[ind,]))
y = rejection_status1

cvK = num_folds  # number of CV folds
cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
cv_acc_knn = cv_acc_svm = cv_acc_rf = c()

n_sim = num_repeats ## number of repeats
for (i in 1:n_sim) {
  
  cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
  cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
  
  for (j in 1:cvK) {
    test_id = cvSets$subsets[cvSets$which == j]
    X_test = X[test_id, ]
    X_train = X[-test_id, ]
    y_test = y[test_id]
    y_train = y[-test_id]
    
    ## KNN
    fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = k_val)
    cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
  }
  cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
} 

  return(cv_50acc5_knn)
}


get_svm = function(num_folds, num_repeats){

largevar = apply(gse, 1, var)
ind = which(largevar > quantile(largevar, 0.9))

X = as.matrix(t(gse[ind,]))
y = rejection_status1

cvK = num_folds  # number of CV folds
cv_50acc5_svm = c()
cv_acc_svm = c()

n_sim = num_repeats ## number of repeats
for (i in 1:n_sim) {
  
  cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
  cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
  
  for (j in 1:cvK) {
    test_id = cvSets$subsets[cvSets$which == j]
    X_test = X[test_id, ]
    X_train = X[-test_id, ]
    y_test = y[test_id]
    y_train = y[-test_id]
    
    
    ## SVM
    svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
    fit <- predict(svm_res, X_test)
    cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
    
  }
  
  cv_50acc5_svm <- append(cv_50acc5_svm, mean(cv_acc_svm))
  
} 

return(cv_50acc5_svm)


}

get_rf = function(num_folds, num_repeats){
  
  largevar = apply(gse, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  X = as.matrix(t(gse[ind,]))
  y = rejection_status1
  
  cvK = num_folds  
  cv_50acc5_svm = c()
  cv_acc_svm = c()
  cv_50acc5_rf = c()
  cv_acc_rf = c()
  
  n_sim = num_repeats 
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK) 
    cv_acc_rf = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
      fit <- predict(rf_res, X_test)
      cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
      
    }
    
    
    cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
    
  } 
  
  return(cv_50acc5_rf)
  
  
}

get_knn_accuracy = function(num_folds, k_val, num_repeats){
  
  largevar = apply(gse, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  X = as.matrix(t(gse[ind,]))
  y = rejection_status1
  
  cvK = num_folds 
  cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
  total_cv_acc_knn = c()
  
  n_sim = num_repeats
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK)
    cv_acc_knn = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      ## KNN
      fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = k_val)
      cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
    }
    cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
    total_cv_acc_knn[i] <- mean(cv_acc_knn)
  } 

  return(total_cv_acc_knn)
}

get_svm_acc = function(num_folds, num_repeats){
  
  largevar = apply(gse, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  X = as.matrix(t(gse[ind,]))
  y = rejection_status1
  
  cvK = num_folds
  total_cv_acc_svm = c()
  
  n_sim = num_repeats
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK)
    cv_acc_svm = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      
      ## SVM
      svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
      fit <- predict(svm_res, X_test)
      cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
      
    }
    
    total_cv_acc_svm[i] <- mean(cv_acc_svm)
    
  } 
  
  return(total_cv_acc_svm)
  
  
}

plot_svm_folds1 = function(repeats){
  svm1 = mean(get_svm_acc(2,repeats))
  svm2 = mean(get_svm_acc(3,repeats))
  svm3 = mean(get_svm_acc(4,repeats))
  svm4 = mean(get_svm_acc(5,repeats))
  svm5 = mean(get_svm_acc(6,repeats))
  svm6 = mean(get_svm_acc(7,repeats))
  svm7 = mean(get_svm_acc(8,repeats))
  svm8 = mean(get_svm_acc(9,repeats))
  svm9 = mean(get_svm_acc(10,repeats))
  results_vec = c(svm1,svm2,svm3,svm4,svm5,svm6,svm7,svm8,svm9)
  return(results_vec)
}


get_all_knn_measurements = function(num_folds, k_val, num_repeats){
  
  
  largevar = apply(gse, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  X = as.matrix(t(gse[ind,]))
  y = rejection_status1
  
  cvK = num_folds 
  
  acc_vector = c()
  fn_vector = c()
  fp_vector = c()
  f1_vector = c()
  
  n_sim = num_repeats ## number of repeats
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
    #cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
    cv_acc_knn = c()
    cv_fn = c()
    cv_fp = c()
    cv_f1 = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      k_acc = c()
      k_fn = c()
      k_fp = c()
      k_f1 = c()
      
      ## KNN
      fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = k_val)
      cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
      k_acc = c(k_acc, table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test)))
      k_fn = c(k_fn, false_neg(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
      k_f1 = c(k_f1, F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
      
      k_fp = c(k_fp, false_pos(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
      cv_acc_knn = rbind(cv_acc_knn, k_acc)
      cv_fn = rbind(cv_fn, k_fn)
      cv_fp = rbind(cv_fp, k_fp)
      cv_f1 = rbind(cv_f1, k_f1)
    }
    
    acc_vector <- rbind(acc_vector, colMeans(cv_acc_knn))
    
    
    fn_vector  <- rbind(fn_vector, colMeans(cv_fn))
    
    fp_vector  <- rbind(fp_vector, colMeans(cv_fp))
    
    
    f1_vector <- rbind(f1_vector, colMeans(cv_f1))
  }
  return (list(acc_vector=colMeans(acc_vector), acc_vector_2d=acc_vector, fn_vector=colMeans(fn_vector), fn_vector_2d=fn_vector, f1_vector=colMeans(f1_vector), f1_vector_2d=f1_vector, fp_vector_2d=fp_vector))
}

get_all_svm_measurements = function(num_folds, k_val, num_repeats){
  
  largevar = apply(gse, 1, var)
  ind = which(largevar > quantile(largevar, 0.9))
  
  X = as.matrix(t(gse[ind,]))
  y = rejection_status1
  
  cvK = num_folds
  fp_vector = c()
  f1_vector = c()
  
  n_sim = num_repeats
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK)
    cv_fp = c()
    cv_f1 = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
    
      k_fp = c()
      k_f1 = c()
      
      
      svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
      fit5 <- predict(svm_res, X_test)
      k_f1 = c(k_f1, F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))) 
      k_fp = c(k_fp, false_pos(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
      cv_fp = rbind(cv_fp, k_fp)
      cv_f1 = rbind(cv_f1, k_f1)
    }
   
    fp_vector  <- rbind(fp_vector, colMeans(cv_fp))
    f1_vector <- rbind(f1_vector, colMeans(cv_f1))
  }
  return (list(f1_vector=colMeans(f1_vector), f1_vector_2d=f1_vector, fp_vector_2d=fp_vector))
}

knn1 = mean(get_knn_accuracy(2,5, 10))
knn2 = mean(get_knn_accuracy(3,5, 10))
knn3 = mean(get_knn_accuracy(4,5, 10))
knn4 = mean(get_knn_accuracy(5,5, 10))
knn5 = mean(get_knn_accuracy(6,5, 10))
knn6 = mean(get_knn_accuracy(7,5, 10))
knn7 = mean(get_knn_accuracy(8,5, 10))
knn8 = mean(get_knn_accuracy(9,5, 10))
knn9 = mean(get_knn_accuracy(10,5, 10))
knn_acc_folds = c(knn1,knn2,knn3,knn4,knn5,knn6,knn7,knn8,knn9)

knn11 = mean(get_knn_accuracy(5,1, 10))
knn12 = mean(get_knn_accuracy(5,2, 10))
knn13 = mean(get_knn_accuracy(5,3, 10))
knn14 = mean(get_knn_accuracy(5,4, 10))
knn15 = mean(get_knn_accuracy(5,5, 10))
knn16 = mean(get_knn_accuracy(5,6, 10))
knn17 = mean(get_knn_accuracy(5,7, 10))
knn18 = mean(get_knn_accuracy(5,8, 10))
knn19 = mean(get_knn_accuracy(5,9, 10))
knn_acc_neighbours = c(knn11,knn12,knn13,knn14,knn15,knn16,knn17,knn18,knn19)

knn21 = mean(get_knn_accuracy(5,5, 5))
knn22 = mean(get_knn_accuracy(5,5, 10))
knn23 = mean(get_knn_accuracy(5,5, 15))
knn24 = mean(get_knn_accuracy(5,5, 20))
knn_acc_sims = c(knn21,knn22,knn23,knn24)

svm21 = mean(get_svm_acc(5, 5))
svm22 = mean(get_svm_acc(5, 10))
svm23 = mean(get_svm_acc(5, 15))
svm24 = mean(get_svm_acc(5, 20))
svm_acc_sims = c(svm21,svm22,svm23,svm24)

svm1 = mean(get_svm_acc(2,5))
svm2 = mean(get_svm_acc(3,5))
svm3 = mean(get_svm_acc(4,5))
svm4 = mean(get_svm_acc(5,5))
svm5 = mean(get_svm_acc(6,5))
svm6 = mean(get_svm_acc(7,5))
svm7 = mean(get_svm_acc(8,5))
svm8 = mean(get_svm_acc(9,5))
svm9 = mean(get_svm_acc(10,5))
svm_acc_folds = c(svm1,svm2,svm3,svm4,svm5,svm6,svm7,svm8,svm9)


