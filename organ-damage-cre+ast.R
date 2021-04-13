shhh <- suppressPackageStartupMessages

shhh(library(mice))
shhh(library(ROCR))
shhh(library(boot))
shhh(library(rms))
shhh(library(ResourceSelection))
shhh(library(LogisticDx))
shhh(library(VIM))

src_data_file <- '../data/2011-15/data.csv'
src_dict_file <- '../data/2011-15/dictionary.csv'

variables <- c('OUT', 'AGE', 'SCNS', 'BLDING', 'AST', 'CRE')
yvar <- 'OUT'
model_string <- 'OUT~AGE+SCNS+BLDING+AST+CRE'
model_name <- 'organ-damage-cre+ast' 
model_folder <- paste0('models/', model_name)

# Number of Multiple imputations
num_imp <- 100
# Number of bootstrap sample
num_boot <- 1000

src_data <- read.table(src_data_file, sep=",", header=TRUE, na.strings="\\N")
src_data <- src_data[variables]

imp_data <- mice(src_data, m=num_imp)
var_drop <- c(".imp", ".id")
imp_data_files <- character(0)
for (iter in 1:num_imp) {
    comp_data <- complete(imp_data, action=iter)  
    comp_data <- comp_data[,!(names(comp_data) %in% var_drop)]
    fn <- paste0(model_folder, '/imp', "/imputation-", iter, ".csv")
    write.csv(comp_data, file=fn, row.names=FALSE)
    imp_data_files <- c(imp_data_files, fn)
}

imp_models <- with(imp_data, glm(family="binomial", 
                                 formula=OUT~AGE+SCNS+BLDING+AST+CRE))

poolmod <- pool(imp_models)
print(summary(poolmod))

sink(paste0(model_folder, "/model.txt"), append=FALSE, split=FALSE)
print(summary(poolmod))
sink()

# Use bootstrap for internal validation

# From Chapter 5 of Analysis of Categorical Data with R:
# http://www.chrisbilder.com/categorical/Chapter5/AllGOFTests.R 
stukel.test <- function(obj) {
    # first, check to see if we fed in the right kind of object
    stopifnot(family(obj)$family == "binomial" && family(obj)$link == "logit")
    high.prob <- (obj$fitted.values >= 0.5) 
    logit2 <- obj$linear.predictors^2
    z1 = 0.5*logit2*high.prob
    z2 = 0.5*logit2*(1-high.prob)
    mf <- obj$model
    trials = rep(1, times = nrow(mf))
    if(any(colnames(mf) == "(weights)")) 
        trials <- mf[[ncol(mf)]]
    prop = mf[[1]]
    # the double bracket (above) gets the index of items within an object
    if (is.factor(prop)) 
        prop = (as.numeric(prop) == 2)  # Converts 1-2 factor levels to logical 0/1 values
    pi.hat = obj$fitted.values 
    y <- trials*prop
    exclude <- which(colnames(mf) == "(weights)")
    vars <- data.frame(z1, z2, y, mf[,-c(1,exclude)])
    full <- glm(formula = y/trials ~ ., family = binomial(link = logit), weights = trials, data = vars)
    null <- glm(formula = y/trials ~ ., family = binomial(link = logit), weights = trials, data = vars[,-c(1,2)])
    LRT <- anova(null,full)
    p.value <- 1 - pchisq(LRT$Deviance[[2]], LRT$Df[[2]])
    return(p.value)
}

# (Adjusted) McFadden R2
# In the future could use this library for calculation
# http://www.inside-r.org/packages/cran/bayloredpsych/docs/PseudoR2        
adjr2 <- function(obj) {
    # For the time being, just get numer of dofs in model (including intercept) 
    # using LogLik: http://stats.stackexchange.com/a/5580
    ll <- logLik(obj)
    K <- attr(ll, "df")
    r2 <- 1 - ((obj$deviance - K) / obj$null.deviance)
    return(r2)
}
        
calib <- function(probs,outcome,nbins=10) {
    c <- 0.0

    # Construct bins
    judgement_bins <- seq(0, nbins)/nbins

    # Which bin is each prediction in?
    bin_num <- .bincode(probs, judgement_bins, TRUE)

    for (j_bin in sort(unique(bin_num))) {
        # Is event in bin
        in_bin <- bin_num == j_bin
        
        # Predicted probability taken as average of preds in bin
        predicted_prob <- mean(probs[in_bin])
        
        # How often did events in this bin actually happen?
        true_bin_prob <- mean(outcome[in_bin])
        
        # Squared distance between predicted and true times num of obs
        c <- c + sum(in_bin) * (predicted_prob - true_bin_prob)^2
    } 
    cal <- c / length(probs)
    return(cal)
}
        
brier <- function(probs,outcome) {
    res <- mean((probs - outcome)^2)
    return(res)
}
 
accu <- function(probs,outcome) {
    preds = 0.5 <= probs
    res <- 1 - mean(abs(preds - outcome))
    return(res)
}        
       
# Transform Z-scores back to score, and calculates CI at 95%
# https://stats.idre.ucla.edu/stata/faq/how-can-i-estimate-r-squared-for-a-model-estimated-with-multiply-imputed-data/
zinv <- function(values, N, M) {
    # Fist, we need the inter-imputation variance
    B <- sum((values - mean(values))^2) / (M - 1)

    # Now, we get the MI estimate of the variance of z
    V <- 1/(N-3) + B/(M+1)
     
    # The confidence interval, using the confidence level for 95%    
    Q <- mean(values)  
    ci_min <- tanh(Q - 1.959964*sqrt(V*Q))^2
    ci_max <- tanh(Q + 1.959964*sqrt(V*Q))^2
    val_mean <- tanh(Q)^2
    
    res <- c(val_mean, ci_min, ci_max)
    return(res)
}
     
optim <- function(src_dat, boot_idx) {
    src_idx <- 1:nrow(src_dat)
    boot_idx <- sample(src_idx, replace=TRUE)
    boot_dat <- src_dat[boot_idx,]

    boot_y <- as.matrix(boot_dat[,1])
    boot_x <- as.matrix(boot_dat[,2:ncol(boot_dat)])
  
    boot_mod <- glm(family="binomial", formula=mod_formula, data=boot_dat)

    # Get the indices of the rows not used in the bootstrap sample (the .632 method)
    rem_idx <- setdiff(src_idx, boot_idx)
    rem_dat <- train_data[rem_idx,] 
    rem_x <- as.matrix(rem_dat[,2:ncol(rem_dat)])
    rem_y <- as.matrix(rem_dat[,1])
    
    boot_prob <- predict(boot_mod, boot_dat, type="response")
    boot_pred <- prediction(boot_prob, boot_y)
    boot_auc <- performance(boot_pred, measure = "auc")

    rem_prob <- predict(boot_mod, rem_dat, type="response")
    rem_pred <- prediction(rem_prob, rem_y)
    rem_auc <- performance(rem_pred, measure = "auc")    
    rem_bri <- brier(rem_prob, rem_y)
    rem_cal <- calib(rem_prob, rem_y)
    rem_acc <- accu(rem_prob, rem_y)
    
    # All values are returned as Z-scores using the method from 
    # https://stats.idre.ucla.edu/stata/faq/how-can-i-estimate-r-squared-for-a-model-estimated-with-multiply-imputed-data/
    auc_value <- atanh(sqrt(rem_auc@y.values[[1]]))
    bri_value <- atanh(sqrt(rem_bri))
    cal_value <- atanh(sqrt(rem_cal))
    acc_value <- atanh(sqrt(rem_acc))    
    r2_value <- atanh(sqrt(adjr2(boot_mod)))    
    
    res <- c(auc_value, bri_value, cal_value, acc_value, r2_value)
    return(res)     
}

auc_app_values <- vector(mode="numeric", length=length(imp_data_files))
auc_values <- vector(mode="numeric", length=length(imp_data_files))
bri_values <- vector(mode="numeric", length=length(imp_data_files)) 
cal_values <- vector(mode="numeric", length=length(imp_data_files))
acc_values <- vector(mode="numeric", length=length(imp_data_files))
r2_values <- vector(mode="numeric", length=length(imp_data_files))     

N <- 0
M <- length(imp_data_files) 
imp_iter <- 0
for (fn in imp_data_files) {
    imp_iter <- imp_iter + 1
    train_data <- read.table(fn, sep=",", header=TRUE)    
    N <- nrow(train_data)
    yvalues <- train_data[yvar]
    mod_formula <- as.formula(model_string)
    model <- glm(family="binomial", formula=mod_formula, data=train_data)

    prob <- predict(model, train_data)
    pred <- prediction(prob, yvalues)
    auc <- performance(pred, measure = "auc")
    auc_app <- auc@y.values[[1]]
    
    bootres <- boot(train_data, optim, R=num_boot, parallel="multicore", ncpus=4)
        
    auc_app_values[imp_iter] <- atanh(sqrt(auc_app))
    auc_values[imp_iter] <- bootres$t[,1]    
    bri_values[imp_iter] <- bootres$t[,2]  
    cal_values[imp_iter] <- bootres$t[,3]
    acc_values[imp_iter] <- bootres$t[,4]    
    r2_values[imp_iter] <- bootres$t[,5]
}
     
auc_app_mean <- zinv(auc_app_values, N, M)
auc_mean <- zinv(auc_values, N, M)
bri_mean <- zinv(bri_values, N, M)
cal_mean <- zinv(cal_values, N, M)
acc_mean <- zinv(acc_values, N, M)
r2_mean <- zinv(r2_values, N, M)

print(sprintf("Apparent AUC : %0.3f 95 CI: %0.3f, %0.3f", auc_app_mean[1], auc_app_mean[2], auc_app_mean[3]))
print(sprintf("Corrected AUC: %0.3f 95 CI: %0.3f, %0.3f", auc_mean[1], auc_mean[2], auc_mean[3]))
print(sprintf("Brier score  : %0.3f 95 CI: %0.3f, %0.3f", bri_mean[1], bri_mean[2], bri_mean[3]))
print(sprintf("Calibration  : %0.3f 95 CI: %0.3f, %0.3f", cal_mean[1], cal_mean[2], cal_mean[3]))
print(sprintf("Accuracy     : %0.3f 95 CI: %0.3f, %0.3f", acc_mean[1], acc_mean[2], acc_mean[3]))        
print(sprintf("Adjusted R2  : %0.3f 95 CI: %0.3f, %0.3f", r2_mean[1], r2_mean[2], r2_mean[3]))
     
sink(paste0(model_folder, "/boot.txt"), append=FALSE, split=FALSE)
print(sprintf("Apparent AUC : %0.3f 95 CI: %0.3f, %0.3f", auc_app_mean[1], auc_app_mean[2], auc_app_mean[3]))
print(sprintf("Corrected AUC: %0.3f 95 CI: %0.3f, %0.3f", auc_mean[1], auc_mean[2], auc_mean[3]))
print(sprintf("Brier score  : %0.3f 95 CI: %0.3f, %0.3f", bri_mean[1], bri_mean[2], bri_mean[3]))
print(sprintf("Calibration  : %0.3f 95 CI: %0.3f, %0.3f", cal_mean[1], cal_mean[2], cal_mean[3]))
print(sprintf("Accuracy     : %0.3f 95 CI: %0.3f, %0.3f", acc_mean[1], acc_mean[2], acc_mean[3]))        
print(sprintf("Adjusted R2  : %0.3f 95 CI: %0.3f, %0.3f", r2_mean[1], r2_mean[2], r2_mean[3]))
sink()