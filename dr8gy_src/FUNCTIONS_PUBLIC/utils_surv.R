generate_sub_model <- function(DATA,X_outcome,X_treatment,X_censure,Y.grid,tau,type_of_model="nonparametric"){
  # The models
  Y.diff <- diff(c(0, Y.grid))
  Y_all <- DATA$T_obs
  Y.grid_extended <- sort(unique(Y_all))[1:(length(Y.grid)+1)]
  Y.grid_extended[1:length(Y.grid)] <- Y.grid
  
  
  #We are going in this function to compute the different models of the nuisances parameters: the outcome, the censoring and the treatment. 
  # Here a summary of the different models we are trying to estimate:
  
  # Outcome model:
  ## S_hat
  ## h_hat
  ## Q.Y.hat
  ## Q.t.hat
  ## S_hat0
  ## S_hat1
#All of the models are derived from the models S_hat computed first with the function cph or the functioon survival_forest from the grf library
  
  # Censoring model:
  ## S_C_hat
  ## h_C_hat
  ## S_C_hat0
  ## S_C_hat1
#All of the models are derived from the models S_C_hat computed first with the function cph or the functioon survival_forest from the grf library

  
  # Treatment model:
  ## e_hat 
# This is estimated using either a logistic regression or the function regression_forest from the grf library

  
  # model S:
  if (type_of_model == "parametric"){
  
    # ## Cox:
    outcome <- 'Surv(T_obs,status)'
    
    DATA0 <- DATA %>% filter(A == 0)
    DATA1 <- DATA %>% filter(A == 1)
    
    f <- as.formula(paste(outcome, paste(c(X_outcome,'A'), collapse = " + "), sep = " ~ "))
    #f <- as.formula(paste(outcome, paste(c(X_outcome), collapse = " + "), sep = " ~ "))
    fitS <- cph(f,data=DATA,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE,surv=TRUE)
    # fitS0 <- cph(f,data=DATA0,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE)
    # fitS1 <- cph(f,data=DATA1,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE)
    fit.pred <- predictCox(fitS, newdata=DATA, times=Y.grid , type = "survival")
    S_hat <- fit.pred$survival
    
    DATA.1 <- DATA
    DATA.1$A <- 1
    DATA.0 <- DATA
    DATA.0$A <- 0
    
    
    fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
    fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
    S_hat1 <- fit.pred1$survival
    S_hat0 <- fit.pred0$survival
    S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
  }
  
  #Forest:
  if (type_of_model == "nonparametric"){
    categorical_name <- names(which(sapply(subset(DATA, select = c(X_outcome)), class) == "factor"))
  numerical_name <- setdiff(X_outcome,categorical_name)
  options(na.action='na.pass')
  X_one_hot <- model.matrix(~ 0 + ., DATA[categorical_name], na.action = "na.pass")
  categorical_name_one_hot <- names(as.data.frame(X_one_hot))
  replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
  categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
  DATA[categorical_name_one_hot] <- X_one_hot

  X_outcome_one_hot <- union(categorical_name_one_hot,numerical_name)
  #
  #
  # ## Forest T learner:
  X <- as.matrix(DATA[X_outcome_one_hot])
  DATA_grf_complete0 <- DATA %>% filter(A == 0)
  X_0 <- as.matrix(DATA_grf_complete0[X_outcome_one_hot])
  Y_0 <- as.matrix(DATA_grf_complete0[,'T_obs'])
  D_0 <- as.matrix(DATA_grf_complete0[,'status'])


  DATA_grf_complete1 <- DATA %>% filter(A == 1)
  X_1 <- as.matrix(DATA_grf_complete1[X_outcome_one_hot])
  Y_1 <- as.matrix(DATA_grf_complete1[,'T_obs'])
  D_1 <- as.matrix(DATA_grf_complete1[,'status'])

  s.forest0 <- survival_forest(X_0, Y_0, D_0,Y.grid)
  s.forest1 <- survival_forest(X_1, Y_1, D_1,Y.grid)

  S_hat0 <- predict(s.forest0, X,failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  S_hat0 <- S_hat0[,c(1:length(Y.grid))]
  S_hat1 <- predict(s.forest1, X,failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  S_hat1 <- S_hat1[,c(1:length(Y.grid))]
  S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
  }
  
  # ## Forest S learner:
  # X_A <- as.matrix(subset(DATA, select = c(X_outcome_one_hot,"A")))
  # X <- as.matrix(subset(DATA, select = c(X_outcome_one_hot)))
  # Y <- as.matrix(DATA[,'T_obs'])
  # D <- as.matrix(DATA[,'status'])
  # s.forest <- survival_forest(X_A, Y, D,Y.grid)
  # 
  # #X0 <- as.matrix(cbind(DATA[c(X_outcome)],0))
  # #X1 <- as.matrix(cbind(DATA[c(X_outcome)],1))
  # 
  # #S_hat <- predict(s.forest, failure.times = Y.grid)$predictions
  # S_hat <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_hat <- S_hat[,c(1:length(Y.grid))]
  # 
  # #X.orig <- s.forest[["X.orig"]]
  # s.forest[["X.orig"]] <- cbind(X, 0)
  # #S_hat0 <- predict(s.forest)$predictions
  # S_hat0 <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_hat0 <- S_hat0[,c(1:length(Y.grid))]
  # #S_hat0 <- predict(s.forest,newdata= X0, failure.times = Y.grid)$predictions
  # 
  # #X.orig <- s.forest[["X.orig"]]
  # s.forest[["X.orig"]] <- cbind(X, 1)
  # #S_hat1 <- predict(s.forest)$predictions
  # S_hat1 <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_hat1 <- S_hat1[,c(1:length(Y.grid))]
  # #S_hat1 <- predict(s.forest,newdata= X1, failure.times = Y.grid)$predictions
  
  #E(T|X,W) is computed through an integration of S_hat
  
  # model h:
  
  
  #fit.pred <- predictCox(fitS, newdata=DATA, times=Y.grid , type = "hazard")
  #h_hat <- fit.pred$survival
  
  
  grid.length <-  ncol(S_hat)
  log.surv <- -log(base::cbind(1, S_hat))
  h_hat <- log.surv[, 2:(grid.length + 1)] - log.surv[, 1:grid.length]
  h_hat <- sweep(h_hat, 2, Y.diff, "/")
  
  
  # model S_c:
  

  X_confoudingC = intersect(X_outcome,X_censure)
  DATA_C <- DATA
  DATA_C$status <- 1- DATA_C$status
  
  if (type_of_model == "parametric"){
    
    
    #DATA_C$diagnostic_code <- fct_lump_min(DATA_C$diagnostic_code,min=30,other_level = "Other")
    #DATA_C$initial_admission_unit <- fct_lump_min(DATA_C$initial_admission_unit,min=30,other_level = "Other")
    #DATA_C$ncm_chronic_disease <- fct_lump_min(DATA_C$ncm_chronic_disease,min=30,other_level = "Other") 
    
  ## Cox
    DATA_C0 <- DATA_C %>% filter(A == 0)
    DATA_C1 <- DATA_C %>% filter(A == 1)
    DATA_C0 <- droplevels(DATA_C0)
    DATA_C1 <- droplevels(DATA_C1)
    
    
    outcome <- 'Surv(T_obs,status)'
    #f <- as.formula(paste(outcome, paste(c(X_censure), collapse = " + "), sep = " ~ "))
    f <- as.formula(paste(outcome, paste(c(X_censure,'A'), collapse = " + "), sep = " ~ "))
    fitC <- cph(f,data=DATA_C,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE)
    #fitC0 <- cph(f,data=DATA_C0,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE)
    #fitC1 <- cph(f,data=DATA_C1,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE)
    fitSC.pred <- predictCox(fitC, newdata=DATA_C, times=Y.grid , type = "survival")
    S_C_hat <- fitSC.pred$survival
    
    DATA_C.1 <- DATA_C
    DATA_C.1$A <- 1
    DATA_C.0 <- DATA_C
    DATA_C.0$A <- 0
    
    fit.pred1 <- predictCox(fitC, newdata=DATA_C.1, times=Y.grid , type = "survival")
    fit.pred0 <- predictCox(fitC, newdata=DATA_C.0, times=Y.grid , type = "survival")
    S_C_hat1 <- fit.pred1$survival
    S_C_hat0 <- fit.pred0$survival
    S_C_hat <-S_C_hat1*(DATA_C$A) + (1-DATA_C$A)*S_C_hat0
  }
  
  if (type_of_model == "nonparametric"){
  ## Forest T learner:
  X <- as.matrix(DATA_C[X_outcome_one_hot])
  DATA_C_grf_complete0 <- DATA_C %>% filter(A == 0)
  X_0 <- as.matrix(DATA_C_grf_complete0[X_outcome_one_hot])
  Y_0 <- as.matrix(DATA_C_grf_complete0[,'T_obs'])
  D_0 <- as.matrix(DATA_C_grf_complete0[,'status'])


  DATA_C_grf_complete1 <- DATA_C %>% filter(A == 1)
  X_1 <- as.matrix(DATA_C_grf_complete1[X_outcome_one_hot])
  Y_1 <- as.matrix(DATA_C_grf_complete1[,'T_obs'])
  D_1 <- as.matrix(DATA_C_grf_complete1[,'status'])

  s.forest0 <- survival_forest(X_0, Y_0, D_0,Y.grid)
  s.forest1 <- survival_forest(X_1, Y_1, D_1,Y.grid)

  S_C_hat0 <- predict(s.forest0, X,failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  S_C_hat0 <- S_C_hat0[,c(1:length(Y.grid))]

  S_C_hat1 <- predict(s.forest1, X,failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  S_C_hat1 <- S_C_hat1[,c(1:length(Y.grid))]
  S_C_hat <- S_C_hat1*DATA_C$A + (1-DATA_C$A)*S_C_hat0
}
  
  
  # ##Forest S learner
  # 
  # X_A <- as.matrix(subset(DATA_C, select = c(X_outcome_one_hot,"A")))
  # X <- as.matrix(subset(DATA_C, select = c(X_outcome_one_hot)))
  # Y <- as.matrix(DATA_C[,'T_obs'])
  # D <- as.matrix(DATA_C[,'status'])
  # s.forest <- survival_forest(X_A, Y, D,Y.grid)
  # S_C_hat <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_C_hat <- S_C_hat[,c(1:length(Y.grid))]
  # 
  # 
  # #X.orig <- s.forest[["X.orig"]]
  # s.forest[["X.orig"]] <- cbind(X, 0)
  # #S_C_hat0 <- predict(s.forest)$predictions
  # S_C_hat0 <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_C_hat0 <- S_C_hat0[,c(1:length(Y.grid))]
  # #S_hat0 <- predict(s.forest,newdata= X0, failure.times = Y.grid)$predictions
  # 
  # #X.orig <- s.forest[["X.orig"]]
  # s.forest[["X.orig"]] <- cbind(X, 1)
  # #S_C_hat1 <- predict(s.forest)$predictions
  # S_C_hat1 <- predict(s.forest, failure.times = c(Y.grid,max(Y.grid)+1))$predictions
  # S_C_hat1 <- S_C_hat1[,c(1:length(Y.grid))]
  # #S_hat1 <- predict(s.forest,newdata= X1, failure.times = Y.grid)$predictions
  
  
  
  # model h_C:
  ## cox
  #fithC.pred <- predictCox(fitC, newdata=DATA_c, times=Y.grid , type = "hazard")
  #h_C_hat <- fithC.pred$hazard
  
  ## Diff
  # The conditional hazard function lambda.C.hat = -d/dt log(C.hat(t, x, w))
  # This simple forward difference approximation works reasonably well.
  
  grid.length <-  ncol(S_C_hat)
  log.surv.C <- -log(base::cbind(1, S_C_hat))
  h_C_hat <- log.surv.C[, 2:(grid.length + 1)] - log.surv.C[, 1:grid.length]
  h_C_hat <- sweep(h_C_hat, 2, Y.diff, "/")
  
  
  
  # model e:
  X_confouding = intersect(X_outcome,X_treatment)
  
  # ## GLM
  if (type_of_model == "parametric"){
  outcome <- 'A'
  f <- as.formula(paste(outcome, paste(c(X_confouding), collapse = " + "), sep = " ~ "))
  fitA <- glm(f,data = DATA,family = binomial(link="logit"))
  e_hat <- predict(fitA,newdata=DATA,type="response")
  
  
  ## regularized GLM
  # DATA <- as.data.frame(DATA)
  # 
  # categorical_name <- names(which(sapply(subset(DATA, select = c(X_confouding)), class) == "factor"))
  # numerical_name <- setdiff(X_confouding,categorical_name)
  # options(na.action='na.pass')
  # X_one_hot <- model.matrix(~ 0 + ., DATA[categorical_name], na.action = "na.pass")
  # categorical_name_one_hot <- names(as.data.frame(X_one_hot))
  # replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
  # categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
  # DATA[categorical_name_one_hot] <- X_one_hot
  # 
  # X_confouding_one_hot <- union(categorical_name_one_hot,numerical_name)
  # 
  # 
  # #xtest <- c("csabad_sofa","initial_admission_unit")
  # x = as.matrix(DATA[X_confouding_one_hot])
  # y = DATA$A
  # 
  # #fitA <- glmnet(x=X, y=y, family = "binomial")
  # model <-glmnet(x, y, family = "binomial", alpha = 1)
  # e_hat <-  model %>% predict(newx = x)
  # e_hat<- e_hat[,1]
  # e_hat_norm <- DATA$A*(e_hat*mean(DATA$A/e_hat)) + (1-DATA$A)*(1-(1-e_hat)*(mean((1-DATA$A)/(1-e_hat))))
  }
  
  
  
  if (type_of_model == "nonparametric"){
  ## Forest
  DATA <- as.data.frame(DATA)
  categorical_name <- names(which(sapply(subset(DATA, select = c(X_confouding)), class) == "factor"))
  numerical_name <- setdiff(X_confouding,categorical_name)
  options(na.action='na.pass')
  X_one_hot <- model.matrix(~ 0 + ., subset(DATA, select = c(categorical_name)), na.action = "na.pass")
  categorical_name_one_hot <- names(as.data.frame(X_one_hot))
  replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
  categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
  DATA[categorical_name_one_hot] <- X_one_hot

  X_confouding_one_hot <- union(categorical_name_one_hot,numerical_name)

  Xipw <- as.matrix(DATA[X_confouding_one_hot])
  Wipw <- as.matrix(DATA$A)
  forest.W <- regression_forest(Xipw, Wipw, tune.parameters = "all")
  e_hat <- predict(forest.W)$predictions
  }
  
  # model Q.t.hat
  ## E(T_i | X_i,A_i, T_i>y) = y + integrale (S(t+y)/S(y))
  
  
  grid.length <- length(Y.grid)
  num.samples <- nrow(DATA) 
  Y.relabeled <- findInterval(DATA$T_obs_tau, Y.grid)
  
  Q.t.hat <- matrix(0, num.samples, grid.length)
  dot.products <- sweep(S_hat[, 1:(grid.length - 1)], 2, Y.diff[2:grid.length], "*")
  Q.t.hat[, 1] <- rowSums(dot.products)
  for (i in 2:(grid.length - 1)) {
    Q.t.hat[, i] <- Q.t.hat[, i - 1] - dot.products[, i - 1]
  }
  Q.t.hat <- Q.t.hat / S_hat
  Q.t.hat[is.infinite(Q.t.hat)] <- 0 # The points where S.hat = 0
  Q.t.hat <- sweep(Q.t.hat, 2, Y.grid, "+") # Add back t
  Q.t.hat[, grid.length] <- max(Y.grid)
  
  # Get Q.Y.hat = Q(Y, X) = E[T | X, W, T >= T_obs]
  
  
  Q.Y.hat <- Q.t.hat[base::cbind(1:num.samples, Y.relabeled)]
  
  return(list("S_hat" =S_hat,"h_hat"= h_hat,"S_C_hat"=S_C_hat,"h_C_hat"=h_C_hat,"e_hat"=e_hat,"Q.Y.hat"=Q.Y.hat,"Q.t.hat"=Q.t.hat,"S_hat0"=S_hat0,"S_hat1"=S_hat1,"S_C_hat0" =S_C_hat0,"S_C_hat1"=S_C_hat1))
}


## Methods:


#CUT:

## CUT_outcome


CUT_outcome <- function(DATA,Y.grid,Q.Y.hat,h_C_hat,S_C_hat){
  T_pseudo <- DATA$status_tau*DATA$T_obs_tau + (1-DATA$status_tau)* Q.Y.hat
  return(T_pseudo)
}


## CUT_IPCW


CUT_IPCW <- function(DATA,Y.grid,S_C_hat){
  index <- matrix(c(c(1:length(DATA$T_obs_tau))  ,match(DATA$T_obs_tau,Y.grid)), nrow = n , ncol = 2)
  S_C_hat_T_obs_tau <- S_C_hat[index] # here we take for each individual i, S_c(T_obs_i)
  
  T_pseudo <- ((DATA$T_obs_tau * DATA$status_tau)/ S_C_hat_T_obs_tau)
  return(T_pseudo)
}


## CUT_AIPCW


CUT_AIPCW <- function(DATA,Y.grid,Q.t.hat,Q.Y.hat,h_C_hat,S_C_hat){
  
  n <- length(DATA$T_obs_tau)
  index <- matrix(c(c(1:length(DATA$T_obs_tau))  ,match(DATA$T_obs_tau,Y.grid)), nrow = n , ncol = 2)
  S_C_hat_T_obs_tau <- S_C_hat[index]
  
  first_term <- (DATA$T_obs_tau * DATA$status_tau)/ S_C_hat_T_obs_tau
  
  second_term <- (Q.Y.hat * (1-DATA$status_tau))/ S_C_hat_T_obs_tau
  
  #third_term <- expected_survival(h_C_hat / S_C_hat *Q.t.hat,Y.grid)
  num.samples <- nrow(DATA) 
  
  integrand <- sweep( ( (h_C_hat) / S_C_hat )* (Q.t.hat), 2, Y.diff, "*")
  third_term <- rep(0, num.samples)
  Y.relabeled <- findInterval(DATA$T_obs_tau, Y.grid)
  
  for (sample in 1:num.samples) {
    Y.index <- Y.relabeled[sample]
    third_term[sample] <- sum(integrand[sample, 1:Y.index]) 
  }
  
  pseudo_T <- first_term+second_term - third_term
  return(pseudo_T)
}


#Classic methods:

## KM


Unadjusted_Kaplan_Meier <- function(dataset) {
  #dataset$status_tau<- as.numeric((dataset$T_obs>=tau) | (dataset$T_obs<tau &  dataset$status == 1 ))
  #dataset$T_obs_tau <- pmin(dataset$T_obs,tau)
  fit <- survfit(Surv(T_obs_tau, status_tau) ~ A, data = dataset)
  # Summarize
  res.sum <- surv_summary(fit, data = dataset)
  res_KM <- attr(res.sum, "table")
  
  mean_naive <- (res_KM$`*rmean`[2] - res_KM$`*rmean`[1])
  return(mean_naive)
}



Adjusted_Kaplan_Meier<- function(dataset,X_confouding,e_hat) {
  # dataset$status_tau<- as.numeric((dataset$T_obs>=tau) | (dataset$T_obs<tau &  dataset$status == 1 ))
  # dataset$T_obs_tau <- pmin(dataset$T_obs,tau)
  # outcome <- 'A'
  # f <- as.formula(
  #   paste(outcome, 
  #         paste(c(X_confouding), collapse = " + "), 
  #         sep = " ~ "))
  # 
  # fit1 <- glm(f, family = binomial(link="logit"), data=dataset)
  # Pr1<-predict(fit1,type="response",newdata=dataset)
  
  W <- (dataset$A == 1) * (1/e_hat) + (dataset$A==0) * (1)/(1-e_hat)
  fitipw <- survfit(Surv(T_obs_tau, status_tau) ~ A, data = dataset,weights =W )
  # Summarize
  resipw.sum <- surv_summary(fitipw, data = dataset)
  res_ipw <- attr(resipw.sum, "table")
  
  
  mean_ipw <- (res_ipw$`*rmean`[2] - res_ipw$`*rmean`[1])
  return(mean_ipw)
}


function_IPCW <- function(DATA,Y.grid,S_C_hat,e_hat){
  
  status <- DATA$status
  ind_status1 <- outer(DATA$T_obs, Y.grid, '>')*status
  ind_status0 <- outer(DATA$T_obs, Y.grid, '>=')*(1-status)
  
  ind<- ind_status1+ind_status0
  A <- DATA$A
  #ind <- outer(DATA$T_obs, Y.grid, '>')*1
  S_1 <- colMeans(ind*(1/S_C_hat)*(DATA$A/e_hat))
  S_0 <- colMeans(ind*(1/S_C_hat)*((1-DATA$A)/(1-e_hat)))
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  E1.hat <- c(c(1, S_1) %*% grid.diff)
  E0.hat <- c(c(1, S_0) %*% grid.diff)
  return(E1.hat- E0.hat)
}




IPTW_IPCW_function <- function(DATA,Y.grid,S_C_hat,e_hat){
  n <- length(DATA$T_obs_tau)

  index <- matrix(c(c(1:length(DATA$T_obs_tau))  ,match(DATA$T_obs_tau,Y.grid)), nrow = n , ncol = 2)
  S_C_hat_T_obs_tau <- S_C_hat[index] # here we take for each individual i, S_c(T_obs_i)
  
  
  # DATA_C <- DATA
  # DATA_C$status <- 1- DATA_C$status
  # fitKMC <- survfit(Surv(T_obs, status) ~ A, data = DATA_C)
  # res.sum <- surv_summary(fitKMC, data = DATA_C)
  # 
  # time1 <- res.sum[res.sum$A == 1,]$time
  # surv1 <- res.sum[res.sum$A == 1,]$surv
  # 
  # time0 <- res.sum[res.sum$A == 0,]$time
  # surv0 <- res.sum[res.sum$A == 0,]$surv
  # S_C_hat_T_obs_tau <-surv1[findInterval(DATA$T_obs_tau,time1)+1]*DATA$A + surv0[findInterval(DATA$T_obs_tau,time0)+1]*(1-DATA$A)
  
  theta_IPTW_IPCW <- mean(((DATA$T_obs_tau * DATA$status_tau)/ S_C_hat_T_obs_tau) * (DATA$A/e_hat - (1-DATA$A)/(1-e_hat)))
  return(theta_IPTW_IPCW)
}




IPTW_outcomeS_function <- function(DATA,Q.Y.hat,e_hat){
  T_pseudo_G_formula_C <- DATA$status_tau*DATA$T_obs_tau + (1-DATA$status_tau)* Q.Y.hat
  
  theta_IPTW_G_formula_C <- mean(T_pseudo_G_formula_C * (DATA$A/e_hat - (1-DATA$A)/(1-e_hat)))
  return(theta_IPTW_G_formula_C)
}



g_formula_cox_S_learner <- function(DATA,X_outcome,Y.grid) {
  
  
  # model S:
  outcome <- 'Surv(T_obs,status)'
  f <- as.formula(paste(outcome, paste(c(X_outcome,'A'), collapse = " + "), sep = " ~ "))
  fitS <- cph(f,data=DATA,y=TRUE,x=TRUE,times = Y.grid,singular.ok=TRUE,surv=TRUE)
  fit.pred <- predictCox(fitS, newdata=DATA, times=Y.grid , type = "survival")
  #S_hat <- fit.pred$survival
  
  
  DATA.1 <- DATA
  DATA.1$A <- 1
  DATA.0 <- DATA
  DATA.0$A <- 0
  
  fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
  fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
  S_hat1 <- fit.pred1$survival
  S_hat0 <- fit.pred0$survival
  
  E_hat1 <- expected_survival(S_hat1,Y.grid)
  E_hat0 <- expected_survival(S_hat0,Y.grid)
  
  
  theta_g_formula <- mean(E_hat1 - E_hat0)
  return(theta_g_formula)
  
}




g_formula_cox_T_learner <- function(DATA,X_outcome,Y.grid) {
  
  outcome <- 'Surv(T_obs,status)'
  DATA0 <- DATA %>% filter(A == 0)
  DATA1 <- DATA %>% filter(A == 1)
  
  #f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  f <- as.formula(paste(outcome, paste(c(X_outcome), collapse = " + "), sep = " ~ "))
  fitS0 <- cph(f,data=DATA0,y=TRUE,x=TRUE,times = Y.grid)
  fitS1 <- cph(f,data=DATA1,y=TRUE,x=TRUE,times = Y.grid)
  #fit.pred <- predictCox(fitS, newdata=DATA, times=Y.grid , type = "survival")
  #S_hat <- fit.pred$survival
  
  DATA.1 <- DATA
  DATA.1$A <- 1
  DATA.0 <- DATA
  DATA.0$A <- 0
  
  
  fit.pred1 <- predictCox(fitS1, newdata=DATA.1, times=Y.grid , type = "survival")
  fit.pred0 <- predictCox(fitS0, newdata=DATA.0, times=Y.grid , type = "survival")
  S_hat1 <- fit.pred1$survival
  S_hat0 <- fit.pred0$survival
  #S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
  
  # outcome <- 'Surv(T_obs,status)'
  # f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  # fitS <- cph(f,data=DATA,y=TRUE,x=TRUE)
  # fit.pred <- predictCox(fitS, newdata=DATA, times=Y.grid , type = "survival")
  # S_hat <- fit.pred$survival
  # 
  # 
  # DATA.1 <- DATA
  # DATA.1$A <- 1
  # DATA.0 <- DATA
  # DATA.0$A <- 0
  # 
  # 
  # fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
  # fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
  # S_hat1 <- fit.pred1$survival
  # S_hat0 <- fit.pred0$survival
  
  E_hat1 <- expected_survival(S_hat1,Y.grid)
  E_hat0 <- expected_survival(S_hat0,Y.grid)
  
  
  theta_g_formula <- mean(E_hat1 - E_hat0)
  return(theta_g_formula)
  
}



G_formula_foret_T_learner <- function(DATA,X_outcome){
  
  DATA <- as.data.frame(DATA)
  
  categorical_name <- names(which(sapply(subset(DATA, select = c(X_outcome)), class) == "factor"))
  numerical_name <- setdiff(X_outcome,categorical_name)
  options(na.action='na.pass')
  X_one_hot <- model.matrix(~ 0 + ., DATA[categorical_name], na.action = "na.pass")
  categorical_name_one_hot <- names(as.data.frame(X_one_hot))
  replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
  categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
  DATA[categorical_name_one_hot] <- X_one_hot
  
  X_outcome_one_hot <- union(categorical_name_one_hot,numerical_name)
  
  
  X <- as.matrix(subset(DATA, select = c(X_outcome_one_hot)))
  
  DATA_grf_complete0 <- DATA %>% filter(A == 0)
  
  X_0 <- as.matrix(subset(DATA_grf_complete0, select = c(X_outcome_one_hot)))
  Y_0 <- as.matrix(DATA_grf_complete0[,'T_obs_tau'])
  D_0 <- as.matrix(DATA_grf_complete0[,'status_tau'])
  
  
  DATA_grf_complete1 <- DATA %>% filter(A == 1)
  X_1 <- as.matrix(subset(DATA_grf_complete1, select = c(X_outcome_one_hot)))
  Y_1 <- as.matrix(DATA_grf_complete1[,'T_obs_tau'])
  D_1 <- as.matrix(DATA_grf_complete1[,'status_tau'])
  
  s.forest0 <- survival_forest(X_0, Y_0, D_0)
  s.forest1 <- survival_forest(X_1, Y_1, D_1)
  
  s.pred0 <- predict(s.forest0, X)
  s.pred1 <- predict(s.forest1, X)
  
  
  times0 <- s.pred0$failure.times
  survrate0<- apply(s.pred0$predictions,2,mean)
  
  times1 <- s.pred1$failure.times
  survrate1<- apply(s.pred1$predictions,2,mean)
  
  
  #The mean survival time in ECD recipients followed-up to 10 years
  rmst0 = rmst(times=times0, surv.rates=survrate0, max.time=tau)
  
  #The mean survival time in SCD recipients followed-up to 10 years
  rmst1 = rmst(times=times1,surv.rates=survrate1, max.time=tau)
  
  return(rmst1 - rmst0)
}



AIPW_G_formula_C_function <- function(DATA,Q.Y.hat,e_hat,S_hat1,S_hat0){
  
  # model S:
  # outcome <- 'Surv(T_obs,status)'
  # f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  # fitS <- cph(f,data=DATA,y=TRUE,x=TRUE)
  # 
  # 
  # 
  # DATA.1 <- DATA
  # DATA.1$A <- 1
  # DATA.0 <- DATA
  # DATA.0$A <- 0
  # 
  # fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
  # fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
  # S_hat1 <- fit.pred1$survival
  # S_hat0 <- fit.pred0$survival
  
  E_hat1 <- expected_survival(S_hat1,Y.grid)
  E_hat0 <- expected_survival(S_hat0,Y.grid)
  
  T_pseudo_G_formula_C <- DATA$status_tau*DATA$T_obs_tau + (1-DATA$status_tau)* Q.Y.hat
  
  
  theta_AIPTW_G_formula_C <- mean(E_hat1 -E_hat0 + DATA$A *( (T_pseudo_G_formula_C - E_hat1 )/ e_hat ) - (1 - DATA$A) * ( (T_pseudo_G_formula_C - E_hat0 )/(1 - e_hat)  ))
  return(theta_AIPTW_G_formula_C)
}




AIPTW_IPCW_function <- function(DATA,Q.Y.hat,e_hat,S_C_hat,S_hat1,S_hat0){
  
  # model S:
  # outcome <- 'Surv(T_obs,status)'
  # f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  # fitS <- cph(f,data=DATA,y=TRUE,x=TRUE)
  # 
  # 
  # 
  # DATA.1 <- DATA
  # DATA.1$A <- 1
  # DATA.0 <- DATA
  # DATA.0$A <- 0
  # 
  # fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
  # fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
  # S_hat1 <- fit.pred1$survival
  # S_hat0 <- fit.pred0$survival
  
  E_hat1 <- expected_survival(S_hat1,Y.grid)
  E_hat0 <- expected_survival(S_hat0,Y.grid)
  
  n <- length(DATA$T_obs_tau)
  index <- matrix(c(c(1:length(DATA$T_obs_tau))  ,match(DATA$T_obs_tau,Y.grid)), nrow = n , ncol = 2)
  S_C_hat_T_obs_tau <- S_C_hat[index]
  
  theta_AIPTW_IPCW <- mean(E_hat1 -E_hat0 + DATA$A *( (((DATA$T_obs_tau * DATA$status_tau)/ S_C_hat_T_obs_tau) - E_hat1 )/ e_hat ) - (1 -  DATA$A) * ( (((DATA$T_obs_tau * DATA$status_tau)/ S_C_hat_T_obs_tau) - E_hat0 )/(1 - e_hat)  ))
  return(theta_AIPTW_IPCW)
}



IPTW_AIPCW_function <- function(DATA,Y.grid,Q.t.hat,Q.Y.hat,h_C_hat,S_C_hat,e_hat){
  
  pseudo_T <- CUT_AIPCW(DATA,Y.grid,Q.t.hat,Q.Y.hat,h_C_hat,S_C_hat)
  
  
  theta_IPTW_AIPCW <- mean(pseudo_T*(DATA$A/e_hat - (1-DATA$A)/(1-e_hat)))
  return(theta_IPTW_AIPCW)
}




AIPTW_AIPCW_function <- function(DATA,Y.grid,Q.t.hat,Q.Y.hat,h_C_hat,S_C_hat,e_hat,S_hat1,S_hat0){
  
  # model S:
  # outcome <- 'Surv(T_obs,status)'
  # f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  # fitS <- cph(f,data=DATA,y=TRUE,x=TRUE)
  # 
  # DATA.1 <- DATA
  # DATA.1$A <- 1
  # DATA.0 <- DATA
  # DATA.0$A <- 0
  # 
  # fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
  # fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
  # S_hat1 <- fit.pred1$survival
  # S_hat0 <- fit.pred0$survival
  
  E_hat1 <- expected_survival(S_hat1,Y.grid)
  E_hat0 <- expected_survival(S_hat0,Y.grid)
  
  
  
  
  pseudo_T <- CUT_AIPCW(DATA,Y.grid,Q.t.hat,Q.Y.hat,h_C_hat,S_C_hat)
  
  theta_AIPTW_AIPCW <- mean(E_hat1 -E_hat0 + DATA$A *( (pseudo_T - E_hat1 )/ e_hat ) - (1 -  DATA$A) * ( (pseudo_T - E_hat0 )/(1 - e_hat)  ))
  theta_AIPTW_AIPCW
  
  
  return(theta_AIPTW_AIPCW)
}


causal_survival_forest_pred <- function(DATA,X_confouding,e_hat){
  #dataset$status_tau<- as.numeric((dataset$T_obs>=tau) | (dataset$T_obs<tau &  dataset$status == 1 ))
  #dataset$T_obs_tau <- pmin(dataset$T_obs_,tau)
  
  DATA <- as.data.frame(DATA)
  
  categorical_name <- names(which(sapply(subset(DATA, select = c(X_confouding)), class) == "factor"))
  
  numerical_name <- setdiff(X_confouding,categorical_name)
  
  options(na.action='na.pass')    
  
  X_one_hot <- model.matrix(~ 0 + ., subset(DATA, select = c(categorical_name)), na.action = "na.pass")
  
  categorical_name_one_hot <- names(as.data.frame(X_one_hot))
  
  replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
  
  categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
  
  DATA[categorical_name_one_hot] <- X_one_hot
  
  X_confouding_one_hot <- union(categorical_name_one_hot,numerical_name)
  
  
  X <- as.matrix(subset(DATA, select = c(X_confouding_one_hot)))
  W <- as.matrix(DATA[, 'A'])
  Y <- as.matrix(DATA[, 'T_obs_tau'])
  D <- as.matrix(DATA[, 'status_tau'])
  cs.forest <- causal_survival_forest(X, Y, W, D,num.trees=1000)
  cs.pred <- predict(cs.forest, X,W.hat= e_hat,estimate.variance = FALSE)
  return(best_linear_projection(cs.forest)[1])
  #return(mean(cs.pred$predictions))
}




ate_riskregression <- function(DATA,X_outcome,X_treatment,X_censure,Y.grid,method = "IPTW"){
  DATA_input_ATE <- DATA
  DATA_input_ATE$A <- as.factor(DATA_input_ATE$A)
  outcome <- 'Surv(T_obs,status)'
  f <- as.formula(paste(outcome, paste(c(X_outcome,'strat(A)'), collapse = " + "), sep = " ~ "))
  fitS <- cph(f,data=DATA_input_ATE,y=TRUE,x=TRUE,times=Y.grid)
  
  
  X_confouding = intersect(X_outcome,X_treatment)
  outcome <- 'A'
  f <- as.formula(paste(outcome, paste(c(X_confouding), collapse = " + "), sep = " ~ "))
  fitA <- glm(f,data = DATA,family = binomial())
  
  
  outcome <- 'Surv(T_obs,status)'
  f <- as.formula(paste(outcome, paste(c(X_censure,'strat(A)'), collapse = " + "), sep = " ~ "))
  DATA_c <- DATA
  DATA_c$A <- as.factor(DATA_c$A)
  DATA_c$status <- 1- DATA_c$status
  fitC <- cph(f,data=DATA_c,y=TRUE,x=TRUE,times=Y.grid)
  
  
  if (method == "G-formula"){
    ateFit1a <- ate(fitS, data = DATA_input_ATE, treatment = "A",cause=1, times = Y.grid,se = FALSE,verbose= FALSE)
  }
  else{
    ateFit1a <- ate(fitS, data = DATA_input_ATE,censor=fitC, treatment = fitA,cause=1,estimator=method, times = Y.grid,se = FALSE,verbose= FALSE)
  }
  
  
  S_hat1 <- ateFit1a$diffRisk$estimate.A 
  S_hat0 <- ateFit1a$diffRisk$estimate.B
  
  #The mean survival time in ECD recipients followed-up to 10 years
  rmst0 = rmst(times=Y.grid, surv.rates=S_hat0, max.time=tau)
  
  #The mean survival time in SCD recipients followed-up to 10 years
  rmst1 = rmst(times=Y.grid,surv.rates=S_hat1, max.time=tau)
  return(rmst1 - rmst0)
}

# this function give the integral of the survival curve given by S.hat ont the time.grid: Y.grid
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  
  c(base::cbind(1, S.hat) %*% grid.diff)
}


threshol_list <- function(l,threshold){
  
  l_threshold = c()
  for (x in (l)){
    if (is.na(x)){
      return(NA)
    }
    else{
      if (x < - threshold){
        l_threshold <- c(l_threshold,- threshold)
      }
      else{
        if(x > threshold){
          l_threshold <- c(l_threshold,threshold)
        }
        else{
          l_threshold <- c(l_threshold, x)
        }
      }
    }
    
    
  }
  return(l_threshold)
}







