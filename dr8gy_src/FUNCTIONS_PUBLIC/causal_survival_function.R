estimate_outcome_models <- function(DATA,outcome_name,Y.grid,tau,type_of_model="cox",learner='Tlearner',extrapolation_term=TRUE){
  outcome_survival_function <- estimate_survival_function(DATA,outcome_name,Y.grid,tau,type_of_model=type_of_model,learner=learner)
  S_hat <- outcome_survival_function$S_hat
  S_hat1 <- outcome_survival_function$S_hat1
  S_hat0 <- outcome_survival_function$S_hat0
  # h_hat <- estimate_hazard_function(S_hat,Y.grid)
  # if (extrapolation_term){
  #   Q.t.hat <- estimate_remaining_conditinal_outcome(S_hat,Y.grid)
  #   Y.relabeled <- findInterval(DATA$T_obs, Y.grid)
  #   ## E(T_i | X_i,A_i, T_i>y) for y equal to the restricted observed time
  #   Y.relabeled <- findInterval(pmin(DATA$T_obs,tau), Y.grid)
  #   Q.Y.hat <- Q.t.hat[base::cbind(1:nrow(DATA) , Y.relabeled)]
  # }
  # else{
  #   Q.t.hat=NULL
  #   Q.Y.hat=NULL
  # }

  
 # return(list("S_hat"=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0,"h_hat"=h_hat,"Q.t.hat"=Q.t.hat,"Q.Y.hat"=Q.Y.hat))
  return(list("S_hat"=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0))
}

estimate_censoring_models <- function(DATA,censoring_name,Y.grid,tau,type_of_model="cox",learner='Tlearner'){
  DATA$status = 1 -DATA$status
  censoring_survival_function <- estimate_survival_function(DATA,censoring_name,Y.grid,tau,type_of_model=type_of_model,learner=learner)
  S_C_hat <- censoring_survival_function$S_hat
  S_C_hat1 <- censoring_survival_function$S_hat1
  S_C_hat0 <- censoring_survival_function$S_hat0
 # h_C_hat <- estimate_hazard_function(S_C_hat,Y.grid)
  #return(list("S_C_hat"=S_C_hat,"S_C_hat1"=S_C_hat1,"S_C_hat0"=S_C_hat0,"h_C_hat"=h_C_hat))
  return(list("S_C_hat"=S_C_hat,"S_C_hat1"=S_C_hat1,"S_C_hat0"=S_C_hat0))
  
}


estimate_propensity_score <- function(DATA,treatment_name,type_of_model="reglog"){
  # ## GLM
  if (type_of_model == "reglog"){
    outcome <- 'A'
    f <- as.formula(paste(outcome, paste(c(treatment_name), collapse = " + "), sep = " ~ "))
    fitA <- glm(f,data = DATA,family = binomial(link="logit"))
    e_hat <- predict(fitA,newdata=DATA,type="response")
    return(e_hat)
  }
  if (type_of_model == "forest"){
    
    # onehot encode factor variables
    categorical_name <- names(which(sapply(subset(DATA, select = c(treatment_name)), class) == "factor"))
    if (length(categorical_name) >0){
      numerical_name <- setdiff(treatment_name,categorical_name)
      na.action <- options()$na.action
      options(na.action='na.pass')
      X_one_hot <- model.matrix(~ 0 + ., DATA[categorical_name], na.action = "na.pass")
      categorical_name_one_hot <- names(as.data.frame(X_one_hot))
      replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
      categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
      DATA[categorical_name_one_hot] <- X_one_hot
      treatment_name<- union(categorical_name_one_hot,numerical_name)
      options(na.action = na.action)
    }
    ## Forest
    
    Xipw <- as.matrix(DATA[treatment_name])
    Wipw <- as.matrix(DATA$A)
    forest.W <- regression_forest(Xipw, Wipw,honesty = FALSE)
    e_hat <- predict(forest.W)$predictions
    return(e_hat)
  }
}

estimate_survival_kaplan_meier <- function(T_obs,status,weights=NULL,tau=NULL,Y.grid=NULL){
  n <- length(T_obs)
  if (is.null(Y.grid)){
    Y <- pmin(T_obs,tau)
    Y.grid <- sort(unique(Y))
  }
  
  if (is.null(weights)){
    weights <- matrix(rep(1/n,each = n, times = length(Y.grid)), ncol =  length(Y.grid))
  }
  if (length(weights)==n){
    weights <- matrix(rep(weights,times = length(Y.grid)), ncol =  length(Y.grid))
  }
  di <- Yi <- rep(NA, length(Y.grid))
  for (i in 1:length(Y.grid))
  {
    di[i] <- sum((T_obs==Y.grid[i] & status==1)*weights[,i])
    Yi[i] <- sum((T_obs >= Y.grid[i] )*weights[,i]) 
  }
  surv <- cumprod(1 - di/Yi)
  
  return(surv)
}


estimate_hazard_function_bis<- function(S,Y.grid){
  Y.diff <- diff(c(0, Y.grid))
  grid.length <-  length(S)
  log.surv.C <- -log(c(1, S))
  h_hat <- log.surv.C[ 2:(grid.length + 1)] - log.surv.C[ 1:grid.length]
  h_hat <- sweep(h_hat, 2, Y.diff, "/")
  return(h_hat)
}


estimate_hazard_function<- function(S_hat,Y.grid){
  Y.diff <- diff(c(0, Y.grid))
  grid.length <-  ncol(S_hat)
  log.surv.C <- -log(base::cbind(1, S_hat))
  h_hat <- log.surv.C[, 2:(grid.length + 1)] - log.surv.C[, 1:grid.length]
  h_hat <- sweep(h_hat, 2, Y.diff, "/")
  return(h_hat)
}

estimate_remaining_conditinal_outcome <- function(S_hat,Y.grid,outcome_type="RMST",t){
  # model Q.t.hat
  ## E(T_i | X_i,A_i, T_i>y) for y in Y.grid
  if (outcome_type=="RMST"){
    # here we have supposed that tau <- Y.grid[length(Y.grid)]
    Y.grid[1:findInterval(t,Y.grid)]
    S_hat[,1:findInterval(t,Y.grid)]
    Y.diff <- diff(c(0, Y.grid))
    n <- nrow(S_hat)
    grid.length <-  ncol(S_hat)
    Q.t.hat <- matrix(0, n, grid.length)
    dot.products <- sweep(S_hat[, 1:(grid.length - 1)], 2, Y.diff[2:grid.length], "*")
    Q.t.hat[, 1] <- rowSums(dot.products)
    for (i in 2:(grid.length - 1)) {
      Q.t.hat[, i] <- Q.t.hat[, i - 1] - dot.products[, i - 1]
    }
    Q.t.hat <- Q.t.hat / S_hat
    Q.t.hat[is.infinite(Q.t.hat)] <- 0 # The points where S.hat = 0
    Q.t.hat <- sweep(Q.t.hat, 2, Y.grid, "+") # Add back t
    Q.t.hat[, grid.length] <- max(Y.grid)
  }
  if (outcome_type=="survival"){
    # estimate_Q_S <- function(Y.grid,S_hat,t){
    #   n_sample <- length(S_hat[,1])
    #   
    #   T_obs_t <- pmin(DATA$T_obs,t)
    #   
    #   index_denum <- matrix(c(c(1:length(DATA$T_obs_t))  ,findInterval(DATA$T_obs_t,Y.grid)), nrow = n_sample , ncol = 2)
    #   S_C_hat_T_obs_denum <- S_hat[index_denum] 
    #   
    #   index_num <- matrix(c(c(1:length(DATA$T_obs_t))  ,findInterval(pmax(DATA$T_obs_t,t),Y.grid)), nrow = n_sample , ncol = 2)
    #   S_C_hat_T_obs_num <- S_hat[index_num] 
    #   
    #   
    #   return(1-S_C_hat_T_obs_num/S_C_hat_T_obs_denum)  
    #   
    # }
    Q.t.hat <- (1 - (S_hat[,findInterval(t,Y.grid)])/(S_hat))
    
    #Q.t.hat <-(sapply(Y.grid, function(t){ estimate_Q_S(Y.grid,S_hat,t) }))
    
  }
  return(Q.t.hat)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

estimate_survival_function <- function(DATA,X_name,Y.grid,tau,type_of_model="cox",learner='Tlearner'){
  
  if (type_of_model == "cox"){
    if (learner == "Slearner" ){
      outcome <- 'Surv(T_obs,status)'
      f <- as.formula(paste(outcome, paste(c(X_name,'A'), collapse = " + "), sep = " ~ "))
      fitS <- suppressWarnings(coxph(f, data=DATA, x=TRUE))
      fitS$coefficients[is.na(fitS$coefficients)] <- 0
      DATA.1 <- DATA
      DATA.1$A <- 1
      DATA.0 <- DATA
      DATA.0$A <- 0
      
      fit.pred1 <- predictCox(fitS, newdata=DATA.1, times=Y.grid , type = "survival")
      fit.pred0 <- predictCox(fitS, newdata=DATA.0, times=Y.grid , type = "survival")
      S_hat1 <- fit.pred1$survival
      S_hat0 <- fit.pred0$survival
      
      if ( length(Y.grid) > 1){
        S_hat0[,length(Y.grid)] <- S_hat0[,length(Y.grid)-1]
        S_hat1[,length(Y.grid)] <- S_hat1[,length(Y.grid)-1]
      }
      S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
      return(list('S_hat'=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0))
    }
    
    if (learner == "Tlearner" ){
      outcome <- 'Surv(T_obs,status)'
      DATA0 <- DATA[DATA$A == 0,]
      DATA1 <- DATA[DATA$A == 1,]
      f <- as.formula(paste(outcome, paste(c(X_name), collapse = " + "), sep = " ~ "))
      fitS0 <- suppressWarnings(coxph(f, data=DATA0, x=TRUE))
      fitS0$coefficients[is.na(fitS0$coefficients)] <- 0
      
      fitS1 <- suppressWarnings(coxph(f, data=DATA1, x=TRUE))
      fitS1$coefficients[is.na(fitS1$coefficients)] <- 0
      fit.pred1 <- predictCox(fitS1, newdata=DATA, times=Y.grid , type = "survival")
      fit.pred0 <- predictCox(fitS0, newdata=DATA, times=Y.grid , type = "survival")
      S_hat1 <- fit.pred1$survival
      S_hat0 <- fit.pred0$survival
      
      if ( length(Y.grid) > 1){
        S_hat0[,length(Y.grid)] <- S_hat0[,length(Y.grid)-1]
        S_hat1[,length(Y.grid)] <- S_hat1[,length(Y.grid)-1]
      }
      
      S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
      return(list('S_hat'=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0))
    }
    
  }
  
  #Forest:
  if (type_of_model == "forest"){
    
    
    # onehot encode factor variables
    categorical_name <- names(which(sapply(subset(DATA, select = c(X_name)), class) == "factor"))
    if (length(categorical_name) >0){
      numerical_name <- setdiff(X_name,categorical_name)
      na.action <- options()$na.action
      options(na.action='na.pass')
      X_one_hot <- model.matrix(~ 0 + ., DATA[categorical_name], na.action = "na.pass")
      categorical_name_one_hot <- names(as.data.frame(X_one_hot))
      replace_string <- function(string) {return(str_replace_all(string, ' ', '_'))}
      categorical_name_one_hot <- sapply(categorical_name_one_hot,replace_string)
      DATA[categorical_name_one_hot] <- X_one_hot
      X_name<- union(categorical_name_one_hot,numerical_name)
      options(na.action = na.action)
    }
    
    if (learner == "Slearner" ){
      ## Forest S learner:
      X <- data.matrix(DATA[,c(X_name,"A")])
      if (max(X[,"A"])==2){ X[,"A"] <- X[,"A"]-1}
      Y <- as.matrix(DATA[,'T_obs'])
      D <- as.matrix(DATA[,'status'])
      s.forest <- survival_forest(X, Y, D, compute.oob.predictions = FALSE, prediction.type = "Nelson-Aalen")
      
      s.forest[["X.orig"]][, ncol(X) ] <- rep(0, nrow(X))
      S_hat0 <- predict(s.forest, failure.times = Y.grid)$predictions
      #S_hat0 <- S_hat0[,c(1:length(Y.grid))]
      
      s.forest[["X.orig"]][, ncol(X) ] <- rep(1, nrow(X))
      S_hat1 <- predict(s.forest,failure.times = Y.grid)$predictions
      #S_hat1 <- S_hat1[,c(1:length(Y.grid))]
      
      if ( length(Y.grid) > 1){
        S_hat0[,length(Y.grid)] <- S_hat0[,length(Y.grid)-1]
        S_hat1[,length(Y.grid)] <- S_hat1[,length(Y.grid)-1]
      }
    
      S_hat <- S_hat1*X[,"A"] + (1-X[,"A"])*S_hat0
      return(list('S_hat'=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0))
      
    }
    
    if (learner == "Tlearner" ){
      ## Forest T learner:
      DATA_grf_complete0 <-  DATA[DATA$A == 0,]
      X_0 <- as.matrix(subset(DATA_grf_complete0, select = c(X_name)))
      Y_0 <- as.matrix(subset(DATA_grf_complete0, select = c('T_obs')))
      D_0 <- as.matrix(subset(DATA_grf_complete0, select = c('status')))
      
      
      DATA_grf_complete1 <- DATA[DATA$A == 1,]
      X_1 <- as.matrix(subset(DATA_grf_complete1, select = c(X_name)))
      Y_1 <- as.matrix(subset(DATA_grf_complete1, select = c('T_obs')))
      D_1 <- as.matrix(subset(DATA_grf_complete1, select = c('status')))
      
      s.forest0 <- survival_forest(X_0, Y_0, D_0,honesty=FALSE)
      s.forest1 <- survival_forest(X_1, Y_1, D_1,honesty=FALSE)
      
      X <- as.matrix(subset(DATA, select = c(X_name)))
      S_hat0 <- predict(s.forest0, X,failure.times = Y.grid)$predictions
      S_hat1 <- predict(s.forest1, X,failure.times = Y.grid)$predictions
      
      if ( length(Y.grid) > 1){
        S_hat0[,length(Y.grid)] <- S_hat0[,length(Y.grid)-1]
        S_hat1[,length(Y.grid)] <- S_hat1[,length(Y.grid)-1]
      }
      
      S_hat <- S_hat1*DATA$A + (1-DATA$A)*S_hat0
      return(list('S_hat'=S_hat,"S_hat1"=S_hat1,"S_hat0"=S_hat0))
    }
  }
  
}

  
  
censoring_unbias_transformation_withoout_survival <- function(observed_outcome,status_tau,Y.grid,method_censoring="IPCW",S_C_hat,Q.t.hat=NULL,Q.Y.hat=NULL,h_C_hat=NULL){
  ### Attention marche seulement pour RMST
  #status_tau <- DATA$status_tau
  n <- length(observed_outcome)
  index <- matrix(c(c(1:length(observed_outcome))  ,findInterval(observed_outcome,Y.grid)), nrow = n , ncol = 2)
  S_C_hat_T_obs_tau <- S_C_hat[index] # here we take for each individual i, S_c(T_obs_i)
  if (method_censoring=="BJ"){
    pseudo_outcome <-  status_tau*observed_outcome + (1-status_tau)* Q.Y.hat
  }
  
  if (method_censoring=="IPCW"){
    pseudo_outcome <- (observed_outcome * status_tau)/ S_C_hat_T_obs_tau
  }
  
  if (method_censoring=="AIPCW"){
    first_term <- (observed_outcome* status_tau)/ S_C_hat_T_obs_tau
    
    second_term <- (Q.Y.hat * (1-status_tau))/ S_C_hat_T_obs_tau
    
    #third_term <- expected_survival(h_C_hat / S_C_hat *Q.t.hat,Y.grid)
    Y.diff <- diff(c(0, Y.grid))
    integrand <- sweep( ( (h_C_hat) / S_C_hat )* (Q.t.hat), 2, Y.diff, "*")
    third_term <- rep(0, n)
    Y.relabeled <- findInterval(observed_outcome, Y.grid)
    
    for (sample in 1:n) {
      Y.index <- Y.relabeled[sample]
      third_term[sample] <- sum(integrand[sample, 1:Y.index]) 
    }
    pseudo_outcome <- first_term+second_term - third_term
  }
  return(pseudo_outcome)
}
  
observe_risk <- function(times,t){as.numeric(times <= t)}
observe_survival <- function(times,t){as.numeric(times > t)}
observe_RMST <- function(times,t){pmin(times,t)}

observe_outcome <- function(times,t,outcome_type="RMST"){
  if (outcome_type=="RMST"){
    return(observe_RMST(times,t))
  }
  
  if (outcome_type=="survival"){
    return(observe_risk(times,t))
  }
}

estimate_ATE <- function(DATA,tau,Y.grid=NULL,method="AIPTW-AIPCW",method_treatment=NULL,method_censoring=NULL,outcome_type = "RMST",estimate_individual=TRUE,outcome_name,censoring_name,treatment_name,type_of_outcome_model="cox",outcome_learner="Tlearner",type_of_censoring_model ="cox",censoring_learner="Tlearner",type_of_treatment_model="reglog",S_hat=NULL,S_hat1=NULL,S_hat0=NULL,h_hat=NULL,Q.hat = NULL,Q.t.hat=NULL,Q.Y.hat=NULL,S_C_hat=NULL,S_C_hat1=NULL,S_C_hat0=NULL,h_C_hat=NULL,e_hat=NULL,max_time_points=1000,trimming_value = NULL,trimming_value_upper = NULL,trimming_value_lower = NULL,stabilized_IPTW=FALSE){
## outcome_type:
  ### survival 
  ### RMST
## method:
  ### Kaplan_meier
  ### IPTW Kaplan_meier
  ### IPCW Kaplan_meier
  ### IPTW-IPCW Kaplan_meier
  ### g_formula
  ### IPTW-IPCW
  ### AIPTW-IPCW
  ### IPTW-BJ
  ### AIPTW-BJ
  ### IPTW-AIPCW
  ### AIPTW-AIPCW
  ### IPW
## method_treatment:
  ### g_formula
  ### IPTW
  ### AIPTW
## method_censoring:
  ### BJ (for Buckley james)
  ### IPCW
  ### AIPCW
## nuisances function
  ###
  
if (is.factor(DATA$A)){
  DATA$A <- as.numeric(DATA$A)-1
}
  
  
# trimming:
  
if (!is.null(trimming_value)){
  trimming_value_upper <- trimming_value
  trimming_value_lower <- trimming_value
}
  

  
if (!is.null(trimming_value_upper) | !is.null(trimming_value_lower) ) {
  if (is.null(trimming_value_upper)){
    quantile_upper <- 1
  }else{
    quantile_upper <- (100-trimming_value_upper)/100
  }
  if (is.null(trimming_value_lower)){
    quantile_lower <- 0
  }else{
    quantile_lower <- trimming_value_lower/100
  }
  
  e_hat_trimming<-estimate_propensity_score(DATA,treatment_name,type_of_model=type_of_treatment_model)
  quantiles <- quantile(e_hat_trimming,c(quantile_lower,quantile_upper))
  
  e_hat_upper <- quantiles[2]
  e_hat_lower <- quantiles[1]
  
  DATA <- DATA[e_hat_trimming <= e_hat_upper & e_hat_trimming >= e_hat_lower,]
}
  
  
  
  
n <- nrow(DATA)
method_Kaplan_meier<- FALSE

if (!is.null(method)) {
  if (method=="Kaplan_meier"|method=="IPTW Kaplan_meier"|method=="IPCW Kaplan_meier"|method=="IPTW-IPCW Kaplan_meier"){
    method_Kaplan_meier <- TRUE
    extrapolation_term <- FALSE
  }
  else{
    method_Kaplan_meier <- FALSE
  }
  if (method=="g_formula"){
    method_treatment <- "g_formula"
    method_censoring <- NULL
    extrapolation_term <- FALSE
  }
  if (method=="IPTW-IPCW"){
    method_treatment <- "IPTW"
    method_censoring <- "IPCW"
    extrapolation_term <- FALSE
  }
  if (method=="AIPTW-IPCW"){
    method_treatment <- "AIPTW"
    method_censoring <- "IPCW"
    extrapolation_term <- FALSE
  }
  if (method=="IPTW-BJ"){
    method_treatment <- "IPTW"
    method_censoring <- "BJ"
    extrapolation_term <- TRUE
  }
  if (method=="AIPTW-BJ"){
    method_treatment <- "AIPTW"
    method_censoring <- "BJ"
    extrapolation_term <- TRUE
  }
  if (method=="IPTW-AIPCW"){
    method_treatment <- "IPTW"
    method_censoring <- "AIPCW"
    extrapolation_term <- TRUE
  }
  if (method=="AIPTW-AIPCW"){
    method_treatment <- "AIPTW"
    method_censoring <- "AIPCW"
    extrapolation_term <- TRUE
  }
}
n <- nrow(DATA)


if (is.null(Y.grid)){
  Y <- pmin(DATA$T_obs,tau)
  Y.grid <- sort(unique(Y))
  
  if (length(Y.grid) > max_time_points){
    Y.grid <- c(Y.grid[1],sort(sample(x = Y.grid[2:(length(Y.grid)-1)],size = (max_time_points-2),replace = FALSE)),Y.grid[(length(Y.grid))])
  }
}



if (outcome_type == "RMST"){
  times_of_interest <- tau
  #observed_outcome <- pmin(DATA$T_obs,tau)
  observed_outcome <- sapply(times_of_interest,function(t) observe_RMST(DATA$T_obs,t))
  #DATA$status_tau <- as.numeric((DATA$T_obs>=tau) | (DATA$T_obs<tau &  DATA$status == 1 ))
  Delta_tau <- sapply(times_of_interest,function(t) as.numeric((DATA$T_obs>=t) | (DATA$T_obs<t &  DATA$status == 1 ))) # matrix with status tau foor each t of Y.grid (a column is one specific t)
  #Delta_tau <- sapply(tau,function(t) as.numeric((DATA$T_obs>=t) | (DATA$T_obs<t &  DATA$status == 1 )))
  last_observed_time_of_interest <- sapply(times_of_interest,function(t) pmin(DATA$T_obs,t))
  reverse = FALSE
}
if (outcome_type == "survival"){
  times_of_interest <- Y.grid
  # we DATA$T_obs <= t to cmpute the risk but we have to do withdraw 1 to the probability t retrieve the survival probability
  observed_outcome <- sapply(times_of_interest,function(t) observe_risk(DATA$T_obs,t))
  
  Delta_tau <- sapply(times_of_interest,function(t) as.numeric((DATA$T_obs>t) | (DATA$T_obs<=t &  DATA$status == 1 ))) # matrix with status tau foor each t of Y.grid (a column is one specific t)
  last_observed_time_of_interest <- sapply(times_of_interest,function(t) pmin(DATA$T_obs,t))
  reverse = TRUE
}



if (is.null(S_hat) | is.null(S_hat1)|is.null(S_hat0)|is.null(h_hat)){
  outcome_models <- estimate_outcome_models(DATA,outcome_name,Y.grid,tau,type_of_model=type_of_outcome_model,learner=outcome_learner,extrapolation_term =extrapolation_term)
  S_hat <- outcome_models$S_hat
  S_hat1 <- outcome_models$S_hat1
  S_hat0 <- outcome_models$S_hat0
}

 
  
  
if (is.null(S_C_hat) | is.null(S_C_hat1)|is.null(S_C_hat0)|is.null(h_C_hat)){
  censoring_models <- estimate_censoring_models(DATA,censoring_name,Y.grid,tau,type_of_model=type_of_censoring_model,learner=censoring_learner)

  S_C_hat <- censoring_models$S_C_hat
  S_C_hat1 <- censoring_models$S_C_hat1
  S_C_hat0 <- censoring_models$S_C_hat0
  
}

  
if  (is.null(e_hat)){
  e_hat<-estimate_propensity_score(DATA,treatment_name,type_of_model=type_of_treatment_model)
}

if (method_Kaplan_meier){
  if (method=="Kaplan_meier"){
    DATA_treated <- DATA[DATA$A == 1,]
    DATA_not_treated <- DATA[DATA$A == 0,]
    
    S1 <- estimate_survival_kaplan_meier(DATA_treated$T_obs,DATA_treated$status,Y.grid=Y.grid)
    S0 <- estimate_survival_kaplan_meier(DATA_not_treated$T_obs,DATA_not_treated$status,Y.grid=Y.grid)
  }
  if (method=="IPTW Kaplan_meier"){
    DATA_treated <- DATA[DATA$A == 1,]
    DATA_not_treated <- DATA[DATA$A == 0,]
    
    
    weights <- DATA$A * (1/e_hat) + (1-DATA$A)/(1-e_hat)
    weights_treated <- weights[DATA$A == 1]
    weights_not_treated <- weights[DATA$A == 0]
    
    S1 <- estimate_survival_kaplan_meier(DATA_treated$T_obs,DATA_treated$status,Y.grid=Y.grid,weights = weights_treated)
    S0 <- estimate_survival_kaplan_meier(DATA_not_treated$T_obs,DATA_not_treated$status,Y.grid=Y.grid,weights = weights_not_treated)
  }
  if (method=="IPCW Kaplan_meier"){
    DATA_treated <- DATA[DATA$A == 1,]
    DATA_not_treated <- DATA[DATA$A == 0,]
    status_tau_t <- sapply(Y.grid,function(t) as.numeric((DATA$T_obs>=t) | (DATA$T_obs<t &  DATA$status == 1 )))
    
    weights <- status_tau_t/S_C_hat
    weights_treated <- weights[DATA$A == 1,]
    weights_not_treated <- weights[DATA$A == 0,]
    
    S1 <- estimate_survival_kaplan_meier(DATA_treated$T_obs,DATA_treated$status,Y.grid=Y.grid,weights = weights_treated)
    S0 <- estimate_survival_kaplan_meier(DATA_not_treated$T_obs,DATA_not_treated$status,Y.grid=Y.grid,weights = weights_not_treated)
  }
  if (method=="IPTW-IPCW Kaplan_meier"){
    DATA_treated <- DATA[DATA$A == 1,]
    DATA_not_treated <- DATA[DATA$A == 0,]
    status_tau_t <- sapply(Y.grid,function(t) as.numeric((DATA$T_obs>=t) | (DATA$T_obs<t &  DATA$status == 1 )))
    
    weights <- status_tau_t/S_C_hat*DATA$A * (1/e_hat) + (1-DATA$A)/(1-e_hat)
    weights_treated <- weights[DATA$A == 1,]
    weights_not_treated <- weights[DATA$A == 0,]
    
    S1 <- estimate_survival_kaplan_meier(DATA_treated$T_obs,DATA_treated$status,Y.grid=Y.grid,weights = weights_treated)
    S0 <- estimate_survival_kaplan_meier(DATA_not_treated$T_obs,DATA_not_treated$status,Y.grid=Y.grid,weights = weights_not_treated)
  }
  
  if (outcome_type == "RMST"){
    grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
    E_hat1 <- c(c(1, S1) %*% grid.diff)
    E_hat0 <- c(c(1, S0) %*% grid.diff)
    theta_t <- S1-S0
    ATE_treated <- E_hat1
    ATE_not_treated <- E_hat0
    ATE <- E_hat1 - E_hat0
    
    
    
  }
  if (outcome_type == "survival"){
    ATE_treated <- S1
    ATE_not_treated <- S0
    ATE <- S1 - S0
  }

}
  
if (!is.null(method_treatment)){
  if (method_treatment=="g_formula"){
    if (outcome_type == "RMST"){
      
      extrapotaled_outcome_1 <-  expected_survival(S_hat1,Y.grid) #(E_hat1) for treated
      extrapotaled_outcome_1 <-matrix(extrapotaled_outcome_1, nrow = length(extrapotaled_outcome_1) , ncol = 1)
      extrapotaled_outcome_0 <- expected_survival(S_hat0,Y.grid) #(E_hat0) for nt treated
      extrapotaled_outcome_0 <-matrix(extrapotaled_outcome_0, nrow = length(extrapotaled_outcome_0) , ncol = 1)
    }
    if (outcome_type == "survival"){
      extrapotaled_outcome_1 <-  S_hat1 #(E_hat1) for treated
      extrapotaled_outcome_0 <- S_hat0 #(E_hat0) for nt treated
    } 
    ATE_treated_i<- extrapotaled_outcome_1 
    ATE_not_treated_i <- extrapotaled_outcome_0
    ATE_i <- ATE_treated_i - ATE_not_treated_i
  }
  else{
    pseudo_outcome <- sapply(times_of_interest,function(t) censoring_unbias_transformation(DATA,Y.grid,t,method_censoring=method_censoring,outcome_type=outcome_type,S_hat=S_hat,S_C_hat=S_C_hat,Q.t.hat=Q.t.hat,Q.Y.hat=Q.Y.hat,h_C_hat=h_C_hat))
    
    if (!is.null(Q.hat) & outcome_type == "survival"){
      #Q.t.hat <- Q.hat[[findInterval(t,times_of_interest)]]
      pseudo_outcome <- sapply(times_of_interest,function(t) censoring_unbias_transformation(DATA,Y.grid,t,method_censoring=method_censoring,outcome_type=outcome_type,S_hat=S_hat,S_C_hat=S_C_hat,Q.t.hat=Q.hat[[findInterval(t,times_of_interest)]],Q.Y.hat=NULL,h_C_hat=h_C_hat))
      
    }
    if (reverse){
      pseudo_outcome <- 1- pseudo_outcome
    }
  }
  
  if (method_treatment=="IPTW"){
    
    if (stabilized_IPTW){
      ATE_treated_i<- pseudo_outcome* DATA$A/e_hat*(mean(DATA$A/e_hat ))
      ATE_not_treated_i <- pseudo_outcome* ((1-DATA$A)/(1-e_hat))*mean(((1-DATA$A)/(1-e_hat)))
      ATE_i <- ATE_treated_i - ATE_not_treated_i
      
    }else{
      ATE_treated_i<- pseudo_outcome* DATA$A/e_hat 
      
      ATE_not_treated_i <- pseudo_outcome* ((1-DATA$A)/(1-e_hat))
      ATE_i <- ATE_treated_i - ATE_not_treated_i
    }
    #theta_i <- pseudo_outcome* (DATA$A/e_hat - (1-DATA$A)/(1-e_hat))
  }
  
  if (method_treatment=="AIPTW"){
    if (outcome_type == "RMST"){
      extrapotaled_outcome_1 <-  expected_survival(S_hat1,Y.grid) #(E_hat1) for treated
      extrapotaled_outcome_1 <-matrix(extrapotaled_outcome_1, nrow = length(extrapotaled_outcome_1) , ncol = 1)
      extrapotaled_outcome_0 <- expected_survival(S_hat0,Y.grid) #(E_hat0) for nt treated
      extrapotaled_outcome_0 <-matrix(extrapotaled_outcome_0, nrow = length(extrapotaled_outcome_0) , ncol = 1)
    }
    if (outcome_type == "survival"){
      extrapotaled_outcome_1 <-  S_hat1 #(E_hat1) for treated
      extrapotaled_outcome_0 <- S_hat0 #(E_hat0) for not treated
    } 
    
    ATE_treated_i<-extrapotaled_outcome_1 + DATA$A *( (pseudo_outcome - extrapotaled_outcome_1 )/ e_hat ) 
    ATE_not_treated_i <- extrapotaled_outcome_0  + (1 -  DATA$A) * ( (pseudo_outcome - extrapotaled_outcome_0 )/(1 - e_hat)  )
    ATE_i <- ATE_treated_i - ATE_not_treated_i
  }
  
  if (method == "IPW"){
    status <- DATA$status
    ind_status1 <- outer(DATA$T_obs, Y.grid, '>')*status
    ind_status0 <- outer(DATA$T_obs, Y.grid, '>=')*(1-status)
    
    ind<- ind_status1+ind_status0
    A <- DATA$A
    #ind <- outer(DATA$T_obs, Y.grid, '>')*1
    
    S_1_i <- ind*(1/S_C_hat)*(DATA$A/e_hat)
    S_0_i <- ind*(1/S_C_hat)*((1-DATA$A)/(1-e_hat))
    
    if (outcome_type == "RMST"){
      
      extrapotaled_outcome_1 <-  expected_survival(S_1_i,Y.grid) #(E_hat1) for treated
      extrapotaled_outcome_1 <-matrix(extrapotaled_outcome_1, nrow = length(extrapotaled_outcome_1) , ncol = 1)
      extrapotaled_outcome_0 <- expected_survival(S_0_i,Y.grid) #(E_hat0) for nt treated
      extrapotaled_outcome_0 <-matrix(extrapotaled_outcome_0, nrow = length(extrapotaled_outcome_0) , ncol = 1)
    }
    if (outcome_type == "survival"){
      extrapotaled_outcome_1 <-  S_1_i #(E_hat1) for treated
      extrapotaled_outcome_0 <- S_0_i #(E_hat0) for nt treated
    } 
    ATE_treated_i<- extrapotaled_outcome_1 
    ATE_not_treated_i <- extrapotaled_outcome_0
    ATE_i <- ATE_treated_i - ATE_not_treated_i
    print(ATE_treated_i)
  }
  
  if (estimate_individual){
    ATE_treated <- ATE_treated_i
    ATE_not_treated <- ATE_not_treated_i
    ATE <- ATE_i
  }
  else{
    ATE_treated <- colMeans(ATE_treated_i)
    ATE_not_treated <- colMeans(ATE_not_treated_i)
    ATE <- colMeans(ATE_i)
    
  }  
  }
  
if (method == "IPW"){
  status <- DATA$status
  ind_status1 <- outer(DATA$T_obs, Y.grid, '>')*status
  ind_status0 <- outer(DATA$T_obs, Y.grid, '>=')*(1-status)
  
  ind<- ind_status1+ind_status0
  A <- DATA$A
  #ind <- outer(DATA$T_obs, Y.grid, '>')*1
  
  
  if (stabilized_IPTW){
    S_1_i <- ind*(1/S_C_hat)*(DATA$A/e_hat)*(mean(DATA$A/e_hat ))
    S_0_i <- ind*(1/S_C_hat)*((1-DATA$A)/(1-e_hat))*mean(((1-DATA$A)/(1-e_hat)))
    
  }else{
    S_1_i <- ind*(1/S_C_hat)*(DATA$A/e_hat)
    S_0_i <- ind*(1/S_C_hat)*((1-DATA$A)/(1-e_hat))
  }

  
  if (outcome_type == "RMST"){
    
    extrapotaled_outcome_1 <-  expected_survival(S_1_i,Y.grid) #(E_hat1) for treated
    extrapotaled_outcome_1 <-matrix(extrapotaled_outcome_1, nrow = length(extrapotaled_outcome_1) , ncol = 1)
    extrapotaled_outcome_0 <- expected_survival(S_0_i,Y.grid) #(E_hat0) for nt treated
    extrapotaled_outcome_0 <-matrix(extrapotaled_outcome_0, nrow = length(extrapotaled_outcome_0) , ncol = 1)
  }
  if (outcome_type == "survival"){
    extrapotaled_outcome_1 <-  S_1_i #(E_hat1) for treated
    extrapotaled_outcome_0 <- S_0_i #(E_hat0) for nt treated
  } 
  ATE_treated_i<- extrapotaled_outcome_1 
  ATE_not_treated_i <- extrapotaled_outcome_0
  ATE_i <- ATE_treated_i - ATE_not_treated_i
  if (estimate_individual){
    ATE_treated <- ATE_treated_i
    ATE_not_treated <- ATE_not_treated_i
    ATE <- ATE_i
  }
  else{
    ATE_treated <- colMeans(ATE_treated_i)
    ATE_not_treated <- colMeans(ATE_not_treated_i)
    ATE <- colMeans(ATE_i)
    
  }  
}
 
  # if (reverse){
  #   ATE_treated_i <- 1- ATE_treated_i
  #   ATE_not_treated_i <- 1- ATE_not_treated_i
  #   ATE_i <- - ATE_i
  # }
  return(list("ATE_treated"=ATE_treated,"ATE_not_treated"=ATE_not_treated,"ATE"=ATE,times=Y.grid))
  
}



#methods=c("Kaplan_meier","IPTW Kaplan_meier","IPCW Kaplan_meier","g_formula","IPTW-IPCW","AIPTW-IPCW","IPTW-BJ","AIPTW-BJ","IPTW-AIPCW","AIPTW-AIPCW")
estimate_ATE_all <- function(DATA,tau,methods=c("g_formula","IPTW-IPCW","AIPTW-IPCW","IPTW-BJ","AIPTW-BJ","IPTW-AIPCW","AIPTW-AIPCW")
,Y.grid=NULL,method_treatment=NULL,method_censoring=NULL,outcome_type = "RMST",outcome_name,censoring_name,treatment_name,type_of_outcome_model="cox",outcome_learner="Tlearner",type_of_censoring_model ="cox",censoring_learner="Tlearner",type_of_treatment_model="reglog",S_hat=NULL,S_hat1=NULL,S_hat0=NULL,h_hat=NULL,Q.t.hat=NULL,Q.Y.hat=NULL,S_C_hat=NULL,S_C_hat1=NULL,S_C_hat0=NULL,h_C_hat=NULL,e_hat=NULL,max_time_points=1000,trimming_value = NULL,trimming_value_upper = NULL,trimming_value_lower = NULL,stabilized_IPTW=FALSE){
  
  
  if (is.factor(DATA$A)){
    DATA$A <- as.numeric(DATA$A)-1
  }
  
  
  # trimming:
  
  if (!is.null(trimming_value)){
    trimming_value_upper <-  trimming_value
    trimming_value_lower <- trimming_value
  }
  
  
  
  if (!is.null(trimming_value_upper) | !is.null(trimming_value_lower) ) {
    if (is.null(trimming_value_upper)){
      quantile_upper <- 1
    }
    else{
      quantile_upper <- (100-trimming_value_upper)/100
    }
    if (is.null(trimming_value_lower)){
      quantile_lower <- 0
    }
    else{
      quantile_lower <- trimming_value_lower/100
    }
    
    e_hat_trimming<-estimate_propensity_score(DATA,treatment_name,type_of_model=type_of_treatment_model)
    quantiles <- quantile(e_hat_trimming,c(quantile_lower,quantile_upper))
    
    e_hat_upper <- quantiles[2]
    e_hat_lower <- quantiles[1]
    
    DATA <- DATA[e_hat_trimming <= e_hat_upper & e_hat_trimming >= e_hat_lower,]
  }
  n <- nrow(DATA)
  
  if (is.null(Y.grid)){
    Y <- pmin(DATA$T_obs,tau)
    Y.grid <- sort(unique(Y))
    
    if (length(Y.grid) > max_time_points){
      Y.grid <- c(Y.grid[1],sort(sample(x = Y.grid[2:(length(Y.grid)-1)],size = (max_time_points-2),replace = FALSE)),Y.grid[(length(Y.grid))])
    }
  }
  
  if (is.null(S_hat) | is.null(S_hat1)|is.null(S_hat0)){
    outcome_models <- estimate_outcome_models(DATA,outcome_name,Y.grid,tau,type_of_model=type_of_outcome_model,learner=outcome_learner,extrapolation_term =TRUE)
    S_hat <- outcome_models$S_hat
    S_hat1 <- outcome_models$S_hat1
    S_hat0 <- outcome_models$S_hat0
  }
  
  
  
  
  if (is.null(S_C_hat) | is.null(S_C_hat1)|is.null(S_C_hat0)){
    censoring_models <- estimate_censoring_models(DATA,censoring_name,Y.grid,tau,type_of_model=type_of_censoring_model,learner=censoring_learner)
    
    S_C_hat <- censoring_models$S_C_hat
    S_C_hat1 <- censoring_models$S_C_hat1
    S_C_hat0 <- censoring_models$S_C_hat0
  }
  
  if  (is.null(e_hat)){
    e_hat<-estimate_propensity_score(DATA,treatment_name,type_of_model=type_of_treatment_model)
  }
  if (outcome_type == "RMST"){
    if (is.null(Q.t.hat)){
      Q.t.hat <- estimate_remaining_conditinal_outcome(S_hat,Y.grid,outcome_type,tau)
    }
    if (is.null(Q.Y.hat)){
      # attention is only Q.Y.hat is not null 
      T_obs_t <- pmin(DATA$T_obs,tau)
      Q.Y.hat <-  apply_time_to_matrix(Q.t.hat,Y.grid,T_obs_t)
    }
    Q.hat <- NULL
  }
  if (outcome_type == "survival"){
    Q.t.hat <- NULL
    Q.Y.hat <<- NULL
    Q.hat <-  lapply(Y.grid,function(t) estimate_remaining_conditinal_outcome(S_hat,Y.grid,outcome_type="survival",t))
  }

  list_ate = list()
  for (method in methods) {
    result_ATE <- estimate_ATE(DATA=DATA,tau=tau,Y.grid=Y.grid,estimate_individual=FALSE,method=method,method_treatment=method_treatment,method_censoring=method_censoring,outcome_type = outcome_type,outcome_name=outcome_name,censoring_name=censoring_name,treatment_name=treatment_name,type_of_outcome_model=type_of_outcome_model,outcome_learner=outcome_learner,type_of_censoring_model =type_of_censoring_model,censoring_learner=censoring_learner,type_of_treatment_model=type_of_treatment_model,S_hat=S_hat,S_hat1=S_hat1,S_hat0=S_hat0,h_hat=h_hat,Q.hat=Q.hat,Q.t.hat=Q.t.hat,Q.Y.hat=Q.Y.hat,S_C_hat=S_C_hat,S_C_hat1=S_C_hat1,S_C_hat0=S_C_hat0,h_C_hat=h_C_hat,e_hat=e_hat,stabilized_IPTW=stabilized_IPTW)$ATE
    if (outcome_type == "RMST"){
      theta <- result_ATE
    }
    if (outcome_type == "survival"){
      grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
      theta <- mean(c(c(1, result_ATE) %*% grid.diff))
    }
    list_ate[[method]] = theta
  
}
  return(list_ate)
}


simulate_DATA_bis <- function(n,R,mu,tau,parsS,parsA,parsC,parsSampling,pars_interaction_A_X,misoutcome ='none',miscensoring ='none',mistreatment ='none',misRCT='none',misrule='none'){
  # Covariable
  #R <- matrix(c(1,0,0,0,
  #              0,1,0,0,
  #              0,0,1,0,
  #              0,0,0,1), 
  #            nrow = 4, ncol = 4)
  #mu <- c(X_1 = 2,  X_2 = 1,X_3 = 1,X_4 = -1)
  # X <- mvrnorm(n, mu = mu, Sigma = R)
  # X<- as.data.frame(X)
  
  X_1 <-  runif(n,-1,1)
  X_2 <- runif(n,-1,1)
  X_3 <- runif(n,-1,1)
  X_4 <- runif(n,-1,1)
  X<- cbind(X_1,X_2,X_3,X_4)
  X<- as.data.frame(X)
  names(X) = c("X_1" ,  "X_2" ,"X_3","X_4")
  
  
  
  #Treatment
  #uniform 
  #A <- rbinom(n, 1, 0.5)
  
  # logistic regression model
  expit <- function(x){ return(1/(1+exp(-x))) }
  #parsA <- c(X_1 = 1,  X_2 = 0.5,X_3 = 0,X_4 = 0.3)
  
  X_treatment <- transform_covariate(X,mistreatment)
  prop_scores<- rowSums(as.matrix(X_treatment) %*% diag(parsA))
  prop_scores <- expit(prop_scores)
  #prop_scores <- 0.01 + (pmin(0.98,prop_scores) - min(prop_scores))/(max(prop_scores)-min(prop_scores))
  A <- sapply(prop_scores, FUN=function(p) rbinom(n=1, size=1, prob=p))
  
  
  
  
  
  Xinteraction<- X

  
  X_outcome <- as.matrix(transform_covariate(X,misoutcome))
  Xinteraction <- as.matrix(transform_covariate(Xinteraction,misrule))
  
  
  epsilon <- runif(n, min = 0.00000001, max = 1)
  
  T1 <- -log(epsilon)/(0.1*exp(  X_outcome%*% parsS + cbind(1,Xinteraction)%*% pars_interaction_A_X))
  T0 <- -log(epsilon)/(0.1*exp(  X_outcome%*% parsS ))
  
  
  # epsilon <- log(rexp(n,0.01))
  # T1 <- log(1+exp( - X_outcome%*% parsS - cbind(1,Xinteraction)%*% pars_interaction_A_X +epsilon)) 
  # T0<- log(1+exp(-X_outcome%*% parsS+epsilon)) 
  
  T_true <- A*T1 + (1-A)*T0
  
  #####Censoring
  epsilon2 <- runif(n, min = 0.00000001, max = 1)
  X_censoring<- as.matrix(transform_covariate(X,miscensoring))
  C <- -log(epsilon2)/(0.1*exp(rowSums(X_censoring %*% diag(parsC))))
  
  
  # epsilon2 <- log(rexp(n,0.01))
  # X_censoring<- as.matrix(transform_covariate(X,miscensoring))
  # C <- log(1+exp(- X_censoring%*% parsC +epsilon2)) 
  
  #C <- runif(n,0,10)
  
  
  # Observed data
  T_obs <- pmin(T_true,C)
  status <- (T_true < C)  
  
  
  
  
  # Restricted survival time
  T_obs_tau <- pmin(T_obs,tau)
  status_tau <- as.numeric((T_obs>tau) | (T_obs<=tau &  status == 1 ))
  RMST_threshold_hat_true <- mean(pmin(T1,tau) - pmin(T0,tau))
  #RMST_threshold_hat_true <- mean(pmin(T + delta,tau) - pmin(T,tau))
  
  
  # Dataset
  #T_obs <- T_obs_threshold
  #status <- status_threshold
  DATA_target_population <- data.frame(T_obs,T_obs_tau, status,status_tau,X,A)
  DATA_target_population$status <- as.numeric(DATA_target_population$status)
  
  
  # logistic regression model for sampling
  X_RCT<- as.matrix(transform_covariate(X,misRCT))
  prop_scores<- rowSums(as.matrix(X_RCT) %*% diag(parsSampling))
  prop_scores <- expit(prop_scores)
  #prop_scores <- 0.01 + (pmin(0.98,prop_scores) - min(prop_scores))/(max(prop_scores)-min(prop_scores))
  sigma <- sapply(prop_scores, FUN=function(p) rbinom(n=1, size=1, prob=p))
  DATA_RCT <-DATA_target_population[sigma==1,]
  DATA_target_population$sigma <- sigma
  
  l_DATA <- list("DATA_RCT" = DATA_RCT,"DATA_target_population"=DATA_target_population,"sigma"=sigma, "RMST_threshold_hat_true" = RMST_threshold_hat_true,"T1"= T1,"T0"=T0,"C"=C,"T_true"=T_true)
  return(l_DATA)
}

## misspecification_functin
### X 4 columns
transform_covariate <- function(X,mispecification){
  
  if (mispecification == 'quadratic'){
    X_transformed <- X
    X_transformed$X_1 <- X$X_1**2
    X_transformed$X_2 <- X$X_2**2
    X_transformed$X_3 <- X$X_3**2
    X_transformed$X_4 <- X$X_4**2
    return(X_transformed)
  }
  
  
  if (mispecification == 'interaction'){
    X_transformed <- X
    X_transformed$X_1 <- X$X_1*X$X_2
    X_transformed$X_2 <- X$X_2*X$X_3
    X_transformed$X_3 <- X$X_3*X$X_4
    X_transformed$X_4 <- X$X_4*X$X_1
    return(X_transformed)
  }
  
  if (mispecification == 'cosinus'){
    X_transformed <- X
    X_transformed$X_1 <- cos(X$X_1)
    X_transformed$X_2 <- cos(X$X_2)
    X_transformed$X_3 <- cos(X$X_3)
    X_transformed$X_4 <- cos(X$X_4)
    return(X_transformed)
  }
  
  if (mispecification == 'none'){
    return(X)
  }
  
}


expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  
  c(base::cbind(1, S.hat) %*% grid.diff)
}





apply_time_to_matrix <- function(S_C_hat,Y.grid,times){index <- matrix(c(c(1:length(times))  ,findInterval(times,Y.grid)), nrow = length(times) , ncol = 2)
return(S_C_hat[index])}

integrate <- function(integrand,Y.grid,times){
  filter <- sapply(c(1:length(Y.grid)),function(i){return(as.numeric(i <=findInterval(times,Y.grid)))})
  integrand_filtered <- filter * integrand
  integrated_value <- rowSums(integrand_filtered)
  return(integrated_value)}

censoring_unbias_transformation <- function(DATA,Y.grid,t,method_censoring="IPCW",outcome_type="RMST",S_hat,S_C_hat,Q.t.hat=NULL,Q.Y.hat=NULL,h_C_hat=NULL,method_aipw=1){
  ### Attention marche seulement pour RMST
  #status_tau <- DATA$status_tau
  Delta_t <-  as.numeric((DATA$T_obs>t) | (DATA$T_obs<=t &  DATA$status == 1 )) # matrix with status tau foor each t of Y.grid (a column is one specific t)
  T_obs_t <- pmin(DATA$T_obs,t)
  observed_outcome_t <- observe_outcome(DATA$T_obs,t,outcome_type=outcome_type)
  # compute intermediate results depending on the method used
  if (method_censoring=="IPCW" | method_censoring=="AIPCW"){
    S_C_hat_T_obs_tau <- apply_time_to_matrix(S_C_hat,Y.grid,T_obs_t)
  }
  if (method_censoring=="BJ" | method_censoring=="AIPCW"){
    if (is.null(Q.t.hat)){
      Q.t.hat <- estimate_remaining_conditinal_outcome(S_hat,Y.grid,outcome_type,t)
    }
    if (is.null(Q.Y.hat)){
      # attention is only Q.Y.hat is not null 
      Q.Y.hat <-  apply_time_to_matrix(Q.t.hat,Y.grid,T_obs_t)
    }
  }
  
  if ( method_censoring=="AIPCW"){
    if (is.null(h_C_hat)){
      h_C_hat <- estimate_hazard_function(S_C_hat,Y.grid)
    }
  }
  
  ### run method
  if (method_censoring=="IPCW"){
    pseudo_outcome <- (observed_outcome_t * Delta_t)/ S_C_hat_T_obs_tau
  }
  if (method_censoring=="BJ"){
    pseudo_outcome <-  observed_outcome_t*Delta_t + (1-Delta_t)* Q.Y.hat
  }
  
  if (method_censoring=="AIPCW"){
    
    #first method
    if (method_aipw == 1){
      first_term <- (observed_outcome_t* Delta_t)/ S_C_hat_T_obs_tau
      second_term <- (Q.Y.hat * (1-Delta_t))/ S_C_hat_T_obs_tau
      #third_term <- expected_survival(h_C_hat / S_C_hat *Q.t.hat,Y.grid)
      Y.diff <- diff(c(0, Y.grid))
      integrand <- sweep( ( (h_C_hat) / S_C_hat )* (Q.t.hat), 2, Y.diff, "*")
      third_term <-  integrate(integrand,Y.grid,T_obs_t)
      pseudo_outcome <- first_term+second_term - third_term
    }
    
    #second methood
    if (method_aipw == 2){
      first_term <- (observed_outcome_t* Delta_t)/ S_C_hat_T_obs_tau
      num.samples <- length(observed_outcome_t)
      N_C <- sapply(Y.grid,function(t){return(as.numeric(DATA$T_obs <= t & DATA$status==0 ))})
      M.diff <-  t(diff(rbind(0, t(N_C+log( S_C_hat)))))
      integrand <- Q.t.hat/ S_C_hat* (M.diff)
      augmentation_term <- integrate(integrand,Y.grid,T_obs_t)
      pseudo_outcome <- first_term+augmentation_term
    }
   
    
  }
  
  
  return(pseudo_outcome)
}



censoring_augmentation_term <- function(DATA,Y.grid,t,method_censoring="IPCW",outcome_type="RMST",S_hat,S_C_hat,Q.t.hat=NULL,Q.Y.hat=NULL,h_C_hat=NULL,method_aipw=2){
  ### Attention marche seulement pour RMST
  #status_tau <- DATA$status_tau
  Delta_t <-  as.numeric((DATA$T_obs>t) | (DATA$T_obs<=t &  DATA$status == 1 )) # matrix with status tau foor each t of Y.grid (a column is one specific t)
  T_obs_t <- pmin(DATA$T_obs,t)
  observed_outcome_t <- observe_outcome(DATA$T_obs,t,outcome_type=outcome_type)
  
  
  # compute intermediate results depending on the method used
  if (method_censoring=="IPCW" | method_censoring=="AIPCW"){
    S_C_hat_T_obs_tau <- apply_time_to_matrix(S_C_hat,Y.grid,T_obs_t)
  }
  
  if (method_censoring=="BJ" | method_censoring=="AIPCW"){
    if (is.null(Q.t.hat)){
      Q.t.hat <- estimate_remaining_conditinal_outcome(S_hat,Y.grid,outcome_type,t)
    }
    if (is.null(Q.Y.hat)){
      # attention is only Q.Y.hat is not null 
      Q.Y.hat <-  apply_time_to_matrix(Q.t.hat,Y.grid,T_obs_t)
    }
  }
  if ( method_censoring=="AIPCW"){
    if (is.null(h_C_hat)){
      h_C_hat <- estimate_hazard_function(S_C_hat,Y.grid)
    }
  }
  
  ### run method
  if (method_censoring=="IPCW"){
    pseudo_outcome <- (observed_outcome_t * Delta_t)/ S_C_hat_T_obs_tau
  }
  if (method_censoring=="BJ"){
    pseudo_outcome <-  observed_outcome_t*Delta_t + (1-Delta_t)* Q.Y.hat
  }
  
  if (method_censoring=="AIPCW"){
    
    #first method
    if (method_aipw == 1){
      first_term <- (observed_outcome_t* Delta_t)/ S_C_hat_T_obs_tau
      second_term <- (Q.Y.hat * (1-Delta_t))/ S_C_hat_T_obs_tau
      #third_term <- expected_survival(h_C_hat / S_C_hat *Q.t.hat,Y.grid)
      Y.diff <- diff(c(0, Y.grid))
      integrand <- sweep( ( (h_C_hat) / S_C_hat )* (Q.t.hat), 2, Y.diff, "*")
      third_term <-  integrate(integrand,Y.grid,T_obs_t)
      augmentation_term <- second_term - third_term
      pseudo_outcome <- first_term+augmentation_term
    }
    
    #second methood
    if (method_aipw == 2){
      first_term <- (observed_outcome_t* Delta_t)/ S_C_hat_T_obs_tau
      num.samples <- length(observed_outcome_t)
      N_C <- sapply(Y.grid,function(t){return(as.numeric(DATA$T_obs <= t & DATA$status==0 ))})
      M.diff <-  t(diff(rbind(0, t(N_C+log( S_C_hat)))))
      integrand <- Q.t.hat/ S_C_hat* (M.diff)
      augmentation_term <- integrate(integrand,Y.grid,T_obs_t)
      pseudo_outcome <- first_term+augmentation_term
    }
    
    
  }
  
  
  return(augmentation_term)
}

estimate_average_hazard_ratio <- function(S1,S2,gamma=0.5){
  return(mean(diff(c(1,S1**gamma))*S1**gamma)/mean(diff(c(1,S2**gamma))*S2**gamma))
}

