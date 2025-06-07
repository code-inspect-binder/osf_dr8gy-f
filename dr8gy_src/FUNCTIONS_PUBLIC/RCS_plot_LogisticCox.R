# put.pvalue ####
# Function to put the p-values accoridng to the journals guidelines 
put.pvalue=function(p){
  sapply(p, function(x){
    if (is.na(x)){pc="NA"} else
      if (x<0.0001){pc="<0.0001"} else
        if (0.0001<=x & x<0.001){pc <- sprintf("%.4f", x)} else
          if (0.049<=x & x<0.05) {pc <- sprintf("%.3f", floor(x*1000)/1000)} else
            if (0.001<=x & x<0.1){pc <- sprintf("%.3f", x)} else 
              if (0.1<=x){pc <- sprintf("%.2f", x)} 
    return(pc)
  })
}

# wald.test ####
# Function to perform Wald test in logistic or Cox model
wald.test <- function(model, ntest){
  d <- length(ntest)
  a <- t(coefficients(model)[ntest])
  mat <- vcov(model)[ntest,ntest]
  value <- as.numeric(a%*%solve(mat)%*%t(a))
  pvalue <- pchisq(value, df=d, lower.tail=F)
  return(list(value=value, df=d, pvalue=pvalue))
}

# knots.harrell ####
# Function to calculate the knots used for the restricted cubic splines according to the Harrell's recommandations
knots.harrell <- function(var,data,nk=3){
  var <- data[[var]]
  if (nk<3 | nk>7) {stop(paste0("'nk' = ", nk, ". Please use an integer between 3 and 7."))}
  if (nk==3){knots <- quantile(var,probs=c(0.1,0.5,0.9),type=2,na.rm=T)}
  if (nk==4){knots <- quantile(var,probs=c(0.05,0.35,0.65,0.95),type=2,na.rm=T)}
  if (nk==5){knots <- quantile(var,probs=c(0.05,0.275,0.5,0.725,0.95),type=2,na.rm=T)}
  if (nk==6){knots <- quantile(var,probs=c(0.05,0.23,0.41,0.59,0.77,0.95),type=2,na.rm=T)}
  if (nk==7){knots <- quantile(var,probs=c(0.025,0.1833,0.3417,0.5,0.6583,0.8167,0.975),type=2,na.rm=T)}
  return(knots)
}

# is.binary ####
# Function to check if a vector is a binary vector (only two possible values 0 and 1)
is.binary <- function(x){
  if (!is.numeric(x)){FALSE} else {if_else(length(unique(x))==2 & min(x)==0 & max(x)==1, TRUE, FALSE)}
}

# is.color ####
# Function to check if each element of a vector appears to be a color
is.color <- function(x){
  sapply(x, function(y){if_else(!is.character(y), FALSE,
                                (y %in%colors()) | ((nchar(y) %in% c(7, 9)) & (substr(y, 1, 1) == "#")))})
}

# rcs_effect ####
rcs_effect <- function(t2event=NULL, event, model = c("cox", "logistic"),
                       var, ref_value = NULL, min = NULL, max = NULL, length = 100, 
                       use_rcs=T, knots=NULL, nk=3, adj=NULL, data){
  # Load packages
  require(tidyverse)
  require(survival)
  require(Hmisc)
  
  # Check if variables t2event, event, factor,  cont (and variables in vector adj if it is not null) exist in the dataset.
  if (model == "cox"){
    if (!t2event%in%names(data)){stop(paste0("Variable ", t2event, " doesn't exist."))}
    if (!is.numeric(data[[t2event]])){stop(paste0("Variable ", t2event, " is not a numeric variable."))}
  } else {t2event <- NULL}
  
  if (!event%in%names(data)){stop(paste0("Variable ", event, " doesn't exist."))}
  if (!is.binary(data[[event]])){stop(paste0("Variable ", event, " is not a binary variable. Convert it into a binary."))}
  
  if (!var%in%names(data)){stop(paste0("Variable ", var, " doesn't exist."))}
  if (!is.numeric(data[[var]])){ stop(paste0("Variable '", var, "' is not a numeric variable."))}
  
  if (!is.null(adj)){ 
    adj <- setdiff(adj, var)
    if (sum(!adj%in%names(data))>0) {
      stop(paste0("The following variables don't exist: ", 
                  paste0("'",adj[which(!adj%in%names(data))],"'", collapse=", "), "."))
    }
  } 
  # Remove observations with mising values for variables required for the analysis
  if (model == "cox"){data <- data %>% dplyr::select(t2event, event, var, adj) %>% drop_na()}
  if (model == "logistic"){data <- data %>% dplyr::select(event, var, adj) %>% drop_na()}
  
  # Define min and max value for the continuous variable
  if (is.null(min)){
    min <- quantile(data[[var]], type=2, probs=0.025, na.rm=T)
  } else {
    if (!is.numeric(min) | length(min)>1){
      stop("'min' should be a numeric vector of length 1.")
    }
  }
  if (is.null(max)){
    max <- quantile(data[[var]], type=2, probs=0.975, na.rm=T)
  } else {
    if (!is.numeric(max) | length(max)>1){
      stop("'max' should be a numeric vector of length 1.")
    }
  }
  
  # Define the reference value for the continuous variable
  if (is.null(ref_value)){
    ref_value = median(data[[var]], na.rm=T)
  } else {
    if (!is.numeric(ref_value) | length(ref_value)>1){
      stop("'ref_value' should be a numeric vector of length 1.")
    }
  }
  
  # Define the knots for RC splines if use_rcs == T
  if (use_rcs==T){
    if(is.null(knots)){ knots <- knots.harrell(var, data, nk=nk) }
    mat <- model.matrix(~ rcspline.eval(x, knots=knots, inclx=T), data.frame(x=data[[var]]))[,-1]
    X <- model.matrix(~ rcspline.eval(x, knots=knots, inclx=T), data.frame(x = seq(min, max, length=length)))[,-1]
    Xref <- model.matrix(~ rcspline.eval(x, knots=knots, inclx=T), data.frame(x=ref_value))[,-1]
  } else {
    knots <- NULL
    mat <- matrix(data[[var]], ncol=1)
    X <- matrix(seq(min, max, length=length), ncol=1)
    Xref <- ref_value
  }
  
  # Initialisation of the results list
  results <- list(model = model, t2event = t2event, event = event, 
                  var = var, ref_value = ref_value, method = if_else(use_rcs==T, "RC splines", "linear"),
                  knots = knots, adj = adj, effect = NULL, global_test = NULL, linearity_test = NULL)
  
  if (model == "cox"){
    if (use_rcs==T){
      fml <- formula(paste0("Surv(", t2event, "," , event, ")~", paste0(c(paste0("I(rcspline.eval(",var,", knots = knots, inclx=T))"), adj), collapse=" + ")))
      index <- seq(1, length(knots)-1,1)
    } else {
      fml <- formula(paste0("Surv(", t2event, "," , event, ")~", paste0(c(var,adj), collapse=" + ")))
      index <- 1
    }
    fit <- coxph(fml, data = data, ties = "breslow")
    mat_cov <- vcov(fit)[index,index]
    coeff <- fit$coefficients[index]
    var <- NULL
    for (j in 1:nrow(X)){
      var <- c(var, t(X[j,]-Xref)%*%mat_cov%*%t(t(X[j,]-Xref)))
    }
    pred <- as.vector(X%*%coeff)
    predref <- as.vector(Xref%*%coeff)
    exp_xbeta <- exp(pred-predref)
    results$effect <- data.frame(value = X[,1],
                                 hr = exp_xbeta, 
                                 lower_ci = exp(log(exp_xbeta)-qnorm(0.975)*sqrt(var)),
                                 upper_ci = exp(log(exp_xbeta)+qnorm(0.975)*sqrt(var)))
  }
  if (model == "logistic"){
    if (use_rcs==T){
      fml <- formula(paste0(event, "~", paste0(c(paste0("I(rcspline.eval(",var,", knots = knots, inclx=T))"), adj), collapse=" + ")))
      index <- seq(2, length(knots),1)
    } else {
      fml <- formula(paste0(event, "~", paste0(c(var,adj), collapse=" + ")))
      index <- 2
    }
    fit <- glm(fml, data = data, family = "binomial")
    mat_cov <- vcov(fit)[index,index]
    coeff <- fit$coefficients[index]
    var <- NULL
    for (j in 1:nrow(X)){
      var <- c(var, t(X[j,]-Xref)%*%mat_cov%*%t(t(X[j,]-Xref)))
    }
    pred <- as.vector(X%*%coeff)
    predref <- as.vector(Xref%*%coeff)
    exp_xbeta <- exp(pred-predref)
    results$effect <- data.frame(value = X[,1],
                                 or = exp_xbeta, 
                                 lower_ci = exp(log(exp_xbeta)-qnorm(0.975)*sqrt(var)),
                                 upper_ci = exp(log(exp_xbeta)+qnorm(0.975)*sqrt(var)))
  }
  # Global effect
  
  # Global effect + Test of linearity
  global_test <- wald.test(fit, index)
  results$global_test <- list(value = global_test$value, df = global_test$df, pvalue = global_test$pvalue)
  if (use_rcs==T){
    linearity_test <- wald.test(fit, index[-1])
    results$linearity_test <- list(value = linearity_test$value, df = linearity_test$df, pvalue = linearity_test$pvalue)
  }
  return(results)
}

# plot_rcs_effect ####
plot_rcs_effect <- function(results, ylim=NULL, xlim=NULL, label_var = NULL,
                            main = NULL, xlab=NULL, ylab=NULL, ylog=F,
                            lwd=2, cex_lab=1.2, cex_axis=1.2, polygon_ci=T,
                            col1, col2, xaxis=NULL, yaxis=NULL, 
                            hline=T, col_hline = "black", lwd_hline=1, 
                            add_pglobal=F, name_pglobal=NULL, xpglobal = NULL, ypglobal=NULL, adj_pglobal=0, cex_pglobal=1,
                            add_plin=F, name_plin=NULL, xplin = NULL, yplin=NULL, adj_plin=0, cex_plin=1){
  # results: object from 'rcs_effect' function
  # ylim: the y limits of the plot.
  # xlim: the x limits of the plot.
  # label_var: the label used for the continuous variable.
  # main: the title for the plot.
  # xlab: a title for the x axis.
  # ylab: a title for the y axis.
  # ylog: logical value to represent hazard ration on the log-scale. By default is FALSE.
  # lwd: the line width for plot hazard ratio and confidence interval.
  # cex_lab: a numerical value giving the amount by which plotting text and symbols should be magnified relative to the default
  # cex_axis: the magnification to be used for axis annotation relative to the current setting of cex.
  # polygon_ci: logical value for draws a polygon for the confidence interval.
  # col1: vector of colors for draw the hazard ratio.
  # col2: vector of colors for draw the confidence interval of the odds/hazard ratio.
  # xaxis: 
  # yaxis:
  # hline: logical value if you want to dra the reference line OR/HR=1.
  # col_hline: colour used for the reference line OR/HR=1
  # lwd_hline: lwd for the reference line OR/HR=1
  # add_pglobal: logical value to add the global test p-value.
  # name_pglobal: text in front of the global test p-value. 
  # xpglobal: x position for the global test p-value.
  # ypglobal: y position for the global test p-value.
  # adj_pglobal: adj for the global test p-value.
  # cex_pglobal: cex for the global test p-value.
  # add_plin: logical value to add the linearity test p-value.
  # name_plin: text in front of the linearity test p-value. 
  # xplin: x position for the linearity test p-value.
  # yplin: y position for the linearity test p-value.
  # adj_plin: adj for the linearity test p-value.
  # cex_plin: cex for the linearity test p-value.
  
  # Load packages
  require(scales) # used for the x-axis and y-axis ticks 
  
  # ylim
  if (is.null(ylim)) {
    ylim <- range(c(results$effect$lower_ci,results$effect$upper_ci))
  } else {
    if (!length(ylim)==2 | !is.numeric(ylim)){
      stop(paste0("'ylim' should be a numeric vector of length 2."))}
  }
  
  # xlim
  if (is.null(xlim)){xlim <- range(results$effect$value)}
  
  # Label for the continuous variable
  if (is.null(label_var)){
    label_var <- results$var
  } else if (!is.character(label_var) | length(label_var)>1){
    stop("'label_var' should be a string character.")
  }
  
  # Main
  if (is.null(main)){
    if (results$method == "RC splines"){main <- paste0("Effect of ", label_var, " (RC splines with ", length(results$knots), " knots)")}
    if (results$method == "linear"){main <- paste0("Effect of ", label_var, " (linear)")}
  }
  # xlab
  if (is.null(xlab)){
    xlab <- label_var
  } else {
    if (!is.character(xlab) | length(xlab)>1){
      stop("'xlab' should be a character vector of length 1.")
    }
  }
  # ylab
  if (is.null(ylab)){
    if (results$model == "cox") {ylab <- "Hazard Ratio"}
    if (results$model == "logistic") {ylab <- "Odds Ratio"} 
  } else {
    if (!is.character(ylab) | length(ylab)>1){
      stop("'xlab' should be a character vector of length 1.")
    }
  }
  
  if (ylog==F){
    plot(results$effect$value, results$effect[[2]], ylim = ylim, xlim=xlim,
         col = "white", type = "p", main = main, xlab = xlab, ylab = ylab,
         lwd = lwd, cex.lab = cex_lab, cex.axis=cex_axis, frame.plot = F, xaxt = "n", yaxt = "n",
         panel.first = grid())
  } else {
    plot(results$effect$value, results$effect[[2]], ylim = ylim, xlim=xlim,
         col = "white", type = "p", main = main, xlab = xlab, ylab = ylab,
         lwd = lwd, cex.lab = cex_lab, cex.axis=cex_axis, frame.plot = F, xaxt = "n", yaxt = "n", log="y",
         panel.first = grid())
  }
  # Check if col1 and col2 are color vectors of same length of names
  if (is.null(col1) | !length(col1)==length(names)){
    stop(paste0("'col1' must be a color vector of length ", length(names),"."))
  } else {
    check_col1 <- is.color(col1)
    if (sum(!check_col1)>0){
      stop(paste0("The following colors in the 'col1' vector does'nt exist: ",
                  paste0("'",names(check_col1)[which(!check_col1)],"'", collapse=", "),"."))  
    }
  }
  if (is.null(col2) | !length(col2)==length(names)){
    stop(paste0("'col2' must be a color vector of length ", length(names),"."))
  } else {
    check_col2 <- is.color(col2)
    if (sum(!check_col2)>0){
      stop(paste0("The following colors in the 'col2' vector does'nt exist: ",
                  paste0("'",names(check_col2)[which(!check_col2)],"'", collapse=", "),"."))  
    }
  }
  
  if (hline==T){abline(h = 1, col=col_hline , lwd=lwd_hline)}
  
  if (polygon_ci==T){ 
    polygon(x=c(results$effect$value, rev(results$effect$value)),
            y=c(results$effect$lower_ci, rev(results$effect$upper_ci)),
            col=col2,border=NA)
  }
  lines(results$effect$value, results$effect[[2]], type="l", col=col1, lwd=2)
  lines(results$effect$value, results$effect$lower_ci, type="l", col=col1, lwd=2, lty=2)
  lines(results$effect$value, results$effect$upper_ci, type="l", col=col1, lwd=2, lty=2)
  
  # Axis
  if (is.null(xaxis)){
    xaxis <- pretty_breaks(n=5)(xlim)
    xaxis <- c(head(xaxis, 1) - diff(xaxis[1:2]), xaxis, tail(xaxis, 1) + diff(xaxis[1:2]))
  }
  if (is.null(yaxis)){
    yaxis <- pretty_breaks(n=5)(ylim)
    yaxis <- c(head(yaxis, 1) - diff(yaxis[1:2]), yaxis, tail(yaxis, 1) + diff(yaxis[1:2]))
  }
  axis(1, at=xaxis, lwd=2, cex.axis=1.2)
  axis(2, at=yaxis, lwd=2, cex.axis=1.2)
  
  # Global test p-value
  if (add_pglobal){
    if (is.null(name_pglobal)){
      if (results$method == "RC splines"){name_pglobal <- paste0("Overall effect (RCS with ",length(results$knots), " knots): ")}
      if (results$method == "linear"){name_pglobal <- "Overall effect (linear): "}
    }
    if (is.null(xpglobal)){xpglobal <- min(xlim)}
    if (is.null(ypglobal)){ypglobal <- min(ylim)+diff(range(ylim))/20}
    text(x=xpglobal, y=ypglobal, paste0(name_pglobal,
                                  if_else(results$global_test$pvalue<0.0001, "p<0.0001", 
                                          paste0("p=", put.pvalue(results$global_test$pvalue)))), 
         adj=adj_pglobal, cex = cex_pglobal)
  }
  
  # Linearity test p-value
  if (add_plin & !is.null(results$linearity_test)){
    if (is.null(name_plin)){name_plin <- paste0("Linearity test (RCS with ",length(results$knots), " knots): ")}
    if (is.null(xplin)){xplin <- min(xlim)}
    if (is.null(yplin)){yplin <- min(ylim)}
    text(x=xplin, y=yplin, paste0(name_plin,
                                  if_else(results$linearity_test$pvalue<0.0001, "p<0.0001", 
                                          paste0("p=", put.pvalue(results$linearity_test$pvalue)))), 
         adj=adj_plin, cex = cex_plin)
  }
}

# linearity_test ####
linearity_test <- function(t2event=NULL, event, model = c("cox", "logistic"), var, adj=NULL, data){
  require(tidyverse)
  results <- data.frame(nk = double(), plin = double(), AIC = double())
  if (model == "cox"){
    data <- data %>% dplyr::select(t2event, event, var, adj) %>% drop_na()
    fml <- formula(paste0("Surv(", t2event, "," , event, ")~",paste0(c(var, adj), collapse=" + ")))
    fit <- coxph(fml, data = data, ties = "breslow")
    results <- rbind(results, data.frame(nk=NA, plin = NA, AIC = AIC(fit)))
    for (k in 3:7){
      knots <- knots.harrell(var=var,data=data,nk=k)
      fml <- formula(paste0("Surv(", t2event, "," , event, ")~", 
                            paste0(c(paste0("I(rcspline.eval(",var,", knots = knots, inclx=T))"), adj), collapse=" + ")))
      fit <- coxph(fml, data = data, ties = "breslow")
      results <- rbind(results, data.frame(nk=k, plin = wald.test(fit, 2:(length(knots)-1))$pvalue, AIC = AIC(fit)))
    }
  }
  
  if (model == "logistic"){
    data <- data %>% dplyr::select(event, var, adj) %>% drop_na()
    fml <- formula(paste0(event , "~",paste0(c(var, adj), collapse=" + ")))
    fit <- glm(fml, data = data, family = "binomial")
    results <- rbind(results, data.frame(nk=NA, plin = NA, AIC = AIC(fit)))
    for (k in 3:7){
      knots <- knots.harrell(var=var,data=data,nk=k)
      fml <- formula(paste0(event, "~", paste0(c(paste0("I(rcspline.eval(",var,", knots = knots, inclx=T))"), adj), collapse=" + ")))
      fit <- glm(fml, data = data, family = "binomial")
      results <- rbind(results, data.frame(nk=k, plin = wald.test(fit, 3:(length(knots)))$pvalue, AIC = AIC(fit)))
    }
  }
  return(results)
}

