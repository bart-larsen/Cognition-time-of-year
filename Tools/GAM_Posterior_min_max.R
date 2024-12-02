#POSTERIOR DISTRIBUTION SMOOTH FUNCTIONS
##Function to simulate the posterior distribution from a fitted GAM, calculate smooths for individual posterior draws, and return smooth max and min values + 95% credible intervals

gam.posterior.smooths <- function(gam.model,smooth_var,factor_var=NULL, set_fx = FALSE, draws=1000, increments=200, return_draws = TRUE,UNCONDITIONAL = FALSE,make_plots=FALSE,newdata=NULL){
  
  # Check the model
  if (any(class(gam.model)=="gam")) {
    gam.model <- gam.model
  } else if (class(gam.model$gam)=="gam") {
    gam.model <- gam.model$gam
  } else {
    stop("Can't find a gam object")
  }
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior smooths estimated
  np <- as.numeric(increments) #number of smooth_var increments to predict fit at from posterior model coefficients
  EPS <- 1e-07 #finite differences
  #UNCONDITIONAL should we account for uncertainty when estimating smoothness parameters? default = FALSE
  set.seed(42)

  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  all_terms <- attr(terms(df),"term.labels")
  
  if (is.null(newdata)) {
    if (is.null(factor_var)) {
      newdata <- with(df,expand.grid(x=seq(from=min(df[,smooth_var]),to=max(df[,smooth_var]),length.out=increments)))
      colnames(newdata)=smooth_var
    } else {
      newdata <- with(df,expand.grid(x=seq(from=min(df[,smooth_var]),to=max(df[,smooth_var]),length.out=increments),f=unique(df[,factor_var])))
      colnames(newdata)=c(smooth_var,factor_var)
    }
    covariates <- all_terms[!(all_terms%in%c(smooth_var,factor_var))]
    if (!is_empty(covariates)) {
      covariate_values <- df %>% summarize(across(.cols = all_of(covariates),.fns = function(x){ifelse(is.numeric(x),yes = mean(x,na.rm=T),no = levels(x)[2])}))
      for (cov in covariates) {
        newdata[,cov]=covariate_values[,cov]
      }
    }
  }
 
  
  #Estimate posterior smooth functions (fitted values) from simulated GAM posterior distribution  
  ##Each of the posterior draws has a fitted spline (+ intercept + covariate coefficients) that includes the uncertainty in the estimated model coefficients  
  ##Akin to fitted_samples from gratia  
  Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
  sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance. 
  X0 <- predict(gam.model, newdata = newdata, type = "lpmatrix") #get matrix of linear predictors that maps model parameters to the smooth fit (outcome measure scale)
  predicted.smooth.values <- X0 %*% t(sims) #generate posterior smooths (fitted y for each set of posterior draw model parameters)
  colnames(predicted.smooth.values) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
  predicted.smooth.values <- newdata %>% bind_cols(predicted.smooth.values) #add smooth_var increments from pred df to first column
  
  #Smooth minimum/maximum values and credible intervals
  #Smooth max + 95% credible interval

if (!is.null(factor_var)) {
  max.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    group_by_at(factor_var)%>% #if there is one
    summarise(across(contains("draw"),
                     .fns = function(x){
                       predicted.smooth.values[,smooth_var][which.max(x)]
                     }))
  max.y.draws<-max.y.range%>%pivot_longer(contains("draw"),names_to = "draw",values_to = "max.ys")
  max.y.out <- max.y.draws %>% group_by_at(factor_var)%>%
    summarise(max.y = median(max.ys),max.y.CI = list(quantile(max.ys, probs = c(0.025, 0.975))))%>%
    unnest_wider(max.y.CI)%>%
    rename(lower=`2.5%`,upper=`97.5%`)

  if (make_plots==TRUE) {
    p<-ggplot(max.y.draws,aes_string(x="max.ys",fill=factor_var))+geom_histogram(position = "identity")+
      ggtitle("Max")
    print(p)
  }
  #Smooth min + 95% credible interval
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    group_by_at(factor_var)%>% #group output by factor level
    summarise(across(contains("draw"),
                     .fns = function(x){
                       predicted.smooth.values[,smooth_var][which.max(x)]
                     }))
  min.y.draws<-min.y.range%>%pivot_longer(contains("draw"),names_to = "draw",values_to = "min.ys")
  min.y.out <- min.y.draws %>% group_by_at(factor_var)%>%
    summarise(min.y = median(min.ys),min.y.CI = list(quantile(min.ys, probs = c(0.025, 0.975))))%>%
    unnest_wider(max.y.CI)%>%
    rename(lower=`2.5%`,upper=`97.5%`)
  
  if (make_plots==TRUE) {
    p<-ggplot(min.y.draws,aes_string(x="min.ys",fill=factor_var))+geom_histogram(position = "identity")+
      ggtitle("Max")
    print(p)
  }
} else {
  max.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       as.numeric(predicted.smooth.values[which.max(x),smooth_var])
                     }))
  max.y.draws<-max.y.range%>%pivot_longer(contains("draw"),names_to = "draw",values_to = "max.ys")

  max.y <- median(max.y.draws$max.ys) #median value 
  max.y.CI <- quantile(max.y.draws$max.ys, probs = c(0.025, 0.975)) #credible interval
  max.y.out <- data.frame(max.y = max.y,lower=as.numeric(max.y.CI[1]),upper=as.numeric(max.y.CI[2]))
  if (make_plots==TRUE) {
    p<-ggplot(max.y.draws,aes(x=max.ys))+geom_histogram()+
      geom_text(label=sprintf("Median = %1.3f; 95%% CI [%1.3f, %1.3f]",max.y,max.y.CI[1],max.y.CI[2]),
                x=max.y,
                y=1,color="red")+
      ggtitle("max value")
    print(p)
  }
  #Smooth min + 95% credible interval
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is lowest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       as.numeric(predicted.smooth.values[which.min(x),smooth_var])
                     }))  
  min.y.range<-min.y.range%>%pivot_longer(contains("draw"),names_to = "draw",values_to = "min.y.range")
  min.y <- median(min.y.range$min.y.range) #median value
  min.y.CI <- quantile(min.y.range$min.y.range, probs = c(0.025, 0.975)) #credible interval
  min.y.out <- data.frame(min.y = min.y,lower=as.numeric(min.y.CI[1]),upper=as.numeric(min.y.CI[2]))
  if (make_plots==TRUE) {
    
    p<-ggplot(min.y.range,aes(x=min.y.range))+geom_histogram()+
      geom_text(label=sprintf("Median = %1.3f; 95%% CI [%1.3f, %1.3f]",min.y,min.y.CI[1],min.y.CI[2]),
                x=min.y,
                y=3,color="red")+
    ggtitle("min value")
    print(p)
  }
}

  if(return_draws == TRUE)
    return(predicted.smooth.values)
    # return(min.y.range)
  if(return_draws == FALSE)
    smooth.features <- list(max.y.out, min.y.out)
  names(smooth.features) <- c(sprintf("%s at max y and CI", smooth_var),sprintf("%s at min y and CI", smooth_var))
  return(smooth.features)
  }
