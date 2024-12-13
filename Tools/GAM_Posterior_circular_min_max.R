#POSTERIOR DISTRIBUTION SMOOTH FUNCTIONS
##Function to simulate the posterior distribution from a fitted GAM, calculate smooths for individual posterior draws, and return smooth max and min values + 95% credible intervals
library(circular)

# Function to calculate the circular median and confidence interval for doy data
circular_median_with_ci <- function(data) {
  radian.data <- (data - 1) * (2 * pi) / 365.25
  
  # Calculate circular median
  # cos_median <- median(cos(radian.data))
  # sin_median <- median(sin(radian.data))
  # median_radian <- atan2(sin_median, cos_median)
  # median_day <- (median_radian * 365.25 / (2 * pi)) + 1
  radian.data <- circular(radian.data,units = "radians")
  median.radian <- median.circular(radian.data)
  if (median.radian<0) {
    # median.radian=attr(median.radian,which = "medians")
    median.day <- ((median.radian* 365.25 / (2 * pi))+1)+365.25
  }else {
    median.day <- (median.radian* 365.25 / (2 * pi))+1
  }

  # Calculate circular confidence interval
  alpha <- 0.025
  # lower_radian <- quantile(radian.data, alpha)
  # upper_radian <- quantile(radian.data, 1 - alpha)
  lower.radian <- quantile.circular(radian.data,probs = alpha)
  upper.radian <- quantile.circular(radian.data,probs = 1-alpha)
  
  # Map confidence intervals back to day of the year
  # lower_day <- (lower_radian * 365.25 / (2 * pi)) + 1
  # upper_day <- (upper_radian * 365.25 / (2 * pi)) + 1
  lower.day <- (lower.radian* 365.25 / (2 * pi))+1
  upper.day <- (upper.radian* 365.25 / (2 * pi))+1
  
  return(list(median = median.day, ci.lower = lower.day, ci.upper = upper.day))
}

gam.posterior.smooths.circular <- function(gam.model,smooth_var,factor_var=NULL, set_fx = FALSE, draws=1000, increments=200, return_draws = TRUE,UNCONDITIONAL = FALSE,make_plots=FALSE,newdata=NULL,cov_list = NULL){
  
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
      covariate_values <- df %>% select(all_of(covariates)) %>% summarise_all(.funs = function(x){ifelse(is.numeric(x),yes = mean(x,na.rm=T),no = levels(factor(x))[1])})
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
  # Get the column indices corresponding to the s(doy) term
  # We actually only want the predicted values for the term of interest, not including contributions of covariates.
  term_labels <- attr(X0, "dimnames")[[2]]
  s_doy_cols <- grep("s\\(doy\\)", term_labels)
  
  # Zero out all columns except for those corresponding to s(doy) or other listed vars
  if (is.null(cov_list)) {
    X0_s_doy <- X0
    X0_s_doy[, -s_doy_cols] <- 0
  } else {
    X0_s_doy <- X0
    cov_cols <- grep(paste(cov_list, collapse = "|"), term_labels)
    total_cols <- unique(c(cov_cols, s_doy_cols))
    X0_s_doy[, -total_cols] <- 0
  }
  
  predicted.smooth.values <- X0_s_doy %*% t(sims) #generate posterior smooths (fitted y for each set of posterior draw model parameters)
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
  
  mci <- circular_median_with_ci(max.y.draws$max.ys)

  max.y <-as.numeric(mci$median) #median value 
  max.y.out <- data.frame(max.y = max.y,lower=as.numeric(mci$ci.lower),upper=as.numeric(mci$ci.upper))
  if (make_plots==TRUE) {
    p<-ggplot(max.y.draws,aes(x=max.ys))+geom_histogram()+
      geom_vline(xintercept = mci$median,show.legend = T)+
      geom_vline(xintercept = mci$ci.lower,color="red",linetype=2,show.legend = T)+
      geom_vline(xintercept = mci$ci.upper,color="red",linetype=2)+
      coord_polar(start = -pi/12)+xlim(0,365.25)+ggtitle("max")
    print(p)
    
    if (mci$ci.upper<mci$ci.lower) {
      lower.coord.max<-365.25
      upper.coord.min<-0
      p<-ggplot(max.y.draws,aes(x=max.ys))+
        geom_rect(aes(xmin=mci$ci.lower,xmax=lower.coord.max,ymin=0,ymax=1,fill="CI"))+
        geom_rect(aes(xmin=upper.coord.min,xmax=mci$ci.upper,ymin=0,ymax=1,fill="CI"))+
        coord_polar(start = -pi/12)+
        geom_vline(aes(xintercept = mci$median,color="median"),show.legend = T)+
        scale_color_manual(name = "statistics", values = c(median = "blue",CI="black"))+
        scale_fill_manual(name = "statistics", values = c(CI = "lightgray"))+theme(legend.title = element_blank())
      print(p)
    }
  }
  #Smooth min + 95% credible interval
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is lowest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       as.numeric(predicted.smooth.values[which.min(x),smooth_var])
                     }))  
  min.y.range<-min.y.range%>%pivot_longer(contains("draw"),names_to = "draw",values_to = "min.y.range")%>%
    mutate(min.y.range.day=lubridate::parse_date_time(as.character(floor(min.y.range)),orders="j"))
  mci <- circular_median_with_ci(min.y.range$min.y.range)
  mci.day <- data.frame(mci)%>%
    mutate(across(.cols=where(is.numeric),.fns=~lubridate::parse_date_time(as.character(floor(.)),orders="j")))
  mci.numeric <- data.frame(mci)%>%
    mutate(across(.cols=where(is.circular),.fns=as.numeric))
  min.y <-as.numeric(mci$median) #median value 
  min.y.out <- data.frame(min.y = min.y,lower=as.numeric(mci$ci.lower),upper=as.numeric(mci$ci.upper))

  if (make_plots==TRUE) {
    p<-ggplot(min.y.range,aes(x=min.y.range))+geom_histogram(alpha=.5)+
      geom_vline(xintercept = mci.numeric$median,show.legend = T)+
      geom_vline(xintercept = mci.numeric$ci.lower,color="red",linetype=2,show.legend = T)+
      geom_vline(xintercept = mci.numeric$ci.upper,color="red",linetype=2)+
      coord_polar(start = -pi/12)+xlim(0,365.25)+
      # scale_x_datetime(date_labels = "%B",limits = lubridate::parse_date_time(c("1","365"),orders="j"))+
      ggtitle("Min")
      
    print(p)
    
    if (mci$ci.upper<mci$ci.lower) {
      lower.coord.max<-365.25
      upper.coord.min<-0
      p<-ggplot(min.y.range,aes(x=min.y.range.day))+
        geom_rect(aes(xmin=mci.day$ci.lower,xmax=lower.coord.max,ymin=0,ymax=1,fill="CI"))+
        geom_rect(aes(xmin=upper.coord.min,xmax=mci.day$ci.upper,ymin=0,ymax=1,fill="CI"))+
        coord_polar(start = -pi/12)+
        scale_x_datetime(date_labels = "%B",limits = lubridate::parse_date_time(c("1","365"),orders="j"))+
        geom_vline(aes(xintercept = mci.numeric$median,color="median"),show.legend = T)+
        scale_color_manual(name = "statistics", values = c(median = "blue",CI="black"))+
        scale_fill_manual(name = "statistics", values = c(CI = "lightgray"))+theme(legend.title = element_blank())
      print(p)
    } else{
      p<-ggplot(min.y.range,aes(x=min.y.range))+
        geom_rect(aes(xmin=mci.numeric$ci.lower,xmax=mci.numeric$ci.upper,ymin=0,ymax=1,fill="CI"))+
        coord_polar(start = -pi/12)+
        # scale_x_datetime(date_labels = "%B",limits = lubridate::parse_date_time(c("1","365"),orders="j"))+
        geom_vline(aes(xintercept = mci.numeric$median,color="median"),show.legend = T)+
        scale_color_manual(name = "statistics", values = c(median = "blue",CI="black"))+
        scale_fill_manual(name = "statistics", values = c(CI = "lightgray"))+theme(legend.title = element_blank())
      print(p)
    }
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
