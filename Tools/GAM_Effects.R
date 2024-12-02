library(mgcv)

partialR2 <- function(full_mod,reduced_mod){
  if (any(class(full_mod)=="gam")) {
    is_GAMM=FALSE
    y <- full_mod$y
    y_hat <- full_mod$fitted.values
    y_reduced <- reduced_mod$y
    y_hat_reduced <- reduced_mod$fitted.values
  } else if (class(full_mod$gam)=="gam") {
    is_GAMM=TRUE
    y <- full_mod$mer@resp$y
    y_hat <- fitted(full_mod$mer)
    y_reduced <- reduced_mod$mer@resp$y
    y_hat_reduced <- fitted(reduced_mod$mer)
  }
  if (is_GAMM==TRUE) {
    
  }

  SSres = sum((y-y_hat)^2)
  SStot = sum((y-mean(y))^2)
  SSreg = SStot - SSres
  R2 = 1-(SSres/SStot)
  
  SSres_reduced = sum((y_reduced-y_hat_reduced)^2)
  SStot_reduced = sum((y_reduced-mean(y_reduced))^2)
  SSreg_reduced = SStot_reduced - SSres_reduced
  R2_reduced = 1-(SSres_reduced/SStot_reduced)
  
  partial_r2 <- (SSreg - SSreg_reduced)/SSres_reduced
  
  # partial_r2 <- (SSres_reduced - SSres)/SSres_reduced #these are both the same
  return(partial_r2)
}

partialf2 <- function(full_mod,reduced_mod){
  SSres = sum((full_mod$y-full_mod$fitted.values)^2)
  SStot = sum((full_mod$y-mean(full_mod$y))^2)
  SSreg = SStot - SSres
  R2 = 1-(SSres/SStot)
  
  SSres_reduced = sum((reduced_mod$y-reduced_mod$fitted.values)^2)
  SStot_reduced = sum((reduced_mod$y-mean(reduced_mod$y))^2)
  SSreg_reduced = SStot_reduced - SSres_reduced
  R2_reduced =1-(SSres_reduced/SStot_reduced)
  
  #cohens f2
  partial_f2 <- (R2-R2_reduced)/(1-R2)
  return(partial_f2)
}

model_effects <- function(modobj,return_penalized=FALSE,modobj2=NULL,random_effects=FALSE,random_slope=FALSE,method="partial_R2"){
  if (any(class(modobj)=="gam")) {
    modobj <- modobj
  } else if (class(modobj$gam)=="gam") {
    modobj <- modobj$gam
  } else {
    stop("Can't find a gam object for modobj")
  }

  model_formula <- modobj$formula
  response_variable <- all.vars(model_formula)[1]
  penalized_formula <- as.formula(paste(gsub(x = deparse(modobj$formula),pattern="fx = T",replacement =  "fx = F")%>%gsub(x=.,pattern = "ti\\(",replacement = "t2("),collapse=""))

  if (!is.null(modobj2)) {
    if (any(class(modobj2)=="gam")) {
      modobj2 <- modobj2
    } else if (class(modobj2$gam)=="gam") {
      modobj2 <- modobj2$gam
    } else {
      stop("Can't find a gam object in modobj2")
    }
  } else { #  # if no second model was given, we automatically compare the given model to a model that drops the last term.

    thisResp <- as.character(modobj$terms[[2]])
    theseVars <- attr(terms(model_formula),"term.labels")
    new_formula <- reformulate(theseVars[0:(length(theseVars)-1)],response = thisResp)
    
    if (random_effects == TRUE) {
      if (random_slope == TRUE) {
        g <- gamm4::gamm4(as.formula(new_formula),
                          data=modobj$model,
                          REML=TRUE,random = ~(BLOOD_AGE|BBL_ID))
        g <-g$gam
        
        penalized_gam<-gamm4::gamm4(as.formula(penalized_formula),
                                    data=modobj$model,
                                    REML=TRUE,random = ~(BLOOD_AGE|BBL_ID))
      } else {
        g <-gamm4::gamm4(as.formula(new_formula),
                         data=modobj$model,REML=TRUE,random = ~(1|BBL_ID))
        g <-g$gam
        penalized_gam <-gamm4::gamm4(as.formula(penalized_formula),
                                     data=modobj$model,REML=TRUE,random = ~(1|BBL_ID))
      }
      
    } else {
      g <-gam(as.formula(new_formula),
              data=modobj$model)
      penalized_gam<-gam(as.formula(penalized_formula),
                         data=modobj$model)
    }
    modobj2 <- g

  }
  
  # get effect size now
  if (method %in% c("partial_f2","partial_F2")) {
    partial_f2 <- partialf2(full_mod = modobj,reduced_mod = modobj2)
    # cat(sprintf("\nCohen's f2 = %1.4f\n",
    #             partial_f2))
    effect_size = partial_f2
  } else {
    partial_r2 <- partialR2(full_mod = modobj,reduced_mod = modobj2)
    # cat(sprintf("\nPartial R2 = %1.4f\n",
    #             partial_r2))
    effect_size = partial_r2
  }

  
  if (return_penalized==TRUE) {
    result <- penalized_gam
    return(result)
  } else { return(effect_size)}
  
}
