resid_plot <- function(modobj,term,int.term=NULL,group.term=NULL,add.intercept=FALSE,xlabel=NULL,ylabel=NULL,plot.title = NULL){
  library(gratia)
  terms=c(term,int.term)
  model_terms <- term_names(modobj)
  # Check the model
  if (any(class(modobj)=="gam")) {
    modobj <- modobj
  } else if (class(modobj$gam)=="gam") {
    modobj <- modobj$gam
  } else {
    stop("Can't find a gam object")
  }
  
  
  df <- modobj$model%>%select(model_terms)
  mod.intercept <- modobj$coefficients["(Intercept)"]
  pterms <- predict(modobj,type = "terms",se.fit = TRUE)

  if (!is.null(int.term)) {
    intcol=df %>% select(int.term)
    int_type <- ifelse(class(intcol[[1]])=="numeric",yes ="numeric",no = "factor")
  }
  
  
  if (add.intercept==TRUE) {
    pterms.fit <- pterms$fit+mod.intercept
  } else{
    pterms.fit <- pterms$fit
  }
  pterms.sefit <- pterms$se.fit

  colnames(pterms.fit) <- gsub(x = colnames(pterms.fit),pattern = "s\\(",replacement = "")%>%
    gsub(pattern = "\\)",replacement = "")
  colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit),pattern = "s\\(",replacement = "")%>%
    gsub(pattern = "\\)",replacement = "")
  
  pterms.df <- data.frame(pterms.fit)%>%select(matches(terms))%>%
    mutate(fit=rowSums(across(where(is.numeric))))%>%
    select(fit)%>%
    # rename(fit:=!!term)%>%
    # cbind(data.frame(pterms.sefit)%>%select(matches(term))%>%rename(se.fit:=!!term))%>%
    cbind(data.frame(pterms.sefit)%>%select(matches(terms))%>%mutate(se.fit=rowSums(across(where(is.numeric)))))%>%
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)%>%
    mutate(partial.residuals=fit+resid(modobj))%>%
    select(fit,se.fit,lwr,upr,partial.residuals)
  
  partial.residuals.df <- cbind(pterms.df,df%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))

  # partial.residuals.df<-data.frame(pterms.fit)%>%
  #   mutate(across(.cols = everything(),.fns = function(x){x+resid(modobj)}))%>%
  #   cbind(df%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))
  
  # plot.df <- cbind(partial.residuals,pterms.df)

  if (is.null(int.term)) {
    plot.obj <- ggplot(partial.residuals.df,aes_string(x=paste0(term,".raw"),y="partial.residuals"))+
      geom_hex()+
      # geom_point(alpha=.75,color="white",pch=21,fill="#444444")+
      geom_ribbon(aes(y=fit,ymin=lwr,ymax=upr),fill="gray",alpha=.25,color=NA)+
      geom_line(aes(y=fit),color="black",linewidth=1.5)+
      xlab(term)+ylab("partial.residuals")
  } else{
    plot.obj <- ggplot(partial.residuals.df,aes_string(x=paste0(term,".raw"),
                              y="partial.residuals",
                              color=paste0(int.term,".raw"),
                              fill=paste0(int.term,".raw")))+
      geom_point(alpha=.75,color="white",pch=21)+
      geom_ribbon(aes(y=fit,ymin=lwr,ymax=upr),alpha=.55,color=NA)+
      geom_line(aes(y=fit),linewidth=1.5)+
      scale_color_brewer(type = "qual",palette = "Set2")+scale_fill_brewer(type = "qual",palette = "Set2")+
      xlab(term)+ylab("partial.residuals")
  }
  if (!is.null(xlabel)) {
    plot.obj <- plot.obj + xlab(xlabel)+ylab(ylabel)
  }
  if (!is.null(plot.title)){
    plot.obj <- plot.obj + ggtitle(plot.title)
  }
  return(plot.obj)

  }
    
  
