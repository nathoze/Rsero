##' @title Plot the force of infection
##'
##' @description Plot the mean and credible interval of the yearly force of infection, given an object of the class \code{FOIfit}. The function \code{plot} uses the package \code{ggplot2}. It returns a graphical object that can be further modified.
##' 
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##'  
##' @param FOIfit The \code{FOIfit} object to be plotted.
##' 
##' @param mean_only Boolean. If \code{FALSE}, plot the mean force of infection as well as an envelope of the 95\% credible interval. Default =  \code{FALSE}. 
##'  
##' @param individual_samples  Integer. Number of individual samples to be plotted additionally to the mean and the credible interval of the force of infection. The \code{individual_samples} samples are randomly chosen in the chains. Default = 0.
##' 
##' @param YLIM Upper limit of the y-axis. Default = 1. The lower limit is set to 0. 
##' 
##' @return A list of \code{ggplot2} objects. 
##' 
##' @examples
##' data <- simulate_SeroData( max_age = 50, epidemic_years = c(1976,1992), foi = c(0.2,0.3))
##' model <- FOImodel(type='outbreak', background=1, K = 2)
##' Fit <- fit(model = model, data = data)
##' p <- plot(Fit)
##' p+ylim(0,1)
##' 
##' 
##' @export
##' @rdname plot
##' @importFrom graphics lines

plot.FOIfit <- function(FOIfit,                   
                        mean_only = FALSE,
                        individual_samples = 0,
                        YLIM=1,
                        maxYears = NULL,
                        ...){
  
  chains <- rstan::extract(FOIfit$fit)
  L1 <- chains$lambda
  plots  <- NULL
  
  
  if(FOIfit$model$cat_bg == 0 && FOIfit$model$cat_lambda==0){
    Ncat = 1
  }  else{
    Ncat = FOIfit$data$Ncategory
  }
  
  # 3/09/2018  ??
  if(FOIfit$model$cat_lambda==0){
    Ncat = 1
  }  else{
    Ncat = FOIfit$data$Ncategory
  }
  # Ne pas repeter les plot si cat_bg = TRUE, seulement si cat_lambda = TRUE?

  
  #### HEREEHEREHEREREERH HERE
  # QU EST CE QUON MONTRE DANS LES FIT? TOUTES LES CATEGORIES??
  # copier plot_seropositive
  
  d= FOIfit$data$category.position.in.table
  
  if(FOIfit$model$cat_lambda & dim(d)[1]>0){ 
    
    Ncat = dim(d)[1]>0
    for(i in seq(1,dim(d)[1])){
      
      name <- paste0('FOI of category ', d[i,]$predictor, " relative to " ,  d[i,]$relative_to)
      params <- add.quantiles.text(params,
                                   variable=chains$Flambda[,d[i,]$index],
                                   name = name,
                                   quants= quants, 
                                   quantilestext=quantilestext )
    }
  } 
  #### HEREEHEREHEREREERH HERE
  
  
  
  
  
  for(k in seq(1,Ncat)){
    
    L =  chains$Flambda[,k]*L1
    
    par_out <- apply(L, 2, function(x)c(mean(x), quantile(x, probs=c(0.025, 0.975))))
    latest_sampling_year <- max(FOIfit$data$sampling_year)
    
    if(is.null(maxYears)  ){
      maxYears =FOIfit$data$A
    } 
    if(maxYears >FOIfit$data$A ){
      maxYears =FOIfit$data$A
    } 
    
    yrs <- latest_sampling_year-seq(1,maxYears)+1
    par_out= par_out[ , seq(1,length(yrs))]
    par_out[which(par_out>YLIM)]=YLIM
    meanData <- data.frame(x = yrs, y = par_out[1, seq(1,length(yrs))])
    p <- ggplot2::ggplot()
    
    if(!mean_only){
      xpoly <- (c(rev(yrs), (yrs)))
      ypoly <-  c(rev(par_out[3, ]), par_out[2, ])
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      p <- p + ggplot2::geom_polygon(data=DataEnvelope, ggplot2::aes(x, y), fill="#d7def3", alpha=0.9)
    }
    
    if(individual_samples>0){
      Index_samples <- sample(nrow(L), individual_samples)
      for (i in Index_samples){
        ind_foi <-  data.frame(x = yrs,y = L[i, ])
        p <- p + ggplot2::geom_line(data = ind_foi, ggplot2::aes(x = x, y = y), size = 0.8, colour = "#bbbbbb", alpha = 0.6)
      }
    } 
    
    p <- p+theme_classic()
    
    p <- p+theme(axis.text.x = element_text(size=22),
                 axis.text.y = element_text(size=22),
                 text=element_text(size=16))
    
    p <- p + ggplot2::geom_line(data = meanData, ggplot2::aes(x = x, y = y), size = 1,color  ="#5e6b91",alpha=1)
    
    p <- p + ggplot2::xlab("Year") + ggplot2::ylab("Force of infection")+ylim(0,YLIM)
    plots[[k]] <- p 
  }
  
  return(plots)
}

