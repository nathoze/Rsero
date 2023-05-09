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
##' model <- FOImodel(type='outbreak', se_sp=1, K = 2)
##' Fit <- fit(model = model, data = data)
##' p <- plot(Fit)
##' p+ylim(0,1)
##' p[[1]]$category # the name of the category plotted 
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
  
  index.plot=0
  for(cat in FOIfit$data$unique.categories){ 
    
    index.plot=index.plot+1
    
    
    w = which( FOIfit$data$category==cat, arr.ind = TRUE)[,1]
    
    
    d = FOIfit$data$categoryindex[w]
    p1=proportions.index(d)
    
    M=dim(chains$P)[1] 
    L=matrix(0, nrow = M, ncol=FOIfit$data$A)
    # weighted average of the force of infection for each subcategory
    
    for(i in 1:length(p1$index)){
      L =  L+  p1$prop[i]*chains$Flambda[, p1$index[i]]*L1 
    }
    
    
    # compute mean and 95% credible interval  
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
    
    plots[[index.plot]] <- p 
    plots[[index.plot]]$category <- cat 
    
    
  }
  
  return(plots)
}
