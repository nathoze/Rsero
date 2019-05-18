##' @title Plot the fit of the seroprevalence vs. age
##'
##' @description Plot the mean and credible interval of the posterior of the seroprevalence, given an object of the class \code{FOIfit}. See function \code{\link{seroprevalence}} for plotting the seroprevalence from the data. 
##'  The function returns graphical objects that can be further modified.
##' 
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##'  
##' @param FOIfit The \code{FOIfit} object to be plotted.
##'  
##' @param individual_samples  Integer. Number of individual samples to be plotted additionally to the mean and the credible interval of the force of infection. The \code{indivual_samples} samples are randomly chosen in the chains. Default = 0.
##' 
##' @param show_data Boolean. If \code{TRUE}, shows the fraction of seropositive as well. By default, ages are associated in groups of \code{age_cats} years.  Default = \code{TRUE}. 
##' 
##' @param age_cats Integer. Length of the age groups (in years). Used if \code{show_data = TRUE}. Default = 10.
##' 
##' @param YLIM Upper limit of the y-axis. Default = 1. The lower limit is set to 0.
##' 
##' @return A list of \code{ggplot2} objects. 
##' 
##' @examples
##' data <- simulate_SeroData( max_age = 50, epidemic_years = c(1976,1992), foi = c(0.2,0.3))
##' model <- FOImodel(type='outbreak', background=TRUE, K = 2)
##' 
##' Fit <- fit(model = model, data = data)
##' 
##' # plot the fraction of seropositive, with 2 individual samples
##' p <- seroprevalence.fit(Fit,  individual_samples = 2,YLIM=0.2)
##' 
##' 
##' @export
##' @importFrom graphics lines

seroprevalence.fit<- function(FOIfit,                   
                              individual_samples = 0,
                              age_cats = 10,
                              YLIM=1,
                              ...){
  
  plots  <- NULL
  data= FOIfit$data
  
  chains <- rstan::extract(FOIfit$fit)
  bg = chains$bg # background infection probability
  A <- FOIfit$data$A
  latest_sampling_year <- max(FOIfit$data$sampling_year)
  years <- seq(1,A)
  
  
  if(FOIfit$model$cat_bg == 0 && FOIfit$model$cat_lambda==0){
    Ncat = 1
  }  else{
    Ncat = FOIfit$data$Ncategory
  }
  index.plot=0
  
  for(sampling_year in sort(unique(FOIfit$data$sampling_year)) ){
    for(k in 1:Ncat){ 
      
      index.plot=index.plot+1
      age_group = data$age_group[which(data$sampling_year ==  sampling_year)][1]
      w = which(data$sampling_year ==  sampling_year & data$categoryindex==k)
      subdat = subset(data,sub = w)
      
      
      years.plotted =  seq(latest_sampling_year-sampling_year+1, dim(chains$P)[2])
      yrs= years[seq(1, dim(chains$P)[2]-latest_sampling_year+sampling_year )]
      Pinf <- (1-(1-bg[,k])*chains$P[,years.plotted,age_group,k]) # infection probability
      
      par_out <- apply(Pinf, 2, function(x)c(mean(x), quantile(x, probs=c(0.025, 0.975))))
      par_out[par_out>YLIM]= YLIM # set to the upper limits for plotting
      meanFit <- data.frame(x = yrs, y = par_out[1, ])
      
      
      # create the envelope
      xpoly <- (c(yrs, rev(yrs)))
      ypoly <-  c(par_out[3, ], rev(par_out[2, ]))
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      
      # data
      histdata <- sero.age.groups(dat = subdat,age_cats = age_cats,YLIM=YLIM)
      
      
      
      # plot the mean and 95% credible interval of the seroprevalence
      p <- ggplot2::ggplot()
      p <- p + ggplot2::geom_polygon(data=DataEnvelope, ggplot2::aes(x, y), fill="#d7def3")
      p <- p + ggplot2::geom_line(data = meanFit, ggplot2::aes(x = x, y = y), size = 1, color ='#5e6b91')
      
      
      # plot individual runs of the chain
      if(individual_samples>0){
        Index_samples <- sample(nrow(S), individual_samples)
        for (i in Index_samples){
          ind_foi <-  data.frame(x = yrs,y = S[i, ])
          p <- p + ggplot2::geom_line(data = ind_foi, ggplot2::aes(x = x, y = y), size = 0.8, colour = "#bbbbbb", alpha = 0.6)
        }
      }
      
      p <- p  +  scale_x_continuous(breaks=histdata$age,labels=levels(histdata$labels))
      p <- p +  geom_point(data = histdata, aes(x=age, y=mean)) 
      p <- p + geom_segment(data=histdata, aes(x=age,y=lower, xend= age, yend=upper))
      
      p <- p + ggplot2::xlab("Age (years)") + ggplot2::ylab("Proportion seropositive")
      p <- p +theme_classic()
      p <- p +theme(axis.text.x = element_text(size=12),
                    axis.text.y = element_text(size=14),
                    text=element_text(size=14))
      p <-  p + ylim(0,YLIM)
      
      
      plots[[index.plot]] <- p 
      
    }
  }
  return(plots)
  
}



