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
##' @param show_data Boolean. If \code{TRUE}, shows the fraction of seropositive as well. By default, ages are associated in groups of \code{age_class} years.  Default = \code{TRUE}. 
##' 
##' @param age_class Integer. Length of the age groups (in years). Used if \code{show_data = TRUE}. Default = 10.
##' 
##' @param YLIM Upper limit of the y-axis. Default = 1. The lower limit is set to 0.
##' 
##' @param fill.color Ribbon color. By default, fill.color  = "#d7def3"
##' 
##' @param line.color Line color. By default, line.color  = "#5e6b91"
##' 
##' @param mid.age.plot Boolean. Whether the data is plotted at the mid-point of the age category (TRUE) or at the mean age of the individuals in the data within this category (FALSE). Default = \code{TRUE}. 
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
                              age_class = 10,
                              YLIM=1,
                              fill.color  = "#d7def3", 
                              line.color  = "#5e6b91", 
                              mid.age.plot = TRUE,
                              ...){
  
  plots  <- NULL
  data= FOIfit$data
  
  chains <- rstan::extract(FOIfit$fit)
  se = FOIfit$model$se 
  sp = FOIfit$model$sp 
  A <- FOIfit$data$A
  latest_sampling_year <- max(FOIfit$data$sampling_year)
  years <- seq(1,A)
  
  index.plot=0
  unique.categories = data$unique.categories
  sorted.year = sort.int(unique(FOIfit$data$sampling_year),index.return = TRUE)
  Y=0
  for(sampling_year in sorted.year$x ){
    Y=Y+1
    for(cat in unique.categories){ 
      
      if(length(unique.categories)==1){
        #  title =  paste0("Sampling year: ", sampling_year)
        title =  ""
      }
      if(length(unique.categories)>1){
        # title= paste0('Category: ',cat," Sampling year: ", sampling_year)
        title= paste0('Category: ',cat)
      }
      
      index.plot=index.plot+1
      age_group = data$age_group[which(data$sampling_year ==  sampling_year)][1]
      w = which(data$sampling_year ==  sampling_year & data$category==cat, arr.ind = TRUE)[,1]
      subdat = subset(data,sub = w)
      
      # compute the proportion of seropositive
      P=chains$P[,,sorted.year$ix[Y], 1]
      d = data$categoryindex[w]
      p1=proportions.index(d)
      
      M=dim(chains$P)[1] 
      Pinf=matrix(0, nrow = M, ncol=FOIfit$data$A)
      
      # infection probability weighted on the categories      
      for(i in 1:length(p1$index)){
        Pinf =  Pinf+  p1$prop[i]*( se-(se+sp-1)*chains$P[,,sorted.year$ix[Y],p1$index[i]] ) 
      }
      
      par_out <- apply(Pinf, 2, function(x)c(mean(x), quantile(x, probs=c(0.025, 0.975))))
      par_out[par_out>YLIM]= YLIM # set to the upper limits for plotting
      
      # X axis
      years.plotted =  seq(latest_sampling_year-sampling_year+1, dim(chains$P)[2])
      years.plotted.normal= years.plotted-min(years.plotted)+1
      meanFit <- data.frame(x = years.plotted.normal, y = par_out[1,years.plotted ])
      
      # create the envelope
      xpoly <- (c(years.plotted.normal, rev(years.plotted.normal)))
      ypoly <-  c(par_out[3,years.plotted ], rev(par_out[2,years.plotted ]))
      DataEnvelope = data.frame(x = xpoly, y = ypoly)
      
      # histogram  of data
      histdata <- sero.age.groups(dat = subdat,age_class = age_class,YLIM=YLIM, mid.age.plot=mid.age.plot) 
      histdata$labels_text <- as.character(histdata$labels)
      
      last_row_index <- tail(which(complete.cases(histdata)), 1)
      if(length(last_row_index)>0){
        histdata = histdata[1:last_row_index,]
      }
      
      max.age= max(histdata$age)
      
      DataEnvelope =  subset(DataEnvelope, x<=max.age)
      meanFit = subset(meanFit, x<=max.age )
      
      # plot the mean and 95% credible interval of the seroprevalence
      p <- ggplot2::ggplot() + 
        ggplot2::geom_polygon(data=DataEnvelope, ggplot2::aes(x, y), fill=fill.color) + 
        ggplot2::geom_line(data = meanFit, ggplot2::aes(x = x, y = y), linewidth = 1, color =line.color)
      
      # if plot individual runs of the chain
      if(individual_samples>0){
        Index_samples <- sample(nrow(Pinf), individual_samples)
        for (i in Index_samples){
          ind_foi <-  data.frame(x = years.plotted.normal,y = Pinf[i, years.plotted])
          p <- p + ggplot2::geom_line(data = ind_foi, ggplot2::aes(x = x, y = y), linewidth = 0.8, colour = "#bbbbbb", alpha = 0.6)
        }
      }
       p <- p  +
        scale_x_continuous(breaks=histdata$age-round(age_class/2),labels=histdata$labels_text)+  
#        scale_x_continuous(breaks=histdata$age,labels=histdata$labels_text)+  
        geom_point(data = histdata, aes(x=age, y=mean))  +
        geom_segment(data=histdata, aes(x=age,y=lower, xend= age, yend=upper))+
        ggplot2::xlab("Age (years)") + 
        ggplot2::ylab("Seroprevalence") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, size=14,vjust = 0.5),
              axis.text.y = element_text(size=14),
              text=element_text(size=14)) +
        ylim(0,YLIM)+
        ggtitle(title)
      
      
      plots[[index.plot]] <- p 
      plots[[index.plot]]$category <- cat 
      plots[[index.plot]]$year <- sampling_year 
      
    }
  }
  return(plots)
  
}

#' @export
proportions.index <- function(d){
  b.x=c()
  b.y=c()
  ii=0
  for(i in unique(d)){
    ii=ii+1
    b.x[ii]  = i  
    b.y[ii]  = sum(d==i)/length(d)
  }
  return(list(index=b.x,prop = b.y ))
  
}
