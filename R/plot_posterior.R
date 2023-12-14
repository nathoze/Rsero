##' @title Plot posterior distribution
##'
##' @description  Plot the posterior distribution of the parameters estimated in a \code{FOIfit}. 
##'
##' @param FOIfit An object of the class \code{FOIfit}.
##' 
##' @return A list with the posterior plots. These plots are generated with the \code{ggplot2} package and can be modified.  
##' 
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @examples 
##' 
##' # Here describe how to access the plots 
##' P <- plot_posterior(FOIfit)
##' P[[1]]+ylim(0,200)
##' 
##' @export
plot_posterior<- function(FOIfit) {
  
  chains <- rstan::extract(FOIfit$fit)
  chainsout <- chains
  plots  <- NULL
  plotindex <- 0
  
  name <- FOIfit$model$type
  
  if(name %in% model.list('All models')){
    
    if(name %in% model.list('Constant models')){# constant or constantoutbreak
      
      gT <- plot_histogram(distribution = data.frame(X= chains$annual_foi), 
                           title = 'Histogram for Annual FOI',
                           xlabel = "FOI",
                           ylabel = 'Count') 
      
      plotindex <- plotindex+1
      plots[[plotindex]]  = gT        
    }
    if(name %in% model.list('I models')){

      C<- chains$Time
      K=FOIfit$model$K
      
      YearMax <- max(FOIfit$data$sampling_year)
      Years = matrix(0,ncol(C),nrow(C))
      
      Years  <- YearMax - C+1
      
      for(i in 1:K){
        if(name != 'constant'){
          chainsout$T[,i] <- Years[,i]
          gT <- plot_histogram(distribution = data.frame(X= chainsout$T[,i] ), 
                               title = paste('Histogram for Time',i),
                               xlabel = "Year",
                               ylabel = 'Count') 
          
          plotindex <- plotindex+1
          plots[[plotindex]]  = gT
        }
      }
      
      for(i in 1:K){
        gT <- plot_histogram(distribution = data.frame(X= chains$annual_foi[,i] ), 
                             title = paste('Histogram for annual FOI',i),
                             xlabel = "FOI",
                             ylabel = 'Count') 
        
        
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT        
      }
    }
    
    #if(name =='outbreak'){
    if(name %in% model.list('Outbreak models')){
      
      C<- chains$T
      S <- chains$S
      Torder <- data.frame(C, t(apply(-C, 1, rank, ties.method='min')))
      K=FOIfit$model$K
      Ranks <- matrix(0,ncol(C),nrow(C))
      
      YearMax <- max(FOIfit$data$sampling_year)
      Years = matrix(0,ncol(C),nrow(C))
      
      for (i in 1:K) { 
        Ranks[i, ] <- which(apply(C,1,function(x) rank(x)) == i);
        Years[i, ] <- YearMax - t(C)[Ranks[i,]]+1
      }
      
      for(i in 1:K){
        if(FOIfit$model$cat_lambda){
          Ncat = FOIfit$data$Ncategory
        }else{
          Ncat = 1
        }
        chainsout$T[,i] <- t(Years[i, ]) 
        chainsout$alpha[, i] <- t(chains$alpha)[Ranks[i,]] 
        chainsout$beta[, i] <- t(chains$beta)[Ranks[i,]]
        
        gT <- plot_histogram(distribution = data.frame(X= chainsout$T[,i] ), 
                             title = paste('Histogram for Time',i),
                             xlabel = "Year",
                             ylabel = 'Count') 
        
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
        for(k in 1:Ncat){
          if(Ncat==1){
            Title = 'Histogram for alpha'
          }
          else{
            Title = paste0('Histogram for alpha, category ',k)
          }
          
          gT <- plot_histogram(distribution = data.frame(X= chains$Flambda[,k]*chainsout$alpha[,i]), 
                               title = paste(Title,i),
                               xlabel = "Value",
                               ylabel = 'Count') 
          
          plotindex <- plotindex+1
          plots[[plotindex]]  = gT
        }
        
        for(k in 1:Ncat){
          if(Ncat==1){
            Title = 'Histogram for total force of infection'
          }
          else{
            Title = paste0('Histogram for total force of infection, category ',k)
          }
          
          gT <- plot_histogram(distribution = data.frame(X= 1-exp(-chains$Flambda[,k]*chainsout$alpha[,i])), 
                               title = paste(Title,i),
                               xlabel = "Value",
                               ylabel = 'Count') 
          plotindex <- plotindex+1
          plots[[plotindex]]  = gT
        }
        
        gT <- plot_histogram(distribution = data.frame(X = chainsout$beta[,i]), 
                             title = paste('Histogram for beta',i),
                             xlabel = "Value",
                             ylabel = 'Count') 
        
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
      }
    }
    
    if(name =='independent' | name=='independent_group'){
      
      # L = chains$lambda
      # for(k in seq(1,dim(L)[2])){
      #   de<-data.frame(t(quantile(L[,k], probs = quants)), mean = mean(L[,k]))
      #   names(de) <- quantilestext
      #   rownames(de)[1] =paste0('Y_',FOIfit$data$sampling_year-k)
      #   params<- rbind(params,de)  
      # }
      # 
     }
    
    if(FOIfit$model$seroreversion){
      
      gT <- plot_histogram(distribution = data.frame(X = chainsout$rho), 
                           title = 'Histogram for Seroreversion',
                           xlabel = "Value",
                           ylabel = 'Count') 
       
      plotindex <- plotindex+1
      plots[[plotindex]]  = gT      
    }
    
    d= FOIfit$data$category.position.in.table
    
    if(FOIfit$model$cat_lambda & dim(d)[1]>0){ 
      for(i in seq(1,dim(d)[1])){
        
        Title <- paste0('Histogram of FOI of category ', d[i,]$predictor, " relative to " ,  d[i,]$relative_to)
        
        #   var=chains$Flambda[,d[i,]$index]
        # distribution<-data.frame(X = var)
        # gT <- ggplot2::ggplot(distribution, aes(X)) +
        #   geom_histogram(bins = 16, fill = "red", col="red",alpha=.4) +
        #   labs(title = Title, x='value', y='Count')
        
        gT <- plot_histogram(distribution = data.frame(X =chains$Flambda[,d[i,]$index]), 
                             title = Title,
                             xlabel = "Value",
                             ylabel = 'Count') 
        
        
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
      }
    } 
  }
  return(plots)
  
}

plot_histogram <- function(distribution, title, xlabel,  ylabel='Count'){
  
  gT <- ggplot(distribution, aes(X)) +
    geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
    labs(title = title, x=xlabel, y=ylabel)+
    theme_bw()+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))  
  
  
  return(gT)
}
