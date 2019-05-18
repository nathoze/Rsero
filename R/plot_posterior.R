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
    
    if(name %in% model.list('I models')){
      
      C<- chains$Time
      K=FOIfit$model$K
      
      YearMax <- max(FOIfit$data$sampling_year)
      Years = matrix(0,ncol(C),nrow(C))
      
      Years  <- YearMax - C+1
      
      
      for(i in 1:K){
        chainsout$T[,i] <- Years[,i]
        
        distribution<-data.frame(Time = chainsout$T[,i])
        gT <- ggplot(distribution, aes(Time)) +
          geom_histogram(bins = 16, fill = "red", col="red",alpha=.4) +
          labs(title = paste('Histogram for Time',i), x='Year', y='Count')
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
      }
      
      for(i in 1:K){
        distribution<-data.frame(X= chains$foi[,i])
        gT <- ggplot(distribution, aes(X)) +
          geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
          labs(title = paste('Histogram for FOI',i), x='FOI', y='Count')
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT        
      }
    }
    
    if(name =='outbreak'){
      
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
        #chainsout$alpha[, i] <- t(chains$alpha[Ranks[i,]]/chains$beta[Ranks[i, ]])
        chainsout$alpha[, i] <- t(chains$alpha)[Ranks[i,]] 
        chainsout$beta[, i] <- t(chains$beta)[Ranks[i,]]
        
        distribution<-data.frame(Time = chainsout$T[,i])
        gT <- ggplot(distribution, aes(Time)) +
          geom_histogram(bins = 16, fill = "red", col="red",alpha=.4) +
          labs(title = paste('Histogram for Time',i), x='Year', y='Count')
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
        for(k in 1:Ncat){
          if(Ncat==1){
            Title = 'Histogram for alpha'
          }
          else{
            Title = paste0('Histogram for alpha, category ',k)
          }
          
          
          distribution<-data.frame(alpha = chains$Flambda[,k]*chainsout$alpha[,i])
          gT <- ggplot(distribution, aes(alpha))+
            geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
            labs(title = paste(Title,i), x='Value', y='Count')
          plotindex <- plotindex+1
          plots[[plotindex]]  = gT
        }
        
        
        distribution<-data.frame(beta = chainsout$beta[,i])
        gT <- ggplot(distribution, aes(beta)) +
          geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
          labs(title = paste('Histogram for beta',i), x='Value', y='Count')
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
      }
    }
    
    if(name =='constant'){
      
      L = chains$lambda[,1]
      
      distribution <- data.frame(L=L) 
      gT <- ggplot(distribution, aes(L))+
        geom_histogram(bins = 12, fill = "red", col="red",alpha=.4)+
        labs(title = 'Histogram for constant force of infection', x='Value', y='Count')
      
      plotindex <- plotindex+1
      plots[[plotindex]]  = gT      
      
    }
    
    if(name =='independent'){
      
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
      distribution<-data.frame(rho = chainsout$rho)
      gT <- ggplot2::ggplot(distribution, aes(rho)) +
        geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
        labs(title = 'Histogram for Seroreversion', x='Value', y='Count')
      
      plotindex <- plotindex+1
      plots[[plotindex]]  = gT      
    }
    
    if(FOIfit$model$background){
      bg = chainsout$bg
      if(FOIfit$model$cat_bg){
        Ncat = FOIfit$data$Ncategory
      }else{
        Ncat = 1
      }
      
      for(k in 1:Ncat){
        if(Ncat==1){
          Title = 'Histogram for background infection probability'
        }
        else{
          Title = paste0('Histogram for background infection probability, category ',k)
        }
        
        distribution<-data.frame(bg = bg[,k])
        gT <- ggplot2::ggplot(distribution, aes(bg)) +
          geom_histogram(bins = 12, fill = "red", col="red",alpha=.4) +
          labs(title = Title, x='value', y='Count')
        
        plotindex <- plotindex+1
        plots[[plotindex]]  = gT
        
      }
    }
  }
  
  return(plots)
  
}
