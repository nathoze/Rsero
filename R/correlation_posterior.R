##' @title Plot the correlation matrix of the output
##'
##' @description Compute and plot the correlation matrix of the parameters obtained from a MCMC fit.  It returns a n-by-n matrix, where n is the number of parameters.
##' 
##'
##' @param FOIfit An object of the class \code{FOIfit}.
##' 
##' @param show_coefficients Boolean. If  \code{TRUE}, shows the value of the correlation coefficient on the plotted matrix. Default  = \code{TRUE}.
##' 
##' @return The correlation matrix.
##' 
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @examples 
##' data('two_peaks_simulation')
##' model <- FOImodel('outbreak',K=2, background = 1)
##' F <- fit(data=data, model=model)
##' M <- correlation_posterior(F)
##' @export
correlation_posterior <- function(FOIfit,
                                  show_coefficients = TRUE) {
  
  chains <- rstan::extract(FOIfit$fit)
  chainsout= chains
  
  newdf <- NULL
  J=0
  if(FOIfit$model$background){
    
    if(FOIfit$model$cat_bg){ 
      for(k in seq(1,FOIfit$data$Ncategory)){  
        J=J+1
        newdf <- cbind(newdf, chainsout$bg[,k] )
        colnames(newdf)[J] <- paste0('bg cat',k)
      }
    }else{
      J=J+1
      newdf <- cbind(newdf, chains$bg[,1])
      colnames(newdf)[J] <- 'bg'
    }
    
  }
  
  
  if(FOIfit$model$seroreversion){
    J=J+1
    newdf <- cbind(newdf, chainsout$rho )
    colnames(newdf)[J] <- 'rho'  
  }
  
  
  if(FOIfit$model$type =='outbreak'){
    
    C<- chains$T
    Torder <- data.frame(C, t(apply(-C, 1, rank, ties.method='min')))
    K=FOIfit$model$K
    Ranks <- matrix(0,ncol(C),nrow(C))
    
    YearMax <- max(FOIfit$data$sampling_year)
    Years = matrix(0,ncol(C),nrow(C))
    
    
    for (i in 1:K) { 
      Ranks[i, ] <- which(apply(C,1,function(x) rank(x)) == i);
      Years[i, ] <- YearMax - t(C)[Ranks[i,]]
    }
    for(i in 1:K){
      
      chainsout$T[,i] <- t(Years[i, ]) 
      chainsout$alpha[, i] <- t(chains$alpha)[Ranks[i,]]
      chainsout$beta[, i] <- t(chains$beta)[Ranks[i,]]
      
      newdf <- cbind(newdf, chainsout$T[,i])
      J=J+1
      colnames(newdf)[J]=paste0('T_',i)
      
      
      
      if(FOIfit$model$cat_lambda){ 
        for(k in seq(1,FOIfit$data$Ncategory)){  
          J=J+1
          newdf <- cbind(newdf,  chainsout$alpha[,i]*chainsout$Flambda[,k] )
          colnames(newdf)[J] <- paste0('alpha ',i, ' cat ',k)
        }
      }else{
        J=J+1
        newdf <- cbind(newdf, chainsout$alpha[,i])
        colnames(newdf)[J] <- paste0('alpha ',i)
      }
      
      
      newdf <- cbind(newdf, chainsout$beta[,i])
      J=J+1
      colnames(newdf)[J]=paste0('beta_',i)
      
    }
  }
  
  
  if(FOIfit$model$type =='constant'){
    
    L = chains$lambda
    
    newdf <- cbind(newdf,  L[,dim(L)[2]])
    J=J+1
    colnames(newdf)[J] = 'C'
    
  }
  
  
  if(FOIfit$model$type =='independent'){
    L = chains$lambda
    for(k in seq(1,dim(L)[2])){
      
      newdf <- cbind(newdf, L[,k])
      J=J+1
      colnames(newdf)[J] = paste0('Y_', FOIfit$data$sampling_year-k)
      
    }
  }
  
  if(FOIfit$model$type %in% model.list('I models')){
    K=FOIfit$model$K
    YearMax <- max(FOIfit$data$sampling_year)
    C<- chains$Time
    YearMax <- max(FOIfit$data$sampling_year)
    Years  <- YearMax - C+1
    
    for(i in 1:K){
      if (i>1){ 
        newdf <- cbind(newdf, Years[,i])
        J=J+1
        colnames(newdf)[J] = paste0('Year_', i)
      }
      foi <- chains$foi[,i]
      newdf <- cbind(newdf, foi)
      J=J+1
      colnames(newdf)[J] = paste0('FOI_', i)
      
    }
  }
  
  
  M<-cor(newdf)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  
  if(show_coefficients){ 
    corrplot::corrplot(M, method="color", col=col(200),
                       addCoef.col = "black",
                       tl.col="black",
                       tl.srt=45)
  }else{
    corrplot::corrplot(M, method="color", col=col(200),
                       tl.col="black",
                       tl.srt=45)
  }
  
  return(M)
  
  
}
