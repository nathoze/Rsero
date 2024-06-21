##' @title Analyse a FOIfit
##'
##' @description Determines the mean and quantiles of the parameters estimated in the fit. By default, it returns the mean, the median and 95 \% credible interval of the parameters.
##'
##' @param FOIfit An object of the class \code{FOIfit}.
##' 
##' @param quants Numeric. Contains the list of the estimated quantiles given between 0 and 1. Default parameters are the 2.5\%, 50 \% and 97.5 \% quantiles given by \code{quants = c(0.025,0.5,0.975)}.
##' 
##' @return A dataframe containing the 2.5%, 50% and 97.5% quantiles and the mean of the posterior distributions. If several categories are specified, the force of infection is given for the reference category and the relative force of infection is given for the other categories.
##' In the case of an outbreak model (outbreak or outbreak+constant), the total force of infection is given by alpha, and the probability of infection by "Outbreak Prob. Inf.". For models with constant phases, the annual probability of infection is given.      
##'
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @examples   
##'  data("1peakSimulation")
##'  model <- FOImodel(type='outbreak', seroreversion=TRUE, K=1)
##'  options(mc.cores = parallel::detectCores())
##'  Fit <- fit(data = data, model = model)
##'  parameters_credible_intervals(Fit)
##' 
##' # an example where other quantiles are specified
##'  data("1peakSimulation")
##'  model <- FOImodel(type='outbreak', seroreversion=TRUE, K=1)
##'  options(mc.cores = parallel::detectCores())
##'  Fit <- fit(data = data, model = model)
##'  parameters_credible_intervals(Fit, quants=c(0.25,0.5,0.75))
##' 
##' @export
parameters_credible_intervals <- function(FOIfit,
                                          quants  = c(0.025,0.5,0.975) ) {
  
  chains <- rstan::extract(FOIfit$fit)
  chainsout= chains
  quantilestext  = paste0(quants*100,'%')
  quantilestext[4]  = 'mean'
  
  params <- NULL
  
  if(FOIfit$model$type %in% model.list('All models')){ 
    
    if(FOIfit$model$type  %in% model.list('Outbreak models')){
      
      C<- chains$T
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
        chainsout$T[,i] <- t(Years[i, ]) 
        chainsout$alpha[, i] <- t(chains$alpha)[Ranks[i, ]]
        # chainsout$beta[, i] <- t(chains$beta)[Ranks[i, ]]
        
        params <- add.quantiles.text(params,
                                     variable=chainsout$T[,i],
                                     name = paste('T',i),
                                     quants= quants,
                                     quantilestext=quantilestext)
        
        name <- paste('alpha',i)
        
        params <- add.quantiles.text(params,
                                     variable=chainsout$alpha[,i],
                                     name = name,
                                     quants= quants,
                                     quantilestext=quantilestext )
         
        # params <- add.quantiles.text(params,
        #                              variable  = chainsout$beta[,i],
        #                              name = paste('beta',i),
        #                              quants= quants,
        #                              quantilestext=quantilestext )      
        # 
        params <- add.quantiles.text(params,
                                     variable  = 1-exp(-chainsout$alpha[,i]),
                                     name = paste('Outbreak Prob. Inf. ',i),
                                     quants= quants,
                                     quantilestext=quantilestext )
      }
    }
    
    if(FOIfit$model$type %in% model.list('I models')){
      K=FOIfit$model$K
      C<- chains$Time
      YearMax <- max(FOIfit$data$sampling_year)
      Years  <- YearMax - C+1
      
      for(i in 1:K){
        if (i>1){ 
          params <- add.quantiles.text(params, variable  = Years[,i],
                                       name = paste('Year_',i),
                                       quants= quants,
                                       quantilestext=quantilestext )
          
        }

        LL <- 100*(1-exp(-chains$annual_foi[,i]))
        params <- add.quantiles.text(params,
                                     variable  = LL,
                                     name = paste0('Annual Prob. Infection (in %)_',i),
                                     quants= quants,
                                     quantilestext=quantilestext )
      }
    }
    
    if(FOIfit$model$type=='constantoutbreak'){
      LL <- 100*(1-exp(-chains$annual_foi))
      params <- add.quantiles.text(params,
                                   variable  = LL,
                                   name = paste0('Annual Prob. Infection (in %)_',i),
                                   quants= quants,
                                   quantilestext=quantilestext )
      
    }
    if(FOIfit$model$type=='independent' | FOIfit$model$type=='independent_group'){
      
      L = chains$lambda
      for(k in seq(1,dim(L)[2])){
        params <- add.quantiles.text(params,
                                     variable=L[,k],
                                     name = paste0('Y_',max(FOIfit$data$sampling_year)-k),
                                     quants=quants,
                                     quantilestext=quantilestext)
        
      }
    }
    
    if(FOIfit$model$seroreversion){
      rho = chains$rho
      params <- add.quantiles.text(params,
                                   variable=rho,
                                   name = paste('rho'),
                                   quants = quants,
                                   quantilestext=quantilestext )
    }
  }
  
  
  d= FOIfit$data$category.position.in.table
  
  if(FOIfit$model$cat_lambda & dim(d)[1]>0){ 
    for(i in seq(1,dim(d)[1])){
      
      name <- paste0('FOI of category ', d[i,]$predictor, " relative to " ,  d[i,]$relative_to)
      params <- add.quantiles.text(params,
                                   variable=chainsout$Flambda[,d[i,]$index],
                                   name = name,
                                   quants= quants, 
                                   quantilestext=quantilestext )
      
    }
    
  } 
  names(params) <- quantilestext
  return(params)
  
}


add.quantiles.text <- function(params, variable, name, quants = quants, quantilestext=quantilestext){
  
  P <- data.frame(t(quantile(variable, probs = quants)), mean = mean(variable))
  names( P) <- quantilestext
  rownames( P)[1] = name
  
  return(rbind(params,P))
  
}

# 
# measure.foi.Gmodel <- function(alpha, beta, Time, years){
#   
#   # To measure the total foi corresponding to a peak, set the parameter alpha of the other peaks to 0 and sum the lambdas
#   J=0
#   S <- matrix(0, ncol= length(alpha), nrow = length(years))
#   
#   for(y in years){
#     J=J+1 
#     S[J, ]=  alpha*exp(-(Time-y)^2/beta^2)  
#     
#   }
#   FOI  = colSums(S )  
#   
#   return(FOI)
#   
# }
# 
