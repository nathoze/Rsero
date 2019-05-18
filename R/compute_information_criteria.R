#' @title Compute the AIC, DIC and WAIC 
#'
#' @description This function returns the Aikake Information Criterion (AIC), the Deviance Information Criterion (DIC) and the
#'  Watanabe Aikake Information Criterion (WAIC) from a fit of the class  \code{\link{FOIfit}}.
#'  
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#'  
#' @param FOIfit A \code{FOIfit} object
#' 
#' @return A list with the class \code{information_criteria}, which contains the
##' following items:
##'
##' \itemize{
##'
##' \item AIC: The Aikake Information Criterion.
##' 
##' \item DIC: The Deviance Information Criterion. 
##' 
##' \item WAIC: The Wakanabe-Aikake Information Criterion. 
##' 
##' \item MLE: The Maximum-likelihood, estimated as the maximal value of the likelihood in the drawn samples, used in the AIC.
##' 
##' \item k: Number of parameters, used in the AIC.
##' 
##' \item Dbar: Mean deviance, used in the DIC.
##' 
##' \item pD: Effective number of parameters of the model, used in the DIC.
##' 
##' \item pwaic: Estimated effective number of parameters, used in the WAIC.
##' 
##' \item lpd: log pointwise predictive density, used in the WAIC.
##' }
##' 
#' @examples
#' data('one_peak_simulation')
#' model <- FOImodel(type='outbreak', K=1, background=1)
#' F1  = fit(model = model, data = data )
#' compute_information_criteria(FOIfit  = F1)
#' 
#' @references See Gelman et al. Stat Comput (2014) 24:997-1016
#' DOI 10.1007/s11222-013-9416-2
#' @export
compute_information_criteria <- function(FOIfit,...){
  
  estimated_parameters <- FOIfit$model$estimated_parameters
  chains <- rstan::extract(FOIfit$fit)
  FOIs <- chains$lambda
  bg <- chains$bg 
  
  M <- nrow(FOIs)
  N <- FOIfit$data$N
  A <- FOIfit$data$A 
  Ncategory <- FOIfit$data$Ncategory
  NAgeGroups <- FOIfit$data$NAgeGroups
  
  
  LogLikelihoods <- matrix(0, nrow = M, ncol = N) # as many elements as there are lambdas
  Y <- FOIfit$data$Y
  age <- FOIfit$data$age
  category <-FOIfit$data$categoryindex
  
  for (i in seq(1,M)){
    
    lk  = chains$Like[i,]
    for (j in seq(1,N)){
      
      if( Y[j] == FALSE){ # if the individual is seronegative
        L = log(1-lk[j])   
      }  else{
        L=log(lk[j])
      }
      
      LogLikelihoods[i,j] <- L
    }
  }
  
  
  # log-likelihood on the mean lambdas 
  
  # posterior mean
  
  pc=colMeans(bg)
  LogLikelihoodMean <- 0
  P <- (colMeans(chains$P))
  
  for (j in seq(1,N) ){
    
    age <- FOIfit$data$age[j]
    age_group <- FOIfit$data$age_group[j]
    cat <- category[j]
    
    p <- P[age,age_group, cat]
    
    if(Y[j] == TRUE){
      LogLikelihoodMean <- LogLikelihoodMean + log(1-(1-pc[cat])*p)#log(1-pc[category[j]])-cumfoi[age[j],category[j]]
    }else{
      LogLikelihoodMean <- LogLikelihoodMean + log((1-pc[cat])*p)# log(1 - (1-pc[category[j]])*exp( -cumfoi[age[j],category[j]] ))
    }
  }
  
  LP <- rowSums(LogLikelihoods)
  # Compute the AIC
  # Assumes the maximum likelihood is reached 
  AIC <- -2*max(LP) +2*estimated_parameters
  
  
  # Compute the DIC
  
  Dbar = -2*mean(LP)
  Dthetabar = -2*LogLikelihoodMean
  pD = Dbar-Dthetabar
  
  DIC = pD+Dbar 
  
  
  # Compute the WAIC
  #variance along the column. Each individual has its own variance measured over all sampled parameters 
  
  V = ColVar(exp(LogLikelihoods))
  pwaic = sum(V)
  lpd = sum(log( colSums(exp(LogLikelihoods))/M))
  WAIC <- -2*(lpd-pwaic)
  
  
  information_criteria <- list(AIC = AIC,
                               DIC = DIC,
                               pD = pD,
                               Dbar=Dbar,
                               WAIC = WAIC,
                               pwaic = pwaic,
                               lpd = lpd,
                               k = estimated_parameters,
                               MLE= max(LP))
  
  
  class(information_criteria) <- "information_criteria"
  
  return(information_criteria)
  
}




#' @export
print.information_criteria <- function(x,...){
  
  cat(sprintf('AIC:  %f, MLE: %f, k:  %f\n' , x$AIC, x$MLE, x$k))
  cat(sprintf('DIC:  %f, Dbar:  %f, pD:  %f\n' , x$DIC, x$Dbar, x$pD))
  cat(sprintf('WAIC:  %f, pwaic:  %f, lpd: %f \n' , x$WAIC, x$pwaic, x$lpd))
  
}


RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

ColVar <- function(x) {
  colSums((x - colMeans(x))^2)/(dim(x)[1] - 1)
}

