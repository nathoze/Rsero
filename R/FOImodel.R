##' @title Model of the force of infection
##'
##' @description This function creates an object of the class \code{FOImodel}. Inputs are the type of the model (required), additional parameters if required by the model, and hyperparameters for the prior distributions (optional). The models to be given as an input are predefined. 
##'  More details are given in the vignette  \code{models}.
##'
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @param type A character with the name of model. The options are: 
##' \itemize{
##'   \item 'constant': Constant force of infection
##'   \item 'outbreak': Series of outbreak modeled with gaussians 
##'   \item 'independent': Annual independent values of the force of infection
##'   \item 'piecewise': Piecewise constant force of infection. The number of phases with a constant level is given by the variable \code{K}
##'   \item 'constantoutbreak': A combination of \code{K} outbreaks with a constant yearly force of infection
##'   \item 'independent_group': Similar to the independent model, but with piecewise constant values of the force annual force of infection in time periods of length \code{group_size} years.
##' }
##' 
##' @param K integer. An additional parameter used in the "outbreak", "constantoutbreak" and "piecewise" stan models. In the case of the outbreaks, this parameter is the number of Gaussians used. In the case of the piecewise constant model it is the number of constant phases. Default = 1.
##' 
##' @param group_size integer. An additional parameter used in the independent_group models. The force of infection is averaged over \code{group_size} year period. By default \code{group_size} = 1, which is equivalent ot the independent model.
##' 
##' @param se numeric, between 0 and 1. If \code{se=1} the assay has a perfect sensitivity. Default = 1. 
##' 
##' @param sp numeric, between 0 and 1. If \code{sp=1} the assay has a perfect specificity. Default = 1. 
##' 
##' @param seroreversion integer, equal to 0 or 1. If \code{seroreversion=0} the model includes a rate of seroreversion (waning immunity). See the vignette \code{models} for details. Default = 0. 
##' 
##' @param cat_lambda integer, equal to 0 or 1. If \code{cat_lambda=1} the force of infection varies across the different categories defined in an objet \code{SeroData}.
##'   See the vignette \code{models} for details. Default = 1. 
##' 
##' @param prioralpha1 First parameter of the prior distribution of the parameter alpha, used as the intensity of the force of infection in the outbreak models. If a normal distribution is chosen (prior_distribution_alpha= "normal"), it corresponds to the mean parameter. If an exponential distribution is chosen (prior_distribution_alpha = "exponential"), it is the rate parameter. Default = 0.2.
##' 
##' @param prioralpha2 Second parameter of the prior distribution of the parameter alpha, used as the intensity of the force of infection in the outbreak models. If a normal distribution is chosen (prior_distribution_alpha= "normal"), it corresponds to the sd parameter on the log scale. Default = 0.2.
##' 
##' @param priorT1 First parameter for the prior distribution for the T parameter. The time of the outbreak is defined as the number of years before the survey (outbreak and constant outbreak models). It is the time interval for the change of FOI in the piecewise constant model. If a normal distribution is chosen (prior_distribution_T= "normal"), it corresponds to the mean parameter.  If an exponential distribution is chosen, it is the rate parameter. Default = 20.
##' 
##' @param priorT2 Second parameter for the prior distribution for the T parameter. If a normal distribution is chosen (prior_distribution_T= "normal"), it corresponds to the sd parameter. Default = 10.
##' 
##' @param priorC1 First parameter of the prior distribution for the constant force of infection, used in the constant and piecewise constant models. If a normal distribution is chosen, it corresponds to the mean parameter on the natural scale.  If an exponential distribution is chosen, it is the rate parameter. Default = 0.01.
##' 
##' @param priorC2 Second parameter of the prior distribution for the constant force of infection, used in the constant and piecewise constant models. If a normal distribution is chosen, it corresponds to the sd parameter on the log scale.  Default = 1.
##' 
##' @param priorY1 First parameter of the prior distribution for the annual hazard of infection, used in the independent models. If a normal distribution is chosen, it corresponds to the mean parameter on the natural scale. If an exponential distribution is chosen, it is the rate parameter. Default = 0.01.
##' 
##' @param priorY2 Second parameter of the prior distribution for the annual hazard of infection, used in the independent models. If a normal distribution is chosen, it corresponds to the sd parameter on the log scale. Default = 1.
##' 
##' @param priorRho1 First parameter of the prior distribution for rho, the seroreversion rate used when seroreversion=1. If a lognormal distribution is chosen (prior_distribution_rho = "normal"), it corresponds to the mean parameter on the natural scale. If an exponential distribution is chosen (prior_distribution_rho = "exponential"), it is the rate parameter. Default = 1.
##' 
##' @param priorRho2 First parameter of the prior distribution for rho, the seroreversion rate used when seroreversion=1. If a lognormal distribution is chosen (prior_distribution_rho = "normal"), it corresponds to the sd parameter on the log scale. Default = 1.
##' 
##' @param ... Additional arguments (not used).
##'
##' @return A list with the class \code{FOImodel}, which contains the
##' following items:
##'
##' \itemize{
##'
##' \item type: The type of the model.
##' 
##' \item stanname: The name of the stan file used. 
##'
##' \item estimated_parameters: The number of estimated parameters.
##' 
##' \item priors: a list with the priors.
##' 
##' \item K: the input parameter used in the outbreak and piecewise models, if given.
##' 
##' \item prior_distribution_alpha: The prior distribution of the alpha parameter used in the outbreaks models. "normal" or "exponential". Default = "normal".
##' 
##' \item prior_distribution_T: The prior distribution of the T parameter used in the outbreaks models and the piecewise constant model. "normal" or "exponential". Default = "normal".
##' 
##' \item prior_distribution_constant_foi: The prior distribution of the constant FOI parameter used in the constant and piecewise constant models. "normal" or "exponential". Default = "normal".
##' 
##' \item prior_distribution_independent_foi: The prior distribution of the annual FOI parameter used in the independent model. "normal" or "exponential". Default = "normal".
##' 
##' \item prior_distribution_rho: The prior distribution of the seroreversion parameter. "normal" or "exponential". Default = "normal".
##' }
##'
##'
##' @examples
##' 
##' 
##' ## A gaussian model, with two gaussians, and user-defined priors
##' model <- FOImodel('outbreak',K = 2, prioralpha1 = 0, prioralpha2 = 1)
##'
##' ## A piecewise model, consisting in two constant phases, with seroreversion and user-defined priors
##' model <- FOImodel('piecewise', K=2, seroreversion=1, priorRho1=0.1)
##'
##' 
##'
##' @export
FOImodel <- function(type = 'constant',
                     K = 1,
                     group_size = 1,
                     seroreversion =0,
                     se = 1, 
                     sp = 1, 
                     prioralpha1 = 0.2,
                     prioralpha2 = 0.2,
                     priorT1 = 20,  
                     priorT2 = 10,
                     priorC1 = 0.01,   
                     priorC2 = 1,  
                     priorY1 =  0.01, 
                     priorY2 = 1, 
                     priorRho1 = 1,
                     priorRho2 = 1,
                     cat_lambda = 1,
                     prior_distribution_alpha = "normal", #  "normal" or "exponential"
                     prior_distribution_T = "normal",#  "normal" or "exponential"
                     prior_distribution_constant_foi = "normal",#  "normal"or "exponential"
                     prior_distribution_independent_foi = "normal",#  "normal" or "exponential"
                     prior_distribution_rho = "normal",#  "normal" or "exponential"
                     ...) {
  
  estimated_parameters <- 0
  
  if (!(type %in% model.list('All models'))){
    print("Model is not defined.")
  }
  
  if(seroreversion){
    estimated_parameters <- estimated_parameters +1
  }
  if(type == "constantoutbreak"){
    estimated_parameters <- estimated_parameters +1
  }
  
  if(type  == 'intervention' | type == 'piecewise'){
    stanname= 'intervention'
  }
  
  if(type  == 'constant'){
    stanname = 'constant'
  }
  
  if(type=='outbreak'){
    stanname= 'outbreak'
  }
  
  if(type=='independent'){
    stanname= 'independent'
  }
  
  if(type=='independent_group'){
    stanname= 'independent_group'
  }
  
  if(type=='constantoutbreak'){
    stanname= 'constantoutbreak'
  }
  
  if(type %in% model.list(whichmodels = 'Outbreak models')  ){
    estimated_parameters <- estimated_parameters + K*2
  } 
  if(type  %in% model.list(whichmodels = 'I models') ){
    estimated_parameters <- estimated_parameters + K*2-1 
  } 
  if(type == 'constant'){
    estimated_parameters <- estimated_parameters + 1
  } 
  # check the necessary parameters are correctly in the input
  if(type %in% model.list('K models') & is.null(K) ){
    print("K not defined.")
  } 
  if(type %in% model.list('K models')  ){
    if(length(priorT1) != length(priorT2)){
      print("priorT1 and priorT2 must have the same length.")
    }
    if(length(priorT1) == 1 ){
      priorT1 = as.array(rep(priorT1,K))
    }
    if(length(priorT2) == 1){
      priorT2 = as.array(rep(priorT2,K))
    }
    if(length(priorT1) != K ){
      print("priorT1 must be of length 1 or K.")
    }
    if(length(priorT2) != K){
      print("priorT2 must be of length 1 or K.")
    }
  } 
  if(type %in% model.list('Outbreak models')){
    if(length(prioralpha1) != length(prioralpha2)){
      print("prioralpha1 and prioralpha2 must have the same length.")
    }
    if(length(prioralpha1) == 1 ){
      prioralpha1 = as.array(rep(prioralpha1,K))
    }
    if(length(prioralpha2) == 1){
      prioralpha2 = as.array(rep(prioralpha2,K))
    }
    if(length(prioralpha1) != K ){
      print("prioralpha1 must be of length 1 or K.")
    }
    if(length(prioralpha2) != K){
      print("prioralpha2 must be of length 1 or K.")
    }
  } 
  
  priors <- list(prioralpha1 = prioralpha1,
                 prioralpha2 = prioralpha2,
                 priorT1 = priorT1,
                 priorT2 = priorT2,
                 priorC1 = priorC1,
                 priorC2 = priorC2,
                 priorY1 = priorY1,
                 priorY2 = priorY2,
                 priorRho1 = priorRho1,
                 priorRho2 = priorRho2)
  
  
  model <- list(type = type,
                stanname = stanname,
                K = K,
                group_size=group_size,
                seroreversion=seroreversion,
                se=se,
                sp=sp,
                cat_lambda = cat_lambda,
                estimated_parameters = estimated_parameters,
                priors = priors,
                prior_distribution_alpha = prior_distribution_alpha,
                prior_distribution_T = prior_distribution_T,
                prior_distribution_constant_foi = prior_distribution_constant_foi,
                prior_distribution_independent_foi = prior_distribution_independent_foi,
                prior_distribution_rho = prior_distribution_rho)
  
  class(model) <- "FOImodel"
  
  return(model)
  
}  

##' @export
##' @rdname FOImodel
##' @param x A object of the class \code{FOImodel}.

print.FOImodel <- function(x, ...){
  cat("<FOImodel object>\n")
  
  cat('Model name: ',x$type,"\n")
  
  if(x$type %in% model.list(whichmodels = 'All models')){
    
    if(x$type== 'independent' |x$type== 'independent_group'){
      cat('\t Estimated parameters: ',x$estimated_parameters ,' + number of age classes \n')
    }  
    else{
      cat('\t Estimated parameters: ',x$estimated_parameters ,'\n')
    }
    if(x$cat_lambda){
      cat('Model with categories for the force of infection \n')
    }
    cat('Parameters: \n')
    
    if(x$type  %in% c('intervention', 'piecewise','outbreak') ){
      cat('\t K: ',x$K ,'\n')
    }
    if(x$se <1 | x$sp < 1){
      cat('\t Sensitivity: ', x$se, '\n')
      cat('\t Specificity: ', x$sp, '\n')
    }
    cat('Priors: \n')
    
    # print(x$prior_distribution_alpha )
    # print(x$prior_distribution_T )
    # print(x$prior_distribution_rho)
    
    if(x$type %in% model.list('Outbreak models')){
      for(i in 1:x$K){
        # if(x$prior_distribution_alpha == "uniform"){
        #   cat('\t alpha: Uniform(',x$priors$prioralpha1[i], ', ', x$priors$prioralpha2[i] ,')\n')
        # }
        if(x$prior_distribution_alpha == "normal"){
          cat('\t alpha: Normal(',x$priors$prioralpha1[i], ', ', x$priors$prioralpha2[i] ,')\n')
        }
        if(x$prior_distribution_alpha == "exponential"){
          cat('\t alpha: Exponential(',x$priors$prioralpha1[i] ,')\n')
        }
        
        # if(x$prior_distribution_T == "uniform"){
        #   cat('\t T: Uniform(',x$priors$priorT1[i], ', ', x$priors$priorT2[i],')\n')
        # }
        if(x$prior_distribution_T == "normal"){
          cat('\t T: Normal(',x$priors$priorT1[i], ', ', x$priors$priorT2[i],')\n')
        } 
        if(x$prior_distribution_T == "exponential"){
        cat('\t T: Exponential(',x$priors$priorT1[i] ,')\n')
      }
        
      }
    }
    
    if(x$type%in% model.list('Constant models')){
      # if(x$prior_distribution_constant_foi == "uniform"){
      #   cat('\t Annual FOI: Uniform(',x$priors$priorC1, ', ', x$priors$priorC2 ,')\n')
      # }
      if(x$prior_distribution_constant_foi == "normal"){
        cat('\t Annual FOI: Normal(',x$priors$priorC1, ', ', x$priors$priorC2 ,')\n')
      }
      if(x$prior_distribution_constant_foi == "exponential"){
        cat('\t Annual FOI: Exponential(',x$priors$priorC1, ')\n')
      }
    }
    
    if(x$type=='independent' | x$type=='independent_group'){
      # if(x$prior_distribution_independent_foi == "uniform"){
      #   cat('\t Annual FOI: Uniform(',x$priors$priorY1, ', ', x$priors$priorY2 ,')\n')
      # }
      if(x$prior_distribution_independent_foi == "normal"){
        cat('\t Annual FOI: Normal(',x$priors$priorY1, ', ', x$priors$priorY2 ,')\n')
      }
      if(x$prior_distribution_independent_foi == "exponential"){
        cat('\t Annual FOI: Exponential(',x$priors$priorY1,')\n')
      }
    }
  
    if(x$type %in% model.list('I models')){  
      for(i in 1:x$K){
        # if(x$prior_distribution_T == "uniform"){
        #   cat('\t T: Uniform(',x$priors$priorT1[i], ', ', x$priors$priorT2[i],')\n')
        # }
        if(x$prior_distribution_T == "normal"){
          cat('\t T: Normal(',x$priors$priorT1[i], ', ', x$priors$priorT2[i],')\n')
        }   
        if(x$prior_distribution_T == "exponential"){
          cat('\t T: Exponential(',x$priors$priorT1[i] ,')\n')
        }
      }
    }
    
    if(x$seroreversion){
      # if(x$prior_distribution_rho == "uniform"){
      #   cat('\t rho: Uniform(',x$priors$priorRho1, ', ', x$priors$priorRho2,')\n')
      # }
      if(x$prior_distribution_rho == "normal"){
        cat('\t rho: Normal(',x$priors$priorRho1, ', ', x$priors$priorRho2,')\n')
      }   
      if(x$prior_distribution_rho == "exponential"){
        cat('\t rho: Exponential(',x$priors$priorRho1,')\n')
      }
    }
    
  }else{
    for(i in seq(2, length(names(x)))){
      if(class(x[[i]]) == 'list'){
        LS.df = as.data.frame(do.call(rbind, x[[i]]))        
        print(LS.df)
      }
      else{
        cat(names(x)[i],": ",  x[[i]],'\n')
      } 
    }
  }
}

#' @export
model.list <-function(whichmodels){
  
  out <- NA
  if(whichmodels == 'K models'){
    out <- c('outbreak','intervention','constantoutbreak','piecewise')
  }
  if(whichmodels == 'I models'){
    #out <- c('intervention','constant')
    out <- c('intervention', 'piecewise')
  } 
  if(whichmodels == 'Constant models'){
    out <- c('constantoutbreak','constant')
  } 
  
  if (whichmodels == 'All models'){
    out <-  c('outbreak','independent','constant','intervention', 'constantoutbreak','independent_group','piecewise')
  } 
  
  if (whichmodels == 'Outbreak models'){
    out <-  c('outbreak', 'constantoutbreak')
  } 
  return(out)
}

