##' @title Model of the force of infection
##'
##' @description This function creates an object of the class \code{FOImodel}. Inputs are the type of the model (required), additional parameters if required by the model, and parameters for the prior distributions. The models to be given as an input are predefined. 
##'  More details are given in the vignette  \code{models}.
##'
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @param type A character with the name of model. The options are: 
##' \itemize{
##'   \item 'constant': Constant force of infection
##'   \item 'outbreak': Series of outbreak modeled with gaussians 
##'   \item 'independent': Yearly independent values of the force of infection
##'   \item 'intervention': Piecewise constant force of infection. The number of phases with a constant level is given by the variable \code{K}
##'   \item 'constantoutbreak': A combination of \code{K} outbreaks with a constant yearly force of infection
##' }
##' 
##' @param K integer. An additional parameter used in the outbreak, constantoutbreak and intervention stan models. This parameter is the number of Gaussians used in the model.  Default = \code{NULL}.
##' 
##' @param background integer, equal to 0 or 1. If \code{background=1}  the model includes a background infection probability. See the vignette \code{models} for details.  Default = 0. 
##' 
##' @param seroreversion integer, equal to 0 or 1. If \code{seroreversion=0} the model includes a rate of seroreversion (waning immunity). See the vignette \code{models} for details.  Default = 0. 
##' 
##' @param cat_lambda integer, equal to 0 or 1.  If \code{cat_lambda=1} the force of infection varies accross the different categories defined in an objet \code{SeroData}.  See the vignette \code{models} for details. Default = 0. 
##' 
##' @param cat_bg integer, equal to 0 or 1. If \code{cat_bg=1} a specific value of the background infection probability \code{pB} is attributed to each category (see \code{SeroData}). See the vignette \code{models} for details. Default = 0. 
##' 
##' @param cat_lambda integer, equal to 0 or 1.  If \code{cat_lambda=1} the force of infection varies accross the different categories defined in an objet \code{SeroData}.  See the vignette \code{models} for details. Default = 0. 
##' 
##' @param prioralpha1 First parameter of the uniform prior distribution  of the parameter alpha, used as the intensity of the force of infection in the outbreak and intervention models. Default = 0.
##' 
##' @param prioralpha2 Second parameter of the uniform prior distribution  of the parameter alpha, used as the intensity of the force of infection in the outbreak and intervention models.  Default = 5.
##' 
##' @param priorbeta1  First parameter of the uniform prior distribution  of the parameter beta, used as the spread of the force of infection in the outbreak and intervention models. Default = 0.
##' 
##' @param priorbeta2  Second parameter of the uniform prior distribution  of the parameter beta, used as the spread of the force of infection in the outbreak and intervention models. Default = 1.
##' 
##' @param priorT1 First parameter for the uniform distribution for the T parameter, used as the time of infection in the outbreak and intervention models. T is defined as the number of years between the survey and the outbreak. Default = 0. 
##' 
##' @param priorT2 Second parameter for the uniform distribution for the T parameter, used as the time of infection in the outbreak and intervention models. Default = 70.
##' 
##' @param priorC1  First parameter of the uniform prior distribution for the constant force of infection, used in the constant and intervention models. Default = 0. 
##' 
##' @param priorC2  Second parameter of the uniform prior distribution for the constant force of infection, used in the constant and intervention models. Default = 10. 
##' 
##' @param priorY1 Mean of the normal distribution for the logit of the annual hazard of infection, used in the independent models. Default = 0. 
##' 
##' @param priorY2 Standard-deviation of the normal distribution for the logit of the annual hazard of infection, used in the independent models.  Default = 100.
##' 
##' @param priorbg1 First  parameter of the uniform prior distribution for the \code{pB} parameter (background probability of infection) when background = 1.  Default = 0.
##' 
##' @param priorbg2 Second shape parameter of the beta distribution for the  \code{pB} parameter  (background probability of infection) when background = 1. Default = 1.
##' 
##' @param priorRho1 First parameter of the uniform prior distribution for rho, the seroreversion rate used when seroreversion=1. Default = 0.
##' 
##' @param priorRho2  Second parameter of the uniform prior distribution for rho, the seroreversion rate used when seroreversion=1. Default = 10.
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
##' \item K: the input parameter used in the outbreak and intervention models, if given.
##' 
##' }
##'
##'
##' @examples
##' 
##' ## A gaussian model, with one gaussian, and background infection
##' model <- FOImodel(type='outbreak', background = 1, K = 1)
##' 
##' ## A gaussian model, with two gaussians, background infection, and user-defined hyperparameters
##' model <- FOImodel('outbreak', background=1, K = 2, prioralpha1 = 0, prioralpha2 = 1)
##'
##' ## An intervention model, consisting in two constant phases, with seroreversion and user-defined hyperparameters
##' model <- FOImodel('intervention', K=2, seroreversion=1, priorRho1=0.1)
##'
##' 
##'
##' @export
FOImodel <- function(type = 'constant',
                     K = 1,
                     seroreversion =0,
                     background = 0, 
                     prioralpha1 = 0,
                     prioralpha2 = 5,
                     priorbeta1 = 0,
                     priorbeta2 = 5,
                     priorT1 = 0,  
                     priorT2 = 70,
                     priorC1 = 0, # 0
                     priorC2 = 10, # 100
                     priorY1 = 0, 
                     priorY2 = 100, 
                     priorbg1 = 0,
                     priorbg2 = 1,
                     priorRho1  =0,
                     priorRho2 = 10,
                     cat_bg = 0,  # cat_background
                     cat_lambda = 0,
                     ...) {
  
  
  estimated_parameters <- 0
  
  if (!(type %in% model.list('All models'))){
    print("Model is not defined.")
  }
  
  if(background){
    estimated_parameters <- estimated_parameters +1
  }
  
  if(seroreversion){
    estimated_parameters <- estimated_parameters +1
  }
  if(type == "constantoutbreak"){
    estimated_parameters <- estimated_parameters +1
  }
  
  if(type %in% model.list(whichmodels = 'All models')){
    stanname= 'intervention'
  }
  if(type=='outbreak'){
    stanname= 'outbreak'
  }
  
  if(type=='independent'){
    stanname= 'independent'
  }
  
  if(type=='constantoutbreak'){
    stanname= 'constantoutbreak'
  }

  if(type %in% model.list(whichmodels = 'Outbreak models')  ){
    estimated_parameters <- estimated_parameters + K*3
  } 
  if(type  %in% model.list(whichmodels = 'I models') ){
    estimated_parameters <- estimated_parameters + K*2-1 
  } 
  
  # check the necessary parameters are correclty in the inputs
  if(type %in% model.list('K models') & is.null(K) ){
    print("K not defined.")
  } 
  
  priors <- list(prioralpha1 = prioralpha1,
                 prioralpha2 = prioralpha2,
                 priorbeta1 = priorbeta1,
                 priorbeta2 = priorbeta2,
                 priorT1 = priorT1,
                 priorT2 = priorT2,
                 priorC1 = priorC1,
                 priorC2 = priorC2,
                 priorY1 = priorY1,
                 priorY2 = priorY2,
                 priorbg1 = priorbg1,
                 priorbg2 = priorbg2,
                 priorRho1 = priorRho1,
                 priorRho2 = priorRho2)
  
  
  model <- list(type = type,
                stanname = stanname,
                K = K,
                seroreversion=seroreversion,
                background=background,
                cat_bg = cat_bg, 
                cat_lambda = cat_lambda,
                estimated_parameters = estimated_parameters,
                priors = priors)
  
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
    
    if(x$type== 'independent'){
      cat('\t Estimated parameters: ',x$estimated_parameters ,' + number of age classes \n')
    }
    else{
      cat('\t Estimated parameters: ',x$estimated_parameters ,'\n')
    }
    
    if(x$cat_bg){
      cat('Model with categories for background infection \n')
    }
    if(x$cat_lambda){
      cat('Model with categories for the force of infection \n')
    }
    cat('Parameters: \n')
    
    if(x$type  %in% c('intervention','outbreak') ){
      cat('\t K: ',x$K ,'\n')
    }
    
    cat('Priors: \n')
    if(x$type %in% model.list('Outbreak models')){
      cat('\t alpha1: ',x$priors$prioralpha1 ,'\n')
      cat('\t alpha2: ',x$priors$prioralpha2 ,'\n')
      cat('\t beta1: ',x$priors$priorbeta1 ,'\n')
      cat('\t beta2: ',x$priors$priorbeta2 ,'\n')
      cat('\t T1: ',x$priors$priorT1 ,'\n')
      cat('\t T2: ',x$priors$priorT2 ,'\n')
    }
    
    
    if(x$type=='constant'){
      cat('\t C1: ',x$priors$priorC1 ,'\n')
      cat('\t C2: ',x$priors$priorC2 ,'\n')
    }
    
    if(x$type=='independent'){
      cat('\t Y1: ',x$priors$priorY1 ,'\n')
      cat('\t Y2: ',x$priors$priorY2 ,'\n')
    }
    if(x$type=='intervention'){
      cat('\t T1: ',x$priors$priorT1 ,'\n')
      cat('\t T2: ',x$priors$priorT2 ,'\n')
    }
    if(x$background){
      cat('\t Background 1: ', x$priors$priorbg1, '\n')
      cat('\t Background 2: ', x$priors$priorbg2, '\n')
    }
    if(x$seroreversion){
      cat('\t Seroreversion 1: ', x$priors$priorRho1, '\n')
      cat('\t Seroreversion 2: ', x$priors$priorRho2, '\n')
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
    out <- c('outbreak','intervention','constantoutbreak')
  }
  if(whichmodels == 'I models'){
    out <- c('intervention','constant')
  } 
  if(whichmodels == 'Constant models'){
    out <- c('constantoutbreak','constant')
  } 
  
  if (whichmodels == 'All models'){
    out <-  c('outbreak','independent','constant','intervention', 'constantoutbreak')
  } 
 # if (whichmodels == 'P models'){
#    out <- c('PG','PGB')
 # } 
  if (whichmodels == 'Outbreak models'){
    out <-  c('outbreak', 'constantoutbreak')
  } 
  
  return(out)
}

