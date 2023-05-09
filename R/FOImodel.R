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
##'   \item 'intervention': Piecewise constant force of infection. The number of phases with a constant level is given by the variable \code{K}
##'   \item 'constantoutbreak': A combination of \code{K} outbreaks with a constant yearly force of infection
##'   \item 'independent_group': Similar to the independent model, but with piecewise constant values of the force annual force of infection in time periods of length \code{group_size} years.
##' }
##' 
##' @param K integer. An additional parameter used in the outbreak, constantoutbreak and intervention stan models. This parameter is the number of Gaussians used in the model.  Default = 1.
##' 
##' @param group_size integer. An additional parameter used in the independent_group models. The force of infection is averaged over \code{group_size} year period. By default \code{group_size} = 1, which is equivalent ot the independent model.
##' 
##' @param se_sp integer, equal to 0 or 1. If \code{background=1}  the model includes an imperfect sensitivity and specificity of the assay.
##'  See the vignette \code{models} for details.  Default = 0. 
##' 
##' @param seroreversion integer, equal to 0 or 1. If \code{seroreversion=0} the model includes a rate of seroreversion (waning immunity). See the vignette \code{models} for details.  Default = 0. 
##' 
##' @param cat_lambda integer, equal to 0 or 1.  If \code{cat_lambda=1} the force of infection varies accross the different categories defined in an objet \code{SeroData}.
##'   See the vignette \code{models} for details. Default = 1. 
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
##' @param priorY1 First parameter of the uniform prior distribution for the annual hazard of infection, used in the independent models. Default = 0. 
##' 
##' @param priorY2 Second parameter of the uniform prior distribution for the annual hazard of infection, used in the independent models.  Default = 10.
##' 
##' @param priorse1 First  parameter of the uniform prior distribution for the \code{se} parameter (sensitivity of the assay) when se_sp = 1. Default = 0.
##' 
##' @param priorse2 Second parameter of the uniform prior distribution for the \code{se} parameter (sensitivity of the assay) when se_sp = 1. Default = 1.
##' 
##' @param priorsp1 First  parameter of the uniform prior distribution for the \code{sp} parameter (specificity of the assay) when se_sp = 1. Default = 0.
##' 
##' @param priorsp2 Second parameter of the uniform prior distribution for the \code{sp} parameter (specificity of the assay) when se_sp = 1. Default = 1.
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
                     group_size=1,
                     seroreversion =0,
                     se_sp = 0, 
                     prioralpha1 = 0,
                     prioralpha2 = 5,
                     priorbeta1 = 0,
                     priorbeta2 = 1,
                     priorT1 = 0,  
                     priorT2 = 70,
                     priorC1 = 0, # 0
                     priorC2 = 10, # 100
                     priorY1 = 0, 
                     priorY2 = 10, 
                     priorse1 = 0,
                     priorse2 = 1,
                     priorsp1 = 0,
                     priorsp2 = 1,
                     priorRho1  =0,
                     priorRho2 = 10,
                     cat_lambda = 1,
                     ...) {
  
  
  estimated_parameters <- 0
  
  if (!(type %in% model.list('All models'))){
    print("Model is not defined.")
  }
  
  if(se_sp){
    estimated_parameters <- estimated_parameters +2
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
  
  if(type=='independent_group'){
    stanname= 'independent_group'
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
                 priorse1 = priorse1,
                 priorse2 = priorse2,   
                 priorsp1 = priorsp1,
                 priorsp2 = priorsp2,
                 priorRho1 = priorRho1,
                 priorRho2 = priorRho2)
  
  
  model <- list(type = type,
                stanname = stanname,
                K = K,
                group_size=group_size,
                seroreversion=seroreversion,
                se_sp=se_sp,
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
    
    if(x$type=='independent' | x$type=='independent_group'){
      cat('\t Y1: ',x$priors$priorY1 ,'\n')
      cat('\t Y2: ',x$priors$priorY2 ,'\n')
    }
    if(x$type=='intervention'){
      cat('\t T1: ',x$priors$priorT1 ,'\n')
      cat('\t T2: ',x$priors$priorT2 ,'\n')
    }
    if(x$se_sp){
      cat('\t Sensitivity 1: ', x$priors$priorse1, '\n')
      cat('\t Sensitivity 2: ', x$priors$priorse2, '\n')
      cat('\t Specificity 1: ', x$priors$priorsp1, '\n')
      cat('\t Specificity 2: ', x$priors$priorsp2, '\n')
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
    out <-  c('outbreak','independent','constant','intervention', 'constantoutbreak','independent_group')
  } 
  # if (whichmodels == 'P models'){
  #    out <- c('PG','PGB')
  # } 
  if (whichmodels == 'Outbreak models'){
    out <-  c('outbreak', 'constantoutbreak')
  } 
  return(out)
}

