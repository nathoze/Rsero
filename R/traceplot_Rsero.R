##' @title Markov chain Traceplot
##'
##' @description  Draw the traceplot corresponding to one or more Markov chains, providing a visual way to inspect sampling behavior and assess mixing across chains and convergence.
##'  It is an adaptation of the function \code{\link[rstan]{traceplot}} in the rstan package. 
##'  
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @param FOIfit A \code{FOIfit} object.
##' 
##' @param object	An instance of class stanfit.
##' @param pars	A character vector of parameter names. Defaults to relevant parameters of the model. 
##' @param include	Should the parameters given by the pars argument be included (the default) or excluded from the plot? Only relevant if pars is not missing.
##' @param inc_warmup TRUE or FALSE, indicating whether the warmup sample are included in the trace plot; defaults to FALSE.
##' @param window A vector of length 2. Iterations between window[1] and window[2] will be shown in the plot. The default is to show all iterations if inc_warmup is TRUE and all iterations from the sampling period only if inc_warmup is FALSE. If inc_warmup is FALSE the iterations specified in window should not include iterations from the warmup period.
##' @param unconstrain Should parameters be plotted on the unconstrained space? Defaults to FALSE.
##' @param nrow, ncol	 Passed to facet_wrap.
##' @param ...	Optional arguments to pass to geom_path (e.g. size, linetype, alpha, etc.).
##' 
##' @examples
##' ##'   data <- simulate_SeroData(number_samples = 1000,
##'   age_class = 1,
##'   epidemic_years = c(1976,1992),
##'   foi = c(0.2,0.3))
##' model <- FOImodel('outbreak', K = 2)
##' F1 <- fit(model = model, data = data)
##' traceplot_Rsero(F1)
##' 
##' @export
traceplot_Rsero <- function(FOIfit,
                      pars = NULL, include = TRUE, unconstrain = FALSE, 
                      inc_warmup = FALSE, window = NULL, nrow = NULL, ncol = NULL,...){
  
  name <- FOIfit$model$type
  if(name %in% model.list('All models')){
    
    if(name %in% model.list('Constant models')){# constant or constantoutbreak
      pars = c(pars, "annual_foi")
    }
    if(name %in% model.list('I models')){
      pars = c(pars, "Time")
      pars = c(pars, "annual_foi")
    }
    
    if(name %in% model.list('Outbreak models')){
      pars = c(pars, "T", "alpha")
    }
    
    if(FOIfit$model$seroreversion){
      pars =c(pars, "rho")
    }
    
    d= FOIfit$data$category.position.in.table
    if(FOIfit$model$cat_lambda & dim(d)[1]>0){ 
      pars =c(pars, "Flambda")
    }
  } 
  
  rstan::traceplot(object = FOIfit$fit,
                   pars = pars,
                   include = include, unconstrain = unconstrain, 
                   inc_warmup = inc_warmup, window = window, nrow = nrow, ncol = ncol)
  
  
}
