#' @title Creates a subset of a serological survey 
#'
#' @description This function creates a serological survey of the class \code{SeroData}.
#'
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#' 

#' @param sex An optional character factor defining the sex of the individuals. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}. 
#'
#' @param sub  Used in the \code{subset} function. A list of the indices of the individuals to be considered in the new \code{SeroData} object.  
#'
#' @param loc  Used in the \code{subset} function. A list of the \code{location} of the individuals to be considered in the new  \code{SeroData} object.
#' 
#' @param ... Additional arguments (not used).
#'  
#' 
#' @return A list with the class \code{SeroData}, which contains the
##' following items:
##'
##' \itemize{
##'
##' \item age: A vector of integers containing the age of the sampled individuals at the time of the latest sampling, of length N.
##' 
##' \item age_at_sampling: A vector of integers containing the age of the sampled individuals at the time of sampling, of length N.
##' 
##' \item Y: A vector of seropositivity status, of length N.
##' 
##' \item N: The number of individuals considered.
##'
##' \item A: The maximal age.
##'
##' \item NGroups: The total number of age groups.
##'
##' \item sampling_year: A vector of sampling years, of length N.
##'
##' \item location: A vector of the sampling location, of length N.
##'
##' \item sex: A vector of the sex of each individual, of length N.
##'
##' }
##'
##'
#' @export
#' @examples
#' 
#' ## A very simple example of a serological survey with three individuals: 
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1), max_age = 50, age_cats = 1, sampling_year = 2017)
#' 
#' seroprevalence(data = data)
#' data2 = subset(data,c(1,3))
#' 
#' # add example with subset(data, location)
#' 


subset.SeroData <- function(data,sub =seq(1,data$N), loc = NULL, category = NULL){

  
  if(!is.null(loc)){
    sub1 = which(data$location %in% loc)
    m = match(sub, sub1)
    sub  = m[!is.na(m)]#??
    sub=sub1
  }

  if(!is.null(category)){ 
    sub1 = which(data$category==category, arr.ind = TRUE)[,1]
    m = match(sub, sub1)
    sub  = m[!is.na(m)]#??
    sub=sub1
  }
  
  Y=as.matrix(data$Y)
  
  
  sub.age  = data$age_at_sampling[sub]+max(data$sampling_year[sub])-data$sampling_year[sub]
  
  age.groups <-compute.age.groups(age = sub.age,sampling_year = data$sampling_year[sub])
  
  
  N=length(sub)
  category =  matrix(data$category[sub,],nrow = N) 
  
  param.category = category.parameters(category,N)
  
  
  
  subdata <- list( A = data$A,
                   NGroups = data$NGroups,
                   N = N,
                   Y = Y[sub,],
                   age = sub.age,
                   age_at_sampling = data$age_at_sampling[sub],
                   #ind_by_age = data$ind_by_age,
                   sampling_year = data$sampling_year[sub],
                   location = data$location[sub],
                   sex = data$sex[sub],
                   category = category,
                   categoryindex=param.category$categoryindex,
                   MatrixCategory = param.category$MatrixCategory,
                   Ncategory = param.category$Ncategory,
                   maxNcategory=param.category$maxNcategory,
                   Ncategoryclass=param.category$Ncategoryclass,
                   unique.categories=param.category$unique.categories,
                   Ncat.unique = param.category$Ncat.unique,
                   NAgeGroups = age.groups$NAgeGroups,
                   age_at_init =  as.array(age.groups$age_at_init), 
                   age_group  = age.groups$age_group)
  
  
  class(subdata) <- 'SeroData'
  
  return(subdata)
  
}
