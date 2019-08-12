#' @title Compute the seroprevalence 
#'
#' @description Compute the seroprevalence using an object of the class \code{SeroData}. If multiple \code{category} are defined, it will also compute the seroprevalence for each category.
#'
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#' 
#' @param data An object of the class \code{SeroData}. 
#' 
#' @return The mean and 95  \%  confidence interval of the seroprevalence
#'
#' @export
#' @examples
#' 
#' ## A very simple example of a serological survey with three individuals: 
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1), max_age = 50, age_class = 1, sampling_year = 2017)
#' seroprevalence(data = data)
#' 
#' 
seroprevalence <- function(data){
  if(data$Ncategory>1){
    for(i in  levels(factor(data$category))){
      d=subset(data,category = i) # not necessarily an integer
      
      print(paste0("Category ",i))
      B= binom.confint(x=length(which(d$Y==TRUE)),n = d$N,methods = "exact")
      A=paste0('Mean: ', round(B$mean,2), '    2.5%: ',round(B$lower,2), '    97.5%: ', round(B$upper,2))
      print(A)
    }
    
    print("All: ") 
  }
  B= binom.confint(x=length(which(data$Y==TRUE)),n = data$N,methods = "exact")
  A=paste0('Mean: ', round(B$mean,2), '    2.5%: ',round(B$lower,2), '    97.5%: ', round(B$upper,2))
 # print(A)
  
  return(A)
  
}