#' @title Compute the seroprevalence 
#'
#' @description Compute the seroprevalence using an object of the class \code{SeroData}. If multiple \code{category} are defined, it will also compute the seroprevalence for each category.
#'
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#' 
#' @param serodata An object of the class \code{SeroData}. 
#' 
#' @return The mean and 95  \%  confidence interval of the seroprevalence
#'
#' @export
#' @examples
#' 
#' ## A very simple example of a serological survey with three individuals: 
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1), max_age = 50, age_class = 1, sampling_year = 2017)
#' seroprevalence(serodata = data)
#' 
 


seroprevalence <- function(serodata){
  if(serodata$Ncategory>1){
    for(i in  1:serodata$Ncategoryclass){
      for(j in  unique(serodata$category[,i])){
        
        d=subset(serodata,sub  = which(serodata$category[,i]==j)) # not necessarily an integer
        
        print(paste0("Category ",j))
        B= binom.confint(x=length(which(d$Y==TRUE)),n = d$N,methods = "exact")
        A=paste0('Mean: ', round(B$mean,2), '    2.5%: ',round(B$lower,2), '    97.5%: ', round(B$upper,2))
        print(A)
      }
    }  
      print("All: ") 
  }
  B= binom.confint(x=length(which(serodata$Y==TRUE)),n = serodata$N,methods = "exact")
  A=paste0('Mean: ', round(B$mean,2), '    2.5%: ',round(B$lower,2), '    97.5%: ', round(B$upper,2))
  # print(A)
  
  return(A)
  
}


