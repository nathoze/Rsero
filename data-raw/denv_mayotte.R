#' Serological survey of Dengue virus in Mayotte Island, 2006
#'
#' These data describe a serological survey conducted in Mayotte in 2006.
#'
#' @docType data
#'
#' @format {
#' A \code{SeroData} object
#'
#' }
#' @rdname denv_mayotte
#'
#' @author Data from Sissoko, Daouda, et al. "Seroepidemiology of dengue virus in Mayotte, Indian Ocean, 2006." PLoS One 5.11 (2010): e14141.
#' 
#' 
#' @source Sissoko et al. (2010)
#'
#' @references
#' Sissoko, Daouda, et al. "Seroepidemiology of dengue virus in Mayotte, Indian Ocean, 2006." PLoS One 5.11 (2010): e14141.
#' 
#' @examples
#' denv_mayotte
#'

transform_data_from_histogram_uniform <- function(n, sero, age.range){
  
  v = !is.na(sero)
  positive <- round(n[v]*sero[v]/100)
  negative <- n[v]-positive
  
  S <- rep(rep(c(TRUE, FALSE), length(n[v])), c(rbind(positive, negative))) # create a list which the seropositivity TRUE/FALSE of each individual
  # random uniform ages
  
  age.at.sampling <- c()
  for( i in 1:dim(age.range)[2]){
    #age <- round(runif(n = n[i],min = age.range[1,i], max=age.range[2,i]) )
    age <- round(runif(n = n[i],min = age.range[1,i], max=age.range[1,i]) )
    age.at.sampling <- c(age.at.sampling, age)
    
  } 
  list(S, age.at.sampling) 
}

#age.range <- matrix(c(2,14,15,24, 25,34,35,44,45,54, 55, 60), nrow=2) # age groups (2-15; 15-24; ...)
n <- c(0,302,294,193,169,107,89) # number of individuals in each age group
ninfected <- c(0,93,72,52,56,53,46,41) # number of seropositive individuals in each age group

sero <-c(0,2.2,19.1, 34,35.8,38.8,28.4)

age.range <- matrix(c(1,1, 2,14,15,24, 25,34,35,44,45,54, 55, 60), nrow=2) # age groups (2-15; 15-24; ...)

# define the age groups
class1 = c()
class2 = c()
for(i in 1:7){
  class1=c(class1, rep(age.range[1,i],age.range[2,i]-age.range[1,i]+1))
  class2=c(class2, rep(age.range[2,i],age.range[2,i]-age.range[1,i]+1))
}



transformed_data <- transform_data_from_histogram_uniform(n, sero, age.range = age.range)

denv_mayotte <- SeroData(age_at_sampling =  transformed_data[[2]],
                 Y = transformed_data[[1]],
                 max_age =  60,
                 age_cats = 1,
                 sampling_year = 2006,
                 location = 'Mayotte',
                 class1 =  class1,
                 class2 = class2)

usethis::use_data(denv_mayotte, overwrite = TRUE)
