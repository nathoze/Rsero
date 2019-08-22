##' @title Simulate a serological survey
##' 
##' @description Simulate a serological survey, based on years of epidemics and force of infection.
##' 
##' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
##' 
##' @param number_samples integer. The number of individuals in the simulated serological survey. Default = 500.
##' 
##' @param age_class integer. The length in years of the age classes.  Default = 1.  
##' 
##' @param equal_number boolean. If FALSE, draw random ages with a uniform distribution between 1 and \code{max_age}. if TRUE, generates a SeroData object with an equal distribution of the ages. Default = FALSE.
##' 
##' @param max_age integer. The maximal age. Default = 70.
##' 
##' @param sampling_year integer. Default = 2019. 
##' 
##' @param epidemic_years numeric. A vector with the years where the force of infection is positive. The force of infection at the non specified years is equal to zero. Default = c(2000,1987,1982).  
##' 
##' @param location An optional character factor defining the sampling location. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.
##' 
##' @param sex An optional character factor defining the sex of the individuals. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}. 
##' 
##' @param foi numeric. The force of infection at the years defined in \code{epidemic_years}. It must be the same size as \code{epidemic_years}. Default = c(0.3,0.2,0.4).
##'
##' @param pb numeric, between 0 and 1. Probability for an individual to be seropositive independently of the force of infection (background probability of infection).  Default = 0.
##' 
##' @param rho numeric, 0 or positive. Seroreversion rate.  Default = 0.
##' 
##' @return An object of the class \code{SeroData}.
##' 
##' 
##' @seealso \code{\link{SeroData}} Define the format of the serological data.
##' 
##' @examples 
##' 
##' ## Example 1: Simulates a serological survey taken in a population that experienced
##' ## a series of three outbreaks in 1972, 1988 and 1996. 
##' 
##' data = simulate_SeroData(number_samples = 1000,
##' age_class = 1,
##' sampling_year = 2015,
##' max_age = 50,
##' epidemic_years = c(1996,1988,1972),
##' foi = c(0.2,0.1,0.3),
##' pb = 0)
##' 
##' seroprevalence(data)
##' 
##' 
##' ## Example 2: 500 individuals sampled in 2019, in a population that experienced two epidemics, in 1962
##' ## and 2012.
##' ## Sampled individuals have additionally a probability of being seropositive equal to pb = 0.1. 
##' years =  seq(1962,2012)
##' alpha=c(0.25, 0.1)
##' T=c(1974, 2000)
##' beta=c(1,0.5)
##' FOI = alpha[1]*exp(-(years-T[1])^2/beta[1]^2) + alpha[2]*exp(-(years-T[2])^2/beta[2]^2)
##' data1 <- simulate_SeroData( epidemic_years = years,foi=FOI, pb=0.1)
##' seroprevalence(data1)
##' 
##' ## Adding another survey, sampled in 2005, where ages are categorized in 5 year groups: 
##' data2 <- simulate_SeroData(age_class = 5,
##'  epidemic_years = years,
##'  foi=FOI,
##'  pb=0.1,
##'  sampling_year = 2005)
##' data <- combine_surveys(data1,data2)
##' seroprevalence(data)
##' 
##' 
##' 
##' @export 
simulate_SeroData <- function(number_samples = 500,
                              equal_number = FALSE,
                              age_class = 1,
                              max_age = 70,
                              sampling_year = 2019,
                              epidemic_years =c(2000,1987,1982),
                              foi = c(0.3,0.2,0.4),
                              location = NULL, 
                              sex = NULL, 
                              category = "Category 1",
                              pb=  0,
                              rho=0){
  
  if(length(epidemic_years)!= length(foi) ){
    stop("The number of epidemic years must be the same as the number input foi")
  }
  
  FOI <- rep(0, max_age)  
  
  epidemics_index <- max_age - (sampling_year - epidemic_years) 
  I <- which(epidemics_index>0 & epidemics_index<=max_age)
  FOI[epidemics_index[I]]  =  foi[I] 
  
  if(rho == 0){
    cumulative_FOI <- cumsumfromright(FOI)
    infection_probability <- 1-(1-pb)*exp(-cumulative_FOI)
    
  } else{
    P=rep(0,1,max_age)
    for(J in seq(max_age,1,by=-1)){
      x=rep(0,1,max_age)
      x[J]=1
      for (i in seq(J,1,by =-1)){
        L = rev(FOI)[i]
        x[i-1] = x[i]*exp(-rho-L)+(rho/(L+rho))*(1-exp(-rho-L))
      }
      P[J]=x[1]
    }
    
    infection_probability <- 1-(1-pb)*P#rev(P)
    
  }
  
  # generate random individuals with seropositivity status
  # dataframe with numbers_samples individuals
  if(equal_number == FALSE){ 
    age <- round(runif(number_samples, 1, max_age))
  }else{
    number_samples=number_samples*max_age
    age <- rep(1:max_age, each = number_samples/max_age)
    
  }  
  
  Y <- rbinom(number_samples, 1, infection_probability[age])
  
  sampling_year <- rep(sampling_year, 1, number_samples) 
  
  age <-  ceiling(age_class*floor((age-1)/age_class)+age_class/2)
  
  dat <- SeroData(age_at_sampling =  age,
                  Y = as.logical(as.integer(Y)),
                  sampling_year = sampling_year,
                  max_age = max_age,
                  age_class = age_class,
                  location = location,
                  sex = sex,
                  category =category)
  
  
  return(dat)
}


cumsumfromright <- function(x) cumsum(rev(x))


