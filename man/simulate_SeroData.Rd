% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_SeroData.R
\name{simulate_SeroData}
\alias{simulate_SeroData}
\title{Simulate a serological survey}
\usage{
simulate_SeroData(
  number_samples = 500,
  equal_number = FALSE,
  age_class = 1,
  max_age = 70,
  sampling_year = 2019,
  epidemic_years = c(2000, 1987, 1982),
  foi = c(0.3, 0.2, 0.4),
  location = NULL,
  sex = NULL,
  category = "Category 1",
  se = 1,
  sp = 1,
  rho = 0
)
}
\arguments{
\item{number_samples}{integer. The number of individuals in the simulated serological survey. Default = 500.}

\item{equal_number}{boolean. If FALSE, draw random ages with a uniform distribution between 1 and \code{max_age}. if TRUE, generates a SeroData object with an equal distribution of the ages. Default = FALSE.}

\item{age_class}{integer. The length in years of the age classes.  Default = 1.}

\item{max_age}{integer. The maximal age. Default = 70.}

\item{sampling_year}{integer. Default = 2019.}

\item{epidemic_years}{numeric. A vector with the years where the force of infection is positive. The force of infection at the non specified years is equal to zero. Default = c(2000,1987,1982).}

\item{foi}{numeric. The force of infection at the years defined in \code{epidemic_years}. It must be the same size as \code{epidemic_years}. Default = c(0.3,0.2,0.4).}

\item{location}{An optional character factor defining the sampling location. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.}

\item{sex}{An optional character factor defining the sex of the individuals. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.}

\item{se}{numeric, between 0 and 1. Sensitivity of the assay.  Default = 1.}

\item{sp}{numeric, between 0 and 1. Specificity of the assay.  Default = 1.}

\item{rho}{numeric, 0 or positive. Seroreversion rate.  Default = 0.}
}
\value{
An object of the class \code{SeroData}.
}
\description{
Simulate a serological survey, based on years of epidemics and force of infection.
}
\examples{

## Example 1: Simulates a serological survey taken in a population that experienced
## a series of three outbreaks in 1972, 1988 and 1996. 

data = simulate_SeroData(number_samples = 1000,
age_class = 1,
sampling_year = 2015,
max_age = 50,
epidemic_years = c(1996,1988,1972),
foi = c(0.2,0.1,0.3),
se = 1)

seroprevalence(data)


## Example 2: 500 individuals sampled in 2019, in a population that experienced two epidemics, in 1962
## and 2012.
## The assay has a sensitivity of 0.9. 
years =  seq(1962,2012)
alpha=c(0.25, 0.1)
T=c(1974, 2000)
beta=c(1,0.5)
FOI = alpha[1]*exp(-(years-T[1])^2/beta[1]^2) + alpha[2]*exp(-(years-T[2])^2/beta[2]^2)
data1 <- simulate_SeroData( epidemic_years = years,foi=FOI, se=0.9)
seroprevalence(data1)

## Adding another survey, sampled in 2005, where ages are categorized in 5 year groups: 
data2 <- simulate_SeroData(age_class = 5,
 epidemic_years = years,
 foi=FOI,
 se = 0.9,
 sampling_year = 2005)
data <- combine_surveys(data1,data2)
seroprevalence(data)



}
\seealso{
\code{\link{SeroData}} Define the format of the serological data.
}
\author{
Nathanael Hoze \email{nathanael.hoze@gmail.com}
}
