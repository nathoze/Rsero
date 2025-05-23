% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeroData.R
\name{SeroData}
\alias{SeroData}
\alias{compute.age.groups}
\title{Create a dataset containing the serological data and other information}
\usage{
SeroData(
  age_at_sampling,
  age = NULL,
  Y,
  age_class = 1,
  max_age = NULL,
  sampling_year = NULL,
  location = NULL,
  sex = NULL,
  category = "Category 1",
  reference.category = NULL,
  class1 = NULL,
  class2 = NULL,
  ...
)

compute.age.groups(age, sampling_year)
}
\arguments{
\item{age_at_sampling}{A vector of integers containing the age of the sampled individuals at the time of sampling. Must be equal or greater than 1.}

\item{Y}{A vector containing the seropositivy status of the sampled individuals. It can be in a numeric form (1 or 0) or boolean (\code{TRUE} or \code{FALSE}). This vector must have the same size as \code{age}.}

\item{age_class}{Integer. The length in years of the age classes. Default = 1.}

\item{max_age}{Integer. The maximal age considered for the individuals in the sample. Individuals older than \code{max_age} are set equal to \code{max_age}.}

\item{sampling_year}{Integer. Defines the sampling year.  It can be a single value or a vector of the same size as the number of sampled individuals. If it is a single value, all the sampled individuals have the same sampling year. Default = 2017.}

\item{location}{An optional character factor defining the sampling location. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.}

\item{sex}{An optional character factor defining the sex of the individuals. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.}

\item{category}{Character. An optional element containing the name of the categories and defining the category of the individuals. 
This feature is used when fitting the models assuming  different risks of infection for the different categories. It can be a single character element or a matrix of characters with row size equal to the number of individuals and column size equal to the number of different category classes. Default = "Category 1".}

\item{reference.category}{Character. A vector containing the name of the reference categories. It must have the same length as the number of columns of \code{category}. By default, it will take as reference the most common element in each category.}

\item{...}{Additional arguments (not used).}
}
\value{
A list with the class \code{SeroData}, which contains the
following items:

\itemize{

\item age: A vector of integers containing the age of the sampled individuals at the time of the latest sampling, of length N.

\item age_at_sampling: A vector of integers containing the age of the sampled individuals at the time of sampling, of length N.

\item Y: A vector of seropositivity status, of length N.

\item N: The number of individuals considered.

\item A: The maximal age.

\item NGroups: The total number of age groups.

\item sampling_year: A vector of sampling years, of length N.

\item location: A vector of the sampling location, of length N.

\item sex: A vector of the sex of each individual, of length N.

}
}
\description{
This function contains the definition of the class \code{SeroData}. It requires the age and seropositivity of individuals, and accepts more input parameters.
}
\examples{

## A very simple example of a serological survey with three individuals: 
data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1), max_age = 50, age_class = 1, sampling_year = 2017)

seroprevalence(data = data)
data2 = subset(data,c(1,3))
seroprevalence(data = data2, age_class = 5) 

## A example using categories 
sex= c('male', 'male', 'female')
data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1),  category= sex )

# defining the reference category 'female'

data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1),  category= sex, reference.category='female')

# Grouping individuals by age categories 1-10, 11-20, 21-30, etc. 
data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1),  category= sex, reference.category='female', age_class = 10)
}
\author{
Nathanael Hoze \email{nathanael.hoze@gmail.com}
}
