% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_posterior.R
\name{plot_posterior}
\alias{plot_posterior}
\title{Plot posterior distribution}
\usage{
plot_posterior(FOIfit)
}
\arguments{
\item{FOIfit}{An object of the class \code{FOIfit}.}
}
\value{
A list with the posterior plots. These plots are generated with the \code{ggplot2} package and can be modified.
}
\description{
Plot the posterior distribution of the parameters estimated in a \code{FOIfit}.
}
\examples{

# Here describe how to access the plots 
P <- plot_posterior(FOIfit)
P[[1]]+ylim(0,200)

}
\author{
Nathanael Hoze \email{nathanael.hoze@gmail.com}
}
