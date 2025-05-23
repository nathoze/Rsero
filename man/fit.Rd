% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{fit}
\alias{fit}
\alias{print.FOIfit}
\alias{summary.FOIfit}
\title{Fit force of infection models to serological data}
\usage{
fit(
  model,
  data = data,
  iter = 5000,
  chains = 4,
  warmup = floor(iter/2),
  thin = 1,
  seed = sample.int(.Machine$integer.max, 1),
  init = "random",
  check_data = TRUE,
  sample_file = NULL,
  diagnostic_file = NULL,
  verbose = FALSE,
  algorithm = c("NUTS", "HMC", "Fixed_param"),
  control = NULL,
  include = TRUE,
  cores = getOption("mc.cores", 1L),
  open_progress = interactive() && !isatty(stdout()) && !identical(Sys.getenv("RSTUDIO"),
    "1"),
  show_messages = TRUE,
  ...
)

\method{print}{FOIfit}(x, digits_summary = 2, ...)

\method{summary}{FOIfit}(x, ...)
}
\arguments{
\item{model}{A \code{FOImodel} object, generated by the function \code{\link{FOImodel}}.}

\item{data}{A \code{SeroData} object, generated by the function 
\code{\link{SeroData}}.}

\item{iter}{integer. Number of iterations for each chain, including warmup. Default = 5000.}

\item{chains}{integer. Number of independent Markov chains. Default = 4.}

\item{warmup}{A positive integer specifying the number of warmup (aka burnin) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup  samples should not be used for inference). The number of warmup iterations should not be larger than iter and the default is iter/2.}

\item{thin}{A positive integer specifying the period for saving samples. 
The default is 1, which is usually the recommended value.}

\item{seed}{The seed for random number generation. The default is generated from 1 to the maximum integer supported by \R on the machine. Even if 
multiple chains are used, only one seed is needed, with other chains having seeds derived from that of the first chain to avoid dependent samples.
When a seed is specified by a number, \code{as.integer} will be applied to it. 
If \code{as.integer} produces \code{NA}, the seed is generated randomly. The seed can also be specified as a character string of digits, such as
\code{"12345"}, which is converted to integer.}

\item{init}{Initial values specification. See the detailed documentation for the init argument in \code{\link{stan}}.}

\item{check_data}{Logical, defaulting to \code{TRUE}. If \code{TRUE} 
the data will be preprocessed; otherwise not.  See the Note section in \code{\link{stan}}.}

\item{sample_file}{An optional character string providing the name of a file. If specified the draws for \emph{all} parameters and other saved quantities
will be written to the file. If not provided, files are not created.   When the folder specified is not writable, \code{tempdir()} is used.   When there are multiple chains, an underscore and chain number are appended
to the file name prior to the \code{.csv} extension.}

\item{diagnostic_file}{An optional character string providing the name of a file. If specified the diagnostics data for \emph{all} parameters will be written
to the file. If not provided, files are not created. When the folder specified is not writable, \code{tempdir()} is used. When there are multiple chains, 
an underscore and chain number are appended to the file name prior to the \code{.csv} extension.}

\item{verbose}{\code{TRUE} or \code{FALSE}: flag indicating whether to print intermediate output from Stan on the console, which might  be helpful for model debugging.}

\item{algorithm}{One of sampling algorithms that are implemented in Stan. Current options are \code{"NUTS"} (No-U-Turn sampler, Hoffman and Gelman 2011, Betancourt 2017),   \code{"HMC"} (static HMC), or \code{"Fixed_param"}. The default and preferred algorithm is \code{"NUTS"}.}

\item{control}{A named \code{list} of parameters to control the sampler's behavior. See the details in the documentation for the \code{control} argument    in \code{\link{stan}}.}

\item{include}{Logical scalar defaulting to \code{TRUE} indicating whether to include or exclude the parameters given by the 
\code{pars} argument. If \code{FALSE}, only entire multidimensional parameters can be excluded, rather than particular elements of them.}

\item{cores}{Number of cores to use when executing the chains in parallel, which defaults to 1 but we recommend setting the \code{mc.cores} option 
to be as many processors as the hardware and RAM allow (up to the     number of chains).}

\item{open_progress}{Logical scalar that only takes effect if 
\code{cores > 1} but is recommended to be \code{TRUE} in interactive
use so that the progress of the chains will be redirected to a file
that is automatically opened for inspection. For very short runs, the
user might prefer \code{FALSE}.}

\item{show_messages}{Either a logical scalar (defaulting to \code{TRUE})
indicating whether to print the summary of Informational Messages to
the screen after a chain is finished or a character string naming a path
where the summary is stored. Setting to \code{FALSE} is not recommended
unless you are very sure that the model is correct up to numerical 
error.}
}
\value{
A list with the class \code{FOIfit}, which contains the
following items:

\itemize{

\item fit: The results of the fit, of class \code{stanfit}.

\item data: The input data.

\item model: The input model.
}
}
\description{
Runs a Bayesian MCMC using rstan to draw samples of a model of the class \code{FOImodel}, with data of the class \code{SeroData}.
 It is an adaptation of the function \code{\link[rstan]{sampling}} in the rstan package.
}
\examples{

data <- simulate_SeroData(number_samples = 1000,
  age_class = 1,
  epidemic_years = c(1976,1992),
  foi = c(0.2,0.3))
model <- FOImodel('outbreak', K = 2)
F1 <- fit(model = model, data = data)
seroprevalence.fit(F1)
plot(F1)

}
\seealso{
\code{\link{SeroData}} Define the format of the serological data.

\code{\link{FOImodel}} Define a model.
}
\author{
Nathanael Hoze \email{nathanael.hoze@gmail.com}
}
