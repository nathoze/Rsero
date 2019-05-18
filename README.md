<!-- README.md is generated from README.Rmd. Please edit that file -->
Rsero: Estimate the annual force of infection using serological data
====================================================================

Rsero is R package dedicated to the implementation of serocatalytic models that are used to estimate the force of infection from age-stratified serological surveys.

Estimations requires:
*The age of each individual *Their seropositivity status
*The year of sampling *A model of pathogen circulation

The package provides a standardized framework to store serological data, analyze serological surveys, use a variety of serocatalytic models, run MCMC algorithm to estimate the parameters of the force of infection, and analyse the results.

Installation
------------

To install the package, type

``` r
devtools::install_github("nathoze/Rsero")
```

(This step requires that the package *devtools* is installed). Installation may take a few minutes due to the compilation of the stan files that encode the serocatalytic models.

MCMC simulations are done using the package rstan, that requires Rtools on Windows computers. Rtools can be found here <https://cran.r-project.org/bin/windows/Rtools/>

More details
------------

More details on the available serocatalytic models can be found in the vignette *models*. This vignette includes several examples and shows the different features of the packages

Worked example.
---------------

In this example we show the basic steps to follow to get a complete analysis of a serological survey. This analysis uses a simulated dataset containing samples of 500 individuals.

``` r
library(Rsero)
data('one_peak_simulation')
```

The data is saved under a custom format *SeroData*, that stores information about the serological surveys and details on the participants. Basic information such as the seroprevalence can be obtained simply by typing

``` r
seroprevalence(one_peak_simulation)
```

and a graph of the age profile of seroprevalence is obtained using

``` r
seroprevalence.plot(one_peak_simulation)
```

Age-stratified serological surveys assess immunological markers of past infections and can be used to reconstruct the historical patterns of the circulation of an infectious disease and the force of infection (the per capita rate at which susceptible individuals will be infected a given year). Such inference is performed using serocatalytic models. The Rsero package proposes several models that can be fitted to the data using the command

``` r
FOIfit = fit( data = one_peak_simulation,  model = seromodel)
```

In its basic form a model is defined by specifying the mode of circulation of the pathogen. The simplest type of model assumes a constant force of infection. We define it with the command

``` r
ConstantModel = FOImodel(type = 'constant')
```

We can now fit the defined model to the data:

`r  FOIfit.constant = fit( data = one_peak_simulation,  model = Constantmodel, chains=1)` (we simulate only one MCMC chain here). We now visualize the result of the fit

`r  seroprevalence.fit(FOIfit.constant, YLIM=0.5)` Here the solid line is the mean annual FOI obtained from the MCMC simulations and the envelope is the 95% credible interval. This seems not a very good fit. Indeed younger individuals are all seronegative, which suggests that the pathogen did not circulate in the recent years. Several other models could explain the data. For instance we define a model of one outbreak.

`r  OutbreakModel = FOImodel( type='outbreak', K=1)` We can now fit the defined model to the data:

``` r
FOIfit.outbreak = fit( data = one_peak_simulation,  model = OutbreakModel, chains=1)
```

and we plot the result

`r  seroprevalence.fit(FOIfit.outbreak)` Visually, this a better fit (and it is not surprising given the way the data was generated). We can compare the results of the fit using the deviance information criterion (DIC).

``` r
DIC.constant = compute_information_criteria(FOIfit.constant)
DIC.outbreak = compute_information_criteria(FOIfit.outbreak)
```

The DIC obtained for the outbreak model is lower than that for the constant model. This indicates a better fit of the outbreak model to the data.

More information
----------------

For more details on the models, go to the vignette model.
