one_peak_simulation = Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2)
devtools::use_data(one_peak_simulation)

two_peaks_simulation =  Rsero::simulate_SeroData(epidemic_years = c(1980,2000), foi =c(0.3,0.2) ) 
devtools::use_data(two_peaks_simulation)

intervention_seroreversion_simulation =  Rsero::simulate_SeroData(epidemic_years = seq(1940,2009), foi = c(rep(0.1,50),rep(0.01,20)),rho = 0.01)
devtools::use_data(intervention_seroreversion_simulation)

constant_seroreversion_simulation =  Rsero::simulate_SeroData(epidemic_years = seq(1940,2009), foi =rep(0.1,70),rho = 0.01)
devtools::use_data(constant_seroreversion_simulation)

d1 =  Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2, rho=0.02,sampling_year = 2019)
d2 =  Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2, rho=0.02, sampling_year = 2005)
intervention_seroreversion_two_surveys_simulation =  Rsero::combine_surveys(d1,d2)
devtools::use_data(intervention_seroreversion_two_surveys_simulation)

two_peaks_background_simulation =  Rsero::simulate_SeroData(epidemic_years = c(1980,2000), foi =c(0.2,0.1),pb = 0.05)
devtools::use_data(two_peaks_background_simulation)
