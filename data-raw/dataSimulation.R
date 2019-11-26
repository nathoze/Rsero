one_peak_simulation = Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2)
usethis::use_data(one_peak_simulation, overwrite = TRUE)

two_peaks_simulation =  Rsero::simulate_SeroData(epidemic_years = c(1980,2000), foi =c(0.3,0.2) ) 
usethis::use_data(two_peaks_simulation, overwrite = TRUE)

intervention_seroreversion_simulation =  Rsero::simulate_SeroData(epidemic_years = seq(1940,2009), foi = c(rep(0.1,50),rep(0.01,20)),rho = 0.01)
usethis::use_data(intervention_seroreversion_simulation, overwrite = TRUE)

constant_seroreversion_simulation =  Rsero::simulate_SeroData(epidemic_years = seq(1940,2009), foi =rep(0.1,70),rho = 0.01)
usethis::use_data(constant_seroreversion_simulation, overwrite = TRUE)

d1 =  Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2, rho=0.02,sampling_year = 2019)
d2 =  Rsero::simulate_SeroData(epidemic_years = 2000, foi =0.2, rho=0.02, sampling_year = 2005)
intervention_seroreversion_two_surveys_simulation =  Rsero::combine_surveys(d1,d2)
usethis::use_data(intervention_seroreversion_two_surveys_simulation, overwrite = TRUE)

two_peaks_background_simulation =  Rsero::simulate_SeroData(epidemic_years = c(1980,2000), foi =c(0.2,0.1),pb = 0.05)
usethis::use_data(two_peaks_background_simulation, overwrite = TRUE)
