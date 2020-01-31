#' Serological survey of Chikungunya in Vietnam, 2015
#'
#' These data describe a serological survey conducted on 546 individuals in four regions of Vietnam
#'
#' @docType data
#'
#' @format {
#' An object of the class \code{SeroData} containing
#' \describe{
#'   \item{age}{age of individuals}
#'   \item{age_at_sampling}{age of individuals}
#'   \item{Y}{Seropositive/Seronegative}
#'   \item{sampling_year}{2015}
#'   \item{location}{One of the four regions}
#'   \item{sex}{none}
#'   \item{category}{One of the four regions}
#' }
#'
#' }
#' @rdname chikv_vietnam
#'
#' @author Data from Quan et al., Evidence of previous but not current transmission of chikungunya virus in southern and central Vietnam: Results from a systematic review and a seroprevalence study in four locations, Plos NTD(2018), under CC BY licence.
#'
#' @source Quan et al. (2018)
#'
#' @references
#' Quan, T. M., Phuong, H. T., Vy, N. H. T., Le Thanh, N. T., Lien, N. T. N., Hong, T. T. K., ... & Clapham, H. E. (2018). 
#' Evidence of previous but not current transmission of chikungunya virus in southern and central Vietnam: Results from a systematic review and a seroprevalence study in four locations. PLoS neglected tropical diseases, 12(2), e0006246.
#'
#'
#' @examples
#' chikv_vietnam
#'

sero.dta =  read.csv('CHIKVVietnam.csv', header = TRUE, sep=';')


chikv_vietnam = Rsero::SeroData(age_at_sampling = round(as.numeric(sero.dta$age)), 
               Y = as.numeric(sero.dta$Y)>9,
               sampling_year = as.numeric(sero.dta$sampling_year), 
               category = as.character(sero.dta$Loc),
               location = as.character(sero.dta$Loc))

usethis::use_data(chikv_vietnam, overwrite = TRUE)
