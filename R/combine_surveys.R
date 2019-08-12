#' @title Combine serological surveys
#' 
#' @description The function \code{combine_surveys} creates a new object of the class \code{SeroData} from two \code{SeroData} objects. In particular, it allows combining serological surveys sampled at different dates.
#'
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#'  
#' @param SeroData1	A first \code{SeroData} object.
#'  
#' @param SeroData2	A second \code{SeroData} object.
#' 
#' @return A \code{SeroData} object.
#'   
#' @seealso \code{\link{SeroData}} Define the format of the serological data.
#' @seealso \code{\link{subset}} Extract a smaller subset of a  \code{SeroData} object.
#' 
#'
#' @examples 
#'  # Combine two simulated surveys, one acquired in 2015, with one-year age categories, 
#'  # and another one acquired in 1995 with 10-year age categories. 
#'  Years <- c(1976,1992)
#'  FOI <- c(0.2,0.3)
#'  data1 <- simulate_SeroData(sampling_year = 2015, epidemic_years = Years, foi = FOI)
#'  data2 <- simulate_SeroData(sampling_year = 1995,age_class=10, epidemic_years = Years, foi = FOI)     
#'  data <- combine_surveys(data1,data2)             
#' @export

combine_surveys <- function(SeroData1,SeroData2){
  
  N = SeroData1$N+SeroData2$N
  y1 =max(SeroData1$sampling_year)
  y2= max(SeroData2$sampling_year)
  
  z1 =min(SeroData1$sampling_year)
  z2= min(SeroData2$sampling_year)
  
  Year1 = max(y1,y2)
  #Year2 = min(z1,z2)
  Year2 = min(y1,y2)
  
  if (y1 > y2){
    Offset = Year1-Year2
    
    dat1 = SeroData1
    dat2 = SeroData2
    
  }else{
    Offset = Year1-Year2
    dat1 = SeroData2
    dat2 = SeroData1
  }
  A = max(dat1$A, dat2$A+Offset)
  
  # u1=unique(dat1$sampling_year)
  # u2=unique(dat2$sampling_year)
  # 
  # if( is.na(match(u1,u2)) ) {
  #   
  #   NAgeGroups = dat1$NAgeGroups+dat2$NAgeGroups -length(which(!is.na(match(u1,u2) ) ) ) 
  #   
  # }
  
  
 # Indices =  combine_ind_by_age(dat1,dat2)  # 7/9/2018
      
  
  
  age2  = dat2$age_at_sampling+Offset #?
  
  sampling_year = c(dat1$sampling_year,dat2$sampling_year)

  a = max(sampling_year)-sampling_year +1
  
  
  age_group = rep(0,N)
   u  = unique(a) 
   for(i in 1:length(u)){
    age_group[which(u[i] == a) ] =i
   } 
    age_at_init =u

    NAgeGroups = length(u)
  
  
  data <- list( A = A,#length(Indices),
                NGroups = A,# max(Indices),
                N = N,
                Y = c(dat1$Y,dat2$Y),
                age = c(dat1$age,age2),
                age_at_sampling = c(dat1$age_at_sampling,dat2$age_at_sampling),
                #  New 06/09/2018    ind_by_age = Indices,
                sampling_year = sampling_year,
                location = c(dat1$location,dat2$location),
                sex = c(dat1$sex,dat2$sex),
                category = c(dat1$category,dat2$category),
                Ncategory = length(unique(c(dat1$category,dat2$category))), 
                #Ngroups = Alength(Indices), #  New 06/09/2018
                NAgeGroups = NAgeGroups,
                age_at_init =  as.array(age_at_init), 
                age_group  = age_group
                )
  
  
  class(data) <- 'SeroData'
  return(data)
}


combine_ind_by_age <- function(data1,data2){
  
  # 
  # y1 =max(data1$sampling_year)
  # y2= max(data2$sampling_year)
  # 
  # Year1 = max(y1,y2)
  # #Year2 = min(z1,z2)
  # Year2 = min(y1,y2)
  # 
  # if (y1 > y2){
  #   Offset = Year1-Year2
  # 
  #   
  # }else{
  #   Offset = Year2-Year1 #??
  # 
  # }
  # 

  
  
  y1 =max(data1$sampling_year)
  y2= max(data2$sampling_year)
  
  z1 =min(data1$sampling_year)
  z2= min(data2$sampling_year)
  
  
  
  
  Year1=max(y1,y2)
  Year2=max(z1,z2)
  
 # Year1 = max(data1$sampling_year,data2$sampling_year)
  #Year2 = min(data1$sampling_year,data2$sampling_year)

  if (Year1 == Year2){
    Offset = Year1-Year2

  }else{
    Offset = Year1-Year2
  }
  
  if(data1$A-Offset<data2$A){
    L1 <- data2$A+Offset
  }else{
    L1 <- data1$A
  }
  IND1 <- rep(0,L1)
  IND1[1:data1$A] = data1$ind_by_age
  IND2 <- rep(0,L1)
  IND2[Offset+1:L1] = data2$ind_by_age
  
  IND <- rep(0,L1)
  startIndex=0
  element=0
  for(i in 1:length(IND)){
    if (IND2[i] == 0){
      IND[i] = IND1[i]
      startIndex = i
      element = IND[i]
    }
  }
  m1=0
  m2=0
  while ( is.finite(m1) &&  is.finite(m2)){
    startIndex=startIndex+1
    ind1 = IND1[startIndex]
    ind2 = IND2[startIndex]
    
#    print(IND1[startIndex:length(IND1)])
    #print(IND)
    
    m1 = min(which(IND1[startIndex:length(IND1)]  > ind1))
    m2 = min(which(IND2[startIndex:length(IND2)]  > ind2))
    if  ( is.finite(m1) &&  is.finite(m2)){
      endIndex =  startIndex+max(m1,m2)-2
      IND[startIndex:endIndex ] = element+1
          element  = element+1
    startIndex = endIndex
    }else{
      IND[startIndex:length(IND)]=IND2[startIndex:length(IND)] + Offset
    }
    

  }
  return(IND)
  
} 