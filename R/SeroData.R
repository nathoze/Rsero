#' @title Create a dataset containing the serological data and other information
#'
#' @description This function contains the definition of the class \code{SeroData}. It requires the age and seropositivity of individuals, and accepts more input parameters. 
#'
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#' 
#' @param age_at_sampling A vector of integers containing the age of the sampled individuals at the time of sampling. Must be equal or greater than 1.  
#' 
#' @param Y A vector containing the seropositivy status of the sampled individuals. It can be in a numeric form (1 or 0) or boolean (\code{TRUE} or \code{FALSE}). This vector must have the same size as \code{age}.
#'  
#' @param age_class Integer. The length in years of the age classes. Default = 1. 
#' 
#' @param max_age Integer. The maximal age considered for the individuals in the sample. Individuals older than \code{max_age} are set equal to \code{max_age}. 
#' 
#' @param sampling_year Integer. Defines the sampling year.  It can be a single value or a vector of the same size as the number of sampled individuals. If it is a single value, all the sampled individuals have the same sampling year. Default = 2017. 
#' 
#' @param location An optional character factor defining the sampling location. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}.
#' 
#' @param sex An optional character factor defining the sex of the individuals. It can be a single chain of characters or a vector of the same size as the number of sampled individuals. Default = \code{NULL}. 
#'
#' @param category Character. An optional element containing the name of the categories and defining the category of the individuals. 
#' This feature is used when fitting the models assuming  different risks of infection for the different categories. It can be a single character element or a matrix of characters with row size equal to the number of individuals and column size equal to the number of different category classes. Default = "Category 1".
#'
#' @param reference.category Character. A vector containing the name of the reference categories. It must have the same length as the number of columns of \code{category}. By default, it will take as reference the most common element in each category.
#' 
#' @param ... Additional arguments (not used).
#'  
#' 
#' @return A list with the class \code{SeroData}, which contains the
#' following items:
#'
#' \itemize{
#'
#' \item age: A vector of integers containing the age of the sampled individuals at the time of the latest sampling, of length N.
#' 
#' \item age_at_sampling: A vector of integers containing the age of the sampled individuals at the time of sampling, of length N.
#' 
#' \item Y: A vector of seropositivity status, of length N.
#' 
#' \item N: The number of individuals considered.
#'
#' \item A: The maximal age.
#'
#' \item NGroups: The total number of age groups.
#'
#' \item sampling_year: A vector of sampling years, of length N.
#'
#' \item location: A vector of the sampling location, of length N.
#'
#' \item sex: A vector of the sex of each individual, of length N.
#'
#' }
#'
#' @export
#' @examples
#' 
#' ## A very simple example of a serological survey with three individuals: 
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1), max_age = 50, age_class = 1, sampling_year = 2017)
#' 
#' seroprevalence(data = data)
#' data2 = subset(data,c(1,3))
#' seroprevalence(data = data2, age_class = 5) 
#' 
#' ## A example using categories 
#' sex= c('male', 'male', 'female')
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1),  category= sex )
#' 
#' # defining the reference category 'female'
#'
#' data = SeroData(age_at_sampling = c(10,32,24), Y=c(0,1,1),  category= sex, reference.category='female')


SeroData <- function(age_at_sampling,
                     age=NULL,
                     Y,         
                     age_class = 1,
                     max_age = NULL,
                     sampling_year = NULL,
                     location = NULL,
                     sex = NULL,
                     category = "Category 1",
                     reference.category = NULL, ## add in documentation  +  define default
                     class1 = NULL,
                     class2 = NULL,
                     ...){
  # Error Messages
  #add  Check : Y and age_at_sampling must have the same length also when Y is multidimensional
  if(dim(as.matrix(Y))[1] != dim(as.matrix(age_at_sampling))[1] ){
    stop("Error : Y and age_at_sampling must have the same length") 
  }
  
  if(sum(is.na(Y))){
    stop("Error : remove NA from the seropositivity Y") 
  }
  if(sum(is.na(age_at_sampling))){
    stop("Error : remove NA from the age age_at_sampling") 
    
  }
  
  
  
  if(!testInteger(age_at_sampling)){ 
    stop("age_at_sampling must be given as integer")
  }
  
  if(min(age_at_sampling)<=0){ 
    stop("age_at_sampling must have only positive numbers (>=1)")
  }
  
  if(!typeof(category)=="character"){
    stop("'category' must contain characters")
  }
  
  if (is.null(sampling_year)){
    sampling_year <- rep(2017,1,length(age_at_sampling))
  }else if(length(sampling_year) == 1){
    sampling_year <- rep(sampling_year,1,length(age_at_sampling))
  }
  
  if (is.null(location)){
    location <- rep('NA',1,length(age_at_sampling))
  }else if(length(location) == 1){
    location <- rep(location,1,length(age_at_sampling))
  }
  
  if (is.null(sex)){
    sex <- rep("All", 1, length(age_at_sampling))
  }else if(length(sex) == 1){
    sex <- rep(sex,1,length(age_at_sampling))
  }
  
  if(length(category) == 1){
    category <- as.matrix(rep(category,1,length(age_at_sampling)))
  }
  
  if(is.null(age)){
    age = age_at_sampling+max(sampling_year)-sampling_year
  } 
  
  if (is.null(max_age)){
    max_age <- max(age)
  }  
  
  if(is.null(class1) & is.null(class2) & age_class==1){
 
    class1 = as.matrix(1:max_age)
    class2= class1

  }
  
  if(is.null(class1) & is.null(class2)){
    class1 = as.matrix(rep(seq(1,max_age,age_class), each = age_class) )
    class2 = as.matrix(age_class-1 + rep(seq(1,max_age,age_class), each = age_class) ) 
  }
   
  
  age[which(age>max_age)] <- max_age
  
  age_at_sampling[which(age_at_sampling>max_age)] <- max_age 
  age.groups <-compute.age.groups(age = age,sampling_year = sampling_year)
  
  
  N=length(age)
  
  param.category =  category.parameters(category=category, N=N, reference.category=reference.category)
  
  
  
  
  data <- list( A = max_age,
                N = N,
                Y = Y,
                age = age,
                age_at_sampling = age_at_sampling,
                sampling_year = sampling_year,
                location = location,
                sex = sex,
                category = param.category$category,
                categoryindex=param.category$categoryindex,
                MatrixCategory = param.category$MatrixCategory,
                Ncategory = param.category$Ncategory,
                maxNcategory=param.category$maxNcategory,
                Ncategoryclass=param.category$Ncategoryclass,
                unique.categories=param.category$unique.categories,
                Ncat.unique = param.category$Ncat.unique,
                category.position.in.table=param.category$category.position.in.table,
                reference.category=param.category$reference.category,
                NGroups = max_age, 
                NAgeGroups = age.groups$NAgeGroups, 
                age_at_init =  as.array(age.groups$age_at_init), # CHECK THIS 25/09/2018 
                age_group  = age.groups$age_group,
                class1  = class1,
                class2 = class2)
  
  
  class(data) <- 'SeroData'
  
  return(data)
}


#' @rdname SeroData
#' @export
compute.age.groups <- function(age,sampling_year){
  
  N=length(age)
  a = max(sampling_year)-sampling_year +1
  age_group = rep(0,N)
  u  = unique(a) 
  for(i in 1:length(u)){
    age_group[which(u[i] == a) ] =i
  } 
  age_at_init =u
  NAgeGroups = length(u)
  
  output <- list(age_at_init=age_at_init, NAgeGroups=NAgeGroups, age_group=age_group)
  return(output)
}



#' @export
print.SeroData <- function(data,...){
  
  
  df <-data.frame(age = data$age,
                  age_at_sampling = data$age_at_sampling,
                  Y=data$Y,
                  sampling_year = data$sampling_year,
                  location = data$location,
                  sex = data$sex,
                  category  = data$category )
  cat("<SeroData object>\n")
  cat(sprintf("%i serological samples ", data$N),'\n')
  cat(sprintf("Max age: %i ", data$A),'\n')
  print(df)
  
}


find.category <- function(x, categories){
  return(tail(which(x-categories>0), 1L))
}

testInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}

#' @export
category.parameters <- function(category,N, reference.category){
  
  category =  as.matrix(category)
  
  Ncategoryclass =  dim(category)[2]
  
  if(is.null(reference.category)){
    reference.category = c()
    for(i in 1:Ncategoryclass){
      most.common=  names(sort(summary(as.factor(category[,i])), decreasing=T)[1])
      reference.category=c(reference.category, most.common )  # add the most common element as reference
    }
  }
  
  if(Ncategoryclass != length(reference.category)){
    stop("'reference.category' and 'category' don't have the same number of elements ")
  }
  
  # # if reference.category doesn't exist in the list 
  # 
  # for(i in 1:Ncategoryclass){
  #   if(sum(category[,i]==reference.category[i])==0) {
  #     print(paste0('The reference category ', reference.category[i], ' was not found in the variable category'))
  #   }
  # }
  
  
  
  #Verfier ICI pour plusieurs dimensions
  for(i in 1:Ncategoryclass){
    if(sum(category[,i]==reference.category[i])==0) {
      most.common=  names(sort(summary(as.factor(category[,i])), decreasing=T)[1])
      reference.category[i] =   most.common    # add the most common element as reference
    }
  }
  
  
  A=apply(category, 2, unique)
  if(typeof(A)=="character"){ # needed in the case there is only one type of category
    B=list()
    B[[1]]=A
    A=B
  }
  
  maxNcategory=0
  V=c()
  for(I in 1:Ncategoryclass){
    
    if(maxNcategory<length(A[[I]])){
      maxNcategory=length(A[[I]])
    }
    
    s =  category[,I]
    v1=rep(0,N)
    
    u = unique(s) 
    # reorder with reference category as first element
    ref =  reference.category[I]
    u=unique(c(ref, u)) # put the reference category as the first element of u
    
    for(i in seq(1,length(u))){
      v1[which(s==u[i])] = i
    }
    V=cbind(V,v1)
  }
  
  # list of all combinations
  l=NULL
  for(I in 1:Ncategoryclass){
    l[I]  = list(unique(V[,I]))
    
  }
  Exp = expand.grid(l)
  
  rowProd <-  function(X){
    return(prod(X==a))
  }
  
  categoryindex=c()
  
  for(i in 1:N){
    a=V[i,]
    apply(Exp,1, FUN = rowProd)
    categoryindex[i]= which(apply(Exp,1, FUN = rowProd)==1)
  }
  
  MatrixCategory= Exp
  Ncategory = dim(Exp)[1]
  
  
  unique.categories  = unique(as.vector(category))
  Ncat.unique =length(unique.categories)
  
  
  L1 <- function(x){ 
    return(length(which(x==1))==Ncategoryclass-1)
  }
  L2 <- function(x){ 
    G=0
    if(L1(x) == TRUE){
      G=which(x!=1)
      
    }
    return(G)
  }
  
  
  index = which(apply(X = Exp, 1, FUN = L1))
  w=(apply(X = Exp, 1, FUN = L2))
  
  col.index = w[w>0]
  
  
  # comment
  category.position.in.table= data.frame(predictor = character(),
                                         relative_to = character(),
                                         index = integer())
  I=0
  for(i in unique(col.index)){
    ref =  reference.category[i]
    U=unique(c(ref,  unique(category[,i]))) # put the reference category as the first element of u
    for(j in 2:(length(U)) ){
      I=I+1
      de<-data.frame(U[j],U[1], index[I])
      category.position.in.table=rbind(category.position.in.table, setNames(de, names(category.position.in.table)))
    }
  }
  
  return(  list(category = category,
                categoryindex=categoryindex,
                MatrixCategory = MatrixCategory,
                Ncategory = Ncategory,
                maxNcategory=maxNcategory,
                Ncategoryclass=Ncategoryclass,
                unique.categories=unique.categories,
                Ncat.unique = Ncat.unique,
                category.position.in.table=category.position.in.table,
                reference.category=reference.category))
  
}
