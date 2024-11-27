#' @title Plot the seroprevalence by age 
#'
#' @description Plot the mean and 95 \%  confidence interval of the seroprevalence by age-class. Accepts as an input a \code{SeroData} object. If multiple \code{category} are defined, it will also compute the seroprevalence for each category.  
#' 
#' @author Nathanael Hoze \email{nathanael.hoze@gmail.com}
#' 
#' @param serodata An object of the class \code{SeroData}.
#' 
#' @param age_class Integer. The length in years of the age classes. Default = 10. 
#' 
#' @param YLIM Upper limit of the y-axis. Default = 1. The lower limit is set to 0. 
#' 
#' @param mid.age.plot Boolean. Whether the data is plotted at the mid-point of the age category (TRUE) or at the mean age of the individuals in the data within this category (FALSE). Default = \code{TRUE}. 
#'  
#' @return a list with plots of the seroprevalence for each category, or one plot if only one category is defined. 
#'
#' @export
#' @examples
#' 
#' dat  = data("one_peak_simulation")
#' seroprevalence.plot(serodata = dat)
#' 

seroprevalence.plot<- function(serodata, age_class = 10, YLIM = 1, mid.age.plot = TRUE, ...){
  plots  <- NULL
  index.plot=0
  unique.categories= serodata$unique.categories
  for(sampling_year in sort(unique(serodata$sampling_year))){
    for(cat in unique.categories){
      index.plot <- index.plot+1
      
      w <- which(serodata$sampling_year ==  sampling_year & serodata$category==cat, arr.ind = TRUE)[,1]
      
      if(length(unique.categories)==1){
        title =  paste0("Sampling year: ", sampling_year)
      }
      if(length(unique.categories)>1){
        title= paste0('Category: ',cat," Sampling year: ", sampling_year)
      }
      
      if(length(w)>0){
        subdata <- subset(serodata,sub = w)
        histdata <- sero.age.groups(dat = subdata, age_class = age_class, YLIM=YLIM,  mid.age.plot = mid.age.plot) 
        histdata$labels_text <- as.character(histdata$labels)
        
        # histdata =histdata[-which(is.na(histdata$mean)),]
        
        # Find the index of the last non-NA row
        last_row_index <- tail(which(complete.cases(histdata)), 1)
        if(length(last_row_index)>0){
          histdata =histdata[1:last_row_index,]
        }
        #  print(histdata)
        x.breaks = as.numeric(histdata$age)#-round(age_class/2))
        histdata$age = histdata$age#+round(age_class/2)
        XLIM=max(subdata$age_at_sampling+round(age_class/2),na.rm = TRUE) #+round(age_class/2)
        
        g <- ggplot(histdata, aes(x=age, y=mean)) +
          geom_point() + 
          geom_segment(aes(x=age,y=lower, xend= age,yend=upper))+
          scale_x_continuous(breaks=x.breaks,labels=histdata$labels_text)+  
          theme_classic()+
          theme(axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                text=element_text(size=14))+
          xlab('Age')+
          ylab('Proportion seropositive')+
          ylim(0,YLIM)+ ggtitle(title)
        
        plots[[index.plot]]= g
        plots[[index.plot]]$category  = cat
        
      }
    }
    
  }
  return(plots)
  
}

#' @export
#' 
sero.age.groups <- function(dat,age_class,YLIM, mid.age.plot = TRUE){
  
  if(age_class<= max(dat$age_at_sampling)){
    age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = age_class)
  }else{
    age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = max(dat$age_at_sampling))
  }
  
  age_bin <- sapply(dat$age, function(x) tail(which(x-age_categories >= 0), 1L)) # find the closest element
  S <- as.integer(as.logical(dat$Y)) 
  S1 <- sapply(1:length(age_categories), function(x) length(which(age_bin==x)) )
  S2 <- sapply(1:length(age_categories), function(x) sum(S[which(age_bin==x)] ))
  C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
  
  data.age = data.frame(age = dat$age_at_sampling, age_bin = age_bin) 
  mean.age = data.age %>% group_by(age_bin) %>% summarise(mean.age = mean(age))
  
  df = data.frame(x=age_categories,y=S2/S1)
  
  G=matrix(NA,nrow =  dim(df)[1], ncol=3)
  
  for(j in seq(1,length(S1))){
    if(S1[j]>1){
      B= binom::binom.confint(x=S2[j],n = S1[j],methods = "exact")
      G[j,1]=B$lower
      G[j,2]=B$upper
      G[j,3]=B$mean
    }
  }
  
  G[which(G >YLIM)] =YLIM
  mid.point.age = c( (age_categories[1:length(age_categories)-1] +age_categories[2:length(age_categories)])/2, age_categories[length(age_categories)] ) 
  mean.point.age  = mid.point.age
  for(i in 1:length(mid.point.age)){
    w=which(mean.age$age_bin == i)
    if(length(w)>0)
      mean.point.age[i] = mean.age$mean.age[w]
  }
  
  
  if(mid.age.plot == TRUE){
    age.plot = mid.point.age
  }else
    age.plot = mean.point.age
  
  C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
  
  if(sum(C[1, ] - C[2, ]) == 0 ){ # means that the age categories are each 1 year long
    histo_label <- append(format(C[1, ]), paste(">=", tail(age_categories, n = 1), sep = ""))
  } else{
    histo_label <- append(apply(format(C), 2, paste, collapse = "-"), paste(">=", tail(age_categories, n = 1), sep = ""))
  }
  
  
  histdata <- data.frame(age = age.plot,
                         mean = G[,3],
                         lower = G[,1],
                         upper = G[, 2],
                         labels = factor(histo_label, levels=histo_label))
  
  return(histdata)
  
  
}

# get the seroprevalence (mean and 95%CI) for each age group
# sero.age.groups <- function(dat,age_class,YLIM){
#   
#   #  age_categories <- seq(from = 0, to = min(dat$A, max(dat$age)), by = age_class)
#   
#   if(age_class<= max(dat$age_at_sampling)){
#     age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = age_class)
#   }else{
#     age_categories <- seq(from = 0, to =  max(dat$age_at_sampling), by = max(dat$age_at_sampling))
#   }
# 
#   age_bin <- sapply(dat$age, function(x) tail(which(x-age_categories >= 0), 1L)) # find the closest element
#   S <- as.integer(as.logical(dat$Y)) 
#   S1 <- sapply(1:length(age_categories), function(x) length(which(age_bin==x)) )
#   S2 <- sapply(1:length(age_categories), function(x) sum(S[which(age_bin==x)] ))
#   C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
#   
#   df = data.frame(x=age_categories,y=S2/S1)
#   
#   G=matrix(NA,nrow =  dim(df)[1], ncol=3)
#   
#   for(j in seq(1,length(S1))){
#     if(S1[j]>1){
#       B= binom::binom.confint(x=S2[j],n = S1[j],methods = "exact")
#       G[j,1]=B$lower
#       G[j,2]=B$upper
#       G[j,3]=B$mean
#     }
#   }
#   
#   G[which(G >YLIM)] =YLIM
#   mean_age =  c( (age_categories[1:length(age_categories)-1] +age_categories[2:length(age_categories)])/2, age_categories[length(age_categories)] ) 
#   
#   C <- (rbind((age_categories[1:length(age_categories)-1]), (age_categories[2:length(age_categories)]-1)))
#   
#   if(sum(C[1, ] - C[2, ]) == 0 ){ # means that the age categories are each 1 year long
#     histo_label <- append(format(C[1, ]), paste(">=", tail(age_categories, n = 1), sep = ""))
#   } else{
#     histo_label <- append(apply(format(C), 2, paste, collapse = "-"), paste(">=", tail(age_categories, n = 1), sep = ""))
#   }
#   
#   
#   histdata <- data.frame(age = mean_age,
#                          mean = G[,3],
#                          lower = G[,1],
#                          upper = G[, 2],
#                          labels = factor(histo_label, levels=histo_label))
#   
#   return(histdata)
#   
#   
# }



