data {
  int <lower=0> A; //the number of age classes
  int <lower=1> NAgeGroups ;  
  
  array[A,NAgeGroups] int <lower=1> class1 ; //lower boundary for the age class corresponding to the indexed age
  array[A,NAgeGroups] int <lower=1> class2 ; //lower boundary for the age class corresponding to the indexed age
  int <lower=0> NGroups; //the number of foi groups      
  int <lower=0> N; //the number of individuals  
  array[N] int <lower=0> age; // Age   
  array[N] int <lower=0, upper=1> Y; // Outcome
  int<lower = 0, upper=1> seroreversion; 
  int<lower = 0, upper=1> age_dependent_foi; 
  array[N] int <lower=1> categoryindex ; 
  int<lower= 1> Ncategoryclass; 
  int <lower=1> Ncategory;  
  int<lower=1> maxNcategory;
  array[Ncategory,Ncategoryclass] int<lower=1> MatrixCategory ;
  array[N] int <lower=0> age_at_sampling ; 
  array[N] int <lower=0> sampling_year ;  
  array[N] int <lower=1> age_group ;   
  array[NAgeGroups] int <lower=1> age_at_init; 

  int <lower=1> group_size_length;   
  array[group_size_length] int <lower=1> group_size_array;   

  real <lower = 0> priorY1;  
  real <lower = 0> priorY2;  
  real <lower = 0> priorRho1;
  real <lower = 0> priorRho2;  
  real<lower =0, upper=1> se;  
  real<lower =0, upper=1> sp;  
  int <lower = 0> cat_lambda; // 1 or 0: characterizes whether we distinguish categories by different FOI
  
  //  1 or 2 to specify the prior distributions
  int prior_distribution_alpha;
  int prior_distribution_T;
  int prior_distribution_constant_foi;
  int prior_distribution_independent_foi;
  int prior_distribution_rho;
}

parameters {
  array[NGroups] real lambda_group_raw; 
  real rho_raw;    
  array[maxNcategory,Ncategoryclass] real Flambda2;  
  real age_risk;

}

transformed parameters {
 
  array[A] real x; 
  real L;
  array[A] real<lower = 0.0> lambda;  
  array[A,NAgeGroups,Ncategory] real<lower = 0, upper = 1> P1; 
  array[A,NAgeGroups,Ncategory] real<lower = 0, upper = 1> P;  
  array[Ncategory] real<lower = 0> Flambda;  
  array[N] real<lower = 0, upper = 1> Likelihood;  
  array[N] real log_lik;  
    real<lower = 0> C1;

  real c; 
  array[NGroups] real<lower = 0> lambda_group; 
  real<lower = 0, upper = 20> rho;    
  
  c=0;
  C1=0;

  for (j in 1:NGroups) {   
    if(prior_distribution_independent_foi == 1){
      lambda_group[j] =  priorY1*exp(priorY2*lambda_group_raw[j]);
    }
    if(prior_distribution_independent_foi == 2){
      lambda_group[j] = lambda_group_raw[j]; 
    }
  }
  
  if(prior_distribution_rho == 1){ //normal distribution
  rho = priorRho1*exp( rho_raw*priorRho2);
  }
  if(prior_distribution_rho == 2){ // exponential distribution
  rho =  rho_raw ;
  } 
  
  for (j in 1:A) {
    // group size 
    lambda[j] = lambda_group[group_size_array[j]];
  }
   
  if(!cat_lambda){
    for(i in 1:Ncategory){
      Flambda[i] = 1; 
    }
  }else{
    for(i in 1:Ncategory){
      c = 0;
      for(I in 1:Ncategoryclass){
        if(MatrixCategory[i,I]>1){ // if ==1, no change in the FOI
        c = c+ Flambda2[MatrixCategory[i,I], I];  
        }
      }   
      Flambda[i] =  exp(c);// exp(Flambda2[I,i]);
    }
  }
  
  L=1;
  
   if(seroreversion==0 && age_dependent_foi == 0){
    for(J in 1:NAgeGroups){
      for(i in 1:Ncategory){      
        P1[1,J,i ] = exp(-Flambda[i]*lambda[1]) ;
        for(j in 1:A-1){
          x[j]=1;         
          if(j<age_at_init[J]){
            P1[j+1,J,i] = exp(- Flambda[i]*lambda[j]) ;    
          }else{
            P1[j+1,J,i] = P1[j,J,i]*exp(- Flambda[i]*lambda[j+1]);                 
          } 
        }
        x[A]=1;
      }
    }
  }
  
  if(seroreversion==1 && age_dependent_foi == 0){
    for(J in 1:NAgeGroups){
      for(i in 1:Ncategory){    
        x[A] =1;
        for(j in 1:A){
          x[j] = exp(-Flambda[i]*lambda[1]) ;
        }
        
        for(j in 1:A){
          //  x[j] = exp(-Flambda[i]*lambda[age_at_init[J]]) ;
          L=Flambda[i]*lambda[age_at_init[J]];
          x[j] = rho/(L+rho) +L/(L+rho)*exp(-L) ;
          if(j >1){
            for(k in 2:j){
              L=Flambda[i]*lambda[j-k+2];
              x[j-k+2-1] = x[j-k+2]*exp(-(rho+L)) +rho/(L+rho)*(1- exp(-(rho+L)));
            }
          }         
          P1[j,J,i]  = x[age_at_init[J]];
          
        }
      }
    }
  }

  if(seroreversion==1 && age_dependent_foi == 1){
    for(J in 1:NAgeGroups){
      for(i in 1:Ncategory){    
        x[A] =1;
        for(j in 1:A){
          x[j] = exp(-Flambda[i]*lambda[1]) ;
        }
        
        for(j in 1:A){
          //  x[j] = exp(-Flambda[i]*lambda[age_at_init[J]]) ;
          L=Flambda[i]*lambda[age_at_init[J]];
          C1= exp(age_risk*(age_at_init[J]-1));// equivalent

          x[j] = rho/(L*C1+rho) +L*C1/(L*C1+rho)*exp(-L*C1) ;
          if(j >1){
            for(k in 2:j){
              L=Flambda[i]*lambda[j-k+2];
              C1= C1*exp(age_risk);
              x[j-k+2-1] = x[j-k+2]*exp(-(rho+L*C1)) +rho/(L*C1+rho)*(1- exp(-(rho+L*C1)));
            }
          }         
          P1[j,J,i]  = x[age_at_init[J]];
          
        }
      }
    }
  }

  if(seroreversion==0 && age_dependent_foi == 1){
    for(J in 1:NAgeGroups){
      for(i in 1:Ncategory){    
        x[A] =1;
        for(j in 1:A){
          x[j] = exp(-Flambda[i]*lambda[1]) ;
        }
        
        for(j in 1:A){
          //  x[j] = exp(-Flambda[i]*lambda[age_at_init[J]]) ;
          L=Flambda[i]*lambda[age_at_init[J]];
          C1= exp(age_risk*(age_at_init[J]-1));// equivalent

          x[j] = exp(-L*C1) ;
          if(j >1){
            for(k in 2:j){
              L=Flambda[i]*lambda[j-k+2];
              C1= C1*exp(age_risk);
              x[j-k+2-1] = x[j-k+2]*exp(-L*C1) ;
            }
          }         
          P1[j,J,i]  = x[age_at_init[J]];
          
        }
      }
    }
  }

  for(J in 1:NAgeGroups){
    for(i in 1:Ncategory){        
      for(j in 1:A){
        P[j,J,i]=0;
        for(k in class1[j,J]:class2[j,J]){
          P[j,J,i]  = P1[k,J,i]+P[j,J,i];
        }
        P[j,J,i] = P[j,J,i]/(class2[j,J]-class1[j,J]+1);
      }
    }
  }
  
  for(j in 1:N){
    Likelihood[j] =se-(se+sp-1)*P[age[j],age_group[j],categoryindex[j]]; 
    log_lik[j] = bernoulli_lpmf( Y[j] | Likelihood[j]);
  }   
}

model {
  for (j in 1:NGroups) {    
    if(prior_distribution_independent_foi == 1){
      lambda_group_raw[j] ~ normal(0,1);
    }
    if(prior_distribution_independent_foi == 2){
      lambda_group_raw[j] ~ exponential(priorY1);
    }
  }
  
  for(I in 1:Ncategoryclass){
    for(i in 1:maxNcategory){      
      Flambda2[i,I] ~ normal(0,1.73) ;
    }
  }
  
  if(prior_distribution_rho == 1){
    rho_raw  ~ normal(0,1);
  } 
  if(prior_distribution_rho == 2){
    rho_raw  ~ exponential(priorRho1);
  }

    age_risk ~ normal(0,1);


  for (j in 1:N) {  
    target += bernoulli_lpmf( Y[j] | Likelihood[j]);
  }
}
