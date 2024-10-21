data {
  int <lower=0> A; //the number of age classes
  int <lower=1> NAgeGroups ;   
  int <lower=1> class1[A,NAgeGroups]; //lower boundary for the age class corresponding to the indexed age
  int <lower=1> class2[A,NAgeGroups]; //upper boundary for the age class corresponding to the indexed age  
  int <lower=0> NGroups; //the number of foi groups      
  int <lower=0> N; //the number of individuals  
  int <lower=0> age[N]; // Age   
  int <lower=0, upper=1> Y[N]; // Outcome
  int<lower = 0, upper=1> seroreversion; 
  int <lower=1> categoryindex[N]; 
  int<lower= 1> Ncategoryclass; 
  int <lower=1> Ncategory;  
  int<lower=1> maxNcategory;
  int<lower=1> MatrixCategory[Ncategory,Ncategoryclass];
  int <lower=0> age_at_sampling[N]; 
  int <lower=0> sampling_year[N];  
  int <lower=1> age_group[N] ;   
  int <lower=1> age_at_init[NAgeGroups]; 
  int <lower=1> K; // the number of peaks of epidemics
  real <lower = 0> prioralpha1[K];
  real <lower = 0> prioralpha2[K];
  real priorT1[K];
  real <lower = 0> priorT2[K];
  real<lower =0, upper=1> se;
  real<lower =0, upper=1> sp;
  real <lower = 0> priorRho1;
  real <lower = 0> priorRho2;
  int <lower = 0> cat_lambda; // 1 or 0: characterizes whether we distinguish categories by different FOI
  
  //  1 or 2 to specify the prior distributions
  int prior_distribution_alpha;
  int prior_distribution_T;
  int prior_distribution_constant_foi;
  int prior_distribution_independent_foi;
  int prior_distribution_rho;
}

parameters {
  real T_raw[K];
  real alpha_raw[K];
  real rho_raw;      
  real  Flambda2[maxNcategory,Ncategoryclass];  
}

transformed parameters {
  real x[A]; 
  real L;
  real<lower=0> lambda[A];
  real S[K]; // Normalization constant
  real<lower =0, upper=1> P1[A,NAgeGroups,Ncategory];  
  real<lower =0, upper=1> P[A,NAgeGroups,Ncategory];  
  real<lower =0> Flambda[Ncategory]; 
  real<lower = 0, upper=1> Likelihood[N];  
  real log_lik[N];  
  real T[K];
  real<lower = 0> alpha[K];
  real<lower = 0, upper = 20> rho;      
  real c; 
  
  for(i in 1:K){
    if(prior_distribution_alpha == 1){ //normal distribution
    alpha[i] = prioralpha1[i]*exp(alpha_raw[i]*prioralpha2[i]);
    }
    if(prior_distribution_alpha == 2){ // exponential distribution
    alpha[i] =   alpha_raw[i];
    } 
    if(prior_distribution_T == 1){ //normal distribution
    T[i] = priorT1[i] + T_raw[i]*priorT2[i];
    }
    if(prior_distribution_T == 2){ // exponential distribution
    T[i] = T_raw[i] ;
    } 
  }
  if(prior_distribution_rho == 1){ //normal distribution
  rho = priorRho1*exp( rho_raw*priorRho2);
  }
  if(prior_distribution_rho == 2){ // exponential distribution
  rho =  rho_raw ;
  } 
  //}
  for(i in 1:K){
    S[i] =0;
    for(j in 1:A){
      S[i]  =  S[i]  + exp(-(j-T[i])^2);
    }
  }
  
  for(j in 1:A){
    lambda[j] =0;
    for(i in 1:K){
      lambda[j]  =  lambda[j]  + alpha[i]/S[i]*exp(-(j-T[i])^2);
    }
  }
  
  c=0;
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
  if(seroreversion==0){
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
  
  if(seroreversion==1){
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
  // accounting for age classes
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
  for (i in 1:K){
    if(prior_distribution_alpha == 1){
      alpha_raw[i] ~ normal(0,1);// T[-prioralpha1[i]/prioralpha2[i],];
    }
    if(prior_distribution_alpha == 2){
      alpha_raw[i] ~  exponential(prioralpha1[i]);
    }    
    
    if(prior_distribution_T == 1){
      T_raw[i] ~ normal(0,1) ;
    }
    if(prior_distribution_T == 2){
      T_raw[i] ~  exponential(priorT1[i]);
    }
  }
  
  for(I in 1:Ncategoryclass){
    for(i in 1:maxNcategory){      
      Flambda2[i,I] ~ normal(0,1.73) ; // the prior corresponds to xx
    }
  }
  
  if(prior_distribution_rho == 1){
    rho_raw  ~ normal(0,1);
  } 
  if(prior_distribution_rho == 2){
    rho_raw  ~ exponential(priorRho1);
  }
  for (j in 1:N) {  
    target += bernoulli_lpmf( Y[j] | Likelihood[j]);
  }
}

