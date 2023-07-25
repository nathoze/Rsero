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

    int <lower=1> Ncategory;  

    int<lower= 1> Ncategoryclass; 

    int<lower=1> maxNcategory; 

    int<lower=1> MatrixCategory[Ncategory,Ncategoryclass];
      
    int <lower=0> age_at_sampling[N]; 

    int <lower=0> sampling_year[N];  

    int <lower=1> age_group[N] ; 
  
    int <lower=1> age_at_init[NAgeGroups]; 

    int <lower=1> K; // the number of peaks of epidemics

    real <lower = 0> prioralpha1;

    real <lower = 0> prioralpha2;

    real <lower = 0> priorbeta1;

    real <lower = 0> priorbeta2;

    real <lower = 0> priorC1;

    real <lower = 0> priorC2;

    real priorT1;

    real <lower = 0> priorT2;
    
    real <lower = 0> priorRho1;
    
    real <lower = 0> priorRho2;
    
    real<lower =0, upper=1> se;
    
    real<lower =0, upper=1> sp;
    
    int <lower = 0> cat_lambda; // 1 or 0: characterizes whether we distinguish categories by different FOI
}



parameters {
    real  T[K];
    real<lower=0> alpha[K];
    real<lower=0> beta[K];
    real<lower = 0, upper = 20> rho;    
   real  Flambda2[maxNcategory,Ncategoryclass];  
    real<lower = 0> constant;
 
}

transformed parameters {

    real x[A]; 
    real L;
    real lambda[A];
    real S[K]; // Normalization constant

    real<lower =0, upper=1> P1[A,NAgeGroups,Ncategory]; //14 08 
    real<lower =0, upper=1> P[A,NAgeGroups,Ncategory]; //14 08 
    real<lower =0> Flambda[Ncategory]; //14 08
    real<lower = 0, upper=1> Likelihood[N];  
    real c; 
      
    for(i in 1:K){
        S[i] =0;
        for(j in 1:A){
           S[i]  =  S[i]  + exp(-((j-T[i])^2)/(beta[i])^2);
        }
    }

     for(j in 1:A){
        lambda[j] =0;
        for(i in 1:K){
            lambda[j]  =  constant + lambda[j]  + alpha[i]/S[i]*exp(-((j-T[i])^2)/(beta[i])^2);
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
                P1[1,J,i] = exp(-Flambda[i]*lambda[1]) ;
                for(j in 1:A-1){
                    x[j]=1;         
                    if(j<age_at_init[J]){
                        P1[j+1,J,i] = exp(-Flambda[i]*lambda[j]) ;    
                    }else{
                        P1[j+1,J,i] = P1[j,J,i]*exp(-Flambda[i]*lambda[j+1]);                 
                    } 
                }
                x[A]=1;
            }
        } 
    }

  

/*
    if(seroreversion==1){
        for(J in 1:NAgeGroups){
            for(i in 1:Ncategory){        
                for(j in 1:A){
                    x[j] = 1; 
                    for(k in 2:j){
                        L=Flambda[i]*lambda[j-k+2];
                        x[j-k+2-1] = x[j-k+2]*exp(-(rho+L)) +rho/(L+rho)*(1- exp(-(rho+L)));
                    }
                    P1[j,J,i]  = x[age_at_init[J]];
                }
            }
        }
    }
 */
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
     //   Like[j] =1-(1-bg)*P[age[j],age_group[j],categoryindex[j]];
       //    Likelihood[j] =1-sp+(se+sp-1)*(1-P[age[j],age_group[j],categoryindex[j]]); 
           Likelihood[j] =se-(se+sp-1)*P[age[j],age_group[j],categoryindex[j]]; 
 }

}

model {

  //FOI by group
    for (i in 1:K){
        T[i] ~ uniform(priorT1, priorT2);
        alpha[i] ~ uniform(prioralpha1, prioralpha2);
        beta[i] ~ uniform(priorbeta1, priorbeta2) ; 
    }
    rho  ~ uniform(priorRho1, priorRho2);
    constant ~ uniform(priorC1, priorC2);

    for(I in 1:Ncategoryclass){
        for(i in 1:maxNcategory){      
            Flambda2[i,I] ~ normal(0,1.73) ;
        }
    }
    for (j in 1:N) {  
        target += bernoulli_lpmf( Y[j] | Likelihood[j]);

    }
}

 
