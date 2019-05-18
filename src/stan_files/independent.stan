data {
    int <lower=0> A; //the number of age classes
  
    int <lower=0> NGroups; //the number of foi groups
      
    int <lower=0> N; //the number of individuals
  
    int <lower=0> age[N]; // Age
  
    int <lower=0, upper=1> Y[N]; // Outcome

    int<lower = 0, upper=1> seroreversion; 

    int<lower = 0, upper=1> background; 

    int <lower=1> categoryindex[N];  

    int <lower=1> Ncategory;  
  
    int <lower=0> age_at_sampling[N];  

    int <lower=0> sampling_year[N];  

    int <lower=1> NAgeGroups ; 

    int <lower=1> age_group[N] ;  
  
    int <lower=1> age_at_init[NAgeGroups]; 
    
  //  int <lower=1, upper=NGroups> ind_by_age[A]; // 

    real <lower = 0> priorY1;

    real <lower = 0> priorY2;

    real <lower = 0> priorbg1;

    real <lower = 0> priorbg2;

    real <lower = 0> priorRho1;

    real <lower = 0> priorRho2;

    int <lower = 0> cat_bg;  // 1 or 0: characterizes whether we distinguish categories by different background infection probability

    int <lower = 0> cat_lambda; // 1 or 0: characterizes whether we distinguish categories by different FOI
}


parameters {
    real logitlambda[NGroups]; 
    real<lower = 0, upper = 1> rho;    
    real<lower = 0, upper=1> bg2[Ncategory];
    real Flambda2[Ncategory]; 
}


transformed parameters {
    real x[A]; 
    real L;
    real lambda[A];
    real<lower =0, upper=1> P[A,NAgeGroups,Ncategory];
    real<lower =0> bg[Ncategory];
    real<lower =0> Flambda[Ncategory];
    real<lower = 0, upper=1> Like[N];   

    for (j in 1:A) {
         lambda[j] = inv_logit(logitlambda[j]);
    
    }

    if(!cat_bg){
        for(i in 1:Ncategory){
            bg[i] = bg2[1];
        }
    }else{
        for(i in 1:Ncategory){
            bg[i] = bg2[i];
        }
    }


    if(background==0){
        for(i in 1:Ncategory){
            bg[i] = 0;
        }
    }
    
    if(!cat_lambda){
        for(i in 1:Ncategory){
            Flambda[i] = 1;
        }
    }else{
        for(i in 1:Ncategory){
            Flambda[i] =exp(Flambda2[i]);
        }   
        Flambda[1] = 1;
    }
    
   
    L=1;


    if(seroreversion==0){
        for(J in 1:NAgeGroups){
            for(i in 1:Ncategory){      
                P[1,J,i] = exp(-Flambda[i]*lambda[1]) ;
                for(j in 1:A-1){
                    x[j]=1;         
                    if(j<age_at_init[J]){
                        P[j+1,J,i] = exp(-Flambda[i]*lambda[j]) ;    
                    }else{
                        P[j+1,J,i] = P[j,J,i]*exp(-Flambda[i]*lambda[j+1]);                 
                    } 
                }
                x[A]=1;
            }
        } 
    }

    if(seroreversion==1){
        for(J in 1:NAgeGroups){
            for(i in 1:Ncategory){        
                for(j in 1:A){
                    x[j] = 1; 
                    for(k in 2:j){
                        L=Flambda[i]*lambda[j-k+2];
                        x[j-k+2-1] = x[j-k+2]*exp(-(rho+L)) +rho/(L+rho)*(1- exp(-(rho+L)));
                    }
                    P[j,J,i]  = x[age_at_init[J]];
                }
            }
        }
    }

   for(j in 1:N){
        Like[j] =1-(1-bg[categoryindex[j]])*P[age[j],age_group[j],categoryindex[j]];
    }

}


model {

  //FOI by group
 
    for (j in 1:NGroups) {    
        logitlambda[j] ~ normal(priorY1,priorY2);
    }
 

    for(i in 1:Ncategory){
        bg2[i] ~uniform(priorbg1, priorbg2);   // category background infection. Size = Ncategory 
    }

    for(i in 1:Ncategory){
        Flambda2[i] ~ normal(0,1.73) ;
    }
    rho  ~ uniform(priorRho1, priorRho2);


    for (j in 1:N) {  
        target += bernoulli_lpmf( Y[j] | Like[j]);
    }
}
