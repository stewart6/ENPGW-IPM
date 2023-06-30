//
// This Stan program defines a population model for Gray Whales
// Import data
data {
  int<lower=1> Nyrs ;
  int<lower=1> Nyrm1 ;
  int<lower=1> N_Abund ;
  int<lower=1> N_Calves ;
  int<lower=1> N_Strand ;
  int<lower=1> N_Cond ;
  int<lower=1> N_CondN ; //Added by Josh for NBound Migration Condition
  int<lower=1> NV ;
  array[NV] int<lower=1> N_Vobs ;
  array[NV] int<lower=1> N_Vmis ;
  array[NV,Nyrm1] int yrobs_V ;
  array[NV,Nyrm1] int yrmis_V ;
  array[NV] vector[Nyrm1] V_z ;
  vector[N_Abund] Abund ;
  vector[N_Abund] Abund_se ;
  array[N_Abund] int Abund_year ;
  vector[N_Calves] Calves ;
  vector[N_Calves] Calves_se ;
  array[N_Calves] int Calves_year ;
  array[N_Strand] int Strand ;       // non-anth strandings (counts)
  array[N_Strand] int AnthStrand ;   // anth strandings (counts)
  array[N_Strand] int Strand_year ;
  array[N_Strand] int<lower=1,upper=3> Strand_per ;
  vector<lower=0,upper=1>[N_Cond] Cond ;
  array[N_Cond] int Cond_year ;
  vector<lower=0,upper=1>[N_CondN] CondN ; // Josh Added Northbound Condition
  array[N_CondN] int CondN_year ; // Josh Added
  real r_max ;
  real N_init_pri ;
  real K_med_approx ;
  real ADF ;
  int<lower=0> Nitxn ; // number interactions - set to 0 if no interactions
  array[Nitxn,Nitxn==0 ? 0 : 2] int Intxvars ; 
}
// The parameters accepted by the model. 
parameters {
  // Estimate missing covariates as parameters 
  array[NV] vector[Nyrm1] V_mis ; // missing cov values: no prior = std normal
  // Other estimated parameters
  vector<lower=0>[NV] sig_V ;
  real N_init ;
  real Cond_base_lgt ;
  real Brate_base_lgt ;
  real<lower=0> gammaN_base ;   // base natural hazards
  real gammaA_base ;            // base anthropogenic hazards
  vector[NV+Nitxn] Beta ;       // coefficients, covariate effects on condition   
  real<lower=0,upper=5> phi ;   // density-dependent effects on condition
  real<lower=0,upper=10> alpha ;// density-dependent effects on Calf production
  real<lower=0,upper=10> gamma ;// density-dependent effects on mortality
  real<lower=0> sig_E ;
  real<lower=0> sig_A ;
  vector[Nyrm1] eps ;
  vector[Nyrm1] eta ;
  real<lower=0> tau ;           // precision of beta-distributed body conditions
  real<lower=0> tauN ;           // precision of beta-distributed NORTHBOUND body conditions //Josh Added
  real<lower=0,upper=1> Ppn_dtct_1 ;
  real<lower=0,upper=1> Ppn_dtct_2 ;  
  real<lower=0,upper=1> Ppn_dtct_3 ;
  real<lower=0,upper=1> NBAdjust ; // Decrease in expected condition from Southbound to Northbound whales
}
// Derived params (including latent population dynamics)
transformed parameters {
  array[NV] vector[Nyrm1] V ;
  matrix[Nyrm1, NV+Nitxn] X ;
  vector[Nyrm1] Env_eff ;
  real C_max ;
  real B_max ;
  real M_min ; 
  vector<lower=0>[Nyrm1] M; 
  real r_max_exp ;
  vector[Nyrm1] C ;
  vector[Nyrm1] CN ; //NBound Migration Condition
  vector[Nyrm1] B ;
  vector<lower=0>[Nyrm1] D_N ;
  vector<lower=0>[Nyrm1] D_A ;
  vector[Nyrs] N ;
  vector<lower=0>[3] Ppn_dtct_N ;
  vector<lower=0>[3] Ppn_dtct_A ;
  vector[Nyrm1] K_approx ; 
  // Assemble covariate matrix, including observed & missing (estimated) data
  for(i in 1:NV){
    V[i][yrobs_V[i,1:N_Vobs[i]]] = V_z[i][1:N_Vobs[i]] ;
    V[i][yrmis_V[i,1:N_Vmis[i]]] = V_mis[i][1:N_Vmis[i]] ;
    X[1:Nyrm1,i] = (V[i] - mean(V[i])) ./ sd(V[i]) ; // re-center/scale
  }
  // If approprite add interaction terms between main covariates
  //  (transformed to ensure intxn term approximately centered with unit variance)
  if(Nitxn>0){
    for(i in 1:Nitxn){
      X[1:Nyrm1,NV+i] = 0.1 .* (5 + X[1:Nyrm1,Intxvars[i][1]]) .* (5 + X[1:Nyrm1,Intxvars[i][2]]) - 2.5 ;
    }
  }
  // Some additional derived params
  C_max = inv_logit(Cond_base_lgt - phi * 0.01) / 4 ;
  B_max = inv_logit(Brate_base_lgt - alpha * 0.01) / 5 ;
  M_min = 1 - exp(-exp(-7 + gammaN_base + gamma * 0.01)) ; 
  r_max_exp = B_max - M_min ;
  Ppn_dtct_N[1] = Ppn_dtct_3 * Ppn_dtct_2 *Ppn_dtct_1 ;
  Ppn_dtct_N[2] = Ppn_dtct_3 * Ppn_dtct_2 ;
  Ppn_dtct_N[3] = Ppn_dtct_3 ;
  Ppn_dtct_A[1] = Ppn_dtct_3 * Ppn_dtct_2 * Ppn_dtct_1 * ADF; // ADF adjusts for incr. detection
  Ppn_dtct_A[2] = Ppn_dtct_2 * Ppn_dtct_3 * ADF ;
  Ppn_dtct_A[3] = Ppn_dtct_3 * ADF ;
  // Process model:
  N[1] = N_init ;
  Env_eff = exp((X * Beta) + sig_E * eps) ; // better conditions = higher values
  K_approx = K_med_approx .* Env_eff ;
  for(t in 1:Nyrm1){
    real haz_N ;
    real haz_A ;
    //real M ;
    C[t] = inv_logit(Cond_base_lgt - phi * (N[t] / K_approx[t]) ) / 4 ;
    CN[t] = C[t] * NBAdjust ;
    B[t] = inv_logit(Brate_base_lgt - alpha * (N[t] / K_approx[t])  ) / 5 ;  
    haz_N = exp(-7 + gammaN_base + gamma * (N[t] / K_approx[t]) ) ; 
    haz_A = exp(-7 + gammaA_base + sig_A * eta[t]) ; 
    M[t] = 1 - exp(-(haz_N + haz_A)) ;
    D_N[t] = N[t] * M[t] * (haz_N / (haz_N + haz_A)) ;
    D_A[t] = N[t] * M[t] * (haz_A / (haz_N + haz_A)) ;
    N[t+1] = N[t] * fmax(.01, 1 + B[t] - M[t]) ;
  }
}
// The model to be estimated. 
model {
  // Observed data model
  Abund ~ normal(N[Abund_year], Abund_se) ;
  Calves ~ normal(N[Calves_year] .* B[Calves_year], Calves_se) ;
  Cond ~ beta(C[Cond_year] * tau * 100, (1 - C[Cond_year]) * tau * 100) ;
  CondN ~ beta(CN[CondN_year] * tauN * 100, (1 - CN[CondN_year]) * tauN * 100) ;
  Strand ~ poisson(D_N[Strand_year] .* Ppn_dtct_N[Strand_per]) ;
  AnthStrand ~ poisson(D_A[Strand_year] .* Ppn_dtct_A[Strand_per]) ;
  r_max ~ normal(r_max_exp,0.001) ;
  // AR(1) variables for estimating missing covaiate data:
  for(i in 1:NV){
    V[i][1] ~ normal(0, 1) ;
    V[i][2:Nyrm1] ~ normal(V[i][1:(Nyrm1 - 1)], sig_V[i]) ;
  }
  // Hierarchical random effect for env. stochasticity
  eps ~ normal(0,1) ;
  // Hierarchical random effect for anthropogenic hazards
  eta ~ normal(0,1) ;
  // Priors
  sig_V ~ cauchy(0, 1) ;
  sig_E ~ cauchy(0, 1) ;  
  sig_A ~ cauchy(0, 1) ;  
  N_init ~ gamma(N_init_pri * .0005, 0.0005) ;
  Cond_base_lgt ~ normal(0, 2.5) ;
  Brate_base_lgt ~ normal(0, 2.5) ;
  gammaN_base ~ normal(0, 2.5) ;
  gammaA_base ~ normal(0, 2.5) ;
  Beta ~ normal(0, 2.5) ;
  phi ~ cauchy(0, 1) ;
  alpha ~ cauchy(0, 1) ;
  gamma ~ cauchy(0, 1) ;
  tau ~ cauchy(0, 2.5) ;
  tauN ~ cauchy(0, 2.5) ;
  Ppn_dtct_1 ~ beta(1,1) ;
  Ppn_dtct_2 ~ beta(1,1) ;
  Ppn_dtct_3 ~ beta(1,1) ;
  NBAdjust ~ beta(1,1) ;
}
// derived params
generated quantities{  
  real K_mean ; // numerical projection to estimate K under current conditions
  real K_med ;  // estimate of median value of instantaneous K
  vector[Nyrm1] K_inst ; // vector of realized instantaneous K by year
  vector[Nyrm1] BN ;     // vector of estimated N births
  real Nsum = 0 ;        // nuisance var
  real Nt = N[Nyrs] ;    // nuisance var
  real phi_star ;        // proportional reduction in Condition, 0% --> 100% K
  real alpha_star ;      // proportional reduction in Birth rt., 0% --> 100% K
  real gamma_star ;      // proportional reduction in Survival, 0% --> 100% K
  vector[N_Abund+N_Calves+N_Cond+N_CondN+2*N_Strand] log_lik ;
  phi_star = 1 - (inv_logit(Cond_base_lgt - phi) / 4) / (inv_logit(Cond_base_lgt) / 4) ;  
  alpha_star = 1 - (inv_logit(Brate_base_lgt - alpha) / 5) / (inv_logit(Brate_base_lgt) / 5) ;  
  gamma_star = 1 - exp(-exp(-7 + gammaN_base + gamma)) / exp(-exp(-7 + gammaN_base)) ;
  BN = B[1:Nyrm1] .* N[1:Nyrm1] ;
  // Project for enough years to allow  stabilization around K
  for(t in 1:200){
    real Env_eff_r ;
    real B_r ;
    real haz_r ;
    real M_r ;
    Env_eff_r = exp(normal_rng(mean(log(Env_eff)), sd(log(Env_eff)))) ;
    while(Env_eff_r>max(Env_eff) || Env_eff_r<min(Env_eff)){
      Env_eff_r = exp(normal_rng(mean(log(Env_eff)), sd(log(Env_eff)))) ;
    }
    B_r = inv_logit(Brate_base_lgt - alpha * (Nt/(K_med_approx * Env_eff_r))) / 5 ;  
    haz_r = exp(-7 + gammaN_base + gamma * (Nt/(K_med_approx * Env_eff_r))) ; 
    M_r = 1 - exp(-(haz_r)) ;
    Nt = Nt * (1 + B_r - M_r) ;
    Nsum = t > 100 ? Nsum + Nt : Nsum + 0 ;
  }
  K_mean = Nsum / 100 ;
  // Calculate realized variable annual K values
  for(i in 1:Nyrm1){
    Nt = K_approx[i] ;
    Nsum = 0 ;
    for(t in 1:150){
      real B_r ;
      real haz_r ;
      real M_r ;
      B_r = inv_logit(Brate_base_lgt - alpha * (Nt/K_approx[i] )) / 5 ;  
      haz_r = exp(-7 + gammaN_base + gamma * (Nt/K_approx[i] )) ; 
      M_r = 1 - exp(-(haz_r)) ;
      Nt = Nt * (1 + B_r - M_r) ;
      Nsum = t > 100 ? Nsum + Nt : Nsum + 0 ;
    }
    K_inst[i] = Nsum / 50 ;
  }
  K_med = quantile(K_inst ./ Env_eff, 0.5) ; 
  // Log-likelihood calcs 
  for(i in 1:N_Abund){
    log_lik[i] = normal_lpdf(Abund[i] | N[Abund_year[i]], Abund_se[i]) ;
  }
  for(i in 1:N_Calves){
    log_lik[N_Abund+i] = normal_lpdf(Calves[i] | N[Calves_year[i]] * B[Calves_year[i]], 
                                      Calves_se[i]) ;
  }
  for(i in 1:N_Cond){
    log_lik[N_Abund+N_Calves+i] = beta_lpdf(Cond[i] | 
                   C[Cond_year[i]] * tau * 100, (1 - C[Cond_year[i]]) * tau * 100) ;
  }
  for(i in 1:N_CondN){
    log_lik[N_Abund+N_Calves+N_Cond+i] = beta_lpdf(CondN[i] | 
                   CN[CondN_year[i]] * tauN * 100, (1 - CN[CondN_year[i]]) * tauN * 100) ;
  }
  for(i in 1:N_Strand){
    log_lik[N_Abund+N_Calves+N_Cond+N_CondN+i] = poisson_lpmf(Strand[i] | 
                   D_N[Strand_year[i]] .* Ppn_dtct_N[Strand_per[i]] ) ;
    log_lik[N_Abund+N_Calves+N_Cond+N_CondN+N_Strand+i] = poisson_lpmf(AnthStrand[i] | 
                   D_A[Strand_year[i]] .* Ppn_dtct_A[Strand_per[i]] ) ;
  }
}
