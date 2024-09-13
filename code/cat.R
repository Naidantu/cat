Statement = 20            
Correlation1 = .40      
Correlation2 = .30      
ExposureControl = 0     
Scoring = 2            

Trait = 5            
Obs = 100      
Item = 5       
P=2     	        


ggum <- function(a,d,t,th){
  tmp1=exp(a*(1*(th-d)-t))
  tmp2=exp(a*(2*(th-d)-t))
  tmp3=exp(a*(3*(th-d)))
  prob=(tmp1+tmp2)/(1+tmp1+tmp2+tmp3)
  return(prob)
}
mupp <- function(par){
  th1=par[1]
  a1=par[2]
  d1=par[3]
  t1=par[4]
  th2=par[5]
  a2=par[6]
  d2=par[7]
  t2=par[8]
  p1=ggum(a1,d1,t1,th1)
  p2=ggum(a2,d2,t2,th2)
  q1=1-p1
  q2=1-p2
  pp=(p1*q2)/(p1*q2+q1*p2)
  return(pp)
}
fcirt_nocov <-'
functions {

  real MUPP(int y, real theta1, real theta2, real alpha1, real alpha2, real delta1, real delta2, real tau1, real tau2) {

    vector[2] prob;
    real num01;
    real num11;
    real denominator1;
    real spr1;
    real num02;
    real num12;
    real denominator2;
    real spr2;
    real pst;
    real pts;

    num01 = 1 + exp(alpha1*(3*(theta1-delta1)));

    num11 = exp(alpha1*((theta1-delta1)-tau1)) + exp(alpha1*((2*(theta1-delta1))-tau1));

    denominator1 = num01 + num11;

    spr1 = num11/denominator1;

    num02 = 1 + exp(alpha2*(3*(theta2-delta2)));

    num12 = exp(alpha2*((theta2-delta2)-tau2)) + exp(alpha2*((2*(theta2-delta2))-tau2));

    denominator2 = num02 + num12;

    spr2 = num12/denominator2;

    pst = spr1*(1-spr2); //10 #p=pair specification map for the given pairs

    pts = spr2*(1-spr1); //01

    prob[1] = pst/(pst+pts);

    prob[2] = pts/(pst+pts);

    //return categorical_lpmf(y|prob);

    return categorical_lpmf(y|prob);
  }
}

data {

  int<lower=1> n_student;
  int<lower=1> n_item;
  int<lower=1> n_pair;
  int<lower=1> N;
  int<lower=1> n_dimension;
  int<lower=0> p[n_pair, 2]; //p=pair specification map for the given pairs
  int<lower=1, upper=2> y[N];
  int<lower=1> I1;
  int<lower=1> J1;
  int<lower=1, upper=I1> II[N];
  int<lower=1, upper=J1> JJ[N];

  vector<lower=0, upper=4> [n_item] alpha;
  vector<lower=-5, upper=5> [n_item] delta; 
  vector<lower=-5, upper=5> [n_item] tau;

  // user-defined priors
  //real ma;
  //real va;
  //real md;
  //real vd;
  //real mt;
  //real vt;

  int<lower=1> ind[2*N];  //pairs: 2*N  triplets: 3*N
  vector[n_dimension] theta_mu;

  //vector[n_item] md;
  //vector[n_item] vd;
  //vector<lower=-5, upper=5>[n_item] Delta_lower;              // lower bounds of delta
  //vector<lower=-5, upper=5>[n_item] Delta_upper;              // upper bounds of delta

}

parameters {

  vector[n_dimension] theta[n_student];
  //matrix [n_student, n_dimension] theta;
  //vector<lower=0, upper=4> [n_item] alpha;
  //vector<lower=-5, upper=5> [n_item] delta;
  //vector<lower=-5, upper=0> [n_item] tau;

  //vector[trait] theta[J];
  cholesky_factor_corr[n_dimension] L_Omega;

  //vector<lower=0, upper=1>[n_item] delta_raw;

}

//transformed parameters {
//  vector[n_item] delta = Delta_lower + (Delta_upper - Delta_lower) .* delta_raw;
//}

model {

  //alpha ~ lognormal(ma,va);
  //delta ~ normal(md,vd);
  //tau ~ normal(mt,vt);
  L_Omega  ~ lkj_corr_cholesky(1);

  theta~ multi_normal_cholesky(theta_mu,L_Omega);

  //for (i in 1:n_item){
  //  delta[i] ~ normal(md[i],vd[i]);
  //}

  for (n in 1:N){

    target += MUPP(y[n],theta[JJ[n],ind[2*n-1]],theta[JJ[n],ind[2*n]],alpha[p[II[n],1]],alpha[p[II[n],2]],delta[p[II[n],1]],delta[p[II[n],2]],tau[p[II[n],1]],tau[p[II[n],2]]);

  }
}

generated quantities{
  matrix[n_dimension,n_dimension] Cor;
  vector[N] log_lik;

  Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = MUPP(y[n],theta[JJ[n],ind[2*n-1]],theta[JJ[n],ind[2*n]],alpha[p[II[n],1]],alpha[p[II[n],2]],delta[p[II[n],1]],delta[p[II[n],2]],tau[p[II[n],1]],tau[p[II[n],2]]);

  }
}
'
fcirt_withcov <-'
functions {

  real MUPP(int y, real theta1, real theta2, real alpha1, real alpha2, real delta1, real delta2, real tau1, real tau2) {

    vector[2] prob;
    real num01;
    real num11;
    real denominator1;
    real spr1;
    real num02;
    real num12;
    real denominator2;
    real spr2;
    real pst;
    real pts;

    num01 = 1 + exp(alpha1*(3*(theta1-delta1)));

    num11 = exp(alpha1*((theta1-delta1)-tau1)) + exp(alpha1*((2*(theta1-delta1))-tau1));

    denominator1 = num01 + num11;

    spr1 = num11/denominator1;

    num02 = 1 + exp(alpha2*(3*(theta2-delta2)));

    num12 = exp(alpha2*((theta2-delta2)-tau2)) + exp(alpha2*((2*(theta2-delta2))-tau2));

    denominator2 = num02 + num12;

    spr2 = num12/denominator2;

    pst = spr1*(1-spr2); //10 #p=pair specification map for the given pairs

    pts = spr2*(1-spr1); //01

    prob[1] = pst/(pst+pts);

    prob[2] = pts/(pst+pts);

    //return categorical_lpmf(y|prob);

    return categorical_lpmf(y|prob);
  }
}

data {

  int<lower=1> n_student;
  int<lower=1> n_item;
  int<lower=1> n_pair;
  int<lower=1> N;
  int<lower=1> n_dimension;
  int<lower=0> p[n_pair, 2]; //p=pair specification map for the given pairs
  int<lower=1, upper=2> y[N];
  int<lower=1> I1;
  int<lower=1> J1;
  int<lower=1, upper=I1> II[N];
  int<lower=1, upper=J1> JJ[N];

  vector<lower=0, upper=4> [n_item] alpha;
  vector<lower=-5, upper=5> [n_item] delta; 
  vector<lower=-5, upper=5> [n_item] tau;

  // user-defined priors
  //real ma;
  //real va;
  //real md;
  //real vd;
  //real mt;
  //real vt;

  int<lower=1> ind[2*N];  //pairs: 2*N  triplets: 3*N
  vector[n_dimension] theta_mu;
  int<lower=1> P;   //number of covariates
  matrix[n_dimension*n_student,P] PC1;   //long format person covariates

  //vector[n_item] md;
  //vector[n_item] vd;
  //vector<lower=-5, upper=5>[n_item] Delta_lower;              // lower bounds of delta
  //vector<lower=-5, upper=5>[n_item] Delta_upper;              // upper bounds of delta

}

parameters {

  vector[n_dimension] theta[n_student];
  //matrix [n_student, n_dimension] theta;
  //vector<lower=0, upper=4> [n_item] alpha;
  //vector<lower=-5, upper=5> [n_item] delta;
  //vector<lower=-5, upper=0> [n_item] tau;

  //vector[trait] theta[J];
  cholesky_factor_corr[n_dimension] L_Omega;
  matrix[P,n_dimension] lambda;    // person covariate regression coefficient

  //vector<lower=0, upper=1>[n_item] delta_raw;

}

//transformed parameters {
//  vector[n_item] delta = Delta_lower + (Delta_upper - Delta_lower) .* delta_raw;
//}

model {

  matrix[n_dimension, n_dimension] center_mu;  //PC1*lambda matrix
  row_vector[n_dimension] center_mu0;    //first row of PC1*lambda matrix
  vector[n_dimension] center_mu1;        //change center_mu0 row to column

  //alpha ~ lognormal(ma,va);
  //delta ~ normal(md,vd);
  //tau ~ normal(mt,vt);
  L_Omega  ~ lkj_corr_cholesky(1);

  //theta~ multi_normal_cholesky(theta_mu,L_Omega);
  for (j in 1:n_student){
      center_mu=PC1[(n_dimension*j-(n_dimension-1)):n_dimension*j,]*lambda;
      center_mu0=center_mu[1,];
      center_mu1=to_vector(center_mu0);
      theta[j]~ multi_normal_cholesky(theta_mu + center_mu1, L_Omega);
    }
     to_vector(lambda) ~ normal (0, 1);

  //for (i in 1:n_item){
  //  delta[i] ~ normal(md[i],vd[i]);
  //}

  for (n in 1:N){

    target += MUPP(y[n],theta[JJ[n],ind[2*n-1]],theta[JJ[n],ind[2*n]],alpha[p[II[n],1]],alpha[p[II[n],2]],delta[p[II[n],1]],delta[p[II[n],2]],tau[p[II[n],1]],tau[p[II[n],2]]);

  }
}

generated quantities{
  matrix[n_dimension,n_dimension] Cor;
  vector[N] log_lik;

  Cor=multiply_lower_tri_self_transpose(L_Omega);

  for (n in 1:N) {
    log_lik[n] = MUPP(y[n],theta[JJ[n],ind[2*n-1]],theta[JJ[n],ind[2*n]],alpha[p[II[n],1]],alpha[p[II[n],2]],delta[p[II[n],1]],delta[p[II[n],2]],tau[p[II[n],1]],tau[p[II[n],2]]);

  }
}
'

#Step 1. 1.	Set the initial trait scores for the respondent to 0 for all the dimensions at the beginning of the testing session. 
ThetaEst <- matrix(0, Obs, Trait)
PSD <- matrix(0, Obs, Trait)
ParVals <- read.table('ParVals.txt', sep="",header=FALSE)
ParVals <- cbind(ParVals, 0, 0, matrix(rep(1:Statement, times = Trait, ncol = 1)))
#create all possible item pairting
Spar <- subset(ParVals, ParVals[, 1] == 1)
Tpar <- subset(ParVals, ParVals[, 1] == 2)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix[(m-1) * nrows2 + n, 1] <- m
    pair_matrix[(m-1) * nrows2 + n, 2] <- n
    pair_matrix[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 1)
Tpar <- subset(ParVals, ParVals[, 1] == 3)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix1 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix1[(m-1) * nrows2 + n, 1] <- m
    pair_matrix1[(m-1) * nrows2 + n, 2] <- n
    pair_matrix1[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix1[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 1)
Tpar <- subset(ParVals, ParVals[, 1] == 4)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix2 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix2[(m-1) * nrows2 + n, 1] <- m
    pair_matrix2[(m-1) * nrows2 + n, 2] <- n
    pair_matrix2[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix2[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 1)
Tpar <- subset(ParVals, ParVals[, 1] == 5)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix3 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix3[(m-1) * nrows2 + n, 1] <- m
    pair_matrix3[(m-1) * nrows2 + n, 2] <- n
    pair_matrix3[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix3[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 2)
Tpar <- subset(ParVals, ParVals[, 1] == 3)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix4 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix4[(m-1) * nrows2 + n, 1] <- m
    pair_matrix4[(m-1) * nrows2 + n, 2] <- n
    pair_matrix4[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix4[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 2)
Tpar <- subset(ParVals, ParVals[, 1] == 4)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix5 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix5[(m-1) * nrows2 + n, 1] <- m
    pair_matrix5[(m-1) * nrows2 + n, 2] <- n
    pair_matrix5[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix5[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 2)
Tpar <- subset(ParVals, ParVals[, 1] == 5)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix6 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix6[(m-1) * nrows2 + n, 1] <- m
    pair_matrix6[(m-1) * nrows2 + n, 2] <- n
    pair_matrix6[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix6[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 3)
Tpar <- subset(ParVals, ParVals[, 1] == 4)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix7 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix7[(m-1) * nrows2 + n, 1] <- m
    pair_matrix7[(m-1) * nrows2 + n, 2] <- n
    pair_matrix7[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix7[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 3)
Tpar <- subset(ParVals, ParVals[, 1] == 5)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix8 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix8[(m-1) * nrows2 + n, 1] <- m
    pair_matrix8[(m-1) * nrows2 + n, 2] <- n
    pair_matrix8[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix8[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Spar <- subset(ParVals, ParVals[, 1] == 4)
Tpar <- subset(ParVals, ParVals[, 1] == 5)
Spar <- as.matrix(Spar)
Tpar <- as.matrix(Tpar)

nrows1 <- nrow(Spar)
nrows2 <- nrow(Tpar)

pair_matrix9 <- matrix(ncol = 16, nrow = nrows1 * nrows2)

for (m in 1:nrows1) {
  # Loop through each row of matrix2
  for (n in 1:nrows2) {
    # Store the pairings in the pair_matrix
    pair_matrix9[(m-1) * nrows2 + n, 1] <- m
    pair_matrix9[(m-1) * nrows2 + n, 2] <- n
    pair_matrix9[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
    pair_matrix9[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
  }
}
Pairings <- rbind(pair_matrix, pair_matrix1, pair_matrix2, pair_matrix3, pair_matrix4, pair_matrix5,
                  pair_matrix6, pair_matrix7, pair_matrix8, pair_matrix9)
Pairings <- Pairings[, -c(7,8,14,15)]
Pairings <- cbind(Pairings, 0)

for (i in 1:Obs){
  ParVals[, 5] <- 0
  response <- matrix(nrow = 1, ncol = 0)
  pairmap <- matrix(nrow = 0, ncol = 2)
  ind <- matrix(nrow = 0, ncol = 1)
  ParInits <- matrix(nrow = 0, ncol = 3)
  Itempar <- matrix(nrow = 0, ncol = 3)
  for (j in 1:(Item*Trait)){
    #Step 2. Create all possible MUPP items by pairing single statements in the statement pool following the predetermined content constraints specified in the test blueprint. 
    Spar <- subset(ParVals, ParVals[, 1] == StateDim[2*j-1,1] & ParVals[, 5] == 0)
    Tpar <- subset(ParVals, ParVals[, 1] == StateDim[2*j,1] & ParVals[, 5] == 0)
    Spar <- as.matrix(Spar)
    Tpar <- as.matrix(Tpar)
    
    nrows1 <- nrow(Spar)
    nrows2 <- nrow(Tpar)
    
    pair_matrix <- matrix(ncol = 16, nrow = nrows1 * nrows2)
    
    for (m in 1:nrows1) {
      # Loop through each row of matrix2
      for (n in 1:nrows2) {
        # Store the pairings in the pair_matrix
        pair_matrix[(m-1) * nrows2 + n, 1] <- m
        pair_matrix[(m-1) * nrows2 + n, 2] <- n
        pair_matrix[(m-1) * nrows2 + n, 3:9] <- Spar[m, ]  
        pair_matrix[(m-1) * nrows2 + n, 10:16] <- Tpar[n, ] 
      }
    }
    #Step 3. Calculate item information at the most recent trait scores for the MUPP items.
    Iteminfo <- matrix(NA, nrows1 * nrows2, 1)
    
    for (w in 1:(nrows1 * nrows2)){
      
      pars <- c(ThetaEst[i, StateDim[2*j-1,1]],pair_matrix[w, 4],pair_matrix[w, 5],pair_matrix[w, 6],
                ThetaEst[i, StateDim[2*j,1]],pair_matrix[w, 11],pair_matrix[w, 12],pair_matrix[w, 13])
      
      prob <- mupp(pars)
      dp <- numDeriv::grad(mupp,pars)
      mgrad <- matrix(dp[c(1,5)],2,1)%*%matrix(dp[c(1,5)],1,2)
      info_mat <- diag(mgrad)/(prob*(1-prob))
      info <- sum(diag(mgrad))/(prob*(1-prob))
      Iteminfo[w, 1] <- info
    }
    
    Iteminfo <- cbind(Iteminfo, pair_matrix[,c(9, 16)])
    
    if (ExposureControl == 1){
      #Step 4. Stratify items evenly into 3 strata (i.e., low, moderate, and high information) based on the calculated item information.
      
      sorted_Iteminfo <- Iteminfo[order(Iteminfo[, 1]), ]
      
      stratum_size <- round(nrow(sorted_Iteminfo) / 3)
      stratum_low <- sorted_Iteminfo[1:stratum_size, ]
      stratum_moderate <- sorted_Iteminfo[(stratum_size + 1):(2 * stratum_size), ]
      stratum_high <- sorted_Iteminfo[(2 * stratum_size + 1):nrow(sorted_Iteminfo), ]
      #stratum_low <- as.matrix(stratum_low)
      #stratum_moderate <- as.matrix(stratum_moderate)
      #stratum_high <- as.matrix(stratum_high)
      
      #Step 5. For the first 1/3 of the test, the item is randomly selected from the stratum with low information items; for the second 1/3 of the test, select the most informative item from the stratum with moderate information items; and for the last 1/3 of the test, select the most informative item from the stratum with high information items. 
      if (j < round(Item*Trait / 3)){
        random_item <- stratum_low[sample(nrow(stratum_low), 1), ]
      }
      if (round(Item*Trait / 3) <= j && j < round(Item*Trait / 3)*2){
        random_item <- stratum_moderate[which.max(stratum_moderate[,1]), ]
      }
      if (j >= round(Item*Trait / 3)*2){
        random_item <- stratum_high[which.max(stratum_high[,1]), ]
      }
    }
    if (ExposureControl == 0){
      
      #In conditions where item exposure control is not used, Steps 4-5 will be replaced with one step: select the item with the maximum item information
      
      random_item <- Iteminfo[which.max(Iteminfo[,1]), ]
      
    }
    
    #Freeze the statements within person and keep track of Item and statement exposure rate across person
    selected_item <- pair_matrix[pair_matrix[, 9] == random_item[2] & pair_matrix[, 16] == random_item[3], ]
    ParVals[(StateDim[2*j-1,1]-1)*Statement+random_item[2], 5] <- 1
    ParVals[(StateDim[2*j,1]-1)*Statement+random_item[3], 5] <- 1
    ParVals[(StateDim[2*j-1,1]-1)*Statement+random_item[2], 6] <- ParVals[(StateDim[2*j-1,1]-1)*Statement+random_item[2], 6] + 1
    ParVals[(StateDim[2*j,1]-1)*Statement+random_item[3], 6] <- ParVals[(StateDim[2*j,1]-1)*Statement+random_item[3], 6] + 1
    Pairings[,13] <- ifelse(Pairings[,3] == StateDim[2*j-1,1] & Pairings[,8] == StateDim[2*j,1] & Pairings[,1] == random_item[2] & Pairings[,2] == random_item[3], Pairings[,13] + 1, Pairings[,13])
    
    #Step 6&7. Administer the selected item & Simulate MUPP item response
    pars <- c(Theta[i, StateDim[2*j-1,1]],selected_item[4],selected_item[5],selected_item[6],
              Theta[i, StateDim[2*j,1]],selected_item[11],selected_item[12],selected_item[13])
    prob <- mupp(pars)
    rand <- runif(1)
    
    if (prob>rand){
      response <- cbind(response, 1)
    } else{
      response <- cbind(response, 2)
    }
    
    ind <- rbind(ind, matrix(selected_item[c(3, 10)], ncol = 1))
    Itempar <- rbind(Itempar, selected_item[c(4, 5, 6)])
    Itempar <- rbind(Itempar, selected_item[c(11, 12, 13)])
    
    if (Scoring == 1){
    #Step 8. Estimate the trait scores of the respondent based on the response data and collateral information
    fcirt.Data <- as.matrix(response)
    pairmap <- rbind(pairmap, c((j * 2) - 1, j * 2))
    #ind <- rbind(ind, matrix(selected_item[c(3, 10)], ncol = 1))
    ParInits <- rbind(ParInits, c(1, ifelse(selected_item[5] < 0, -1, 1), -1))
    ParInits <- rbind(ParInits, c(1, ifelse(selected_item[12] < 0, -1, 1), -1))
    #Itempar <- rbind(Itempar, selected_item[c(4, 5, 6)])
    #Itempar <- rbind(Itempar, selected_item[c(11, 12, 13)])
    
    N1 <- nrow(fcirt.Data)
    
    Missing <- matrix(NA,nrow=(ncol(fcirt.Data))*N1,ncol=1)
    MissPattern<-data.frame(Missing=((as.numeric(is.na(t(fcirt.Data)))*-1)+1),ID=seq(1,((ncol(fcirt.Data))*N1),1))
    Miss<-subset(MissPattern,Missing==0)
    if (nrow(Miss)==0){
      ind1<-rep(ind,N1)
    }else{
      ind1<-rep(ind,N1)[-c(Miss$ID*2-1, Miss$ID*2)]
    }
    
    Data<-suppressWarnings(edstan::irt_data(response_matrix =fcirt.Data))
    #sample size
    I2<-dim(fcirt.Data)[1]
    #number of pairs
    J2<-dim(fcirt.Data)[2]
    #S <- J2*2
    S <- max(pairmap)
    #D <- max(ind)
    
    #initial values
    init_fun <- function() {
      list(theta=matrix(0,nrow=I2,ncol=Trait))
    }
    
    rstan::rstan_options(auto_write = TRUE,javascript = FALSE)
    
    if (Correlation2 == 0){
      if (j==1){
        data_list<-list(n_student = I2, n_item = S, n_pair = J2, n_dimension = Trait, ind = ind1, p = pairmap,
                        alpha=Itempar[,1], 
                        delta=Itempar[,2], 
                        tau=Itempar[,3],
                        N=Data$N,
                        II=array(Data$ii, dim = c(1)),
                        JJ=array(Data$jj, dim = c(1)),
                        y=array(Data$y, dim = c(1)),
                        theta_mu=as.array(c(rep(0,Trait))),
                        I1=Data$I,
                        J1=Data$J)
      }else{
        data_list<-list(n_student = I2, n_item = S, n_pair = J2, n_dimension = Trait, ind = ind1, p = pairmap,
                        alpha=Itempar[,1], 
                        delta=Itempar[,2], 
                        tau=Itempar[,3],
                        N=Data$N,
                        II=Data$ii,
                        JJ=Data$jj,
                        y=Data$y,
                        theta_mu=as.array(c(rep(0,Trait))),
                        I1=Data$I,
                        J1=Data$J)
      }
    
    fcirt <- rstan::stan(model_code=fcirt_nocov, data=data_list,
                         iter=5000, chains=2, cores=2, warmup=2500,
                         init=init_fun, thin=1,
                         control=list(adapt_delta=0.9,max_treedepth=15))
    }
    if (Correlation2 == 0.30){
      P <- ncol(Covairates)
      PC1<-t(replicate(5, Covairates[i,]))
      
      if (j==1){
        data_list<-list(n_student = I2, n_item = S, n_pair = J2, n_dimension = Trait, ind = ind1, p = pairmap,
                        alpha=Itempar[,1], 
                        delta=Itempar[,2], 
                        tau=Itempar[,3],
                        N=Data$N,
                        II=array(Data$ii, dim = c(1)),
                        JJ=array(Data$jj, dim = c(1)),
                        y=array(Data$y, dim = c(1)),
                        theta_mu=as.array(c(rep(0,Trait))),
                        I1=Data$I,
                        J1=Data$J,
                        P=P,
                        PC1=PC1)
      }else{
        data_list<-list(n_student = I2, n_item = S, n_pair = J2, n_dimension = Trait, ind = ind1, p = pairmap,
                        alpha=Itempar[,1], 
                        delta=Itempar[,2], 
                        tau=Itempar[,3],
                        N=Data$N,
                        II=Data$ii,
                        JJ=Data$jj,
                        y=Data$y,
                        theta_mu=as.array(c(rep(0,Trait))),
                        I1=Data$I,
                        J1=Data$J,
                        P=P,
                        PC1=PC1)
      }
      
      fcirt <- rstan::stan(model_code=fcirt_withcov, data=data_list,
                    iter=5000, chains=2, cores=2, warmup=2500,
                    init=init_fun, thin=1,
                    control=list(adapt_delta=0.9,max_treedepth=15))
    }
    theta1 <- rstan::summary(fcirt, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary[,1]
    theta1 <- matrix(theta1, nrow=Trait)
    theta1 <- t(theta1)
    ThetaEst[i, ] <- theta1
    }
    if (Scoring == 2){
      
      ##########Transform fc response data to separate ss response data for each dimension
      response1 <- vector(mode = "list", length = Trait)
      
      for (e in 1:Trait){
        response1[[e]] <- matrix(NA, 1, ncol(response))
      }
      
      index <- vector(mode = "list", length = Trait)
      for (m in 1:Trait){
        index[[m]] <- 1
      }
      
      for (j in 1:ncol(response)){
        if (response[1,j] == 1){
          response1[[ind[2*j-1]]][1,index[[ind[2*j-1]]]] <- 1
          response1[[ind[2*j]]][1,index[[ind[2*j]]]] <- 0
          if (index[[ind[2*j-1]]]<ncol(response)){
            index[[ind[2*j-1]]]=index[[ind[2*j-1]]]+1
          }
          if (index[[ind[2*j]]]<ncol(response)){
            index[[ind[2*j]]]=index[[ind[2*j]]]+1
          }
        }
        else{
          response1[[ind[2*j-1]]][1,index[[ind[2*j-1]]]] <- 0
          response1[[ind[2*j]]][1,index[[ind[2*j]]]] <- 1
          if (index[[ind[2*j-1]]]<ncol(response)){
            index[[ind[2*j-1]]]=index[[ind[2*j-1]]]+1
          }
          if (index[[ind[2*j]]]<ncol(response)){
            index[[ind[2*j]]]=index[[ind[2*j]]]+1
          }
        }
      }
      
      #########################Re-organize statement parameters for each dimension for MUPPscore script and output
      parameters <- vector(mode = "list", length = Trait)
      inde <- vector(mode = "list", length = Trait)
      
      for (m in 1:Trait){
        inde[[m]] <- 1
      }
      for (m in 1:Trait){
        parameters[[m]] <- matrix(NA, ncol(response)*2, 3)
      }
      
      for (f in 1:Trait){
        for (j in 1:ncol(response)){
          if (ind[2*j-1]==f){
            parameters[[f]][2*inde[[ind[2*j-1]]]-1, ] <- Itempar[2*j-1,]
            parameters[[f]][2*inde[[ind[2*j-1]]], ] <- Itempar[2*j,]
            if (inde[[ind[2*j-1]]]<ncol(response)*2){
              inde[[ind[2*j-1]]]=inde[[ind[2*j-1]]]+1
            }
          }
          if (ind[2*j]==f){
            parameters[[f]][(2*inde[[ind[2*j]]]-1), ] <- Itempar[2*j,]
            parameters[[f]][2*inde[[ind[2*j]]], ] <- Itempar[2*j-1,]
            if (inde[[ind[2*j]]]<ncol(response)*2){
              inde[[ind[2*j]]]=inde[[ind[2*j]]]+1
            }
          }
        }
      }
      
      ###########################################5. MUPPscore script 
      library(LearnBayes)
      N=1 ## Change sample size
      P=29 ## number of quadrature points (ONLY advanced users are allowed to change)
      PS=841 ## square the P (ONLY advanced users are allowed to change)
      
      for (a in 1:Trait){
        
        Parameters <- as.matrix(parameters[[a]])
        Response <- as.matrix(response1[[a]])
        
        if (any(!is.na(Response))==TRUE){
          Response <- t(as.matrix(Response[, !is.na(Response)]))
          PAIRS=ncol(Response) ##Change number of pairs
          NewResponse=Response[rep(1:nrow(Response),each=PS),]
          
          prior=t(t(seq(-3.5,3.5,.25)))
          PriorSForOne=t(t(prior[rep(1:nrow(prior),each=P),]))
          PriorS=PriorSForOne[rep(1:PS,N)]
          
          PriorT=prior[rep(1:P,N*P),]
          PriorTForOne=prior[rep(1:P,P),]
          
          if (Correlation2 == 0.30){
            mu=c(rowSums(Covairates*Correlation2)[i],rowSums(Covairates*Correlation2)[i])             
            Sigma=matrix(c(1,Correlation1,Correlation1,1),2,2)    
          }
          
          if (Correlation2 == 0){
            mu=c(0,0)
            Sigma=matrix(c(1,0,0,1),2,2)
          }
          
          PriorST=matrix(c(PriorS,PriorT),,2)
          ProbS=dnorm(PriorS, mu, Sigma, log=FALSE)
          ProbT=dnorm(PriorT, mu, Sigma, log=FALSE)
          ProbST=dmnorm(PriorST, mu, Sigma, log=FALSE)
          
          OneItemResponse=matrix(,nrow=(PS*N),ncol=1)
          PickedProb=matrix(,nrow=(PS*N),ncol=PAIRS)
          ThetaNumerator=matrix(,nrow=N,ncol=1)
          ThetaDenominator=matrix(,nrow=N,ncol=1)
          Theta1=matrix(,nrow=N,ncol=1)
          StdErrorNumeratorAll=matrix(,nrow=(PS*N),ncol=1)
          StdErrorNumerator=matrix(,nrow=N,ncol=1)
          StdError=matrix(,nrow=N,ncol=1)
          ProbSgreaterTIRF=matrix(,nrow=PS,ncol=1)
          
          for (b in 1:PAIRS){
            AlphaS=Parameters[(2*b-1),1]
            BetaS=Parameters[(2*b-1),2]
            TauS=Parameters[(2*b-1),3]
            AlphaT=Parameters[(2*b),1]
            BetaT=Parameters[(2*b),2]
            TauT=Parameters[(2*b),3]
            
            GammaS=1+exp(AlphaS*(3*(PriorS-BetaS)))+exp(AlphaS*((PriorS-BetaS)-TauS))+exp(AlphaS*(2*(PriorS-BetaS)-TauS))
            GammaT=1+exp(AlphaT*(3*(PriorT-BetaT)))+exp(AlphaT*((PriorT-BetaT)-TauT))+exp(AlphaT*(2*(PriorT-BetaT)-TauT))
            PS1=(exp(AlphaS*((PriorS-BetaS)-TauS))+exp(AlphaS*(2*(PriorS-BetaS)-TauS)))/GammaS
            PS0=(1+exp(AlphaS*(3*(PriorS-BetaS))))/GammaS
            PT1=(exp(AlphaT*((PriorT-BetaT)-TauT))+exp(AlphaT*(2*(PriorT-BetaT)-TauT)))/GammaT
            PT0=(1+exp(AlphaT*(3*(PriorT-BetaT))))/GammaT
            ProbSgreaterT=(PS1*PT0)/(PS1*PT0+PS0*PT1)
            ProbSlessT=(PS0*PT1)/(PS1*PT0+PS0*PT1)
            
            OneItemResponse=as.matrix(NewResponse)[,b]
            
            for (c in 1:(PS*N))
            {
              if(OneItemResponse[c]==1){
                PickedProb[c,b]=ProbSgreaterT[c]}
              if(OneItemResponse[c]==0){
                PickedProb[c,b]=ProbSlessT[c]}
              else if (OneItemResponse[c]==9){
                PickedProb[c,b]=1}
            }
          }
          
          rowProds=function(X){apply(X,1,FUN="prod")}
          LVfA=rowProds(PickedProb) 
          LVf=LVfA*ProbST  
          Af=ProbST  
          NumeratorRaw=t(t(PriorS*LVf*Af))
          DenominatorRaw=t(t(LVf*Af))
          
          for (k in 1:N){
            ThetaNumerator[k]=sum(NumeratorRaw[(PS*(k-1)+1):(PS*(k-1)+PS),])
            ThetaDenominator[k]=sum(DenominatorRaw[(PS*(k-1)+1):(PS*(k-1)+PS),])
            Theta1[k]=ThetaNumerator[k]/ThetaDenominator[k]
          }
          
          Theta1=t(t(Theta1))  #for 1 trait
          ThetaAll=Theta1[rep(1:nrow(Theta1),each=PS),]
          StdErrorNumeratorAll=t(t(Af*LVf*((PriorS-ThetaAll)^2)))
          
          for (m in 1:N){
            StdErrorNumerator[m]=sum(StdErrorNumeratorAll[(PS*(m-1)+1):(PS*(m-1)+PS),])
            ThetaDenominator[m]=sum(DenominatorRaw[(PS*(m-1)+1):(PS*(m-1)+PS),])
            StdError[m]=sqrt(StdErrorNumerator[m]/ThetaDenominator[m])
          }
          ThetaEst[i, a] <- Theta1
          PSD[i, a] <- StdError
        }
      }
    }
  }
  if (Scoring == 1){
    psd <- rstan::summary(fcirt, pars = c("theta"), probs = c(0.025,0.5,0.975))$summary[,3]
    psd <- matrix(psd, nrow=Trait)
    psd <- t(psd)
    PSD[i, ] <- psd
  }
}

#Item and statement exposure rate
rs <- ParVals[,6]/Obs
ri <- Pairings[, 13]/Obs
var(rs)
  
#Measurement precision
(sd <- mean(PSD))

#Scoring accuracy
(abs <- mean(abs(ThetaEst-Theta))) 
(rmse <- mean(rowMeans((ThetaEst-Theta)^2)))
(cor <- mean(diag(cor(ThetaEst,Theta))))

##################Stan 25 items
#cov + ec
#cor: 0.927704
#abs: 0.3189
#rmse: 0.1812603
#sd:0.5510069

#nocov + ec
#cor: 0.850009
#abs: 0.3916787
#rmse: 0.2812734
#sd:0.5268665

#cov + noec
#cor: 0.9001863
#abs: 0.398274
#rmse: 0.2727045
#sd:0.5772867

#nocov + noec
#cor: 0.8897288
#abs: 0.3681444
#rmse: 0.201461
#sd:0.4506281
  
##################EAP 25 items 100 person
##################Preliminary conclusion: (1) Having exposure control balanced statement exposure rates; 
##################(2) Having cov improved scoring accuracy
##################(3) Having no ec improved scoring accuracy slightly
##################(4) measurement precision?

#expected: higher scoring accuracy, more balanced exposure
#cov + ec  
#cor: 0.8445521
#abs: 0.4078244
#rmse: 0.2943296
#sd:0.437747
#rs
# [1] 0.50 0.39 0.90 0.69 0.85 0.34 0.31 0.91 0.63 0.86 0.55 0.76 0.54 0.37 0.55 0.61 0.39 0.90 0.49 0.46 0.32 0.31 0.75 0.44
# [25] 0.69 0.35 0.26 0.90 0.36 0.90 0.49 0.70 0.34 0.36 0.38 0.48 0.33 0.86 0.35 0.43 0.38 0.42 0.78 0.57 0.81 0.29 0.29 0.81
# [49] 0.39 0.79 0.42 0.71 0.24 0.24 0.49 0.52 0.40 0.81 0.37 0.27 0.36 0.31 0.75 0.47 0.77 0.18 0.21 0.77 0.31 0.82 0.29 0.72
# [73] 0.28 0.26 0.40 0.48 0.29 0.71 0.32 0.30 0.33 0.28 0.85 0.41 0.71 0.17 0.23 0.82 0.39 0.81 0.30 0.67 0.27 0.31 0.40 0.47
# [97] 0.27 0.76 0.25 0.30
#var(rs): 0.04746263

#expected: lowest scoring accuracy, more balanced exposure
#nocov + ec
#cor: 0.7981087
#abs: 0.4740792
#rmse: 0.3934522
#sd:0.4373379
#rs
# [1] 0.39 0.47 0.87 0.67 0.91 0.40 0.38 0.91 0.57 0.87 0.50 0.72 0.48 0.45 0.47 0.72 0.49 0.90 0.50 0.33 0.19 0.23 0.76 0.44
# [25] 0.82 0.30 0.29 0.86 0.46 0.87 0.39 0.72 0.42 0.36 0.40 0.44 0.42 0.84 0.39 0.40 0.27 0.36 0.76 0.57 0.85 0.21 0.26 0.88
# [49] 0.36 0.76 0.34 0.76 0.41 0.29 0.51 0.52 0.42 0.80 0.31 0.36 0.20 0.29 0.79 0.53 0.88 0.25 0.20 0.86 0.37 0.77 0.24 0.53
# [73] 0.30 0.15 0.38 0.55 0.38 0.76 0.32 0.25 0.32 0.26 0.81 0.52 0.81 0.19 0.22 0.87 0.28 0.73 0.34 0.63 0.16 0.23 0.38 0.48
# [97] 0.30 0.81 0.37 0.29
#var(rs): 0.0519899


#expected: highest scoring accuracy, less balanced exposure
#cov + noec
#cor: 0.8535807
#abs: 0.417944
#rmse: 0.2885809
#sd:0.423443
#rs
# [1] 0.84 0.05 0.98 0.97 1.00 0.14 0.54 1.00 0.70 0.96 0.55 0.95 0.00 0.14 0.90 0.94 0.04 0.99 0.00 0.31 0.58 0.00 0.96 0.95 0.98 0.01
# [27] 0.24 1.00 0.51 0.93 0.24 0.87 0.00 0.01 0.84 0.88 0.01 0.95 0.00 0.04 0.58 0.00 0.99 0.96 1.00 0.00 0.24 1.00 0.55 0.87 0.21 0.87
# [53] 0.00 0.00 0.81 0.87 0.00 0.97 0.00 0.08 0.32 0.00 0.98 0.91 1.00 0.00 0.09 0.99 0.42 0.90 0.14 0.86 0.00 0.00 0.63 0.82 0.00 0.92
# [79] 0.00 0.02 0.27 0.00 0.97 0.86 0.98 0.00 0.11 0.97 0.52 0.90 0.16 0.80 0.00 0.00 0.68 0.81 0.00 0.96 0.00 0.01
#var(rs): 0.1751596


#expected: lower scoring accuracy, less balanced exposure
#nocov + noec
#cor: 0.8228474
#abs: 0.4412997
#rmse: 0.3545612
#sd:0.4107955
#rs
# [1] 0.94 0.11 0.99 1.00 0.99 0.05 0.46 1.00 0.63 0.94 0.52 0.98 0.01 0.16 0.92 0.94 0.04 0.98 0.01 0.33 0.51 0.00 1.00 0.94 1.00 0.01
# [27] 0.24 1.00 0.56 0.88 0.25 0.86 0.00 0.01 0.75 0.89 0.01 0.96 0.00 0.13 0.63 0.00 0.99 0.93 0.98 0.00 0.28 1.00 0.48 0.84 0.21 0.91
# [53] 0.00 0.00 0.80 0.90 0.00 0.97 0.00 0.08 0.29 0.00 1.00 0.90 1.00 0.00 0.03 0.95 0.49 0.89 0.12 0.82 0.00 0.00 0.65 0.84 0.00 1.00
# [79] 0.00 0.02 0.26 0.00 0.97 0.91 0.97 0.00 0.12 0.99 0.38 0.82 0.19 0.80 0.00 0.01 0.70 0.89 0.00 0.95 0.00 0.04
#var(rs): 0.1755071
  
  

