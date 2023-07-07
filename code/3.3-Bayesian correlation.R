setwd('YOUR_PATH')


rob.cor.mcmc = function(x, iter = 2000, warmup = 500, chains = 4,treedepth=10) {
  
  library(rstan)
  library(coda)
  library(gridExtra)
  stopifnot(ncol(x) == 2)
  
  # Set up model data
  model.data = list(N=nrow(x), x=x)
  
  # Stan model definition
  stan.model = "
        data {
            int<lower=1> N;  // number of observations
            vector[2] x[N];  // input data: rows are observations, columns are the two variables
        }
        
        parameters {
            vector[2] mu;                 // locations of the marginal t distributions
            real<lower=0> sigma[2];       // scales of the marginal t distributions
            real<lower=1> nu;             // degrees of freedom of the marginal t distributions
            real<lower=-1, upper=1> rho;  // correlation coefficient
        }
        
        transformed parameters {
            // Covariance matrix
            cov_matrix[2] cov = [[      sigma[1] ^ 2       , sigma[1] * sigma[2] * rho],
                                 [sigma[1] * sigma[2] * rho,       sigma[2] ^ 2       ]];
        }
        
        model {
            // Likelihood
            // Bivariate Student's t-distribution instead of normal for robustness
            x ~ multi_student_t(nu, mu, cov);
            
            // Noninformative priors on all parameters
            sigma ~ normal(0, 1000);
            mu ~ normal(0, 1000);
            nu ~ gamma(2, 0.1);
        }
        
        generated quantities {
            // Random samples from the estimated bivariate t-distribution (for assessment of fit)
            vector[2] x_rand;
            x_rand = multi_student_t_rng(nu, mu, cov);
        }"
  
  # Run the model
  stan.cor = stan(model_name="robust_correlation",
                  model_code=stan.model, data=model.data,
                  seed=21, iter=iter, warmup=warmup, chains=chains,
                  control =list(max_treedepth=treedepth))
  
  # Obtain the MCMC samples of rho
  stan.rho = extract(stan.cor, "rho")[[1]]
  hpd95 = HPDinterval(as.mcmc(as.numeric(stan.rho)), prob=0.95)
  hpd99 = HPDinterval(as.mcmc(as.numeric(stan.rho)), prob=0.99)
  
  # Write some descriptive statistics
  cat("POSTERIOR STATISTICS OF RHO\n",
      "Posterior mean and standard deviation:     Mean = ",
      mean(stan.rho), ", SD = ", sd(stan.rho), "\n",
      "Posterior median and MAD:                  Median = ",
      median(stan.rho), ", MAD = ", mad(stan.rho), "\n",
      "Rho values with 99% posterior probability: 99% HPDI = [", 
      hpd99[,"lower"], ", ", hpd99[,"upper"], "]\n",
      "Rho values with 95% posterior probability: 95% HPDI = [", 
      hpd95[,"lower"], ", ", hpd95[,"upper"], "]\n",
      "Posterior probability that rho is ≤0:      P(rho ≤ 0) = ",
      mean(stan.rho <= 0), "\n",
      "Posterior probability that rho is ≥0:      P(rho ≥ 0) = ", 
      mean(stan.rho >= 0), "\n",
      "Posterior probability that rho is weak:    P(-0.1 < rho < 0.1) = ", 
      mean(stan.rho > -0.1 & stan.rho < 0.1), "\n\n",
      sep="")
  l1 <- list(rho=mean(stan.rho),p_overzero=mean(stan.rho >= 0),p_weak=mean(stan.rho > -0.1 & stan.rho < 0.1))
  # Return stanfit object
  return(l1)
  
}

Dat = readRDS('VE_jq.rds')
subinfo = readRDS('VE_info_jq.rds')

result <- data.frame(features=colnames(Dat),rho1=NA,p1_strong=NA,p1_overzero=NA,
                     rho2=NA,p2_strong=NA,p2_overzero=NA)
fea <- result$features

for(i in 1:length(fea)){
  print(paste0('-------------------i = ',i,'---------------------'))
  subinfo$value2 <- Dat[,fea[i]] 
  v1 <- subinfo$value[subinfo$phase2==0]
  v2 <- subinfo$value[subinfo$phase2==1]
  v3 <- subinfo$value2[subinfo$phase2==0]
  v4 <- subinfo$value2[subinfo$phase2==1] #cbw
  
  wp <- rob.cor.mcmc(data.frame(LDL=v1,Strain=v3),iter = 2000,warmup = 1000,chain=4,treedepth = 10)
  result$rho1[i] <- wp$rho
  result$p1_strong[i] <- 1-wp$p_weak
  result$p1_overzero[i] <- wp$p_overzero
  
  wp <- rob.cor.mcmc(data.frame(LDL=v2,Strain=v4),iter = 2000,warmup = 1000,chain=4,treedepth = 10)
  result$rho2[i] <- wp$rho
  result$p2_strong[i] <- 1-wp$p_weak
  result$p2_overzero[i] <- wp$p_overzero
  
}
write.csv(result,'VE_LDL_strain_rho_BC.csv',quote = F,row.names = F)


