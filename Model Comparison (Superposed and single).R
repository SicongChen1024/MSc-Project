#Packages and functions
library(hawkes)
##Poisson process simulation function
sim_pp1 <- function(t, rate) { 
  
  
  
  jumps_time <- rexp(1, rate)
  
  while(jumps_time[length(jumps_time)] < t) {
    
    
    
    jumps_time <- c(jumps_time,
                    jumps_time[length(jumps_time)] + rexp(1, rate))
  }
  
  
  jumps_time[-(length(jumps_time))]
}
##Single Hawkes likelihood function
loglik_pure_hawkes<-function(params,arrivals){
  mu_i <- params[1]
  alpha_i <- params[2]
  beta_i <- params[3]
  
  n<- length(arrivals)
  term_1<- -mu_i*Total_Time_Interval
  term_2<- sum((alpha_i/beta_i)*(exp(-beta_i*(Total_Time_Interval-arrivals))-1))
  A_i <- c(0,sapply(2:n,function(z){
    sum(exp(-beta_i*(arrivals[z]-arrivals[1:(z-1)])))
  }))
  term_3 <- sum(log(mu_i + alpha_i* A_i))
  return(-term_1-term_2-term_3)
}
###Superposed hawkes likelihood function
loglik_mixture_hawkes<-function(params,arrivals,Label){
  mu_i <- params[1]
  alpha_i <- params[2]
  beta_i <- params[3]
  
  n<- length(arrivals)
  term_1<- -mu_i*Total_Time_Interval
  term_2<- sum((alpha_i/beta_i)*(exp(-beta_i*(Total_Time_Interval-arrivals))-1)*(1-Label))
  A_i <- c(0,sapply(2:n,function(z){
    sum(exp(-beta_i*(arrivals[z]-arrivals[1:(z-1)]))*(1-Label[1:(z-1)]))
  }))
  term_3 <- sum(log(mu_i + alpha_i* A_i)*(1-Label))
  return(-term_1-term_2-term_3)
}

##Conditional likelihood function
f_t_given_z<-function(t,z,lambda,mu,alpha,beta,observations){
  if (z==1){
    Poisson_obs<- observations[observations[,2]==1,1]
    obs_less_tp<- Poisson_obs[Poisson_obs<t]
    if (length(obs_less_tp)==0){
      tp<-0
    }else{
      tp<- max(obs_less_tp)
    }
    ft<- lambda*exp(-lambda*(t-tp))
  }else{
    Hawkes_obs<- observations[observations[,2]==0,1]
    obs_less_th<- Hawkes_obs[Hawkes_obs<t]
    if (length(obs_less_th)==0){
      ft<- mu*exp(-mu*t)
    }else{
      th<- max(obs_less_th)
      A_i<- sum(exp(-beta*(t-obs_less_th)))
      term_1<- mu + alpha*A_i
      sum_i<- sum(exp(-beta*(t-obs_less_th))-exp(-beta*(th-obs_less_th)))
      term_2<- exp(-mu*(t-th)+(alpha/beta)*sum_i)
      ft<- term_1 * term_2
    }
    
  }
  ft
}
##Q_theta of Hawkes observation function
Q_theta<-function(params,obs,label_list,w_k_vec,n_sample){
  mu_n<- params[1]
  alpha_n<- params[2]
  beta_n<- params[3]
  
  T_Z_list<-c()
  for (i in 1:n_sample){
    complete_obs_i<- data.frame(t=obs,Label=label_list[i,])
    Hawkes_obs<- complete_obs_i[complete_obs_i[,2]==0,1]
    n<- length(Hawkes_obs)
    term_1<- -mu_n * Total_Time_Interval
    term_2<- sum((alpha_n/beta_n)*(exp(-beta_n*(Total_Time_Interval-Hawkes_obs))-1))
    
    A_i <- c(0,sapply(2:n,function(z){
      if (z==2){
        A_i_vals<<- exp(-beta_n*(Hawkes_obs[z]-Hawkes_obs[z-1]))
      }else{
        A_i_vals<<- exp(-beta_n *(Hawkes_obs[z]-Hawkes_obs[z-1]))*(1+A_i_vals)
      }
      A_i_vals
    }))
    term_to_log<- mu_n + alpha_n* A_i
    
    term_3 <- sum(log(term_to_log))
    sum_prob<- term_1 + term_2 +term_3
    
    T_Z_list<- c(T_Z_list,sum_prob)
    
  }
  Q<- (sum(w_k_vec * T_Z_list))/(sum(w_k_vec))
  return(-Q)
}
##Hawkes intensity function
intensity_hawkes<- function(t,obshawkes,mu,alpha,beta){
  obs_less_t<- obshawkes[obshawkes<t]
  if (length(obs_less_t) == 0){
    return(mu)
  }else{
    term_list<- alpha * exp(-beta * (t-obs_less_t))
    return(mu + sum(term_list))
  }
  
}
##wk lies between Q1 and Q3 function
IQR.not.outliers <- function(x) {
  if(any(is.na(x)))
    stop("x is missing values")
  if(!is.numeric(x))
    stop("x is not numeric")
  Q3<-quantile(x,0.75)
  Q1<-quantile(x,0.25)
  
  x[x>Q1 & x<Q3]
}
##E step function
e_step<- function(observed_t,lambda.n,mu.n,alpha.n,beta.n,pi_vec){
  ###Initialization
  w_k_list<- c()
  sample_label_list<- c()
  loglik_z_list<- c()
  ###Conducting the simulation algorithm 100 times (m=50)
  for (i in 1:100){
    
    
    q_sample_obs<- data.frame(t=c(),Label=c())
    joint_likelihood_z<- 1
    z_label_list<- c()
    loglik_z_val<- 0
    for (j in 1:nm){
      tj<- observed_t[j]
      ###Alternative distribution
      poisson_int<- lambda.n
      hawkes_index<- which(z_label_list == 0)
      hawkes_obs<- observed_t[hawkes_index]
      hawkes_int<- intensity_hawkes(t=tj, obshawkes = hawkes_obs, mu= mu.n,alpha= alpha.n,beta= beta.n)
      denominator_prob<- poisson_int + hawkes_int
      poisson_prob<- poisson_int/denominator_prob
      hawkes_prob<- hawkes_int/denominator_prob
      zj<- sample(c(0,1),1,replace = TRUE, prob = c(hawkes_prob,poisson_prob))
      z_label_list<- c(z_label_list, zj)
      complete_jth_obs<- data.frame(t=tj, Label= zj)
      ###True probability density
      q_sample_obs<- rbind(q_sample_obs, complete_jth_obs)
      ftz_0<- f_t_given_z(t=tj,z=0,lambda = lambda.n,mu=mu.n,alpha=alpha.n,beta=beta.n,observations=q_sample_obs) * pi_vec[2]
      ftz_1<- f_t_given_z(t=tj,z=1,lambda = lambda.n,mu=mu.n,alpha=alpha.n,beta=beta.n,observations=q_sample_obs) * pi_vec[1]
      
      if (zj==1){
        
        ftz_tj<- ftz_1
        q_sample_prob<- poisson_prob
      }else{
        
        ftz_tj<- ftz_0
        q_sample_prob<- hawkes_prob
      }
      
      
      conditional_p_z<- (ftz_tj)/((ftz_1 + ftz_0)*(q_sample_prob))
      loglik_z_val<- loglik_z_val + log((ftz_tj)/(ftz_1 + ftz_0))
      joint_likelihood_z<- joint_likelihood_z * conditional_p_z
      
    }
    sample_label_list<- rbind(sample_label_list,z_label_list)
    w_k_list<- c(w_k_list,joint_likelihood_z)
    loglik_z_list<- c(loglik_z_list, loglik_z_val)
  }
  ###Find samples lie between Q1 and Q3
  wk_without_outliers<- IQR.not.outliers(w_k_list)
  index_without_outliers<- sapply(1:length(wk_without_outliers) ,function(z){
    which(w_k_list == wk_without_outliers[z])
  })
  simulate_q_new<- sample_label_list[index_without_outliers,]
  loglik_zlist_new<- loglik_z_list[index_without_outliers]
  poisson_number_list<-c()
  nwk<- length(wk_without_outliers)
  for(i in 1: nwk){
    number_of_poisson<- sum(simulate_q_new[i,])
    poisson_number_list<- c(poisson_number_list,number_of_poisson)
  }
  hawkes_number_list<- nm - poisson_number_list
  ###Compute total likelihood
  likelihoodk<- log(pi_vec[1]) * poisson_number_list + poisson_number_list * log(lambda.n) - lambda.n * Total_Time_Interval + hawkes_number_list * log(pi_vec[2])
  likelihood_total<- sum(wk_without_outliers * likelihoodk)/sum(wk_without_outliers)
  likelihood_z_total<- sum(wk_without_outliers * loglik_zlist_new)/sum(wk_without_outliers)
  Q_theta_star<- likelihood_total - Q_theta(params = c(mu.n,alpha.n,beta.n),obs = observed_t,label_list = simulate_q_new,w_k_vec = wk_without_outliers,n_sample = nwk) - likelihood_z_total
  list("w_k_vals" = wk_without_outliers,
       "simulate_q_samples" = simulate_q_new,
       "Q" = Q_theta_star)
}
###M step function
m_step<- function(observed_time_vals, w_k_vals, simulate_q_samples){
  poisson_number_list<-c()
  n_w_k<- length(w_k_vals)
  for(i in 1: n_w_k){
    number_of_poisson<- sum(simulate_q_samples[i,])
    poisson_number_list<- c(poisson_number_list,number_of_poisson)
  }
  hawkes_number_list<- nm - poisson_number_list
  denominator_pi<- sum(w_k_vals) * nm
  pi_1_mle<- (sum(w_k_vals * poisson_number_list))/(denominator_pi)
  pi_0_mle<- (sum(w_k_vals * hawkes_number_list))/(denominator_pi)
  lambda_mle<- (sum(w_k_vals * poisson_number_list))/(sum(w_k_vals * 60))
  hawkes_mle<- optim(c(6,3,5), Q_theta, obs= observed_time_vals,label_list=simulate_q_samples,w_k_vec= w_k_vals,n_sample = n_w_k)
  mu_mle<- hawkes_mle$par[1]
  alpha_mle<- hawkes_mle$par[2]
  beta_mle<- hawkes_mle$par[3]
  params<- c(lambda_mle,mu_mle,alpha_mle,beta_mle,pi_1_mle,pi_0_mle)
  
  list("lambda" = lambda_mle,
       "mu" = mu_mle,
       "alpha" = alpha_mle,
       "beta" = beta_mle,
       "pi" = c(pi_1_mle, pi_0_mle),
       "params" = params)
}

#Model Comparison
pi_1_list<- c()
pi_0_list<- c()
pure_poisson_AIClist<- c()
pure_hawkes_AIClist<- c()
lambda_knownlabels_list<- c()
mu_knownlabels_list<- c()
alpha_knownlabels_list<- c()
beta_knownlabels_list<- c()
mixture_knownlabels_AIClist<-c()
lambda_MCEM_list<- c()
mu_MCEM_list<- c()
alpha_MCEM_list<- c()
beta_MCEM_list<- c()
pi1_MCEM_list<- c()
pi0_MCEM_list<- c()
mixture_MCEM_AIClist<- c()
for (z in 1:2){
  Total_Time_Interval<- 60
  True_lambda<- 9
  True_mu<- 6
  True_alpha<- 3
  True_beta<- 5
  
  #Poisson
  Model_Poisson<- sim_pp1(t=Total_Time_Interval,True_lambda)
  np<- length(Model_Poisson)
  labelp<- replicate(np,1)
  Model_PoissonwithLabels<- data.frame(t=Model_Poisson,Label=labelp)
  #Hawkes
  
  Model_Hawkes<- simulateHawkes(lambda0 = True_mu,alpha = True_alpha,beta = True_beta,horizon = Total_Time_Interval)[[1]]
  nh<- length(Model_Hawkes)
  labelh<- replicate(nh,0)
  Model_HawkeswithLabels<- data.frame(t=Model_Hawkes,Label=labelh)
  #Mixture
  Mixture_sim<- rbind(Model_PoissonwithLabels,Model_HawkeswithLabels)
  Mixture_observations<- Mixture_sim[order(Mixture_sim$t),]
  nm<- np+nh
  pi_poisson<- np/nm
  pi_1_list<- c(pi_1_list,pi_poisson)
  pi_hawkes<- nh/nm
  pi_0_list<- c(pi_0_list,pi_hawkes)
  ###Mixture of Poisson and Hawkes with Konwn Labels
  MLE_mixture_hawkes<- optim(c(1.5,4,5),loglik_mixture_hawkes, arrivals=Mixture_observations$t,Label=Mixture_observations$Label)
  mu_knownlabels_list<- c(mu_knownlabels_list,MLE_mixture_hawkes$par[1])
  alpha_knownlabels_list<- c(alpha_knownlabels_list,MLE_mixture_hawkes$par[2])
  beta_knownlabels_list<- c(beta_knownlabels_list,MLE_mixture_hawkes$par[3])
  number_poisson_mixture<- sum(Mixture_observations$Label)
  MLE_mixture_lambda<- number_poisson_mixture/Total_Time_Interval
  lambda_knownlabels_list<- c(lambda_knownlabels_list,MLE_mixture_lambda)
  
  loglik_val<- 0
  for (j in 1:nm){
    tj<- Mixture_observations[j,1]
    ftz_0<- f_t_given_z(t=tj,z=0,lambda = MLE_mixture_lambda,mu=MLE_mixture_hawkes$par[1],alpha=MLE_mixture_hawkes$par[2],beta=MLE_mixture_hawkes$par[3],observations=Mixture_observations) * pi_hawkes
    ftz_1<- f_t_given_z(t=tj,z=1,lambda = MLE_mixture_lambda,mu=MLE_mixture_hawkes$par[1],alpha=MLE_mixture_hawkes$par[2],beta=MLE_mixture_hawkes$par[3],observations=Mixture_observations) * pi_poisson
    if (Mixture_observations[j,2]==1){
      
      ftz_tj<- ftz_1
      
    }else{
      
      ftz_tj<- ftz_0
      
    }
    loglik_val<- loglik_val + log((ftz_tj)/(ftz_0+ftz_1))
    
  }
  loglikelihoodmixture<- number_poisson_mixture * log(MLE_mixture_lambda) - MLE_mixture_lambda * Total_Time_Interval - MLE_mixture_hawkes$value + np * log(pi_poisson) + nh * log(pi_hawkes) - loglik_val
  AIC_mixture<- 2 * 4 - 2 * loglikelihoodmixture
  mixture_knownlabels_AIClist<-c(mixture_knownlabels_AIClist,AIC_mixture)
  ###MCEM
  for (i in 1:30 ){
    if (i ==1){
      e.step<- e_step(observed_t = Mixture_observations[,1], lambda.n = True_lambda, mu.n = True_mu, alpha.n = True_alpha, beta.n = True_beta, pi_vec = c(pi_poisson,pi_hawkes))
      m.step<- m_step(observed_time_vals = Mixture_observations[,1], w_k_vals = e.step[["w_k_vals"]], simulate_q_samples = e.step[["simulate_q_samples"]])
      cur.params<- m.step[["params"]]
    }else{
      e.step<- e_step(observed_t = Mixture_observations[,1], lambda.n = m.step[["lambda"]], mu.n = m.step[["mu"]], alpha.n = m.step[["alpha"]], beta.n = m.step[["beta"]], pi_vec = m.step[["pi"]])
      m.step<- m_step(observed_time_vals = Mixture_observations[,1], w_k_vals = e.step[["w_k_vals"]], simulate_q_samples = e.step[["simulate_q_samples"]])
      params_dif<- cur.params - m.step[["params"]]
      norm_dif<- sqrt(sum(params_dif^(2)))
      if (norm_dif < 0.01){
        break
      }else{
        cur.params<- m.step[["params"]]
      }
    }
  }
  lambda_MCEM_list<- c(lambda_MCEM_list,m.step[["lambda"]])
  mu_MCEM_list<- c(mu_MCEM_list,m.step[["mu"]])
  alpha_MCEM_list<- c(alpha_MCEM_list,m.step[["alpha"]])
  beta_MCEM_list<- c(beta_MCEM_list,m.step[["beta"]])
  pi1_MCEM_list<- c(pi1_MCEM_list,m.step[["pi"]][1])
  pi0_MCEM_list<- c(pi0_MCEM_list,m.step[["pi"]][2])
  e.step<- e_step(observed_t = Mixture_observations[,1], lambda.n = m.step[["lambda"]], mu.n = m.step[["mu"]], alpha.n = m.step[["alpha"]], beta.n = m.step[["beta"]], pi_vec = m.step[["pi"]])
  AIC_MCEM<- 2* 6 - 2 * e.step[["Q"]]
  mixture_MCEM_AIClist<- c(mixture_MCEM_AIClist,AIC_MCEM)
  
  ###Single Poisson
  MLE_pure_lambda<- nm/ Total_Time_Interval
  
  loglikelihood_pure_poisson<- nm * log(MLE_pure_lambda) - MLE_pure_lambda * Total_Time_Interval + nm * log(m.step[["pi"]][1])
  AIC_pure_poisson<- 2 * 1 - 2 * loglikelihood_pure_poisson
  pure_poisson_AIClist<- c(pure_poisson_AIClist,AIC_pure_poisson)
  ###Single Hawkes
  MLE_pure_hawkes<- optim(c(6,3,5),loglik_pure_hawkes, arrivals=Mixture_observations$t)
  AIC_pure_hawkes<- 2 * 3 - 2*(-(MLE_pure_hawkes$value) + nm * log(m.step[["pi"]][2]))
  pure_hawkes_AIClist<- c(pure_hawkes_AIClist,AIC_pure_hawkes)
  
}

output_results<- data.frame(AIC_poisson = pure_poisson_AIClist, AIC_hawkes = pure_hawkes_AIClist, 
                            AIC_Mixture_known = mixture_knownlabels_AIClist, AIC_Mixture_MCEM = mixture_MCEM_AIClist, 
                            true_pi_0 = pi_0_list, true_pi_1 = pi_1_list, 
                            lambda_mle = lambda_knownlabels_list, mu_mle = mu_knownlabels_list, 
                            alpha_mle = alpha_knownlabels_list, beta_mle = beta_knownlabels_list, 
                            lambda_MCEM = lambda_MCEM_list, mu_MCEM = mu_MCEM_list, 
                            alpha_MCEM = alpha_MCEM_list, beta_MCEM = beta_MCEM_list, 
                            pi_0_MCEM = pi0_MCEM_list, pi_1_MCEM = pi1_MCEM_list)

print(output_results)