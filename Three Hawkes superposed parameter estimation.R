#Packages and functions
library(hawkes)
##Conditional likelihood function
f_t_given_z_three_hawkes<-function(t,z,mu0,alpha0,beta0,mu1,alpha1,beta1,mu2,alpha2,beta2,observations){
  if (z==0){
    Hawkes_obs_0<- observations[observations[,2]==0,1]
    obs_less_th0<- Hawkes_obs_0[Hawkes_obs_0<t]
    if (length(obs_less_th0)==0){
      ft<- mu0*exp(-mu0*t)
    }else{
      th0<- max(obs_less_th0)
      A_i0<- sum(exp(-beta0*(t-obs_less_th0)))
      term_10<- mu0 + alpha0*A_i0
      sum_i0<- sum(exp(-beta0*(t-obs_less_th0))-exp(-beta0*(th0-obs_less_th0)))
      term_20<- exp(-mu0*(t-th0)+(alpha0/beta0)*sum_i0)
      ft<- term_10 * term_20
    }
  }else if (z==1){
    Hawkes_obs_1<- observations[observations[,2]==1,1]
    obs_less_th1<- Hawkes_obs_1[Hawkes_obs_1<t]
    if (length(obs_less_th1)==0){
      ft<- mu1*exp(-mu1*t)
    }else{
      th1<- max(obs_less_th1)
      A_i1<- sum(exp(-beta1*(t-obs_less_th1)))
      term_11<- mu1 + alpha1*A_i1
      sum_i1<- sum(exp(-beta1*(t-obs_less_th1))-exp(-beta1*(th1-obs_less_th1)))
      term_21<- exp(-mu1*(t-th1)+(alpha1/beta1)*sum_i1)
      ft<- term_11 * term_21
    }
  }else{
    Hawkes_obs_2<- observations[observations[,2]==2,1]
    obs_less_th2<- Hawkes_obs_2[Hawkes_obs_2<t]
    if (length(obs_less_th2)==0){
      ft<- mu2*exp(-mu2*t)
    }else{
      th2<- max(obs_less_th2)
      A_i2<- sum(exp(-beta2*(t-obs_less_th2)))
      term_12<- mu2 + alpha2*A_i2
      sum_i2<- sum(exp(-beta2*(t-obs_less_th2))-exp(-beta2*(th2-obs_less_th2)))
      term_22<- exp(-mu2*(t-th2)+(alpha2/beta2)*sum_i2)
      ft<- term_12 * term_22
    }
    
  }
  ft
}
##Q_theta function
Q_theta<-function(params,indexh,obs,label_list,w_k_vec,n_sample){
  mu_n<- params[1]
  alpha_n<- params[2]
  beta_n<- params[3]
  
  T_Z_list<-c()
  for (i in 1:n_sample){
    complete_obs_i<- data.frame(t=obs,Label=label_list[i,])
    Hawkes_obs<- complete_obs_i[complete_obs_i[,2]==indexh,1]
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
###E step function
e_step_threehawkes<- function(observed_t,mu.n0,alpha.n0,beta.n0,mu.n1,alpha.n1,beta.n1,mu.n2,alpha.n2,beta.n2,pi_vec){
  ###Initialization
  w_k_list<- c()
  sample_label_list<- c()
  loglik_z_list<- c()
  for (i in 1:100){
    
    
    q_sample_obs<- data.frame(t=c(),Label=c())
    joint_likelihood_z<- 1
    z_label_list<- c()
    loglik_z_val<- 0
    for (j in 1:nm){
      tj<- observed_t[j]
      ###Alternative distribution
      hawkes_index0<- which(z_label_list == 0)
      hawkes_obs0<- observed_t[hawkes_index0]
      hawkes_int0<- intensity_hawkes(t=tj, obshawkes = hawkes_obs0, mu= mu.n0,alpha= alpha.n0,beta= beta.n0)
      hawkes_index1<- which(z_label_list == 1)
      hawkes_obs1<- observed_t[hawkes_index1]
      hawkes_int1<- intensity_hawkes(t=tj, obshawkes = hawkes_obs1, mu= mu.n1,alpha= alpha.n1,beta= beta.n1)
      
      hawkes_index2<- which(z_label_list == 2)
      hawkes_obs2<- observed_t[hawkes_index2]
      hawkes_int2<- intensity_hawkes(t=tj, obshawkes = hawkes_obs2, mu= mu.n2,alpha= alpha.n2,beta= beta.n2)
      
      denominator_prob<- hawkes_int0 + hawkes_int1 + hawkes_int2
      
      hawkes_prob0<- hawkes_int0/denominator_prob
      hawkes_prob1<- hawkes_int1/denominator_prob
      hawkes_prob2<- hawkes_int2/denominator_prob
      
      
      zj<- sample(c(0,1,2),1,replace = TRUE, prob = c(hawkes_prob0,hawkes_prob1,hawkes_prob2))
      z_label_list<- c(z_label_list, zj)
      complete_jth_obs<- data.frame(t=tj, Label= zj)
      ###True probability density
      q_sample_obs<- rbind(q_sample_obs, complete_jth_obs)
      
      ftz_0<- f_t_given_z_three_hawkes(t = tj,z = 0,mu0 = mu.n0,alpha0 = alpha.n0,beta0 = beta.n0,mu1 = mu.n1,alpha1 = alpha.n1,beta1 = beta.n1,mu2 = mu.n2,alpha2 = alpha.n2,beta2 = beta.n2,observations = q_sample_obs) * pi_vec[1]
      ftz_1<- f_t_given_z_three_hawkes(t = tj,z = 1,mu0 = mu.n0,alpha0 = alpha.n0,beta0 = beta.n0,mu1 = mu.n1,alpha1 = alpha.n1,beta1 = beta.n1,mu2 = mu.n2,alpha2 = alpha.n2,beta2 = beta.n2,observations = q_sample_obs) * pi_vec[2]
      ftz_2<- f_t_given_z_three_hawkes(t = tj,z = 2,mu0 = mu.n0,alpha0 = alpha.n0,beta0 = beta.n0,mu1 = mu.n1,alpha1 = alpha.n1,beta1 = beta.n1,mu2 = mu.n2,alpha2 = alpha.n2,beta2 = beta.n2,observations = q_sample_obs) * pi_vec[3]
      
      
      
      
      if (zj==0){
        
        ftz_tj<- ftz_0
        q_sample_prob<- hawkes_prob0
      }else if (zj == 1){
        
        ftz_tj<- ftz_1
        q_sample_prob<- hawkes_prob1
      }else{
        ftz_tj<- ftz_2
        q_sample_prob<- hawkes_prob2
      }
      
      
      conditional_p_z<- (ftz_tj)/((ftz_0 + ftz_1 + ftz_2)*(q_sample_prob))
      loglik_z_val<- loglik_z_val + log((ftz_tj)/(ftz_0 + ftz_1 + ftz_2))
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
  hawkes_number_list0<- c()
  hawkes_number_list1<-c()
  nwk<- length(wk_without_outliers)
  for(i in 1: nwk){
    number_of_hawkes0<- sum(simulate_q_new[i,]==0)
    number_of_hawkes1<- sum(simulate_q_new[i,]==1)
    hawkes_number_list0<- c(hawkes_number_list0, number_of_hawkes0)
    hawkes_number_list1<- c(hawkes_number_list1,number_of_hawkes1)
  }
  hawkes_number_list2<- nm - hawkes_number_list0 - hawkes_number_list1
  ###Compute total likelihood
  likelihoodk<- log(pi_vec[1]) * hawkes_number_list0 +  hawkes_number_list1 * log(pi_vec[2]) + hawkes_number_list2 * log(pi_vec[3])
  likelihood_total<- sum(wk_without_outliers * likelihoodk)/sum(wk_without_outliers)
  likelihood_z_total<- sum(wk_without_outliers * loglik_zlist_new)/sum(wk_without_outliers)
  Q_theta_star<- likelihood_total - Q_theta(params = c(mu.n0,alpha.n0,beta.n0),indexh = 0, obs = observed_t,label_list = simulate_q_new,w_k_vec = wk_without_outliers,n_sample = nwk) - Q_theta(params = c(mu.n1,alpha.n1,beta.n1),indexh = 1, obs = observed_t,label_list = simulate_q_new,w_k_vec = wk_without_outliers,n_sample = nwk) - Q_theta(params = c(mu.n2,alpha.n2,beta.n2),indexh = 2, obs = observed_t,label_list = simulate_q_new,w_k_vec = wk_without_outliers,n_sample = nwk) - likelihood_z_total
  
  list("w_k_vals" = wk_without_outliers,
       "simulate_q_samples" = simulate_q_new,
       "Q" = Q_theta_star)
}
##M step function
m_step_threehawkes<- function(observed_time_vals, w_k_vals, simulate_q_samples){
  hawkes_number_list0<- c()
  hawkes_number_list1<-c()
  n_w_k<- length(w_k_vals)
  for(i in 1: n_w_k){
    number_of_hawkes0<- sum(simulate_q_samples[i,]==0)
    number_of_hawkes1<- sum(simulate_q_samples[i,]==1)
    hawkes_number_list0<- c(hawkes_number_list0, number_of_hawkes0)
    hawkes_number_list1<- c(hawkes_number_list1, number_of_hawkes1)
  }
  hawkes_number_list2<- nm - hawkes_number_list0 - hawkes_number_list1
  
  denominator_pi<- sum(w_k_vals) * nm
  pi_0_mle<- (sum(w_k_vals * hawkes_number_list0))/(denominator_pi)
  pi_1_mle<- (sum(w_k_vals * hawkes_number_list1))/(denominator_pi)
  pi_2_mle<- (sum(w_k_vals * hawkes_number_list2))/(denominator_pi)
  
  
  hawkes_mle0<- optim(c(1.5,5,6), Q_theta, indexh = 0, obs= observed_time_vals,label_list=simulate_q_samples,w_k_vec= w_k_vals,n_sample = n_w_k)
  hawkes_mle1<- optim(c(3,4,6), Q_theta, indexh = 1, obs= observed_time_vals,label_list=simulate_q_samples,w_k_vec= w_k_vals,n_sample = n_w_k)
  hawkes_mle2<- optim(c(5,3,6), Q_theta, indexh = 2, obs= observed_time_vals,label_list=simulate_q_samples,w_k_vec= w_k_vals,n_sample = n_w_k)
  
  mu_mle0<- hawkes_mle0$par[1]
  alpha_mle0<- hawkes_mle0$par[2]
  beta_mle0<- hawkes_mle0$par[3]
  mu_mle1<- hawkes_mle1$par[1]
  alpha_mle1<- hawkes_mle1$par[2]
  beta_mle1<- hawkes_mle1$par[3]
  mu_mle2<- hawkes_mle2$par[1]
  alpha_mle2<- hawkes_mle2$par[2]
  beta_mle2<- hawkes_mle2$par[3]
  params<- c(mu_mle0,alpha_mle0,beta_mle0,mu_mle1,alpha_mle1,beta_mle1, mu_mle2,alpha_mle2,beta_mle2,pi_0_mle,pi_1_mle,pi_2_mle)
  
  list("mu0" = mu_mle0,
       "alpha0" = alpha_mle0,
       "beta0" = beta_mle0,
       "mu1" = mu_mle1,
       "alpha1" = alpha_mle1,
       "beta1" = beta_mle1,
       "mu2" = mu_mle2,
       "alpha2" = alpha_mle2,
       "beta2" = beta_mle2,
       "pi" = c(pi_0_mle, pi_1_mle, pi_2_mle),
       "params" = params)
}

#Simulation
Total_Time_Interval<- 60
True_mu0<- 1.5
True_alpha0<- 5
True_beta0<- 6
True_mu1<- 3
True_alpha1<- 4
True_beta1<- 6
True_mu2<- 5
True_alpha2<- 3
True_beta2<- 6

#Hawkes 0
Model_Hawkes_0<- simulateHawkes(lambda0 = True_mu0,alpha = True_alpha0,beta = True_beta0,horizon = Total_Time_Interval)[[1]]
nh0<- length(Model_Hawkes_0)
labelh0<- replicate(nh0,0)
Model_HawkeswithLabels_0<- data.frame(t=Model_Hawkes_0,Label=labelh0)
#Hawkes 1

Model_Hawkes_1<- simulateHawkes(lambda0 = True_mu1,alpha = True_alpha1,beta = True_beta1,horizon = Total_Time_Interval)[[1]]
nh1<- length(Model_Hawkes_1)
labelh1<- replicate(nh1,1)
Model_HawkeswithLabels_1<- data.frame(t=Model_Hawkes_1,Label=labelh1)
#Hawkes 2

Model_Hawkes_2<- simulateHawkes(lambda0 = True_mu2,alpha = True_alpha2,beta = True_beta2,horizon = Total_Time_Interval)[[1]]
nh2<- length(Model_Hawkes_2)
labelh2<- replicate(nh2,2)
Model_HawkeswithLabels_2<- data.frame(t=Model_Hawkes_2,Label=labelh2)

#Mixture
Mixture_sim<- rbind(Model_HawkeswithLabels_0,Model_HawkeswithLabels_1,Model_HawkeswithLabels_2)
Mixture_observations<- Mixture_sim[order(Mixture_sim$t),]

nm<- nh0 + nh1 + nh2
pi_hawkes_0<- nh0/nm
pi_hawkes_1<- nh1/nm
pi_hawkes_2<- nh2/nm

for (i in 1:10){
  if (i ==1){
    e.step<- e_step_threehawkes(observed_t = Mixture_observations[,1],mu.n0 = True_mu0,alpha.n0 = True_alpha0,beta.n0 = True_beta0,mu.n1 = True_mu1,alpha.n1 = True_alpha1,beta.n1 = True_beta1,mu.n2 = True_mu2,alpha.n2 = True_alpha2,beta.n2 = True_beta2,pi_vec = c(pi_hawkes_0,pi_hawkes_1, pi_hawkes_2))
    m.step<- m_step_threehawkes(observed_time_vals = Mixture_observations[,1], w_k_vals = e.step[["w_k_vals"]], simulate_q_samples = e.step[["simulate_q_samples"]])
    cur.params<- m.step[["params"]]

  }else{
    e.step<- e_step_threehawkes(observed_t = Mixture_observations[,1],mu.n0 = m.step[["mu0"]],alpha.n0 = m.step[["alpha0"]],beta.n0 = m.step[["beta0"]],mu.n1 = m.step[["mu1"]],alpha.n1 = m.step[["alpha1"]],beta.n1 = m.step[["beta1"]],mu.n2 = m.step[["mu2"]],alpha.n2 = m.step[["alpha2"]],beta.n2 = m.step[["beta2"]],pi_vec = m.step[["pi"]])
    m.step<- m_step_threehawkes(observed_time_vals = Mixture_observations[,1], w_k_vals = e.step[["w_k_vals"]], simulate_q_samples = e.step[["simulate_q_samples"]])
    params_dif<- cur.params - m.step[["params"]]
    norm_dif<- sqrt(sum(params_dif^(2)))
    if (norm_dif < 0.01){
      break
    }else{
      cur.params<- m.step[["params"]]
    }
  }
}
print(m.step)