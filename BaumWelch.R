alphabeta<-function(u_1, parameter_vec, gamma_mat, type_){
  if (type_==1){
    #Forward/Backward Algorithm 
    #Alpha1/ BetaT-1
    #Defining a dummy vector which we will replace the values of 
    alpha_1<-matrix(0:0,1,length(parameter_vec)) 
    emission_mat_1<-emissionmat(dataset$obs[1])
    for(j in 1:length(parameter_vec)){
      for(k in 1:length(parameter_vec)){
        alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emission_mat_1[k,j]
      }
    }
    #Logging the result and storing it in a data frame 
    log_alpha_dataframe<-data.frame(log(alpha_1), check.names = TRUE) #First column represents time, each row represents the component of the likelihood
    sum_log<-data.frame(matrix(0:0, length(parameter_vec), length(parameter_vec))) 
    max_alpha<- data.frame(matrix(0:0, 1, length(parameter_vec)),check.names = TRUE)
    emission_mat_2<-emissionmat(dataset$obs[2])
    for (j in c(1:ncol(gamma_mat))){
      for (k in 1:length(parameter_vec)){
        sum_log[j,k]<- log_alpha_dataframe[1,][1,k] + log(gamma_mat[k,j])+log(emission_mat_2[j,j])
      }
      max_alpha[1,j]<- max(sum_log[j,])
    }
    #Defining the rest of the alphas 
    if (time>1){
      for (t in c(2:time)){
        #defining a dummy row for alpha_t which we will fill in 
        log_alpha_dataframe[t,]<-data.frame(matrix(0:0, nrow=1, ncol=ncol(u_1)))
        #defining a dummy row for max_alpha which we will fill it
        max_alpha[t,]<-data.frame(matrix(0:0, nrow=1, ncol=2))
        for (j in c(1:ncol(gamma_mat))){
          for (k in c(1:ncol(u_1))){ 
            #First do the summation 
            log_alpha_dataframe[t,][1,j]<- log_alpha_dataframe[t,][1,j]+exp((sum_log[j,k]-max_alpha[t-1,j]))
            #Work out new value for sum in new row 
          }
          #Then log the summation and add max
          log_alpha_dataframe[t,][1,j]<- max_alpha[t-1,j]+log(log_alpha_dataframe[t,][1,j])
        }
        emission_mat_t_plus_one<-emissionmat(dataset$obs[t+1])
        for (j_1 in c(1:ncol(gamma_mat))){
          for (k_1 in 1:length(parameter_vec)){
            #t+1 so, in the next loop for t, everything is correct (use next obs, so next loop of t uses correct value)
            sum_log[j_1,k_1]<- log_alpha_dataframe[t,][1,k_1] + log(gamma_mat[k_1,j_1])+log(emission_mat_t_plus_one[j_1,j_1])
          }
          max_alpha[t,j_1]<- max(sum_log[j_1,])
        }
      }
    }
    log_alpha_dataframe
  }else{
    if (type_==2){ 
      beta_T_1<-matrix(0:0,1,length(parameter_vec)) 
      #Filling in the dummy vector  
      for(j in 1:length(parameter_vec)){
        for(k in 1:length(parameter_vec)){
          beta_T_1[1,j]<-beta_T_1[1, j]+gamma_mat[j,k]*emissionmat(dataset$obs[time])[k, k]
        }
      }
      
      
      log_beta_T_1<-log(beta_T_1)
      log_beta_dataframe<-data.frame(matrix(c(NA), ncol=length(parameter_vec), nrow=1), check.names = TRUE)
      log_beta_dataframe[time-1,]<-data.frame(log_beta_T_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
      sum_log_beta<-data.frame(matrix(0:0, length(parameter_vec), length(parameter_vec))) 
      max_beta<- data.frame(matrix(0:0, 1, length(parameter_vec)),check.names = TRUE)
      
      for (j in c(1:ncol(gamma_mat))){
        for (k in 1:length(parameter_vec)){
          sum_log_beta[j,k] <- log_beta_dataframe[time-1,][1,k]+ log(gamma_mat[j,k])+log(emissionmat(dataset$obs[(time-1)])[k,k])
        }
        max_beta[1,j] <-max(sum_log_beta[j,])
      }
      
      #Defining the rest of the alphas 
      if (time>1){
        for (t in c(2:time)){
          if (t<time){
            #defining a dummy row for alpha_t which we will fill in 
            log_beta_dataframe[time-t,]<-data.frame(matrix(0:0, nrow=1, ncol=ncol(u_1)))
            #defining a dummy row for max_alpha which we will fill it
            max_beta[t,]<-data.frame(matrix(0:0, nrow=1, ncol=2))
            for (j in c(1:ncol(gamma_mat))){
              for (k in c(1:ncol(u_1))){ 
                #First do the summation 
                log_beta_dataframe[time-t,][1,j]<- log_beta_dataframe[time-t,][1,j]+exp((sum_log_beta[j,k]-max_beta[t-1,j]))
                #Work out new value for sum in new row 
              }
              #Then log the summation and add max
              log_beta_dataframe[time-t,][1,j]<- max_beta[t-1,j]+log(log_beta_dataframe[time-t,][1,j])
            }
            for (j_1 in c(1:ncol(gamma_mat))){
              for (k_1 in 1:length(parameter_vec)){
                #t+1 so, in the next loop for t, everything is correct (use next obs, so next loop of t uses correct value)
                sum_log_beta[j_1,k_1]<- log_beta_dataframe[time-t,][1,k_1] + log(gamma_mat[j_1,k_1])+log(emissionmat(dataset$obs[(time-t)])[k_1,k_1])
              }
              max_beta[t,j_1]<- max(sum_log_beta[j_1,])
            }
          }
          if (t==time){
            back_sum<-matrix(0:0,1,length(parameter_vec)) 
            for (j in c(1:length(parameter_vec))){
              for (k in c(1:length(parameter_vec))){
                #First doing the summation
                back_sum[1,j]<-back_sum[1,j]+exp(log(u_1[1,k])+log(emissionmat(dataset$obs[1])[k,j]))
              }
              #logging the sum and adding the jth element of B_1
              back_sum[1,j]<- log_beta_dataframe[1,j]+log(back_sum[1,j])
            }
          }
        }
      } 
      log_beta_dataframe[time,]<-log(data.frame(matrix(1:1, nrow=1, ncol=length(parameter_vec))))
      log_beta_dataframe
    }
  }
}
#######
tic("EM:")
#While loop 
vals<-data.frame(matrix(0, nrow=1, ncol=9))
forward_prob_df<-alphabeta(u_1, parameter_vec, gamma_mat, 1)
max_loglike<-final.loglike(forward_prob_df[time,])
print(max_loglike)
backward_prob_df<-alphabeta(u_1, parameter_vec, gamma_mat, 2)
b<-2
Q<-0
Q[2]<-as.numeric(max_loglike)


while(abs(Q[b]-Q[b-1])>0.001){
  #log_u_est[t,j] is the probability of being in state j at time t given the T=time obs
  log_u_est<-data.frame(matrix(0:0, nrow=1, ncol=length(parameter_vec)))
  for (t in c(1:time)){
    for (j in c(1:length(parameter_vec))){
      log_u_est[t,j]<-forward_prob_df[t,j]+
        backward_prob_df[t,j]-
        max_loglike
    }
  }
  #log_v_est[j, k, t] is the probability of being in state j at time t-1 and state k at time t given the T=time obs
  log_v_est = array(rep(NA, length(parameter_vec)*length(parameter_vec)*(time)),
                    dim=c(length(parameter_vec), length(parameter_vec), time))
  for (j in c(1:length(parameter_vec))){
    for (k in c(1:length(parameter_vec))){
      for (t in c(2:time)){
        log_v_est[j,k,t]<-as.numeric(forward_prob_df[t-1,j]+log(gamma_mat[j,k])+log(emissionmat(dataset$obs[t])[k,k])+backward_prob_df[t,k]
                                     -max_loglike)
      }
    }
  } #sum(exp(log_v_est[,x,y])) = exp(log_u_est[y,x])
  #sum(exp(log_v_est[x,,y+1])) = exp(log_u_est[y,x])
  
  vals[(b-1),1]<-parameter_vec[1] 
  vals[(b-1),2]<-parameter_vec[2]
  print(parameter_vec)
  vals[(b-1),3]<-gamma_mat[1,1]
  vals[(b-1),4]<-gamma_mat[1,2]
  vals[(b-1),5]<-gamma_mat[2,1]
  vals[(b-1),6]<-gamma_mat[2,2]
  print(gamma_mat)
  vals[(b-1),7]<-u_1[1]
  vals[(b-1),8]<-u_1[2]
  print(u_1)
  vals[(b-1),9]<-max_loglike
  print(max_loglike)
  
  #M-step: Maximise:- 
  #1
  log_u_1<-matrix(0:0, nrow=1, ncol=length(parameter_vec))
  for (j in c(1:length(parameter_vec))){
    log_u_1[j]<-log_u_est[1,j]
  }
  u_1<-exp(log_u_1)
  #2
  f<-matrix(0:0, nrow=length(parameter_vec), ncol=length(parameter_vec))
  for (j in c(1:length(parameter_vec))){
    for (k in c(1:length(parameter_vec))){
      for (t in c(2:time)){
        f[j,k]<-f[j,k]+exp(log_v_est[j,k,t])
      } 
    }
  }
  gamma_mat<-matrix(0:0, nrow=length(parameter_vec), ncol=length(parameter_vec))
  for (j in c(1:length(parameter_vec))){
    for (k in c(1:length(parameter_vec))){
      gamma_mat[j,k]<-f[j,k]/sum(f[j,])
    }
  }
  #3
  for (j in c(1:length(parameter_vec))){
    parameter_vec[j]<-sum(exp(log_u_est[,j])%*%dataset$obs)/sum(exp(log_u_est[,j]))
  }
  
  forward_prob_df<-alphabeta(u_1, parameter_vec, gamma_mat, 1)
  max_loglike<-final.loglike(forward_prob_df[time,]); 
  # vals[(b-1),1]<-parameter_vec[1] 
  # vals[(b-1),2]<-parameter_vec[2]
  # print(parameter_vec)
  # vals[(b-1),3]<-gamma_mat[1,1]
  # vals[(b-1),4]<-gamma_mat[1,2]
  # vals[(b-1),5]<-gamma_mat[2,1]
  # vals[(b-1),6]<-gamma_mat[2,2]
  # print(gamma_mat)
  # vals[(b-1),7]<-u_1[1]
  # vals[(b-1),8]<-u_1[2]
  # print(u_1)
  # vals[(b-1),9]<-max_loglike
  # print(max_loglike)
  backward_prob_df<-alphabeta(u_1, parameter_vec, gamma_mat, 2)
  b<-b+1
  print(b)
  Q[b]<-as.numeric(max_loglike)
}
toc()







ggplot(data.frame(vals), aes(1:30)) + 
  geom_line(aes(y = X1, colour = "Parameter 1")) + 
  geom_line(aes(y = X2, colour = "Parameter 2"))+
  scale_color_manual(values = c("darkred", "darkblue"))+
  ggtitle('Real parameter value (0.8 0.3) ') +
  ylab('Parameter Value')+
  xlab("Iteration of Baum-Welch")+
  geom_hline(yintercept=c(0.8, 0.3), linetype="dotted")
