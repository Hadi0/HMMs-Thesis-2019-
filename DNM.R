#Packages
library(sigmoid)
library(tictoc)
library(Rlab)
###

n<-0  
vals<-data.frame(matrix(0, nrow=1, ncol=6))

#make.NegLogLik initial function which takes as arguments fixed parameters
make.NegLogLik<-function(par){
  n <<- n+1
  vals[n,1:5] <<-par
  print("#####################")
  print(sigmoid(vals[n, 1:2]))
  print(tautogamma(matrix(c(NA, vals[n,3], vals[n,4], NA), nrow=2,ncol=2, byrow=TRUE)))
  print(tautogamma(matrix(c(NA,vals[n, 5]), nrow=1, ncol=2, byrow=TRUE)))
  # Inserting in working parameters
  log_parameter_vec<-matrix(c(as.numeric(par[1]), as.numeric(par[2])), nrow=1, ncol=2, byrow=TRUE)
  taumat<- matrix(c(NA, as.numeric(par[3]), 
                    as.numeric(par[4]), NA), 
                  nrow=2, ncol=2, byrow=TRUE)
  #Working to natural
  gamma_mat<-tautogamma(taumat)
  parameter_vec<-sigmoid(log_parameter_vec)
  
  #Returns stationary distribution 
  # stationary_dist<-function(gamma_mat){
  #   eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))]/sum(eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))])
  # }
  #Initial distribution 
  #only positive eigenvalue 
  u_1<-matrix(c(NA, par[5]), nrow=1, ncol=length(exp(log_parameter_vec)))
  u_1<-tautogamma(u_1)
  ################################################################################
  #Number of time periods in dataset; if header, do -1
  time<-nrow(data.frame(dataset$obs));
  
  ################################################################################
  #FUNCTIONS
  #returns emission matrix with diagonal entries being P(X |State)
  a<-matrix(0L, nrow =length(parameter_vec), ncol=length(parameter_vec))
  emissionmat <- function(obs){
    for (j in c(1:length(parameter_vec))){
      a[j,j]<- dbern(obs, parameter_vec[j])  
    }
    a
  }
  
  #Returns -inf if NaN
  NaNtoInf<-function(difference){
    if (is.nan(difference) == TRUE)
    {
      log(0)
    }
    else{
      sum_log[j,k]-max_alpha[t-1,j]
    }
  }
  ################################################################################
  
  #Forward Algorithm 
  #Alpha1 
  #Defining a dummy vector which we will replace the values of 
  alpha_1<-matrix(0:0,nrow(u_1),length(exp(log_parameter_vec))) 
  #Filling in the dummy vector
  emiss_first<-emissionmat(dataset$obs[1])
  for(j in 1:length(exp(log_parameter_vec)))
    for(k in 1:ncol(u_1))
      alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emiss_first[k,j]
  #Logging the result and storing it in a data frame 
  log_alpha_1<-log(alpha_1)
  log_alpha_dataframe<-data.frame(log_alpha_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
  
  sum_log<-data.frame(matrix(0:0, 2, length(exp(log_parameter_vec)))) 
  max_alpha<- data.frame(matrix(0:0, 1, 2),check.names = TRUE)
  emiss_second<-emissionmat(dataset$obs[2])
  for (j in c(1:ncol(gamma_mat))){
    for (k in 1:length(exp(log_parameter_vec))){
      sum_log[j,k]<- log_alpha_dataframe[1,][1,k] + log(gamma_mat[k,j])+log(emiss_second[j,j])
    }
    max_alpha[1,j]<- max(sum_log[j,])
  }
  #tic("main loop")
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
          log_alpha_dataframe[t,][1,j]<- log_alpha_dataframe[t,][1,j]+exp(NaNtoInf(sum_log[j,k]-max_alpha[t-1,j]))
          #Work out new value for sum in new row 
        }
        #Then log the summation and add max
        log_alpha_dataframe[t,][1,j]<- max_alpha[t-1,j]+log(log_alpha_dataframe[t,][1,j])
      }
      emiss_t_1<-emissionmat(dataset$obs[t+1])
      for (j_1 in c(1:ncol(gamma_mat))){
        for (k_1 in 1:length(exp(log_parameter_vec))){
          #t+1 so, in the next loop for t, everything is correct (use next obs, so next loop of t uses correct value)
          sum_log[j_1,k_1]<- log_alpha_dataframe[t,][1,k_1] + log(gamma_mat[k_1,j_1])+log(emiss_t_1[j_1,j_1])
        }
        max_alpha[t,j_1]<- max(sum_log[j_1,])
      }
    }
    #toc()
  }
  
  #Function to calculate the likelihood of the t=TIME obs 
  final.negloglike<-function(last_log_alpha){
    log_like<-0
    for (i in c(1:length(parameter_vec))){
      a_max<-max(last_log_alpha)
      log_like<- log_like+exp(last_log_alpha[i]-a_max)
    }
    log_like <-a_max+log(log_like)
    log_like
  }
  likelihood__<--as.numeric(final.negloglike(log_alpha_dataframe[time,]))
  vals[n,6] <<-likelihood__
  print(likelihood__)
}

#NLM
tic()
nlm(make.NegLogLik,c(logit(0.4), logit(0.2), 
                     gammatotau(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE))[2], 
                     gammatotau(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE))[3],
                     gammatotau(matrix(stationary_dist(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE)), nrow=1))[2]
))
toc() 

#OPTIM
tic()
optim(par=c(logit(0.4), logit(0.2), 
            gammatotau(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE))[2], 
            gammatotau(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE))[3],
            gammatotau(matrix(stationary_dist(matrix(c(0.2, 0.8, 0.8, 0.2), nrow=2, byrow=TRUE)), nrow=1))[2]
),
fn=make.NegLogLik)
toc()