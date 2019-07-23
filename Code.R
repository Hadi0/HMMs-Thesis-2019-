  #Test data
    #dataset<-read.csv(file="~/desktop/HMM_Dissertation/HMM_project/rainySample.csv")
      dataset<-data.frame(matrix(c(1,1,1,0), colnames("obs"), nrow=4, ncol=1, byrow = FALSE))
      colnames(dataset)<-c("obs")
  ########################################
      
  
  #gamma_mat Prob Matrix 
  gamma_mat<-matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2, byrow=TRUE)
  #gamma_mat<-matrix(c(0.7, 0.2, 0.1, 0.2, 0.5, 0.3, 0.5, 0.1, 0.4), nrow = 3, ncol = 3, byrow=TRUE)
  
  
  #Emission Probabilties (in order from 1 to nth possible emmission)
  emission_A<- matrix(c(0.9, 0, 0, 0.3),nrow = 2, ncol = 2, byrow=TRUE) #Sun
  emission_B<- matrix(c(0.1, 0, 0, 0.7),nrow = 2, ncol = 2, byrow=TRUE) #Rain 
  #emission_A<- matrix(c(0.7, 0, 0, 0, 0.3, 0, 0, 0, 0.4),nrow = 3, ncol = 3, byrow=TRUE)
  #emission_B<- matrix(c(0.1, 0, 0,0, 0.7,0,0,0,0.5),nrow = 3, ncol = 3, byrow=TRUE)
  
  #Initial distribution 
  u_1<-matrix(c(1/1.5, 0.5/1.5), nrow=1, ncol =nrow(emission_A))
  #u_1<-matrix(c(0.1, 0.5, 0.4), nrow=1, ncol =nrow(emission_A))
  
  #Number of time periods in dataset; if header, do -1
  time<-nrow(data.frame(dataset$obs));
  ########################################
  
  
  #Defining necessary functions 
  #returns emission sun or emission rain depending on given observation
  emissionmat <- function(obs){
    if(obs == 1){
      emission_A 
    }  
    else{
      emission_B
    }
  }
  #Returns -inf if NaN
  NaNtoInf<-function(difference){
    if (is.nan(difference) == TRUE)
    {
      log(0)
    }
    else{
      log_alpha_dataframe[t-1,][1,k]+log(gamma_mat[k,j])+log(emissionmat(dataset$obs[t])[j,j])-max_alpha
    }
  }
  ########################################
  
  
  #Forward Algorithm 
  #Alpha1 
  #Defining a dummy vector which we will replace the values of 
   alpha_1<-matrix(0:0,nrow(u_1),ncol(emission_A)) 
  #Filling in the dummy vector  
  for(j in 1:ncol(emission_A))
    for(k in 1:ncol(u_1))
      alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emissionmat(dataset$obs[1])[k,j]
  #Logging the result and storing it in a data frame 
  log_alpha_1<-log(alpha_1)
  log_alpha_dataframe<-data.frame(log_alpha_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
  colnames(log_alpha_dataframe)<-c("alpha(1)", "alpha(2)")
  
  #Defining the rest of the alphas 
  if (time>1){
    for (t in c(2:time)){
      log_alpha_dataframe[t,]<-data.frame(matrix(0:0, nrow=1, ncol=nrow(u_1)))
      for (j in c(1:ncol(gamma_mat))){
        for (k in c(1:ncol(u_1))){ 
          #Defining Max 
          max_alpha<-max(log_alpha_dataframe[t-1,])+log(gamma_mat[which(log_alpha_dataframe[t-1,] == max(log_alpha_dataframe[t-1,])),j])+log(emissionmat(dataset$obs[t])[j,j])
          #Defining the Sum of the logs 
          sum_log<-log_alpha_dataframe[t-1,][1,k]+log(gamma_mat[k,j])+log(emissionmat(dataset$obs[t])[j,j])
          
          #First do the summation 
          log_alpha_dataframe[t,][1,j]<- log_alpha_dataframe[t,][1,j]+exp(NaNtoInf(sum_log-max_alpha))
        }
        #Then log the summation and add max
        log_alpha_dataframe[t,][1,j]<- max_alpha+log(log_alpha_dataframe[t,][1,j])
      }
    }
  }
  ########################################
  
  #Printing off results 
  print(dataset)
  print(u_1)
  print(gamma_mat)
  print(emission_A)
  print(emission_B)
  log_alpha_dataframe
  print(exp(log_alpha_dataframe))
  
