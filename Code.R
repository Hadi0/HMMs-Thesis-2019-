library(Rlab)
#Test data
#dataset<-read.csv(file="~/desktop/HMM_Dissertation/HMM_project/rainySample.csv")
dataset<-data.frame(matrix(c(10,20,1,5,1,2,3,1), nrow=4, ncol=2, byrow = FALSE))
colnames(dataset)<-c("obs")

#gamma_mat Prob Matrix 
gamma_mat<-c()

#Parameter (in order from 1 to nth possible state)
parameter_vec<-c(1, 2)

#Initial distribution 
u_1<-matrix(c(2/3, 1/3), nrow=1, ncol =length(parameter_vec))


#Number of time periods in dataset; if header, do -1
time<-nrow(data.frame(dataset$obs));
################################################################################


#returns emission matrix with diagonal entries being P(X |State)
a<-matrix(0L, nrow =length(parameter_vec), ncol=length(parameter_vec))
emissionmat <- function(obs){
  for (j in c(1:length(parameter_vec))){
    a[j,j]<- dbern(obs, parameter_vec[j])  
    # matrix(c(dpois(obs, lambda_1), 0, 0, dpois(obs, lambda_2)), nrow=2, ncol=2, byrow=FALSE)
    # matrix(c(dbern(obs, 0.9), 0, 0, dbern(obs, 0.3)), nrow=2, ncol=2, byrow=FALSE)
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
########################################


#Forward Algorithm 
#Alpha1 
#Defining a dummy vector which we will replace the values of 
alpha_1<-matrix(0:0,nrow(u_1),length(parameter_vec)) 
#Filling in the dummy vector  
for(j in 1:length(parameter_vec))
  for(k in 1:ncol(u_1))
    alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emissionmat(dataset$obs[1])[k,j]
#Logging the result and storing it in a data frame 
log_alpha_1<-log(alpha_1)
log_alpha_dataframe<-data.frame(log_alpha_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
colnames(log_alpha_dataframe)<-c("alpha(1)", "alpha(2)")

sum_log<-data.frame(matrix(0:0, 2, length(parameter_vec))) 
colnames(sum_log)<-c("k=1", "k=2")
#rownames(sum_log)<-c("j=1","j=2")
max_alpha<- data.frame(matrix(0:0, 1, 2),check.names = TRUE)
colnames(max_alpha)<-c("max j=1", "maxj=2")
#rownames(max_alpha)<-c("t=1")

for (j in c(1:ncol(gamma_mat))){
  # sum_log[j,]<-data.frame(matrix(0:0, nrow=1, length(parameter_vec)))
  # max_alpha[j,]<-data.frame(matrix(0:0, nrow=1, length(parameter_vec)))
  for (k in 1:length(parameter_vec)){
    sum_log[j,k]<- log_alpha_dataframe[1,][1,k] + log(gamma_mat[k,j])+log(emissionmat(dataset$obs[2])[j,j])
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
        log_alpha_dataframe[t,][1,j]<- log_alpha_dataframe[t,][1,j]+exp(NaNtoInf(sum_log[j,k]-max_alpha[t-1,j]))
        #Work out new value for sum in new row 
      }
      #Then log the summation and add max
      log_alpha_dataframe[t,][1,j]<- max_alpha[t-1,j]+log(log_alpha_dataframe[t,][1,j])
    }
    for (j_1 in c(1:ncol(gamma_mat))){
      for (k_1 in 1:length(parameter_vec)){
        #t+1 so, in the next loop for t, everything is correct (use next obs, so next loop is correct)
        sum_log[j_1,k_1]<- log_alpha_dataframe[t,][1,k_1] + log(gamma_mat[k_1,j_1])+log(emissionmat(dataset$obs[t+1])[j_1,j_1])
      }
      max_alpha[t,j_1]<- max(sum_log[j_1,])
    }
  }
}
########################################

#Printing off results 
print(dataset)
print(u_1)
print(gamma_mat)
log_alpha_dataframe
print(exp(log_alpha_dataframe))