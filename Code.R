#Test data
#dataset<-read.csv(file="~/desktop/HMM_Dissertation/HMM_project/rainySample.csv")
#dataset$obs
    dataset<-matrix()
    dataset$obs<-matrix(c(1,1), nrow=2, ncol=1)

####################

#gamma_matition Prob Matrix 
gamma_mat<-matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2, byrow=TRUE)
#gamma_mat<-matrix(c(0.7, 0.2, 0.1, 0.2, 0.5, 0.3, 0.5, 0.1, 0.4), nrow = 3, ncol = 3, byrow=TRUE)
#gamma


#Emission Probabilties (in order from 1 to nth possible emmission)
emission_sun<- matrix(c(1, 0, 0, 0),nrow = 2, ncol = 2, byrow=TRUE)
emission_rain<- matrix(c(0.1, 0, 0, 0.7),nrow = 2, ncol = 2, byrow=TRUE)
#emission_sun<- matrix(c(0.7, 0, 0, 0, 0.3, 0, 0, 0, 0.4),nrow = 3, ncol = 3, byrow=TRUE)
#emission_rain<- matrix(c(0.1, 0, 0,0, 0.7,0,0,0,0.5),nrow = 3, ncol = 3, byrow=TRUE)

#Initial distribution 
u_1<-matrix(c(1/1.5, 0.5/1.5), nrow=1, ncol =nrow(emission_sun))
#u_1<-matrix(c(0.1, 0.5, 0.4), nrow=1, ncol =nrow(emission_sun))

#Number of time periods in dataset  
#If header, do -1
time<-nrow(data.frame(dataset$obs));

#Forward Algorithm 
#returns emission sun or emission rain depending on given observation
emissionmat <- function(obs){
  if(obs == 1){
    emission_sun 
  }  
  else{
    emission_rain
  }
}

#Alpha1 
#Defining a dummy vector which we will replace the values of 
 alpha_1<-matrix(0:0,nrow(u_1),ncol(emission_sun)) 
 for(j in 1:ncol(emission_sun))
   for(k in 1:ncol(u_1))
     alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emissionmat(dataset$obs[1])[k,j]

alpha_1<-matrix(0:0,nrow(u_1),ncol(emission_sun)) 
for(j in 1:ncol(emission_sun)){
  for(k in 1:ncol(u_1)){
    alpha_1[1,j]<-alpha_1[1,j]+u_1[1,k]*emissionmat(dataset$obs[1])[j,k]
  }
}
alpha_dataframe<-data.frame(alpha_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
colnames(alpha_dataframe)<-c("alpha(1)", "alpha(2)")

#Defining the rest of the log alphas
#  if (time>1){
#    for (t in c(2:time)){
#      alpha_dataframe[t,]<-data.frame(matrix(0:0, nrow=1, ncol=nrow(u_1)))
#     for (j in c(1:ncol(gamma_mat))){
#       for (k in c(1:ncol(u_1))){
#         alpha_dataframe[t,][1,j]<- alpha_dataframe[t,][1,j]+alpha_dataframe[t-1,][1,k]*gamma_mat[k,j]*emissionmat(dataset$obs[t])[j,j]
#       }
#     }
#   }
# }

log_alpha_1<-log(alpha_1)
log_alpha_dataframe<-data.frame(log_alpha_1, check.names = TRUE) #First column represents time, each row represents the component of the likelihood
colnames(log_alpha_dataframe)<-c("alpha(1)", "alpha(2)")
if (time>1){
  for (t in c(2:time)){
    log_alpha_dataframe[t,]<-data.frame(matrix(0:0, nrow=1, ncol=nrow(u_1)))
    for (j in c(1:ncol(gamma_mat))){
      for (k in c(1:ncol(u_1))){ 
        #Defining Max 
        max_alpha<-max(log_alpha_dataframe[t-1,])+log(gamma_mat[which(log_alpha_dataframe[t-1,] == max(log_alpha_dataframe[t-1,])),j])+log(emissionmat(dataset$obs[t])[j,j])
  
        #First do the summation 
        log_alpha_dataframe[t,][1,j]<- log_alpha_dataframe[t,][1,j]+exp(log_alpha_dataframe[t-1,][1,k]+log(gamma_mat[k,j])+log(emissionmat(dataset$obs[t])[j,j])
                                                                        - max_alpha)
          if (is.nan(log_alpha_dataframe[t,][1,j]) == TRUE){
          log_alpha_dataframe[t,][1,j] = -1/0
          }
        #if (any(sapply(log_alpha_dataframe[t,][1,j], is.nan)) =TRUE){
          #log_alpha_dataframe[t,][1,j] = -1/0
        }
        }
      #Then log the summation 
      log_alpha_dataframe[t,][1,j]<- max_alpha+log(log_alpha_dataframe[t,][1,j])
    }
  }

#######################
# max_log_alpha<- data.frame(matrix(0:0, nrow=1, ncol=1))
# max_alpha<-max(log_alpha_dataframe[t-1,])+log(gamma_mat[which(log_alpha_dataframe[t-1,] == max(log_alpha_dataframe[t-1,])),j])+log(emissionmat(dataset$obs[t])[j,j])
print(dataset)
print(u_1)
print(gamma_mat)
print(emission_sun)
print(emission_rain)
print(log_alpha_dataframe)

