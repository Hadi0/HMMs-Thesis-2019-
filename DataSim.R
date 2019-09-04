#Size of dataset
size_data<-700

#Vector of parameters 
parameter_vec<-c(0.6,0.3)

#Sampling first  state based on initial distribution 
gamma_mat<-matrix(c(0.9, 0.1, 
                    0.2, 0.8), nrow=2, ncol=2, byrow=TRUE)

##########################################
data_gen<-function(gamma_mat, parameter_vec,size_data){
  #Returns stationary distribution
  stationary_dist<-function(gamma_mat){
    eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))]/sum(eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))])
  }
  #Generating States based on transition probabilities 
  #Function (State generator)
  state_gen<-function(gamma_mat){
    dataset<-data.frame(matrix(0L, nrow=size_data, ncol=2))
    colnames(dataset)<-c("obs", "state")
    dataset$state[1]<-sample(c("State_1","State_2"), 1, prob = abs(matrix(stationary_dist(gamma_mat), nrow=1, ncol=length(parameter_vec))))
    for (i in 2:size_data){
      if (dataset$state[i-1]=="State_1"){
        dataset$state[i]<-sample(c("State_1","State_2"), 1, prob = c(gamma_mat[1,]))
      } 
      else{ 
        if (dataset$state[i-1]=="State_2"){
          dataset$state[i]<-sample(c("State_1","State_2"), 1, prob = c(gamma_mat[2,]))
        }
      }
    }
    dataset
  }
  
  dataset<-state_gen(gamma_mat)
  #Generating observations based on emission probabilities  
  obs_gen<-function(dataset, parameter_vec){
    for (i in 1:size_data){
      if (dataset$state[i]=="State_1"){
        dataset$obs[i]<-rbern(1, parameter_vec[1])
      } else{ if (dataset$state[i]=="State_2"){
        dataset$obs[i]<-rbern(1, parameter_vec[2])
      }
      }
    }
    dataset
  }
  dataset<-obs_gen(dataset, parameter_vec)
  dataset
}
dataset<-data_gen(gamma_mat, parameter_vec,size_data)
setwd(dir="~/Desktop/HMM_Dissertation/HMM_Project/CSV")
write.csv(dataset, "BernSevenG1P5.csv")