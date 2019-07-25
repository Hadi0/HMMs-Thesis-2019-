#Returns stationary distribution 
stationary_dist<-function(gamma_mat){
  eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))]/sum(eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))])
}

#Let nrow = no. of simulated obs 
dataset<-data.frame(matrix(0L, nrow=10, ncol=2))
colnames(dataset)<-c("obs", "state")

#Vector of parameters 
parameter_vec<-c(10,20)

#Sampling first  state based on initial distribution 
gamma_mat<-matrix(c(0.8, 0.2, 0.4, 0.6), nrow=2, ncol=2, byrow=TRUE)
dataset$state[1]<-sample(c("State_1","State_2"), 1, prob = abs(matrix(stationary_dist(gamma_mat), nrow=1, ncol=length(parameter_vec))))


##########################################

#Generating States based on transition probabilities 
for (i in 2:10){
  if (dataset$state[i-1]=="State_1"){
    dataset$state[i]<-sample(c("State_1","State_2"), 1, prob = c(0.8, 0.2))
  } else{
      dataset$state[i]<-sample(c("State_1","State_2"), 1, prob = c(0.4, 0.6))
  }
}

#Generating observations based on emission probabilities  
for (i in 1:10){
  if (dataset$state[i]=="State_1"){
    dataset$obs[i]<-rpois(1, parameter_vec[1])
  } else{
      dataset$obs[i]<-rpois(1, parameter_vec[2])
  }
}
write.csv(dataset, "LambdaTenLambdaTwenty.csv")