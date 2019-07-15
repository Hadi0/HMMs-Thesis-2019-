#Test data
#dataset<-read.csv(file="~/desktop/HMM_Dissertation/rainySample.csv")

dataset<-matrix()
dataset$obs<-matrix(c(1,0,1,1), nrow=4, ncol=1)

#Generating a dataset
#States
#x<-c()
#x<-sample(c("High","Low"), 1, replace=TRUE, prob=c(0.5, 0.5) )
#Sampling states
#for (i in c(1:999)){ #999 since we set 1st entry of x and then add onto at each step. Thus, after step i, x will have i+1 entries 
#  if (x[i]=="High"){
#    x<- rbind(x,sample(c("High","Low"), 1, replace=TRUE, prob=c(0.9, 0.1) ))
#  } else{ 
#    x<- rbind(x,sample(c("High","Low"), 1, replace=TRUE, prob=c(0.2, 0.8) ))

#  }
#}

#y<-c()
#for (i in c(1:1000)){
#  if (x[i]=="High"){
#    y<- rbind(y,sample(c(1,2), 1, replace=TRUE, prob=c(0.9, 0.1) ))

#  }else{
#    y<- rbind(y,sample(c(1,2), 1, replace=TRUE, prob=c(0.3, 0.7) ))

#  }
#}

#dataset<-data.frame(cbind(x,y)) 
#colnames(dataset)<-c("states", "obs")

####################

#gammaition Prob Matrix 
gamma<-matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2, byrow=TRUE)

#Emission Probabilties (in order from 1 to nth possible emmission)
emission_sun<- matrix(c(0.9, 0, 0, 0.3),nrow = 2, ncol = 2, byrow=TRUE)
#Gamma%*%Emiss 
gamma_em_sun<-matrix(0:0, nrow=nrow(gamma), ncol=ncol(emission_sun))
for(i in 1:nrow(gamma))
  for(j in 1:ncol(emission_sun))
    for(k in 1:ncol(gamma))
      gamma_em_sun[i,j]<-gamma_em_sun[i,j]+gamma[i,k]*emission_sun[k,j]

emission_rain<- matrix(c(0.1, 0, 0, 0.7),nrow = 2, ncol = 2, byrow=TRUE)
#Gamma%*%Emiss 
gamma_em_rain<-matrix(0:0, nrow=nrow(gamma), ncol=ncol(emission_rain))
for(i in 1:nrow(gamma))
  for(j in 1:ncol(emission_rain))
    for(k in 1:ncol(gamma))
      gamma_em_rain[i,j]<-gamma_em_rain[i,j]+gamma[i,k]*emission_rain[k,j]

#Initial distribution 
u_1<-matrix(c(1/1.5, 0.5/1.5), nrow=1, ncol=2)

#Number of time periods in dataset  
#If header, do -1
time<-nrow(dataset$obs);

#Forward Algorithm 
#Alpha1 
#Defining a dummy vector which we will replace the values of 
alpha_1<-matrix(0:0,nrow(u_1),ncol(emission_sun)) 
#Defining Alpha1

if (dataset$obs[1] == 1){
  for(i in 1:nrow(u_1))
    for(j in 1:ncol(emission_sun))
      for(k in 1:ncol(u_1))
        alpha_1[i,j]<-alpha_1[i,j]+u_1[i,k]*emission_sun[k,j]
}else{
  for(i in 1:nrow(u_1))
    for(j in 1:ncol(emission_rain))
      for(k in 1:ncol(u_1))
        alpha_1[i,j]<-alpha_1[i,j]+u_1[i,k]*emission_rain[k,j]
      
}

alpha_dataframe<-data.frame(alpha_1, check.names = TRUE)
colnames(alpha_dataframe)<-c("alpha(1)", "alpha(2)")
#Defining the rest of the alphas
count=0; 
#if one step, can just print off alpha calcu
if (time>1){
  for(i in c(2:time)){ 
    alpha_dataframe[i,]<-data.frame(matrix(0:0, nrow=1, ncol=1))
    count=i;
    if ((count<=time)&(dataset$obs[i]==1)){
      for(h in 1:nrow(u_1))
        for(j in 1:ncol(gamma_em_sun))
          for(k in 1:ncol(u_1))
            alpha_dataframe[i,][h,j]<-alpha_dataframe[i,][h,j]+alpha_dataframe[i-1,][h,k]*gamma_em_sun[k,j];
          alpha_dataframe<-alpha_dataframe
    }
    else{
      for(h in 1:nrow(u_1))
        for(j in 1:ncol(gamma_em_rain))
          for(k in 1:ncol(u_1))
            alpha_dataframe[i,][h,j]<-alpha_dataframe[i,][h,j]+alpha_dataframe[i-1,][h,k]*gamma_em_rain[k,j];
          alpha_dataframe<-alpha_dataframe
    }
    
  }
}
alpha_dataframe
likelihood<-(sum(alpha_dataframe[time,]))
