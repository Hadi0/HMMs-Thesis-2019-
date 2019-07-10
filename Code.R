#Test data
dataset<-read.csv(file="~/desktop/HMM_Dissertation/rainySample.csv")

dataset<-matrix()
dataset$obs<-matrix(c(1,2,1), nrow=3, ncol=1)

#Generating a dataset
#States
x<-c()
x<-sample(c("High","Low"), 1, replace=TRUE, prob=c(0.5, 0.5) )
#Sampling states
for (i in c(1:999)){ #999 since we set 1st entry of x and then add onto at each step. Thus, after step i, x will have i+1 entries 
  if (x[i]=="High"){
    x<- rbind(x,sample(c("High","Low"), 1, replace=TRUE, prob=c(0.9, 0.1) ))
  } else{ 
    x<- rbind(x,sample(c("High","Low"), 1, replace=TRUE, prob=c(0.2, 0.8) ))
    
  }
}

y<-c()
for (i in c(1:1000)){
  if (x[i]=="High"){
    y<- rbind(y,sample(c(1,2), 1, replace=TRUE, prob=c(0.9, 0.1) ))
    
  }else{
    y<- rbind(y,sample(c(1,2), 1, replace=TRUE, prob=c(0.3, 0.7) ))
    
  }
}
str(x)
str(y)

dataset<-data.frame(cbind(x,y)) 
colnames(dataset)<-c("states", "obs")

####################

#Transition Prob Matrix 
trans<-matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2, byrow=TRUE)

#Emission Probabilties (in order from 1 to nth possible emmission)
emission_sun<- matrix(c(0.9, 0, 0, 0.3),nrow = 2, ncol = 2, byrow=TRUE)
emission_rain<- matrix(c(0.1, 0, 0, 0.7),nrow = 2, ncol = 2, byrow=TRUE)

#Initial distribution 
u_1<-c(1/1.5, 0.5/1.5)

#Number of time periods in dataset  
#If header, do -1
time<-nrow(dataset)-1;

#Forward Algorithm 
#Alpha1 
#Defining a dummy vector which we will replace the values of 
alpha<-matrix(seq(1,5, length.out = length(u_1)), nrow =1)
#Defining Alpha1
alph<-c()
if (dataset$obs[1]==1){
  alpha[c(1,2)]<- matrix(c(u_1%*%emission_sun), nrow=1)
}else{
  alpha[c(1,2)]<- matrix(c(u_1%*%emission_rain), nrow=1)
}

#Defining the rest of the alphas 
count=0; 
#if one step step, can just print off alpha calcu
if (time>1){
  for(i in c(2:time)){ 
    count=i;
    if ((count<=time)&(dataset$obs[i]==1)){
      alpha[c(1,2)]<-alpha[c(1,2)]%*%trans%*%emission_sun;
    }else{
      alpha[c(1,2)]<-alpha[c(1,2)]%*%trans%*%emission_rain
    }
    
  }
}

likelihood<-sum(alpha)


