#Poisson Examples 
#1
dataset<-data.frame(matrix(c(1,2), nrow=2, ncol=1, byrow = FALSE))
colnames(dataset)<-c("obs")
gamma_mat<-matrix(c(1,0,0,1), nrow = 2, ncol = 2, byrow=TRUE)
parameter_vec<-c(1, 2)
u_1<-matrix(c(2/3, 1/3), nrow=1, ncol =length(parameter_vec))
#Ans: matrix(c(-1.405, -2.405, -3.09, -3.71), nrow=2, ncol=2, byrow=TRUE)

#2 #MrZero 
dataset<-data.frame(matrix(c(1, 10, 3, 0, 100), nrow=5, ncol=1, byrow = FALSE))
colnames(dataset)<-c("obs")
parameter_vec<-c(3, 4, 7, 1, 100)
u_1<-matrix(c(0.2, 0.3, 0.01, 0.4, 0.09), nrow=1, ncol =length(parameter_vec))
gamma_mat<-matrix(c(0.5, 0.4, 0.1, 0, 0, 0, 0.8, 0.2, 0, 0, 0.1, 0.2, 0.5, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0, 0, 0, 0.5, 0.5), nrow = 5, ncol = 5, byrow=TRUE)

#Bernouli Examples 
#1
dataset<-data.frame(matrix(c(1, 1, 1, 0, 1), nrow=5, ncol=1, byrow = FALSE))
colnames(dataset)<-c("obs")
head(dataset)
gamma_mat<-matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, ncol = 2, byrow=TRUE)
parameter_vec<-c(0.9, 0.3)
u_1<-matrix(c(2/3, 1/3), nrow=1, ncol =length(parameter_vec))
#Ans: matrix(c(-0.5108256, -2.302585, -0.6851790, -3.170086,-0.8775509, -3.680911,-3.2721184, -3.141498, -3.2571152, -4.464530), nrow=5, ncol=2, byrow=TRUE)