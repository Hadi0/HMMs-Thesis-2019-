#List of all functions used 
#Returns stationary distribution 
stationary_dist<-function(gamma_mat){
  eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))]/sum(eigen(t(gamma_mat))[[2]][,which.min(abs(1-eigen(gamma_mat)[[1]]))])
}

# # returns emission matrix with diagonal entries being P(X |State): POISSON
# a<-matrix(0L, nrow =length(parameter_vec), ncol=length(parameter_vec))
# emissionmat <- function(obs){
#   for (j in c(1:length(parameter_vec))){
#     a[j,j]<- dpois(obs, parameter_vec[j])
#   }
#   a
# }

#returns emission matrix with diagonal entries being P(X |State): BERN
emissionmat <- function(obs){
  a<-matrix(0L, nrow =length(parameter_vec), ncol=length(parameter_vec))

  for (j in c(1:length(parameter_vec))){
    a[j,j]<- dbern(obs, parameter_vec[j])
  }
  a
}

#Converts (log(a) log(b)) to log(a+b)
final.loglike<-function(last_log_alpha){
  log_like<-0
  for (i in c(1:length(parameter_vec))){
    a_max<-max(last_log_alpha)
    log_like<- log_like+exp(last_log_alpha[i]-a_max)
  }
  log_like <-a_max+log(log_like)
  log_like
}

#Converts logB_1 to log(L_T)
final.loglike.backward<-function(backwardfirst){
  sumsum<-0
  for (k in c(1:length(parameter_vec))){
    for (j in c(1:length(parameter_vec))){
      sumsum<-sumsum+exp(log(u_1[j])+log(emissionmat(dataset$obs[1])[j,k])+backwardfirst[k]-max(backwardfirst))
    }
  }
  max(store)+log(sumsum)
}

#function to convert tau space matrix to gamma space
tautogamma<-function(taumat){
  expmat<-matrix(0L, nrow=nrow(taumat), ncol=ncol(taumat))
  for (i in c(1:nrow(taumat))){
    for (j in c(1:ncol(taumat))){
      if (i!=j){
        expmat[i, j]<- exp(taumat[i, j])/(1+sum(exp(taumat[i,-i])))
      }
      else{
        expmat[i,i]<- 1/(1+sum(exp(taumat[i,-i])))
      }
    }
  }
  expmat
}

#gamma space to tau space 
gammatotau<-function(gammamat){
  logmat<-matrix(0L, nrow=nrow(gammamat), ncol=ncol(gammamat))
  for (i in c(1:nrow(gammamat))){
    for (j in c(1:ncol(gammamat))){
      if (i!=j){
       logmat[i, j]<- log(gammamat[i,j]/gammamat[i,i]) 
      }
      else{
        logmat[i,i]<-NA
      }
    }
  }
  logmat
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