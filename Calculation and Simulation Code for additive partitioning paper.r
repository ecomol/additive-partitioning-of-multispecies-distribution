################################################################################################################
################################################################################################################
#Additive partitioning and numerical simulation of conspecific-encounter index
#Written by Chen's lab (2024)
#2024-03-20
########################################


########################################
## v index, also used in Solow (2000) and Chen et al. (2019) Methods in Ecology and Evolution
## z is a vector of sequential species' label in the sample
v.est = function(z){
  m = length(z)
  TT = sum(z[-m]==z[-1])
  v = TT/(m-1)
  est = v
  est
}

########################################
##estimation of variance of v index
##derived from Chen et al. (2019) Methods in Ecology and Evolution
v.est.se = function(z){
  m = length(z)
  
  intra.pair = (z[-m]==z[-1])
  TT = sum(intra.pair)
  v = TT/(m-1)
  
  intra.3 = (intra.pair[-(m-1)]==1 & intra.pair[-1]==1)
  TTT = sum(intra.3)
  t3 = TTT/(m-2)
  

  if(v<1) d=1/(1-v) else d=0
  est.var = ((m-1)*(v-v^2)+(m-1)*2*d*(t3-v^2)*(t3>v^2))/(m-1)^2
  
  return(sqrt(est.var))
}


########################################
#exact standard error of estimated v
##derived from Chen et al. (2019) Methods in Ecology and Evolution
v.est.EXACT.se = function(z){
  m = length(z)
  
  intra.pair = (z[-m]==z[-1])
  TT = sum(intra.pair)
  v = TT/(m-1)
  
  intra.3 = (intra.pair[-(m-1)]==1 & intra.pair[-1]==1)
  TTT = sum(intra.3)
  t3 = TTT/(m-2)
  
  
  k = m-1
  d=1/(k)*((k)/(1-v)+(v^k-1)/(1-v)^2)
  est.var = ((m-1)*(v-v^2)+(m-1)*2*d*(t3-v^2)*(t3>v^2))/(m-1)^2
  
  return(sqrt(est.var))
}

########################################################
#make the multivariate aggregate transition matrix based on given parameters
#the first element of pars is the positive non-independent parameter pai:
#which has a value ranging from 0 to 1
#the second and third elements are p and q respectively, p/(p+q) is the true detection rate
make.Pai.matrix=function(pars)
{
S=length(pars)-1
pai=pars[1]
p=pars[-1]
p=p/sum(p)
#
Q=matrix(0,nrow=S,ncol=S)
for(i in 1:S)
{
for(j in 1:S)
{
if(i==j)
{
Q[i,j]=(1-pai)*p[i]+pai
}else
{
Q[i,j]=(1-pai)*p[j]
}
}#j
}#i
#
return(Q)
}#
########################################################



#########################################################
#simulate muti-species Pai model (positive non-independent model)
#p is a vector of species relative abundance
#np is the non-independent parameter
#using equal-time or equal-distance scheme
sim.Pai<-function(steps=1000,p,np)
{
S=length(p)
Q=make.Pai.matrix(c(np,p))
#
v=sample(1:S,1,replace=FALSE,prob=p)
p0=rep(0,1,S)
p0[v]=1 #initial probability vector
#####################
for(i in 1:steps)
{
p1=as.vector(p0%*%Q)
cp=cumsum(p1)
r=runif(1)
id=as.integer(rank(c(r,cp)))
###########
#next-step probability vector
p0=rep(0,1,S)
p0[id[1]]=1
v=c(v,id[1])
###########
}#i
#####################
return(v)
}#
#########
#


########################################
##derived from Chen et al. (2019) Methods in Ecology and Evolution
## and Chen et al. (2023) Frontiers in Plant Science
#### pi estimator
## z is a vector of sequential species' label in the sample
pi.est = function(z){
  m = length(z)
  Xi=table(z) ## Xi is the species abundance data
  TT = sum(z[-m]==z[-1])
  v = TT/(m-1)
  b = m*(1-v)
  Y=m^2-sum(Xi^2)
  if((m+1)^2*b^2>8*b*Y) {
    est = 1-((m+1)*b+sqrt((m+1)^2*b^2-8*b*Y))/2/Y
  } else {
    est = v
  }
  if(v==1 || est > v){
    est = v
  }
  if(est<0) est = 0
  est
}







################################################################################################################
################################################################################################################
################################################################################################################
####################################NUMERICAL SIMULATIONS#######################################################
################################################################
################################################################
#TWO-DIMENSIONL SIMULATIONS WITH MANY REPLICATES
#using lognormal distribution to simulate regional abundance variability
#while inertia value pai varies from 0 to 1
delta=seq(0.1,5,len=30)
ppai=seq(0,.95,len=30)
dp=expand.grid(delta,ppai)
simlen=300
reptime=50
Snum=10
################################
ONE<-function(id)
{
pi=rlnorm(Snum, meanlog = 0, sdlog =dp[id,1])
z1=sim.Pai(simlen,pi,dp[id,2])
return(z1)
}#
mat=array(0,dim=c(simlen+1,dim(dp)[1],reptime))
for(i in 1:reptime)
{
m=as.matrix(1:dim(dp)[1])
res=apply(m,1,ONE)
mat[,,i]=res
print(i)
}
################################
################################
total=a=b=vector()
for(i in 1:dim(dp)[1])
{
t0=a0=b0=vector()
for(j in 1:reptime)
{
z1=mat[,i,j]
t0[j]=v.est(z1)
############################
#aggregation partitioning:
a0[j]=pi.est(z1) #spatial autocorrelation effect
b0[j]=v.est(z1)-pi.est(z1) #abundance variability effect
}#j
total=rbind(total,t0)
a=rbind(a,a0)
b=rbind(b,b0)
}#i
#
################################
################################################################
z=z1=z2=matrix(0,nrow=length(delta),ncol=length(ppai))
for(i in 1:length(delta))
{
for(j in 1:length(ppai))
{
id=which(dp[,1]==delta[i] & dp[,2]==ppai[j])
z[i,j]=mean(total[id,])
z1[i,j]=mean(a[id,])
z2[i,j]=mean(b[id,])
}#j
}#i
################################
windows()
par(mfrow=c(2,2))
par(cex.main = 1,cex.lab=1.5)
fields::image.plot(delta,ppai,z,main="A) multi-species distributional aggregation",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z1,main="B) local inertia component",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z2,main="C) regional abundance variability component",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z2/(1-z1),main="D) true regional abundance variability component",
xlab=expression(delta),ylab=expression(omega))
#PERFECT FOR PUBLICATION!!!!!
#PERFECT FOR PUBLICATION!!!!!
################################################################





################################################################
################################################################
#TWO-DIMENSIONL SIMULATIONS WITH MANY REPLICATES
#using uniform distribution along with an extreme vlaue to simulate regional abundance variability
#for characterizing extremely uneven SAD
#while inertia value pai varies from 0 to 1
delta=seq(0.5,100,len=30)
ppai=seq(0,.95,len=30)
dp=expand.grid(delta,ppai)
simlen=300
reptime=50
Snum=10
################################
ONE<-function(id)
{
pi=c(dp[id,1],runif(Snum-1))
z1=sim.Pai(simlen,pi,dp[id,2])
return(z1)
}#
mat=array(0,dim=c(simlen+1,dim(dp)[1],reptime))
for(i in 1:reptime)
{
m=as.matrix(1:dim(dp)[1])
res=apply(m,1,ONE)
mat[,,i]=res
print(i)
}
################################
################################
total=a=b=vector()
for(i in 1:dim(dp)[1])
{
t0=a0=b0=vector()
for(j in 1:reptime)
{
z1=mat[,i,j]
t0[j]=v.est(z1)
############################
#aggregation partitioning:
a0[j]=pi.est(z1) #spatial autocorrelation effect
b0[j]=v.est(z1)-pi.est(z1) #abundance variability effect
}#j
total=rbind(total,t0)
a=rbind(a,a0)
b=rbind(b,b0)
}#i
#
################################
################################################################
z=z1=z2=matrix(0,nrow=length(delta),ncol=length(ppai))
for(i in 1:length(delta))
{
for(j in 1:length(ppai))
{
id=which(dp[,1]==delta[i] & dp[,2]==ppai[j])
z[i,j]=mean(total[id,])
z1[i,j]=mean(a[id,])
z2[i,j]=mean(b[id,])
}#j
}#i
################################
windows()
par(mfrow=c(2,2))
par(cex.main = 1,cex.lab=1.5)
fields::image.plot(delta,ppai,z,main="A) multi-species distributional aggregation",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z1,main="B) local inertia component",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z2,main="C) regional abundance variability component",
xlab=expression(delta),ylab=expression(omega))
fields::image.plot(delta,ppai,z2/(1-z1),main="D) true regional abundance variability component",
xlab=expression(delta),ylab=expression(omega))
#PERFECT FOR PUBLICATION!!!!!
#PERFECT FOR PUBLICATION!!!!!
################################################################





################################################################
################################################################
#TWO-DIMENSIONL SIMULATIONS WITH MANY REPLICATES
#using Symmetric Dirichlet model to simulate regional abundance variability
#for characterizing extremely uneven SAD
#while inertia value pai varies from 0 to 1
library(MCMCpack)
library(dirmult)
##############
#random variate generation
rdirichlet<-function (n = 1, alpha) 
{
    Gam <- matrix(0, n, length(alpha))
    for (i in 1:length(alpha)) Gam[, i] <- rgamma(n, shape = alpha[i])
    Gam/rowSums(Gam)
}
#################
#avoid NAs:
rdirichlet1<-function (n = 1, alpha) 
{
rv=rdirichlet(n,alpha)
if(n==1)
{
rv=t(as.matrix(rv))
}
#
ids=which(is.na(rv[,1]))
if(length(ids)>0)
{
for(i in 1:length(ids))
{
rv[ids[i],]=rep(0,1,length(alpha))
rv[ids[i],sample(1:length(alpha),1)]=1
}#i
}#
#
return(rv)
}
################################################################
################################################################
delta=seq(0.01,1,len=30)
ppai=seq(0,.95,len=30)
dp=expand.grid(delta,ppai)
simlen=50
reptime=50
Snum=10
################################
ONE<-function(id)
{
pi=rdirichlet1(n=1,rep(dp[id,1],1,Snum))
z1=sim.Pai(simlen,pi,dp[id,2])
return(z1)
}#
mat=array(0,dim=c(simlen+1,dim(dp)[1],reptime))
for(i in 1:reptime)
{
m=as.matrix(1:dim(dp)[1])
res=apply(m,1,ONE)
mat[,,i]=res
print(i)
}
################################
################################
total=a=b=vector()
for(i in 1:dim(dp)[1])
{
t0=a0=b0=vector()
for(j in 1:reptime)
{
z1=mat[,i,j]
t0[j]=v.est(z1)
############################
#aggregation partitioning:
a0[j]=pi.est(z1) #spatial autocorrelation effect
b0[j]=v.est(z1)-pi.est(z1) #abundance variability effect
}#j
total=rbind(total,t0)
a=rbind(a,a0)
b=rbind(b,b0)
}#i
#
################################
################################################################
z=z1=z2=matrix(0,nrow=length(delta),ncol=length(ppai))
for(i in 1:length(delta))
{
for(j in 1:length(ppai))
{
id=which(dp[,1]==delta[i] & dp[,2]==ppai[j])
z[i,j]=mean(total[id,])
z1[i,j]=mean(a[id,])
z2[i,j]=mean(b[id,])
}#j
}#i
################################
windows()
par(mfrow=c(2,2))
par(cex.main = 1,cex.lab=1.5)
fields::image.plot(delta,ppai,z,main="A) multi-species distributional aggregation",
xlab=expression(alpha),ylab=expression(omega))
fields::image.plot(delta,ppai,z1,main="B) local inertia component",
xlab=expression(alpha),ylab=expression(omega))
fields::image.plot(delta,ppai,z2,main="C) regional abundance variability component",
xlab=expression(alpha),ylab=expression(omega))
fields::image.plot(delta,ppai,z2/(1-z1),main="D) true regional abundance variability component",
xlab=expression(alpha),ylab=expression(omega))
#PERFECT FOR PUBLICATION!!!!!
#PERFECT FOR PUBLICATION!!!!!
################################################################