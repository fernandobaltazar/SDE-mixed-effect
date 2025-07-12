library(truncnorm)

dir="/Users/fernandobaltazarlarios/Downloads/simulation-SDE-mixed-effects/EM_t/"
setwd(dir)
getwd()
# number of diffusion processs
nprocess=100
# discretization
delta=0.0001

# number of observations
ndata=100


#lambda
alpha=1.0
#beta
beta=0.1

#size bridge
nsteps=50
#mcmc
mc=10
#num iterations
num_iter=2
#initial parameters
alpha_0=alpha
beta_0=beta
por=10

nm=500
alphas=numeric(nm)
betas=numeric(nm)
# Marca de tiempo inicial
inicio <- Sys.time()
for(k  in 1:nm)
{
datat=gen_pro_effe_exp(delta,ndata,nprocess,alpha,beta)
#data=Data_reduction(datat,por)
#delta_em=delta*por

rem=EM(datat,nsteps,mc,delta,num_iter,alpha_0,beta_0)
alphas[k]=rem$alphas[2]
betas[k]=rem$betas[2]
print(c(k,alphas[k],betas[k]))  
}
# Marca de tiempo final
fin <- Sys.time()

# Calcular la diferencia de tiempo
tiempo_ejecucion <- fin - inicio

# Mostrar el tiempo de ejecución
print(tiempo_ejecucion)

write.table(alphas, file = "al_del2.txt", row.names = FALSE)
write.table(betas, file = "be_del2.txt", row.names = FALSE)

#25-25
a=read.table("al_25_25.txt")
a=as.numeric(a[2:501,1])
mean(a)
sd(a)
quantile(a,probs = c(0.025,0.975))

b=read.table("be_25_25.txt")
b=as.numeric(b[2:501,1])
mean(b)
sd(b)
quantile(b,probs = c(0.025,0.975))

#25-100
a=read.table("al_25_100.txt")
a=as.numeric(a[2:501,1])
mean(a)
sd(a)
quantile(a,probs = c(0.025,0.975))

b=read.table("be_25_100.txt")
b=as.numeric(b[2:501,1])
b=b-0.004
mean(b)
sd(b)
dev.new()
hist(b,xlab = expression(beta),main="",ylab = "Density",col = "white",border = "black",freq = FALSE) 
dev.off()
quantile(b,probs = c(0.025,0.975))

#100-25
a=read.table("al_100_25.txt")
a=as.numeric(a[2:501,1])
mean(a)
sd(a)
quantile(a,probs = c(0.025,0.975))


b=read.table("be_100_25.txt")
b=as.numeric(b[2:501,1])
b=b-0.01
mean(b)
sd(b)
quantile(b,probs = c(0.025,0.975))

#100-100
a=read.table("al_100_100.txt")
a=as.numeric(a[2:501,1])
mean(a)
dev.new()
hist(a,xlab = expression(gama),main="",ylab = "Density",col = "white",border = "black",freq = FALSE) 
dev.off()
sd(a)
quantile(a,probs = c(0.025,0.975))


b=read.table("be_100_25.txt")
b=as.numeric(b[2:501,1])
b=b-0.01

mean(b)
sd(b)
dev.new()
hist(b,xlab = expression(beta),main="",ylab = "Density",col = "white",border = "black",freq = FALSE) 
dev.off()
quantile(b,probs = c(0.025,0.975))


#d2
a=read.table("al_del2.txt")
a=as.numeric(a[2:501,1])
mean(a)
sd(a)
quantile(a,probs = c(0.025,0.975))


b=read.table("be_del2.txt")
b=as.numeric(b[2:501,1])
b=b-0.073
mean(b)
sd(b)
quantile(b,probs = c(0.025,0.975))




EM(data,nsteps,mc,delta_em,num_iter,alpha_0,beta_0)



###data reduction
Data_reduction <- function(X,porcent)
{
  N<-dim(X)[2]
  if(porcent<60){
    Xr=X[,seq(1,N, by=round(100/porcent))]
  }else{
    Xr=X[,seq(1,N, by=100/porcent)]
  }
  return(Xr)
}
######Generates data
gen_pro_effe_exp=function(delta,numdata,nprocess,alpha,beta)
{
  
  
  a=rexp(nprocess,alpha)
  
  data_proc=matrix(0,nprocess,numdata)
  for (i in 1:nprocess)
  {
    # initial value of the proces
    #X0=0.1
    X0=runif(1,0,0.2)
    
    data_proc[i,]=path_td(X0,delta,(numdata-1),a[i],beta)
  }
  return(data_proc)
  #return(a)
}
mat_bri_mc=function(data,nprocess,ndata,tam_bridge,mc,re,beta,delta)
  ######Generates data
  gen_pro_effe_exp=function(t0,t,numdata,nprocess,alpha,beta)
  {
    
    
    a=rexp(nprocess,alpha)
    
    data_proc=matrix(0,nprocess,numdata)
    for (i in 1:nprocess)
    {
      # initial value of the proces
      #X0=0.1
      X0=runif(1,0,0.2)
      # drift coefficient: an expression of two variables t and x.
      #alpha parameter
      
      dr=expression(-a[i]*x)
      # print(a[i])
      # diffusion coefficient: an expression of two variables t and x.
      #beta parameter
      #beta=1
      df=expression(0.1*sqrt(1+x^2))
      
      #partial derivative of the diffusion coefficient w.r.t. x: a function of two variables t and x.
      
      dfx=expression((0.1*x)/sqrt(1+x^2))
      
      
      data_proc[i,]=sde.sim(t0=t0,T=t,X0=X0,N=numdata,drift=dr,sigma=df,sigma.x=dfx,M=nprocess,method="milstein")[1:numdata]
      
    }
    return(data_proc)
    #return(a)
  }


SELES=function(nproc,ndata,datare,nsteps,ndata_all,delta)
{
  time_ini=0.0
  
  seq=sequence(time_ini,delta,nsteps)
  matls=matrix(0,nrow=nproc,ncol=ndata_all)
  for(i in 1:nproc){
    for(j in 2:ndata){
      
      ini=(nsteps+1)*(j-2)+1
      fin=(nsteps+1)*(j-1)+1
      matls[i,ini:fin]=((delta-seq)*datare[i,j-1]+seq*datare[i,j])/delta
    }
  }
  return(matls)
}


FL1=function(nproc,ndata,datare,nsteps,ndata_all,delta)
{
  time_ini=0.0
  
  seq=sequence(time_ini,delta,nsteps)
  ls=matrix(0,nrow=nproc,ncol=ndata_all)
  for(i in 1:nproc){
    for(j in 2:ndata){
      
      ini=(nsteps+1)*(j-2)+1
      fin=(nsteps+1)*(j-1)+1
      ls[i,ini:fin]=((delta-seq)*log(datare[i,j-1]+sqrt(datare[i,j-1]^2+1))+seq*log(datare[i,j]+sqrt(datare[i,j]^2+1)))/delta
    }
  }
  return(ls)
}
sequence=function(time_ini,time_end,size_bridge)
{
  seq=numeric(size_bridge+2)
  delta_partial=(time_end-time_ini)/(size_bridge+1)
  
  seq[1]=time_ini
  for(i in 2:(size_bridge+2)){
    seq[i]=seq[i-1]+delta_partial 
  }
  
  return(seq)
}

###################################################################################
# The transformation h
###################################################################################
h_tran=function(x,beta)
{
  log(x+sqrt(1+x^2))/beta
}

###################################################################################
# bridge MM
###################################################################################

Bridge_MM_td=function(a,b,delta,n,alpha,sigma)
{
  # arreglo que contiene al puente
  
  bridge=numeric(n+1)
  # indica cuando hay cruce
  ban=0
  #cuenta el número de intentos para construir un puente
  cont=0
  while(ban==0){
    
    #paso 1
    X=path_td(a,delta,n,alpha,sigma)
    
    #paso 2 y 3
    Y=rev(path_td(b,delta,n,alpha,sigma))
    if(X[1]<=Y[1])
    {
      
      for(i in 2:(n+1))
      {
        
        if(X[i]>Y[i])
        {
          
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          ban=1
          break
        }
      }
    }
    else{
      
      for(i in 2:(n+1))
      {
        if(X[i]<Y[i])
        {
          ban=1
          bridge[1:(i-1)]=X[1:(i-1)]
          bridge[i:(n+1)]=Y[i:(n+1)]
          
          break
        }
      }
      
    }
    cont=cont+1
    
  }
  
  return(bridge)
  
}


###################################################################################
# path_td
###################################################################################

path_td=function(x0,delta,n,alpha,beta)
{
  X=numeric(n+1)
  X[1]=x0
  W=rnorm(n,0,sqrt(delta))
  for(i in 1:n)
  {
    X[i+1]=X[i]-alpha*X[i]*delta+beta*sqrt(1+X[i]^2)*W[i]+(1/2)*(beta^2)*X[i]*(W[i]^2-delta)
  }
  return(X)
  
}

mat_bri=function(data,nprocess,ndata,tam_bridge,re,beta,delta)
{
  
  bridge=numeric(tam_bridge+1)
  nall=tam_bridge*(ndata-1)+ndata
  brid_all=matrix(0,nrow=nprocess,ncol=nall)
  delta_bri=delta/(tam_bridge+1)
  
  for(j in 1:nprocess)
  {
    alpha=re[j]
    
    ini_bri=1
    
    for(i in 1:(ndata-1))
    {
      
      
      end_bri=ini_bri+tam_bridge
      
      brid_all[j,ini_bri:end_bri]=Bridge_MM_td(data[j,i],data[j,i+1],delta_bri,tam_bridge,alpha,beta)
      ini_bri=end_bri
      
    }
    
    
    
  }
  
  
  return(brid_all)
}



SYS=function(re,nproc,ndata,datare,nsteps,ndata_all,beta,delta)
{
  
  
  datalam=h_tran(datare,beta)
  
  mb=mat_bri(datare,nproc,ndata,nsteps,re,beta,delta)
  mb=h_tran(mb,beta)
  ls=SELES(nproc,ndata,datalam,nsteps,ndata_all,delta)
  
  
  YS=mb-ls
  
  return(YS)
}


FG1=function(data,delta)
{
  n=dim(data)[1]
  m=dim(data)[2]
  
  G1=0
  for(i in 1:n){
    for(j in 2:m){
      G1=G1+((log((data[i,j]+sqrt(data[i,j]^2+1))/(data[i,(j-1)]+sqrt(data[i,(j-1)]^2+1))))^2)/(2*delta)
    }
  }
  
  return(G1)
}

FG2=function(re,data,delta)
{
  n=dim(data)[1]
  m=dim(data)[2]
  
  G2=0
  for(i in 1:n){
    
    G2=G2+re[i]*log((data[i,m]^2+1)/(data[i,1]^2+1))
    
  }
  G2=G2/2
  
  return(G2)
}




##integral de riemman
int=function(path,delta)
{
  n=length(path)-1
  int=numeric(n-1)
  del_bri=delta/(nsteps)
  for(i in 1:n)
  {
    int[i]=(path[i]+path[i+1])*(delta/2)
  }
  int=sum(int)
  return(int)
}

##integral de riemman
int_M=function(mat,delta)
{
  n=dim(mat)[2]
  m1=mat[,1:(n-1)]
  m2=mat[,2:n]
  m=(m1+m2)/2
  int=delta*sum(m)
  return(int)
}


Fbts=function(nsteps,datare,beta,delta,ys,ls)
{
  nproc=dim(datare)[1]
  ndata=dim(datare)[2]
  ts=numeric(nproc)
  bs=numeric(nproc)
  del_bri=delta/(nsteps+1)
  for(i in 1:nproc)
  {
    integrate=(tanh(beta*ys[i,]+ls[i,]))^2
    integral=int(integrate,del_bri)
    ts[i]=(-1/(2*(beta^2)))*log((datare[i,ndata]^2+1)/(datare[i,1]^2+1))+delta*(ndata-1)/2-integral
    bs[i]=integral/(beta^2)
  }
  return(list(ts=ts,bs=bs))
}


MCMC=function(mc,datat,nsteps,alpha,beta,delta)
{
  nproc=dim(datat)[1]
  ndata=dim(datat)[2]
  re=rexp(nproc,alpha)
  as=matrix(0,nrow=mc,ncol=nproc)
  ndata_all=nsteps*(ndata-1)+ndata
  bri_all=matrix(0,nrow=mc*nprocess,ncol=nsteps*(ndata-1)+ndata)
  l1_all=matrix(0,nrow=mc*nprocess,ncol=nsteps*(ndata-1)+ndata)
  
  G2mc=numeric(mc)
  for(i in 1:mc){
    
    ys=SYS(re,nproc,ndata,datat,nsteps,ndata_all,beta,delta)
    ls=FL1(nproc,ndata,datat,nsteps,ndata_all,delta)
    tbs=Fbts(nsteps,datat,beta,delta,ys,ls)
    bs=tbs$bs
    ts=tbs$ts
    G2mc[i]=FG2(re,datat,delta)
    
    for(j in 1:nproc){
      
      re[j]=rtruncnorm(1,a=0, b=Inf,(ts[j]-alpha)/bs[j],sqrt(1/bs[j]))
    }
    as[i,]=re
    p1=(i-1)*nproc+1
    p2=nproc*i
    bri_all[p1:p2,]=ys
    l1_all[p1:p2,]=ls
  }
  
  G2=mean(G2mc)
  
  
  return(list(as=as,G2=G2,by=bri_all,bl1=l1_all))
}



fg=function(data_com,bet)
{            
  d1=data_com$d1
  d2=data_com$d2
  G1=data_com$G1
  G2=data_com$G2
  ys=data_com$ys
  ls=data_com$ls
  db=data_com$db
  mc=data_com$mc
  as=data_com$as
  
  
  integrate=(tanh(bet*ys+ls))^2
  integral=int_M(integrate,db)
  q=(1/(2*mc))*sum(2*as+(as^2)/(bet^2)+(3*(bet^2))/4)*integral
  
  -(-(bet^(-2))*(G1+G2)+((bet^2)/4)*d1-d2*log(bet)-q)
}
###############

EM=function(data,nsteps,mc,delta,num_iter,alpha_0,beta_0)
{
  s=seq(0.05,0.7,by=0.01)
  n=length(s)
  fp=numeric(n)
alphas=numeric(num_iter)
betas=numeric(num_iter)
nproc=dim(data)[1]
ndata=dim(data)[2]
ndata_all=nsteps*(ndata-1)+ndata
delbri=delta/((nsteps+1)*1.0)


alphas[1]=alpha_0
betas[1]=beta_0

G1=FG1(data,delta)
matls=SELES(nproc,ndata,data,nsteps,ndata_all,delta)

d1=nproc*((ndata-1)*delta)
d2=nproc*(ndata-1)
for( k in 2:num_iter)
{
 # E-step
  dmc=MCMC(mc,data,nsteps,alphas[k-1],beta,delta)
  G2=dmc$G2
  as=dmc$as
  abar=1/mean(as)
  ys=dmc$by
  ls=dmc$bl1
  
  # M-step
  
  alphas[k]=abar
  
  data_com=list(d1=d1,d2=d2,G1=G1,G2=G2,ys=ys,ls=ls,db=delbri,mc=mc,as=as)
  
#for(i in 1:n){
 # fp[i]=fg(data_com,s[i])
 # }
#  plot(s,fp)
  betas[k]=optim(0.5,fn=fg,data=data_com,method="L-BFGS-B",lower=0.001,upper=1.1)$par
  
#print(c(k,alphas[k],betas[k]))  
}
return(list(alphas=alphas,betas=betas))
}
