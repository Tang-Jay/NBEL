#===============================================
#                 Codes for NBEL 
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')
library('MASS')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')
SARAR_Data <- function(m,indexs,Delta=0,nsim=2000){
  # ==================修改参数====================
  rou1 = 0.85 
  rou2 = 0.15
  n = m*m
  ps = round(3*n^indexs)
  # ==================函数准备====================
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v =c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  # ==================数据准备====================
  Wnb = cell2nb(m,m,type='queen')
  Ws = nb2listw(Wnb)
  Wn = listw2mat(Ws)
  Mn = Wn
  In = diag(n)
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  # 估计方程准备
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  g = diag(Gnn)
  h = diag(Hnn)
  # ==================开始模拟====================
  # 数据保存空间
  Sta = matrix(NA,length(ps),3)
  colnames(Sta)=c('p','EL','NBEL')
  EL = c()
  NBEL = c()
  for(p in ps){
    mu = rep(0,p)
    Sigma = diag(p)
    
    f1 = 0
    f2 = 0
    for(i in 1:nsim){
      Xn = mvrnorm(n,mu,Sigma)
      En = rnorm(n,0,1)
      beta0 = matrix(rep(1,p),p,1)
      Yn = Ani%*%Xn%*%beta0 + Ani%*%Bni%*%En
      
      beta1 = beta0 + Delta
      sigma2 = 1 + Delta
      b = t(Xn)%*%t(Bn)
      s = Bn%*%Wn%*%Ani%*%Xn%*%beta1
      e = c(Bn%*%An%*%Yn - Bn%*%Xn%*%beta1)
      
      # 估计方程赋值
      z = matrix(NA,nrow=n,ncol=p+3)
      z[,1:p] = b*e
      z[,p+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
      z[,p+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
      z[,p+3] = e*e - rep(sigma2, n)
      
      # 计算EL值
      lam = lambdaChen(z)
      el = 2*sum(log(1+t(lam)%*%t(z)))
      if(el>qchisq(0.95,p+3)) f1 = f1 + 1
      
      # 计算NBEL值
      z = z%*%rep(1,p+3)
      lam = lambdaChen(z)                						
      nbel = 2*sum(log(1+t(lam)%*%t(z)))
      if(nbel>qchisq(0.95,1)) f2 = f2 + 1
    }
    cat(n,p,' ',f1/nsim,' ',f2/nsim,'\n')
    # 储存数据
    EL = c(EL,f1/nsim)
    NBEL = c(NBEL,f2/nsim)
  }
  
  Sta[,1] = ps
  Sta[,2] = EL
  Sta[,3] = NBEL
}
#===============================================
#                   赋值运算            
#=============================================== 
ms = c(10, 13)
Deltas = c(0, 0.1, 1, 2)
indexs = c(0, 0.1, 0.3)
for(m in ms){
  for(Delta in Deltas){
    cat('Delta =',Delta,'\n')
    cat('n  ','p ',' EL  ','  NBEL ','\n')
    SARAR_Data(m,indexs,Delta,nsim=2000)
    cat('\n')
  }
}



