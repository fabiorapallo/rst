library(Rmpfr)
require(gmp)
library(zoo)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)

#epidemic sizes through simulation
reedfrost = function(p,I0,S0,n) 
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  ##n number of simulated trajectories
  
  S = rep(S0,n)
  St = S
  I = rep(I0,n)
  It = I
  q = 1-p
  time = 0
  
  while (sum(It)>0)
  {
    It = rbinom(n,St,1-(q^It))
    St = St-It
    I = rbind(I, It)
    S = rbind(S, St)
    time = time+1
  }
  I = as.matrix(I)
  S = as.matrix(S)
  list(I0=I0,S0=S0,p=p,n=n,I=I,S=S)
}

#barplot of simulated epidemic sizes
plotepidemicsizes = function(results, title="Sizes of epidemics") 
{
  ##results output of the reedfrost function
  pop = results$I0+results$S0
  tes = table(factor(pop-results$S[nrow(results$S),],levels=1:pop))
  tes = tes/sum(tes)
  df.tes = data.frame(x=1:pop,est.p=as.vector(tes))
  ggplot(df.tes, aes(x=x, y=est.p)) +
    geom_bar(stat="identity") +
    labs(title=title,x="number of cases",y="fraction of trials") +
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15)
    )
}

#exact distribution of the epidemic sizes
finsize = function(p,I0,S0)
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  
  q = 1-p
  G = rep(NA,S0+1)
  U = q^(I0:(I0+S0))
  
  G[1] = 1
  for (k in 1:S0)
  {
    G[k+1] = 1/as.double(factorialZ(k))
    for (j in 0:(k-1)){
      G[k+1] = G[k+1]-q^((I0+j)*(k-j))/as.double(factorialZ(k-j))*G[j+1] 
    }
  }
  
  distr = rep(NA,S0+1)
  for (k in 0:S0)
  {
    distr[k+1] = as.double(factorialZ(S0)/factorialZ(S0-k))*q^((I0+k)*(S0-k))*G[k+1]
  }
  distr = c(rep(0,I0-1),distr)
  return(distr)
}

#barplot of the exact distribution
plotexactsize = function(p,I0,S0, title="Sizes of epidemics") 
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  
  pop = I0+S0
  exd = finsize(p,I0,S0)
  df.exd = data.frame(x=1:pop,ex.p=exd)
  ggplot(df.exd, aes(x=x, y=ex.p)) +
    geom_bar(stat="identity") +
    labs(title = title,x="number of cases",y="probability") +
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15)
    )
}

#Montecarlo and exact distribution on the same plot
compareepidemicsizes = function(p,I0,S0,n, title="Sizes of epidemics") 
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  ##n number of simulated trajectories
  
  results = reedfrost(p,I0,S0,n)
  ed = finsize(p,I0,S0)
  pop = results$I0+results$S0
  tes = table(factor(pop-results$S[nrow(results$S),],levels=1:pop))
  tes = tes/sum(tes)
  df.both = data.frame(x=1:pop,tes=as.vector(tes),ed=ed)
  ggplot(df.both, aes(x=x)) +
    geom_point(aes(y=tes)) +              
    geom_point(aes(y=ed),color="red") +
    labs(title=title, x="number of cases",y="fraction of trials") +
    theme(axis.title=element_text(size=15),axis.text=element_text(size=15)
    )
}


#Montecarlo with sd estimation and exact distribution on the same plot
compareepidemicsizes2 = function(p,I0,S0,n,B, title="Sizes of epidemics") 
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  ##n number of simulated trajectories
  ##B bootstrap replicates to estimate sd
  
  pop = I0+S0
  R0 = p*S0
  res = matrix(rep(NA,pop*B),nrow=B)
  for (b in 1:B)
  {
    rf.results = reedfrost(p,I0,S0,n)
    res[b,] = table(factor(pop-rf.results$S[nrow(rf.results$S),],levels=1:pop))/n
  }
  res = data.frame(res)
  res.long = res %>%
    pivot_longer(cols=everything(), names_to="Category", values_to="Value")
    res.long$Category = factor(res.long$Category, levels=unique(res.long$Category), labels=1:pop)
  
  summary.res = res.long %>%
    group_by(Category) %>%
    summarise(
      Mean=mean(Value),
      SD=sd(Value)
    )
  summary.res$TrueValue = finsize(p,I0,S0)
  mx = floor(pop/5)*5
    g = ggplot(summary.res, aes(x=Category, y=Mean)) +
    geom_bar(stat="identity", fill="grey") +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2) +
    geom_point(aes(y=TrueValue), color="red", size=2) +
    scale_x_discrete(breaks=c(1, seq(5,mx,by=5),pop)) +
    labs(y="probability", x="number of cases") +
    ggtitle(substitute(paste(R[0],'=',R0,', ',s[0],'=',S0))) +
    theme(plot.title=element_text(hjust=0.5,size=15),axis.title=element_text(size=15),axis.text=element_text(size=15))
    return(g)
}


#only Montecarlo for simulations with large pop size
# (simple version used only for TABLE IV)
compareepidemicsizesMCmod = function(p,I0,S0,n,B, title="Sizes of epidemics") 
{
  ##p the probability of infection
  ##I0 initial number of infectives
  ##S0 initial number of susceptibles
  ##n number of simulated trajectories
  ##B bootstrap replicates to estimate sd
  
  pop = I0+S0
  R0 = p*S0
  res = matrix(rep(NA,pop*B),nrow=B)
  for (b in 1:B)
  {
    rf.results = reedfrost(p,I0,S0,n)
    res[b,] = table(factor(pop-rf.results$S[nrow(rf.results$S),],levels=1:pop))/n
  }
  res = data.frame(res)
  val = colMeans(res)
  csm = cumsum(val)
  
  reslt = list(sum(csm<0.50),sum(csm<0.90),csm[0.05*S0])
  return(reslt)
}

##critical point and probabilities -- based on both exact formula and MC
crpp = function(p,I0,S0,n,B)
{
  pop = I0+S0
  crsize = rep(NA,B)
  estprob = rep(NA,B)
  for (b in 1:B)
  {
    results = reedfrost(p,I0,S0,n)
    rm = round(2*(pop/16)+1)
    ds = rollmean(table(factor(pop-results$S[nrow(results$S),],levels=1:pop))/results$n,rm)
    dds = sign(diff(ds))
    i = 1
    cgsn=c()
    for (i in (1:(length(dds)-1)))
      if (((dds[i]==-1)&(dds[i+1]==1))|(dds[i]==0))
      {cgsn=c(cgsn,i+(rm-1)/2)
      }
    if (length(cgsn)>0){
      crsize[b]=round(median(cgsn))
      estprob[b]=sum(pop-results$S[nrow(results$S),]<=crsize[b])/ncol(results$S)
    }
  }
  
  exresults = finsize(p,I0,S0)
  #plot(exresults)
  ddsex = sign(diff(exresults))
  excrsize = NA
  i = 1
  exprob = NA
  for (i in (I0:(length(ddsex)-1)))
    if (((ddsex[i]==-1)&(ddsex[i+1]==1))|(ddsex[i]==0))
    {excrsize=i}
  exprob = sum(exresults[1:excrsize])
  
  rslts=list(median(crsize,na.rm=T),mean(estprob,na.rm=T),sd(estprob,na.rm=T),excrsize,exprob)
  return(rslts)
}


##critical point and probabilities -- only MC
crppMC = function(p,I0,S0,n,B)
{
  pop = I0+S0
  crsize = rep(NA,B)
  estprob = rep(NA,B)
  for (b in 1:B)
  {
    results = reedfrost(p,I0,S0,n)
    ##rm=round(2*(pop/16)+1)
    ds = density((pop-results$S[nrow(results$S),]))
    dds = sign(diff(ds$y))
    i = 1
    cgsn = c()
    for (i in (1:(length(dds)-1)))
      if (((dds[i]==-1)&(dds[i+1]==1))|(dds[i]==0))
      {cgsn = c(cgsn,i)
      }
    if (length(cgsn)>0){
      crsize[b] = round(ds$x[round(median(cgsn))])
      estprob[b] = sum(pop-results$S[nrow(results$S),]<=crsize[b])/ncol(results$S)
    }
  }
  
  rslts=list(median(crsize,na.rm=T),mean(estprob,na.rm=T),sd(estprob,na.rm=T))
  return(rslts)
}



################################
##run  this for TABLE IV
################################

compareepidemicsizesMCmod(1.1/1000,1,1000,5000,1000)  
compareepidemicsizesMCmod(1.08/1000,1,1000,5000,1000)
compareepidemicsizesMCmod(1.06/1000,1,1000,5000,1000)  
compareepidemicsizesMCmod(1.04/1000,1,1000,5000,1000)  
compareepidemicsizesMCmod(1.02/1000,1,1000,5000,1000)  
compareepidemicsizesMCmod(1.0/1000,1,1000,5000,1000) 
compareepidemicsizesMCmod(0.98/1000,1,1000,5000,1000) 
compareepidemicsizesMCmod(0.96/1000,1,1000,5000,1000) 
compareepidemicsizesMCmod(0.94/1000,1,1000,5000,1000) 
compareepidemicsizesMCmod(0.92/1000,1,1000,5000,1000) 
compareepidemicsizesMCmod(0.9/1000,1,1000,5000,1000) 

compareepidemicsizesMCmod(1.1/10000,1,100000,5000,1000)  
compareepidemicsizesMCmod(1.08/10000,1,10000,5000,1000)
compareepidemicsizesMCmod(1.06/10000,1,10000,5000,1000)  
compareepidemicsizesMCmod(1.04/10000,1,10000,5000,1000)  
compareepidemicsizesMCmod(1.02/10000,1,10000,5000,1000)  
compareepidemicsizesMCmod(1.0/10000,1,10000,5000,1000) 
compareepidemicsizesMCmod(0.98/10000,1,10000,5000,1000) 
compareepidemicsizesMCmod(0.96/10000,1,10000,5000,1000) 
compareepidemicsizesMCmod(0.94/10000,1,10000,5000,1000) 
compareepidemicsizesMCmod(0.92/10000,1,10000,5000,1000) 
compareepidemicsizesMCmod(0.9/10000,1,10000,5000,1000) 



################################
##run  this for FIGURE 4
################################

decs = 4
I0 = 1
init = 4
fin = 47
lmv = rep(NA,fin)

for (s0 in init:fin)
{
  lm = 1
  for (d in 1:decs)
  {
    for (i in 1:10)
    {
      x = finsize((lm+i*10^(-d))/s0,I0,s0)
      is_local_min = c(FALSE, 
                       x[2:(length(x)-1)] < x[1:(length(x)-2)] & x[2:(length(x)-1)] < x[3:length(x)],
                       FALSE)
      if (sum(is_local_min)>0)
      {
        lm = lm+(i-1)*10^(-d)
        break
      }
    }
  }
  lmv[s0] = lm
}

print(lmv)
##function of R0
plot(1:fin,lmv)
##function of p (or tilde p)
plot(1:fin,lmv/(1:fin))


# Create a unique data frame
df = data.frame(
  R0 = 1:fin,
  lmv = lmv,
  p_func = lmv/(1:fin)
)

# Plot 1: R0
plot1 = ggplot(df, aes(x=R0, y=lmv)) +
  geom_point(color="blue",size=2) +
  labs(
    title=expression("Critical "*R[0]*" with 1 initial infected"),
    x=expression(s[0]),
    y=expression("Critical "*R[0])
  ) +
  theme_minimal() +
  theme(
    plot.title=element_text(hjust=0.5, size=16, face="bold"),
    axis.title=element_text(size=14),
    axis.text=element_text(size=12)
  )

# Plot 2: tilde p
plot2 = ggplot(df, aes(x = R0, y = p_func)) +
  geom_point(color="red",size=2) +
  labs(
    title=expression("Critical "*tilde(p) * " with 1 initial infected"),
    x=expression(s[0]),
    y=expression("Critical "*tilde(p))
  ) +
  scale_y_continuous(expand=c(0, 0), limits=c(0,0.30)) +
  theme_minimal() +
  theme(
    plot.title=element_text(hjust= 0.5, size=16, face="bold"),
    axis.title=element_text(size=14),
    axis.text=element_text(size=12)
  )

plot1
plot2
# Arrange the two plots side by side
p = grid.arrange(plot2, plot1, ncol=2)
p



#######################################
##run  this for FIGURE 3
#######################################

g1 = compareepidemicsizes2(1.1/32,1,32,5000,2000)  ##R_0=1.1  s_0=32 
g2 = compareepidemicsizes2(1.5/32,1,32,5000,2000)  ##R_0=1.5  s_0=32
g3 = compareepidemicsizes2(2.5/32,1,32,5000,2000)  ##R_0=2.5  s_0=32  


g4 = compareepidemicsizes2(1.1/64,1,64,5000,2000)  ##R_0=1.1  s_0=64
g4MC = compareepidemicsizesMC(1.1/64,1,64,5000,2000)  ##R_0=1.1  s_0=64 only MC
g5 = compareepidemicsizes2(1.5/64,1,64,5000,2000)  ##R_0=1.5  s_0=64
g6 = compareepidemicsizes2(2.5/64,1,64,5000,2000)  ##R_0=2.5  s_0=64  


grid.arrange(g1,g4MC,g2,g5,g3,g6, ncol=2)



#######################################
##run  this for TABLES II AND III
#######################################

crpp(1.2/16,1,16,5000,2000)
crpp(1.2/32,1,32,5000,2000)
crpp(1.2/64,1,64,5000,2000)   ##unstable exact formula

crpp(1.4/16,1,16,5000,2000)
crpp(1.4/32,1,32,5000,2000)
crpp(1.4/64,1,64,5000,2000)   ##unstable exact formula

crpp(1.6/16,1,16,5000,2000)
crpp(1.6/32,1,32,5000,2000)
crpp(1.6/64,1,64,5000,2000) 

crpp(1.8/16,1,16,5000,2000)
crpp(1.8/32,1,32,5000,2000)
crpp(1.8/64,1,64,5000,2000) 

crpp(2/16,1,16,5000,2000)
crpp(2/32,1,32,5000,2000)
crpp(2/64,1,64,5000,2000)

crpp(2.2/16,1,16,5000,2000)
crpp(2.2/32,1,32,5000,2000)
crpp(2.2/64,1,64,5000,2000)

crpp(2.4/16,1,16,5000,2000)
crpp(2.4/32,1,32,5000,2000)
crpp(2.4/64,1,64,5000,2000)

crpp(2.6/16,1,16,5000,2000)
crpp(2.6/32,1,32,5000,2000)
crpp(2.6/64,1,64,5000,2000)

crpp(2.8/16,1,16,5000,2000)
crpp(2.8/32,1,32,5000,2000)
crpp(2.8/64,1,64,5000,2000)

crpp(3/16,1,16,5000,2000)
crpp(3/32,1,32,5000,2000)
crpp(3/64,1,64,5000,2000)


##with two initial infectives

crpp(1.2/16,2,16,5000,2000)
crpp(1.2/32,2,32,5000,2000)
crpp(1.2/64,2,64,5000,2000)   ##unstable exact formula

crpp(1.4/16,2,16,5000,2000)
crpp(1.4/32,2,32,5000,2000)
crpp(1.4/64,2,64,5000,2000)   ##unstable exact formula

crpp(1.6/16,2,16,5000,2000)
crpp(1.6/32,2,32,5000,2000)
crpp(1.6/64,2,64,5000,2000) 

crpp(1.8/16,2,16,5000,2000)
crpp(1.8/32,2,32,5000,2000)
crpp(1.8/64,2,64,5000,2000) 

crpp(2/16,2,16,5000,2000)
crpp(2/32,2,32,5000,2000)
crpp(2/64,2,64,5000,2000)

crpp(2.2/16,2,16,5000,2000)
crpp(2.2/32,2,32,5000,2000)
crpp(2.2/64,2,64,5000,2000)

crpp(2.4/16,2,16,5000,2000)
crpp(2.4/32,2,32,5000,2000)
crpp(2.4/64,2,64,5000,2000)

crpp(2.6/16,2,16,5000,2000)
crpp(2.6/32,2,32,5000,2000)
crpp(2.6/64,2,64,5000,2000)

crpp(2.8/16,2,16,5000,2000)
crpp(2.8/32,2,32,5000,2000)
crpp(2.8/64,2,64,5000,2000)

crpp(3/16,2,16,5000,2000)
crpp(3/32,2,32,5000,2000)
crpp(3/64,2,64,5000,2000)


##only MonteCarlo

crppMC(1.2/128,1,128,5000,2000)
crppMC(1.2/256,1,256,5000,2000)
crppMC(1.2/512,1,512,5000,2000)

crppMC(1.4/128,1,128,5000,2000)
crppMC(1.4/256,1,256,5000,2000)
crppMC(1.4/512,1,512,5000,2000)

crppMC(1.6/128,1,128,5000,2000)
crppMC(1.6/256,1,256,5000,2000)
crppMC(1.6/512,1,512,5000,2000)

crppMC(1.8/128,1,128,5000,2000)
crppMC(1.8/256,1,256,5000,2000)
crppMC(1.8/512,1,512,5000,2000)

crppMC(2/128,1,128,5000,2000)
crppMC(2/256,1,256,5000,2000)
crppMC(2/512,1,512,5000,2000)

crppMC(2.2/128,1,128,5000,2000)
crppMC(2.2/256,1,256,5000,2000)
crppMC(2.2/512,1,512,5000,2000)

crppMC(2.4/128,1,128,5000,2000)
crppMC(2.4/256,1,256,5000,2000)
crppMC(2.4/512,1,512,5000,2000)

crppMC(2.6/128,1,128,5000,2000)
crppMC(2.6/256,1,256,5000,2000)
crppMC(2.6/512,1,512,5000,2000)

crppMC(2.8/128,1,128,5000,2000)
crppMC(2.8/256,1,256,5000,2000)
crppMC(2.8/512,1,512,5000,2000)

crppMC(3/128,1,128,5000,2000)
crppMC(3/256,1,256,5000,2000)
crppMC(3/512,1,512,5000,2000)


##with two initial infectives

crppMC(1.2/128,2,128,5000,2000)
crppMC(1.2/256,2,256,5000,2000)
crppMC(1.2/512,2,512,5000,2000)

crppMC(1.4/128,2,128,5000,2000)
crppMC(1.4/256,2,256,5000,2000)
crppMC(1.4/512,2,512,5000,2000)

crppMC(1.6/128,2,128,5000,2000)
crppMC(1.6/256,2,256,5000,2000)
crppMC(1.6/512,2,512,5000,2000)

crppMC(1.8/128,2,128,5000,2000)
crppMC(1.8/256,2,256,5000,2000)
crppMC(1.8/512,2,512,5000,2000)

crppMC(2/128,2,128,5000,2000)
crppMC(2/256,2,256,5000,2000)
crppMC(2/512,2,512,5000,2000)

crppMC(2.2/128,2,128,5000,2000)
crppMC(2.2/256,2,256,5000,2000)
crppMC(2.2/512,2,512,5000,2000)

crppMC(2.4/128,2,128,5000,2000)
crppMC(2.4/256,2,256,5000,2000)
crppMC(2.4/512,2,512,5000,2000)

crppMC(2.6/128,2,128,5000,2000)
crppMC(2.6/256,2,256,5000,2000)
crppMC(2.6/512,2,512,5000,2000)

crppMC(2.8/128,2,128,5000,2000)
crppMC(2.8/256,2,256,5000,2000)
crppMC(2.8/512,2,512,5000,2000)

crppMC(3/128,2,128,5000,2000)
crppMC(3/256,2,256,5000,2000)
crppMC(3/512,2,512,5000,2000)


