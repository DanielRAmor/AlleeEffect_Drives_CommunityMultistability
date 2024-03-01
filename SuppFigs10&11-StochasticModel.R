# This code runs stochastic simulations for the dynamics of the 6-species community
#composed of the microbial species Se, Cn, Sa, Mc, Lp and Cp described in the manuscript

library(GillespieSSA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# Parameter values:
parms <- c(r1=0.3, r2=0.3, r3 =0.3, r4=0.3, r5=0.3, r6=0.3,
           k1=10^3, k2=10^3,k3=10^3, k4=10^3, k5=10^3, k6=10^3, 
           a12 = 0.4, a13 = 0.4, a14 = 1.2, a15 = 0.1, a16 = 0.1, 
           a21 = 1.05, a23 = 1.2, a24 = 1.6, a25 = 0.3, a26 = 0.3,
           a31 = 0.4, a32 = 0.4, a34 = 0.8, a35 = 0.1, a36 = 0.1,
           a41 = 0.1, a42 = 0.1, a43 = 0.1, a45 = 0.1, a46 = 0.1,
           a51 = 1.3, a52 = 1.3, a53 = 1.3, a54 = 1.5, a56 = 0.1,
           a61 = 1.3, a62 = 1.3, a63 = 1.3, a64 = 1.5, a65 = 0.1,
           c1=0.3, c3=0.27, c4=0.05)
#where r_i are the growth rates, k_i the carrying capacities, 
#a_ij is the strength of the (competitive) interaction from species j to species i, 
#c1 and c2 capture the strength of a weak Allee effect, and c4 the strength of a strong Allee effect 
#the changes associated to a weak or strong Allee effect involve  
#a change in the sign of the corresponding rates, which is captured in the matrix 'nu' below

#species convention:
#Y1=Se,Y2=Cn, Y3=Sa, Y4=Mc, Y5=Lp, Y6=Cp

  
#stochastic reaction rates for Gillespie (Tau-leap method) simulations under the Allee effect
a <- c("r1*(1-c1)*(Y1^2)/k1","r1*c1*Y1","r1*(Y1^3)/(k1^2)","a12*r1*Y1*Y2/k2","a13*r1*Y1*Y3/k3","a14*r1*Y1*Y4/k4","a15*r1*Y1*Y5/k5","a16*r1*Y1*Y6/k6",
       "r2*Y2","r2*(Y2^2)/k2","a21*r2*Y1*Y2/k1","a23*r2*Y3*Y2/k3","a24*r2*Y4*Y2/k4","a25*r2*Y2*Y5/k5","a26*r2*Y2*Y6/k6",
       "r3*(1-c3)*(Y3^2)/k3","r3*c3*Y3","r3*(Y3^3)/(k3^2)","a31*r3*Y1*Y3/k1","a32*r3*Y2*Y3/k2","a34*r3*Y4*Y3/k4","a35*r3*Y3*Y5/k5","a36*r3*Y3*Y6/k6",
       "r4*(1+c4)*(Y4^2)/k4","r4*c4*Y4","r4*(Y4^3)/(k4^2)","a41*r4*Y1*Y4/k1","a42*r4*Y2*Y4/k2","a43*r4*Y3*Y4/k3","a45*r4*Y4*Y5/k5","a46*r4*Y4*Y6/k6",
       "r5*Y5","r5*(Y5^2)/k5","a51*r5*Y5*Y1/k1","a52*r5*Y5*Y2/k2","a53*r5*Y3*Y5/k3","a54*r5*Y4*Y5/k4","a56*r5*Y6*Y5/k6",
       "r6*Y6","r6*(Y6^2)/k6","a61*r6*Y6*Y1/k1","a62*r6*Y6*Y2/k2","a63*r6*Y3*Y6/k3","a64*r6*Y4*Y6/k4","a65*r6*Y6*Y5/k5")
  
#presence and sign of the reactions for simulations under the Allee effect
nu <- matrix(c(+1,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1),
             nrow=6,byrow=TRUE)
  
 #final simulated time
  Tfinal = 150 

#Loops to generate several simulations for each initial condition in terms of species abundances 
  #initial abundance
  x0 <- c(Y1=40,Y2=40, Y3=40, Y4=40, Y5=40, Y6=40)
  for (i in 1:12) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="EvenInoculum-Allee")
    jpeg(file=paste("Even-Allee-long", i,".jpeg"))
  #  pdf(paste("Even-Allee-long", i,".pdf"), width=6.6, height=4.5)
      print(ssa.plot(out))
    dev.off()
  }  
    
  x0 <- c(Y1=200,Y2=12, Y3=12, Y4=12,Y5=12, Y6=12)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="SeUnevenInoculum-Allee")
    
    jpeg(file=paste("SeUneven-Allee-long", i,".jpeg"))
#    pdf(paste("SeUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  
  x0 <- c(Y1=12,Y2=200, Y3=12, Y4=12,Y5=12, Y6=12)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="CnUnevenInoculum-Allee")
    
    jpeg(file=paste("CnUneven-Allee-long", i,".jpeg"))
    #pdf(paste("CnUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=12,Y2=12, Y3=200, Y4=12,Y5=12, Y6=12)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="SaUnevenInoculum-Allee")
    jpeg(file=paste("SaUneven-Allee-long", i,".jpeg"))
    #pdf(paste("SaUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=12,Y2=12, Y3=12, Y4=200,Y5=12, Y6=12)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="McUnevenInoculum-Allee")
    
    jpeg(file=paste("McUneven-Allee-long", i,".jpeg"))
    #pdf(paste("McUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=12,Y2=12, Y3=12, Y4=12,Y5=200, Y6=12)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="LpUnevenInoculum-Allee")
    
    jpeg(file=paste("LpUneven-Allee-long", i,".jpeg"))
    #pdf(paste("LpUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=12,Y2=12, Y3=12, Y4=12,Y5=12, Y6=200)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="CpUnevenInoculum-Allee")
    
    jpeg(file=paste("CpUneven-Allee-long", i,".jpeg"))
    #pdf(paste("CpUneven-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  
  
  
# Parameter values under reduced Allee effect
parms <- c(r1=0.3, r2=0.3, r3 =0.3, r4=0.3, r5=0.3, r6=0.3,
             k1=10^3, k2=10^3,k3=10^3, k4=10^3, k5=10^3, k6=10^3, 
             a12 = 0.4, a13 = 0.4, a14 = 1.2, a15 = 0.1, a16 = 0.1, 
             a21 = 1.05, a23 = 1.2, a24 = 1.6, a25 = 0.3, a26 = 0.3,
             a31 = 0.4, a32 = 0.4, a34 = 0.8, a35 = 0.1, a36 = 0.1,
             a41 = 0.1, a42 = 0.1, a43 = 0.1, a45 = 0.1, a46 = 0.1,
             a51 = 1.3, a52 = 1.3, a53 = 1.3, a54 = 1.5, a56 = 0.1,
             a61 = 1.3, a62 = 1.3, a63 = 1.3, a64 = 1.5, a65 = 0.1,
             c1=0.3, c3=0.27, c4=0.02)
#stochastic reaction rates for Gillespie (Tau-leap method) simulations under reduced Allee effect conditions
a <- c("r1*Y1","r1*(Y1^2)/k1","a12*r1*Y1*Y2/k2","a13*r1*Y1*Y3/k3","a14*r1*Y1*Y4/k4","a15*r1*Y1*Y5/k5","a16*r1*Y1*Y6/k6",
       "r2*Y2","r2*(Y2^2)/k2","a21*r2*Y1*Y2/k1","a23*r2*Y3*Y2/k3","a24*r2*Y4*Y2/k4","a25*r2*Y2*Y5/k5","a26*r2*Y2*Y6/k6",
       "r3*Y3","r3*(Y3^2)/k3","a31*r3*Y1*Y3/k1","a32*r3*Y2*Y3/k2","a34*r3*Y4*Y3/k4","a35*r3*Y3*Y5/k5","a36*r3*Y3*Y6/k6",
       "r4*(1+c4)*(Y4^2)/k4","r4*c4*Y4","r4*(Y4^3)/(k4^2)","a41*r4*Y1*Y4/k1","a42*r4*Y2*Y4/k2","a43*r4*Y3*Y4/k3","a45*r4*Y4*Y5/k5","a46*r4*Y4*Y6/k4",
       "r5*Y5","r5*(Y5^2)/k5","a51*r5*Y5*Y1/k1","a52*r5*Y5*Y2/k2","a53*r5*Y3*Y5/k3","a54*r5*Y4*Y5/k4","a56*r5*Y6*Y5/k6",
       "r6*Y6","r6*(Y6^2)/k6","a61*r6*Y6*Y1/k1","a62*r6*Y6*Y2/k2","a63*r6*Y3*Y6/k3","a64*r6*Y4*Y6/k4","a65*r6*Y6*Y5/k5")
 
#presence and sign of the reactions for simulations under the Allee effect
nu <- matrix(c(+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,-1,-1,-1,-1,-1),
             nrow=6,byrow=TRUE)

#Loops to generate several simulations for each initial condition in terms of species abundances 
#initial abundance
  x0 <- c(Y1=40,Y2=40, Y3=40, Y4=40, Y5=40, Y6=40)
  for (i in 1:12) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="EvenInoculum-Allee")
    jpeg(file=paste("Even-Allee-long", i,".jpeg"))
    #  pdf(paste("Even-Allee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  
  
  x0 <- c(Y1=200,Y2=15, Y3=15, Y4=15,Y5=15, Y6=15)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="SeUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Se", i,".jpeg"))
    #pdf(paste("NoAllee-Se", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  
  x0 <- c(Y1=15,Y2=200, Y3=15, Y4=15,Y5=15, Y6=15)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="CnUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Cn", i,".jpeg"))
    #pdf(paste("CnUneven-NoAllee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=15,Y2=15, Y3=200, Y4=15,Y5=15, Y6=15)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="SaUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Sa", i,".jpeg"))
    #pdf(paste("SaUneven-NoAllee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=15,Y2=15, Y3=15, Y4=200,Y5=15, Y6=15)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="McUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Mc", i,".jpeg"))
    #pdf(paste("McUneven-NoAllee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=15,Y2=15, Y3=15, Y4=15,Y5=200, Y6=15)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="LpUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Lp", i,".jpeg"))
    #pdf(paste("LpUneven-NoAllee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
  x0 <- c(Y1=15,Y2=15, Y3=15, Y4=15,Y5=15, Y6=200)
  for (i in 1:3) {
    i = i+1
    out <- ssa(x0,a,nu,parms,tf=Tfinal,method=ssa.otl(),simName="CpUnevenInoculum")
    
    jpeg(file=paste("NoAllee-Cp", i,".jpeg"))
    #pdf(paste("CpUneven-NoAllee-long", i,".pdf"), width=6.6, height=4.5)
    print(ssa.plot(out))
    dev.off()
  }  
