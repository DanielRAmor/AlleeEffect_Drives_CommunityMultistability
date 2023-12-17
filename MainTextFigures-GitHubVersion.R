#Code for a 2-species LV model with one or both species subject to an Allee effect
library(deSolve)
library(ggplot2)
library(phaseR)
library(rootSolve)


##The microbial species' color code for the figures:
#Cn	#f9ac95
#Mc	#b5f1b4
#Sa	#ffee93
#Se	#a0ced9
#Lp	#CAA787
#Cp	#F2C7F9

#The color code for community outcomes:
#Bistability	#C12727
#Exclusion	#333333
#Coexistance #3C5CC2

#My Functions
#Time derivatives of species abundances
  # Species 1 under Allee effect:
  LVdyn <- function (t, y, mu) {
    list(c(
      #With Allee
      mu[1]*y[1] * ((mu[3] - y[1])*(y[1]-mu[7]) - mu[5] * y[2] ) / mu[3] - mu[8] * y[1],
      #No Allee:
     # mu[1]*y[1] * ((mu[3] - y[1]) - mu[5] * y[2] ) / mu[3],# - mu[] * y[1], 
      
      mu[2]*y[2] * (mu[4] - y[2] - mu[6] * y[1] ) / mu[4]  - mu[8] * y[2]
      # for quadratic Allee effect use (y[i]-A) * (1-y[i]) / Ki
      
        ))
  }
  
  #No Allee effect
  LVdynNoAllee <- function (t, y, mu) {
    list(c(
      #No Allee:
      mu[1]*y[1] * ((mu[3] - y[1]) - mu[5] * y[2] ) / mu[3],# - mu[] * y[1], 
      
      mu[2]*y[2] * (mu[4] - y[2] - mu[6] * y[1] ) / mu[4]  - mu[8] * y[2]
    ))
  }
  
  
  
  # 2 species under Allee effect
  LV_2Allee <- function (t, y, mu) {
    list(c(
      #With Allee
      mu[1]*y[1] * ((mu[3] - y[1])*(y[1]-mu[7]) - mu[5] * y[2] ) / mu[3] - mu[8] * y[1],
      #No Allee:
      # mu[1]*y[1] * ((mu[3] - y[1]) - mu[5] * y[2] ) / mu[3],# - mu[] * y[1], 
      
      mu[2]*y[2] * ((mu[4] - y[2])*(y[2]-mu[9]) - mu[6] * y[1] ) / mu[4]  - mu[8] * y[2]
      # for quadratic Allee effect use (y[i]-A) * (1-y[i]) / Ki
      
    ))
  }
  
#some parameter values
  pars <- c(0.3, 0.3, #growth rates
            1.0, 1.0, #carrying capacities (here normalized)
            0.4, #interaction towards species 1 
            1.3, #interaction towards species 2
            -0.2, #Allee effect
            0.00) #death rate    
  Tf = 1000
  stepT = 0.1
  yini = c(0.0, 0.04)
#check running a simulation
  LV.simulation <- ode(y = yini, func = LVdyn,
                       times = seq(0,Tf,by = stepT), parms = pars)
  
#Phase diagram 
  #Preparing a background for our plots
    # from line 79 to line 129 we just prepare a data frame and other variables that we may later use for our plots. 
    #data frame to store the results
    a = c(1.0,1.0)
    phaseX <- data.frame(a,a,a,a)
    colnames(phaseX) = c('xini', 'yini', 'xf','yf')
    #Simulations to generate the phase space begin here
    #resolution of the phase space (Nsims * Nsims)
    Nsims = 10
    Tf = 1000
    stepT = 0.1
    yini = c(0.0, 0.0)
    for (i in 1:(Nsims+1)){
      #explore the axis for initial abundance of F
      if (i > 1) {yini[1] = yini[1] + (1.2/Nsims)} #where 4 is the order of magnitudes that will be covered
      
      for (j in 1:(Nsims+1)) {
        if (j < 2) {yini[2] = 0.0
        } else {
          #explore the axis for initial abundance of F
          yini[2] = yini[2] +(1.2/Nsims) #where 4 is the order of magnitudes that will be covered
        }
        #run a simulation
        LV.simulation <- ode(y = yini, func = LVdyn,
                             times = seq(0,Tf,by = stepT), parms = pars)
        #export the result of the simulation at the final time
        phaseX[(i-1)*Nsims + j,1] = yini[1]
        phaseX[(i-1)*Nsims + j,2] = yini[2]
        phaseX[(i-1)*Nsims + j,3] = LV.simulation[Tf/stepT,2]
        phaseX[(i-1)*Nsims + j,4] = LV.simulation[Tf/stepT,3]
        
      }
    }
    
    #compute the fraction of species F in the outcome of the simulations
    phaseX[5] = phaseX[3] / (phaseX[3] + phaseX[4]) 
    colnames(phaseX)[5] = c('FinalSaFraction')
    ploti <- ggplot(data = phaseX, aes(x= xini, y= yini ))
    background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
    titulos<- labs(title="", x=expression(paste("Normalized abundance Se")), y=expression(paste("Normalized abundance Cn")))
    escalaX <- scale_x_continuous(expand= c(0,0), limits=c(0,1.1), breaks=c(0.0,0.5,1.0))
    escalaY <- scale_y_continuous(expand= c(0,0), limits=c(-0,1.1), breaks = c(0.0,0.5,1.0)) 
    couleur <- scale_fill_gradient(low = "#947B64", high = "#DBA43B", limits=c(0.0000, 1.0), guide="colorbar")
    rast <- geom_tile(aes(fill = FinalSaFraction))
    # nullcline with Allee effect
    null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = '#a0ced9', linetype = "dashed", linewidth = 1.3 )
    # nullcline with no Allee effect
    #null1 <- geom_function(fun = function(x) ((1-x)/pars[5]), color = 'deepskyblue4', linetype = "dashed", size = 1.3)
    null2 <- geom_function(fun = function(x) ((1-x*pars[6])-pars[8]/pars[2]), color = '#f9ac95', linetype = "dashed", linewidth = 1.3)
    plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2  
    #plotii
  #The background for our plots is ready
   
    
##CODE FOR FIGURE 3
    #PANEL WITH ALLEE EFFECT
    pars <- c(0.3, 0.3, #growth rates, Cn is faster grower
              1.0, 1.0, #carrying capacities 
              0.8, #interaction towards species 1 (Se). Supernatant data -> positive interaction
              1.1, #interaction towards species 2 (Cn)
              -0.3, #Allee effect
              0.00) #death rate    
    background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
    #nullcline with Allee effect
    null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = '#a0ced9', linetype = "dashed", linewidth = 1.3 )
    #the other nullcline
    null2 <- geom_function(fun = function(x) ((1-x*pars[6])-pars[8]/pars[2]), color = '#f9ac95', linetype = "dashed", linewidth = 1.3)
    
    #Trajectories:
    yini = c(0.1, 0.2)
    LV.traj1 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
    colnames(LV.traj1) = c("time", "xini", "yini")
    traj1 <- geom_path(data=LV.traj1, color = 'black')
    
    yini = c(0.35, 0.1)
    LV.traj2 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
    colnames(LV.traj2) = c("time", "xini", "yini")
    traj2 <- geom_path(data=LV.traj2, color = 'black')
    
    yini = c(0.4, 0.65)
    LV.traj3 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
    colnames(LV.traj3) = c("time", "xini", "yini")
    traj3 <- geom_path(data=LV.traj3, color = 'black')
    
    yini = c(0.6, 0.6)
    LV.traj4 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
    colnames(LV.traj4) = c("time", "xini", "yini")
    traj4 <- geom_path(data=LV.traj4, color = 'black')
    
    #Direction arrows
    dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                         arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                         aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                               xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                               y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                               yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                            ) 
                        )
    dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                         arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                         aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                               xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                               y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                               yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                            ) 
                         )
    dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                         arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                         aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                               xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                               y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                               yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                            ) 
                         )
    dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                         arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                         aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                               xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                               y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                               yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                            ) 
                        )
    
    #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      Xsteady = polyroot(c(c1, b1, a1))
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])-pars[8]/pars[2]
      
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      
      SteadyXaxis <- data.frame (xini = c(1),
                                   yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                   yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      #make the plot
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
              traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 + PtsAxis +
        puntos1 + puntos2 + puntos3 
      plotii
      #export the plot
      pdf(paste("TestFig3D-panel2-white.pdf"), width=4.5, height=4.5)
        print(plotii)
      dev.off()
    

  #PANEL WITHOUT ALLEE EFFECT
      pars <- c(0.3, 0.3, #growth rates, Cn is faster grower
                1.0, 1.0, #carrying capacities 
                0.8, #interaction towards species 1 (Se). Supernatant data -> positive interaction
                1.1, #interaction towards species 2 (Cn)
                -0.3, #Allee effect
                0.00) #death rate    
      
      LVdynNoAllee <- function (t, y, mu) {
        list(c(
          #With Allee
          #mu[1]*y[1] * ((mu[3] - y[1])*(y[1]-mu[7]) - mu[5] * y[2] ) / mu[3] - mu[8] * y[1],
          #No Allee:
           mu[1]*y[1] * ((mu[3] - y[1]) - mu[5] * y[2] ) / mu[3],# - mu[] * y[1], 
          
          mu[2]*y[2] * (mu[4] - y[2] - mu[6] * y[1] ) / mu[4]  - mu[8] * y[2]
          # for quadratic Allee effect use (y[i]-A) * (1-y[i]) / Ki
          
        ))
      }
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      null1 <- geom_function(fun = function(x) (1-x)/pars[5], color = '#a0ced9', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) (1-x*pars[6]), color = '#f9ac95', linetype = "dashed", linewidth = 1.3)
      
      #Steady state points
      Xsteady = (1-pars[5])/(1-pars[5]*pars[6])
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])
      
      #take only the real, positive solutions:
      print(c(Xsteady, Ysteady))
      print(c(Re(Xsteady),Re(Ysteady)))
      Steady1 <- data.frame (xini = c(Re(Xsteady)),
                             yini = c(Re(Ysteady)))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #make the plot
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
                puntos1 + puntos2 + puntos3 + PtsAxis +
                traj1 + traj2 + traj3 + traj4 +
                dir1 + dir2 + dir3 + dir4
      plotii
      #export the plot
      pdf(paste("TestFig3D-panel1-white.pdf"), width=4.5, height=4.5)
        print(plotii)
      dev.off()
      
      
#CODE FOR FIGURE 5 
  ##Panel 1 LEFT  (1 species w Allee)
      titulos<- labs(title="", x=expression(paste("N1")), y=expression(paste("N2")))
      pars <- c(0.3, 0.3, #growth rates, Cn is faster grower
                1.0, 1.0, #carrying capacities 
                0.8, #interaction towards species 1 (Se). Supernatant data -> positive interaction
                1.1, #interaction towards species 2 (Cn)
                -0.3, #Allee effect
                0.00) #death rate    
      #nullcline with Allee effect
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) ((1-x*pars[6])-pars[8]/pars[2]), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      Xsteady = polyroot(c(c1, b1, a1))
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])-pars[8]/pars[2]
      
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
        traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 + PtsAxis +
        puntos1 + puntos2 + puntos3 
      plotii
      
      pdf(paste("Fig5B-panel1-Left.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
 ##PANEL 1 Right
      pars <- c(0.2, 0.3, #growth rates, Cn is faster grower
                1.0, 1.0, #carrying capacities 
                0.8, #interaction towards species 1 (Se). Supernatant data -> positive interaction
                1.1, #interaction towards species 2 (Cn)
                -0.3, #Allee effect
                0.00) #death rate
      
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      
      #nullcline with NO Allee effect and NO death rate
      null1 <- geom_function(fun = function(x) (1-x)/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) (1-x*pars[6]), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      
      #Steady state points
      Xsteady = (1-pars[5])/(1-pars[5]*pars[6])
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])
      #take only the real, positive solutions:
      print(c(Xsteady, Ysteady))
      print(c(Re(Xsteady),Re(Ysteady)))
      Steady1 <- data.frame (xini = c(Re(Xsteady)),
                             yini = c(Re(Ysteady)))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
        puntos1 + puntos2 + puntos3 + PtsAxis +
        traj1 + traj2 + traj3 + traj4 +
        dir1 + dir2 + dir3 + dir4
      plotii
      
      pdf(paste("Fig5B-panel1-right.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
           
 ##Panel 2 Left (Allee effect on 2 species)
      pars <- c(0.2, 0.3, #growth rates, Se is faster grower
                1.0, 1.0, #carrying capacities 
                0.2, #interaction towards species 1 (Mc). Supernatant data -> positive interaction
                1.2, #interaction towards species 2 (Se)
                0.1, #Allee effect
                0.00, # death rate
                -0.3) #Allee effect on species 2
      
      null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) ((1+pars[9] + sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      null3 <- geom_function(fun = function(x) ((1+pars[9] - sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
    # substract null1 to null 2 (or null3) to find the roots (steady states)
      cross1 <- function(x) ((1+pars[9] + sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2) - ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5]
      cross2 <- function(x) ((1+pars[9] - sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2) - ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5]
      uniroot.all(cross1, c(0,1))
      uniroot.all(cross2, c(0,1))
      Xsteady = uniroot.all(cross1, c(0,1))
      #the y value of the steady state is therefore:
      Ysteady = ((1-Xsteady)*(Xsteady-pars[7])-(pars[8]/pars[1]))/pars[5]
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      #make the plot
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 + null3 +
        traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 +
        puntos1 + puntos2 + puntos3 + PtsAxis
      plotii
      
      pdf(paste("Fig5b-Panel2-left.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
 ## Panel 2 Right      
      pars <- c(0.2, 0.3, #growth rates, Se is faster grower
                1.0, 1.0, #carrying capacities 
                0.2, #interaction towards species 1 (Mc). Supernatant data -> positive interaction
                1.2, #interaction towards species 2 (Se)
                0.1, #Allee effect
                0.00) # death rate
      #nullcline with Allee effect
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      #null1 <- geom_function(fun = function(x) ((1-x)/pars[5]), color = 'deepskyblue4', linetype = "dashed", size = 1.3)
      null2 <- geom_function(fun = function(x) ((1-x*pars[6])-pars[8]/pars[2]), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.03),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.03)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.03),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.03)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      
      Xsteady = polyroot(c(c1, b1, a1))
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])-pars[8]/pars[2]
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      PtsAxis <- coord_cartesian(clip = 'off')
      
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
        traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 +
        puntos1 + puntos2 + puntos3 + PtsAxis
      plotii
      
      pdf(paste("Fig5B-panel2-Right.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
  ## Panel 3 right -- Mc vs Se
      pars <- c(0.2, 0.3, #growth rates, Sa is faster grower
                1.0, 1.0, #carrying capacities 
                0.2, #interaction towards species 1 (Mc). Supernatant data -> positive interaction
                0.9, #interaction towards species 2 (Sa)
                0.1, #Allee effect
                0.00) # death rate
      #nullcline with Allee effect
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) ((1-x*pars[6])-pars[8]/pars[2]), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdyn, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.03),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.03)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.03),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.03)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      Xsteady = polyroot(c(c1, b1, a1))
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])-pars[8]/pars[2]
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      
      Steady4 <- data.frame (xini = c(Re(Xsteady[2])),
                             yini = c(Re(Ysteady[2])))
      puntos4 <- geom_point(data = Steady4, fill = "black", shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      
      PtsAxis <- coord_cartesian(clip = 'off')
      
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
        traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 +
        puntos1 + puntos2 + puntos3 + puntos4 + PtsAxis
      plotii
      
      pdf(paste("Fig5B-panel3-Right.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
      
 ## Panel 3 Left -- Sa vs Se
      pars <- c(0.3, 0.3, #growth rates, Se is faster grower
                1.0, 1.0, #carrying capacities 
                0.4, #interaction towards species 1 (Sa). Supernatant data -> positive interaction
                0.4, #interaction towards species 2 (Se)
                -0.3, #Allee effect
                0.00, # death rate
                -0.3) #Allee effect on species 2
      
      null1 <- geom_function(fun = function(x) ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      #null1 <- geom_function(fun = function(x) ((1-x)/pars[5]), color = 'deepskyblue4', linetype = "dashed", size = 1.3)
      null2 <- geom_function(fun = function(x) ((1+pars[9] + sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      null3 <- geom_function(fun = function(x) ((1+pars[9] - sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      # substract null1 to null 2 (or null3) to find the roots (steady states)
      cross1 <- function(x) ((1+pars[9] + sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2) - ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5]
      cross2 <- function(x) ((1+pars[9] - sqrt((1+pars[9])^2-4*(pars[9]+pars[6]*x)))/2) - ((1-x)*(x-pars[7])-(pars[8]/pars[1]))/pars[5]
      uniroot.all(cross1, c(0,1))
      uniroot.all(cross2, c(0,1))
      Xsteady = uniroot.all(cross1, c(0,1))
      #the y value of the steady state is therefore:
      Ysteady = ((1-Xsteady)*(Xsteady-pars[7])-(pars[8]/pars[1]))/pars[5]
      #take only the real, positive solutions:
      print(c(Xsteady[1], Ysteady[1]))
      print(c(Re(Xsteady[2]),Re(Ysteady[2])))
      Steady1 <- data.frame (xini = c(Re(Xsteady[1])),
                             yini = c(Re(Ysteady[1])))
      puntos1 <- geom_point(data = Steady1, shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, fill = "black", shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.4, 0.65)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LV_2Allee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.02),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.02),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.02)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      #Steady state points
      #the two nullclines cross at (x,y) so that x is a root of this polynomial:
      # (-1/pars[5])*x^2 + (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])*x + (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))
      a1 = (-1/pars[5])
      b1 = (1/pars[5] + pars[7]/pars[5] + pars[8]/(pars[5]*pars[1]) + pars[6])
      c1 = (pars[8]/pars[2] -1 -pars[7]/pars[5] - pars[8]/(pars[5]*pars[1]))   
      
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      #C12727
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 + null3 +
        traj1 + traj2 + traj3 + traj4 + dir1 + dir2 + dir3 + dir4 +
        puntos1 + puntos2 + puntos3 + PtsAxis
      plotii
      
      pdf(paste("Fig5b-Panel3-left.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
 ## Panel 4 Right 
      pars <- c(0.3, 0.3, #growth rates, Se is faster grower
                1.0, 1.0, #carrying capacities 
                0.4, #interaction towards species 1 (Sa). Supernatant data -> positive interaction
                0.4, #interaction towards species 2 (Se)
                -0.3, #Allee effect
                0.00, # death rate
                -0.3) #Allee effect on species 2
      background <- theme_bw() + theme (panel.background = element_rect(fill= 'white'), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
      
      #nullcline with NO Allee effect and NO death rate
      null1 <- geom_function(fun = function(x) (1-x)/pars[5], color = 'darkturquoise', linetype = "dashed", linewidth = 1.3 )
      null2 <- geom_function(fun = function(x) (1-x*pars[6]), color = 'darkorange', linetype = "dashed", linewidth = 1.3)
      
      #Steady state points
      Xsteady = (1-pars[5])/(1-pars[5]*pars[6])
      #the y value of the steady state is therefore:
      Ysteady = (1-Xsteady*pars[6])
      #take only the real, positive solutions:
      print(c(Xsteady, Ysteady))
      print(c(Re(Xsteady),Re(Ysteady)))
      Steady1 <- data.frame (xini = c(Re(Xsteady)),
                             yini = c(Re(Ysteady)))
      puntos1 <- geom_point(data = Steady1, fill = "black", shape =21, stat = "identity", position = "identity", size =3, stroke = 1.5)
      SteadyXaxis <- data.frame (xini = c(1),
                                 yini = c(0))
      puntos2 <- geom_point(data = SteadyXaxis, shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      SteadyYaxis <- data.frame (xini = c(0),
                                 yini = c(1))
      puntos3 <- geom_point(data = SteadyYaxis, shape =21, stat = "identity", position = "identity", size = 3, stroke = 1.5)
      PtsAxis <- coord_cartesian(clip = 'off')
      
      #Trajectories:
      yini = c(0.1, 0.2)
      LV.traj1 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj1) = c("time", "xini", "yini")
      traj1 <- geom_path(data=LV.traj1, color = 'black')
      
      yini = c(0.35, 0.1)
      LV.traj2 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj2) = c("time", "xini", "yini")
      traj2 <- geom_path(data=LV.traj2, color = 'black')
      
      yini = c(0.5, 0.9)
      LV.traj3 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj3) = c("time", "xini", "yini")
      traj3 <- geom_path(data=LV.traj3, color = 'black')
      
      yini = c(0.6, 0.6)
      LV.traj4 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj4) = c("time", "xini", "yini")
      traj4 <- geom_path(data=LV.traj4, color = 'black')
      
      yini = c(0.9, 0.5)
      LV.traj5 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj5) = c("time", "xini", "yini")
      traj5 <- geom_path(data=LV.traj5, color = 'black')
      
      yini = c(0.85, 0.85)
      LV.traj6 <- as.data.frame(ode(y = yini, func = LVdynNoAllee, times = seq(0,Tf,by = stepT), parms = pars))
      colnames(LV.traj6) = c("time", "xini", "yini")
      traj6 <- geom_path(data=LV.traj6, color = 'black')
      
      #Direction arrows
      dir1 <- geom_segment(data=LV.traj1, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj1[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj1[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj1[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir2 <- geom_segment(data=LV.traj2, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj2[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj2[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj2[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir3 <- geom_segment(data=LV.traj3, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj3[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj3[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj3[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj3[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir4 <- geom_segment(data=LV.traj4, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj4[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj4[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj4[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      dir5 <- geom_segment(data=LV.traj5, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj5[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj5[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj5[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj5[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      dir6 <- geom_segment(data=LV.traj6, size = 0.3, color = 'black',
                           arrow = arrow(length=unit(.3, 'cm')), #, lwd=3),  #, type = 'closed'),
                           aes(  x = LV.traj6[as.integer(nrow(LV.traj1)*0.01),2], 
                                 xend = LV.traj6[as.integer(nrow(LV.traj1)*0.01)+2,2], 
                                 y = LV.traj6[as.integer(nrow(LV.traj1)*0.01),3], 
                                 yend = LV.traj6[as.integer(nrow(LV.traj1)*0.01)+2,3]
                           ) 
      )
      
      plotii <- ploti + background + titulos + couleur + escalaX + escalaY  + null1 + null2 +
        puntos1 + puntos2 + puntos3 + PtsAxis +
        traj1 + traj2 + traj3 + traj4 + traj5 + traj6 +
        dir1 + dir2 + dir3 + dir4 + dir5 + dir6 
      plotii
      
      pdf(paste("Fig5B-panel4-right.pdf"), width=4.5, height=4.5)
      print(plotii)
      dev.off()
      
      
  