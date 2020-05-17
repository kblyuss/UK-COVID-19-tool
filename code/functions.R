# This file defines all the functions that are used in the simulation

# ----------------------------------------------------------------------------
# SEIR_opt function:
# -----------------------
# Defines the system of differential equations describing the SEIR model 
# INPUT: p - named list of parameter values
# OUTPUT: list of derivatives of each variable

SEIR_opt=function(t,y,p){
  
  with(as.list(p),{
    
    z=rep(0,18*(K1+4*K2+9))
    
    Ntot=rep(0,18)
    
    ind1<-18*seq(0,K1+4*K2+8)

    for (i in 1:18){
      Ntot[i]=sum(y[ind1+i])
    }
    
    EEE=10.8336
    
    yn=matrix(y,nrow=K1+4*K2+9,byrow=TRUE)
    
    Aaux=yn[(K1+2):(K1+K2+1),]
    A=colSums(Aaux)
    
    IRaux=yn[(K1+K2+2):(K1+2*K2+1),]
    IR=colSums(IRaux)
    
    IHaux=yn[(K1+2*K2+2):(K1+3*K2+1),]
    IH=colSums(IHaux)
    
    ICaux=yn[(K1+3*K2+2):(K1+4*K2+1),]
    IC=colSums(ICaux)
    
    Itot=A+IR+IH+IC
    Itots=Itot/Ntot
    
    z[1:18]=-beta*y[1:18]*my_mm(as.matrix(mix_matrix),Itots)/(EEE)
    
    z[19:36]=-z[1:18]-K1*nu*y[19:36] # E1
    
    ind2<-18*seq(2,K1)

    for (i in 1:18){
      z[ind2+i]=K1*nu*y[ind2+i-18]-K1*nu*y[ind2+i] # E2 to E_K1
    }
    
    z[(18*(K1+1)+1):(18*(K1+2))]=pa*K1*nu*y[(18*K1+1):(18*(K1+1))]-K2*g*y[(18*(K1+1)+1):(18*(K1+2))] # A1
    
    ind3<-18*(K1+1)+18*seq(1,K2)
    
    for (i in 1:18){
      z[ind3+i]=K2*g*y[ind3+i-18]-K2*g*y[ind3+i] # A2 to A_K2
    }
    
    z[(18*(K1+K2+1)+1):(18*(K1+K2+2))]=K1*nu*pr*y[(18*K1+1):(18*(K1+1))]-K2*g*y[(18*(K1+K2+1)+1):(18*(K1+K2+2))] # IR1
    
    ind4<-18*(K1+K2+1)+18*seq(1,K2)
    
    for (i in 1:18){
      z[ind4+i]=K2*g*y[ind4+i-18]-K2*g*y[ind4+i] # IR2 to IR_K2
    }
    
    z[(18*(K1+2*K2+1)+1):(18*(K1+2*K2+2))]=K1*nu*ph*y[(18*K1+1):(18*(K1+1))]-K2*g*y[(18*(K1+2*K2+1)+1):(18*(K1+2*K2+2))] #IH1
    
    ind5<-18*(K1+2*K2+1)+18*seq(1,K2)
    
    for (i in 1:18){
      z[ind5+i]=K2*g*y[ind5+i-18]-K2*g*y[ind5+i] # IH2 to IH_K2
    }
    
    z[(18*(K1+3*K2+1)+1):(18*(K1+3*K2+2))]=K1*nu*pc*y[(18*K1+1):(18*(K1+1))]-K2*g*y[(18*(K1+3*K2+1)+1):(18*(K1+3*K2+2))] #IC1
    
    ind6<-18*(K1+3*K2+1)+18*seq(1,K2) 
    
    for (i in 1:18){
      z[ind6+i]=K2*g*y[ind6+i-18]-K2*g*y[ind6+i] # IC2 to IC_K2
    }
    
    z[(18*(K1+4*K2+1)+1):(18*(K1+4*K2+2))]=K2*g*y[(18*(K1+3*K2)+1):(18*(K1+3*K2+1))]-dh*y[(18*(K1+4*K2+1)+1):(18*(K1+4*K2+2))] #Hh
    
    z[(18*(K1+4*K2+2)+1):(18*(K1+4*K2+3))]=dh*y[(18*(K1+4*K2+1)+1):(18*(K1+4*K2+2))] #Rh
    
    z[(18*(K1+4*K2+3)+1):(18*(K1+4*K2+4))]=K2*g*y[(18*(K1+4*K2)+1):(18*(K1+4*K2+1))]-dc*y[(18*(K1+4*K2+3)+1):(18*(K1+4*K2+4))] #Hc
    
    z[(18*(K1+4*K2+4)+1):(18*(K1+4*K2+5))]=dc*y[(18*(K1+4*K2+3)+1):(18*(K1+4*K2+4))]-(xic+mu)*y[(18*(K1+4*K2+4)+1):(18*(K1+4*K2+5))] #Cc
    
    z[(18*(K1+4*K2+5)+1):(18*(K1+4*K2+6))]=xic*y[(18*(K1+4*K2+4)+1):(18*(K1+4*K2+5))] # Rc
    
    z[(18*(K1+4*K2+6)+1):(18*(K1+4*K2+7))]=K2*g*y[(18*(K1+K2)+1):(18*(K1+K2+1))]  #RA
    
    z[(18*(K1+4*K2+7)+1):(18*(K1+4*K2+8))]=K2*g*y[(18*(K1+2*K2)+1):(18*(K1+2*K2+1))] # RR
    
    z[(18*(K1+4*K2+8)+1):(18*(K1+4*K2+9))]=mu*y[(18*(K1+4*K2+4)+1):(18*(K1+4*K2+5))] # DD
    
    return(list(z))
    
  })
}

# ----------------------------------------------------------------------------
# Comix_opt function:
# -----------------------
# Defines the system of differential equations describing the SEIR model 
# INPUT: p - named list of parameter values
# OUTPUT: list of derivatives of each variable

Comix_opt=function(t,y,p){
  
  with(as.list(p),{
    
    z=rep(0,8*(K1+4*K2+9))
    
    Ntot=rep(0,8)
    
      ind1<-8*seq(0,K1+4*K2+8)
    
    for (i in 1:8){
      Ntot[i]=sum(y[ind1+i])
    }
    
    EEE=10.8336
    
    yn=matrix(y,nrow=K1+4*K2+9,byrow=TRUE)
    
    Aaux=yn[(K1+2):(K1+K2+1),]
    A=colSums(Aaux)
    
    IRaux=yn[(K1+K2+2):(K1+2*K2+1),]
    IR=colSums(IRaux)
    
    IHaux=yn[(K1+2*K2+2):(K1+3*K2+1),]
    IH=colSums(IHaux)
    
    ICaux=yn[(K1+3*K2+2):(K1+4*K2+1),]
    IC=colSums(ICaux)
    
    Itot=A+IR+IH+IC
    Itots=Itot/Ntot
    
     z[1:8]=-beta*y[1:8]*my_mm(as.matrix(mix_matrix),Itots)/(EEE)
    
     z[9:16]=-z[1:8]-K1*nu*y[9:16] # E1
    
    ind2<-8*seq(2,K1)
    
    for (i in 1:8){
      z[ind2+i]=K1*nu*y[ind2+i-8]-K1*nu*y[ind2+i] # E2 to E_K1
    }
    
    z[(8*(K1+1)+1):(8*(K1+2))]=pa*K1*nu*y[(8*K1+1):(8*(K1+1))]-K2*g*y[(8*(K1+1)+1):(8*(K1+2))] # A1
    
    ind3<-8*(K1+1)+8*seq(1,K2)
    
    for (i in 1:8){
      z[ind3+i]=K2*g*y[ind3+i-8]-K2*g*y[ind3+i] # A2 to A_K2
    }
    
    z[(8*(K1+K2+1)+1):(8*(K1+K2+2))]=K1*nu*pr*y[(8*K1+1):(8*(K1+1))]-K2*g*y[(8*(K1+K2+1)+1):(8*(K1+K2+2))] # IR1
    
     ind4<-8*(K1+K2+1)+8*seq(1,K2)
    
    for (i in 1:8){
      z[ind4+i]=K2*g*y[ind4+i-8]-K2*g*y[ind4+i] # IR2 to IR_K2
    }
    
     z[(8*(K1+2*K2+1)+1):(8*(K1+2*K2+2))]=K1*nu*ph*y[(8*K1+1):(8*(K1+1))]-K2*g*y[(8*(K1+2*K2+1)+1):(8*(K1+2*K2+2))] #IH1
    
    ind5<-8*(K1+2*K2+1)+8*seq(1,K2)
    
    for (i in 1:8){
      z[ind5+i]=K2*g*y[ind5+i-8]-K2*g*y[ind5+i] # IH2 to IH_K2
    }
    
   z[(8*(K1+3*K2+1)+1):(8*(K1+3*K2+2))]=K1*nu*pc*y[(8*K1+1):(8*(K1+1))]-K2*g*y[(8*(K1+3*K2+1)+1):(8*(K1+3*K2+2))] #IC1
    
    ind6<-8*(K1+3*K2+1)+8*seq(1,K2) 
    
    for (i in 1:8){
      z[ind6+i]=K2*g*y[ind6+i-8]-K2*g*y[ind6+i] # IC2 to IC_K2
    }
    
     z[(8*(K1+4*K2+1)+1):(8*(K1+4*K2+2))]=K2*g*y[(8*(K1+3*K2)+1):(8*(K1+3*K2+1))]-dh*y[(8*(K1+4*K2+1)+1):(8*(K1+4*K2+2))] #Hh
    
    z[(8*(K1+4*K2+2)+1):(8*(K1+4*K2+3))]=dh*y[(8*(K1+4*K2+1)+1):(8*(K1+4*K2+2))] #Rh
    
    z[(8*(K1+4*K2+3)+1):(8*(K1+4*K2+4))]=K2*g*y[(8*(K1+4*K2)+1):(8*(K1+4*K2+1))]-dc*y[(8*(K1+4*K2+3)+1):(8*(K1+4*K2+4))] #Hc
    
    z[(8*(K1+4*K2+4)+1):(8*(K1+4*K2+5))]=dc*y[(8*(K1+4*K2+3)+1):(8*(K1+4*K2+4))]-(xic+mu)*y[(8*(K1+4*K2+4)+1):(8*(K1+4*K2+5))] #Cc
    
    z[(8*(K1+4*K2+5)+1):(8*(K1+4*K2+6))]=xic*y[(8*(K1+4*K2+4)+1):(8*(K1+4*K2+5))] # Rc
    
    z[(8*(K1+4*K2+6)+1):(8*(K1+4*K2+7))]=K2*g*y[(8*(K1+K2)+1):(8*(K1+K2+1))]  #RA
    
    z[(8*(K1+4*K2+7)+1):(8*(K1+4*K2+8))]=K2*g*y[(8*(K1+2*K2)+1):(8*(K1+2*K2+1))] # RR
    
    z[(8*(K1+4*K2+8)+1):(8*(K1+4*K2+9))]=mu*y[(8*(K1+4*K2+4)+1):(8*(K1+4*K2+5))] # DD
    
    return(list(z))
    
  })
}

# ----------------------------------------------------------------------------
# GetModelParams function:
# --------------------
# Function to take the parameters entered by the user and turn them into the rate parameters used by the model
# INPUT: input - structure containing all the user entered information
# OUTPUT: named list consisting of the population size N and another list of the model parameters, pModel

GetModelParams = function(input){
  
  IncubPeriod=input$IncubPeriod  #Incubation period, days
  RecovPeriod=input$RecovPeriod #Duration of recovery period, days
  FracAsym=input$FracAsym/100 #Fraction of infections that are asymptomatic
  HospRecPer=input$HospRecPer #Duration of hospital recovery for non-critical cases, days
  CritProg=input$CritProg #Time to progress to critical care in hospital, days
  CritRec=input$CritRec #Time to recover from critical care, days
  Region=input$Region #Region being studied
  K1=input$K1 # Number of incubation stages
  K2=input$K2 # Number of recovery stages
  beta=input$beta # Transmission rate
  Tstart=7*input$Tstart # Time for introduction of lockdown
  Tfinish=7*input$Tfinish # Time for lifting of lockdown
  Tmax=input$Tmax
  InitInf=input$InitInf

  pClin=c(IncubPeriod=IncubPeriod,FracAsym=FracAsym,RecovPeriod=RecovPeriod,HospRecPer=HospRecPer,CritProg=CritProg,CritRec=CritRec)
  
  # Turn these clinical parameters into the rate constants of the model
  pModel=GetParams_SEIR(pClin)
  
  pModel=c(list("beta"=beta,"Region"=Region,"K1"=K1,"K2"=K2,"Tstart"=Tstart,"Tfinish"=Tfinish,"Tmax"=Tmax,"InitInf"=InitInf),pModel) #
 
  return(list("pModel"=pModel))
  
}

# ----------------------------------------------------------------------------
# GetParams_SEIR function:
# --------------------
# Function to relate the clinical parameters entered by the user into the rate parameters used by the model
# INPUT: pClin - named list of the clinical parameters
# OUTPUT: named list of the model rate parameters, excluding the Betas

GetParams_SEIR = function(pClin){
  
  with(as.list(pClin),{
    
    nu=1/IncubPeriod
    
    g=1/RecovPeriod
    
    dh=1/HospRecPer
    
    dc=1/CritProg
    
    xic=1/CritRec
    
    pa=FracAsym
    
    return(list("nu"=nu,"g"=g,"dh"=dh,"dc"=dc,"xic"=xic,"pa"=pa)) #,"pr"=pr,"ph"=ph,"pc"=pc,"mu"=mu))
  })
  
}

# ----------------------------------------------------------------------------
# GetSpread_SEIR_opt function:
# --------------------
# This function numerically intergrates the system of differential equations for a given set of parameter values, initial conditions, and maximum time
# INPUT: p- named list of parameter values
#        Tmax - max time to integrate for
#        y0 - named list of initial conditions for each variable
# OUTPUT: Dataframe with rows as timepoints and columns as variables

GetSpread_SEIR_opt = function(p,TT,y0){
  
  t = seq(from=0, to=TT, by=1)
  
  fun1=cmpfun(SEIR_opt)
  
  out = ode(y=y0, times=t, func=fun1, parms=p ,method="ode45")
  
  return(list("out"=out))
}

# ----------------------------------------------------------------------------
# GetSpread_Comix_opt function:
# --------------------
# This function numerically intergrates the system of differential equations for a given set of parameter values, initial conditions, and maximum time
# INPUT: p- named list of parameter values
#        Tmax - max time to integrate for
#        y0 - named list of initial conditions for each variable
# OUTPUT: Dataframe with rows as timepoints and columns as variables
GetSpread_Comix_opt = function(p,TT,y0){
  
  t = seq(from=0, to=TT, by=1)
  
  fun1=cmpfun(Comix_opt)
  
  out = ode(y=y0, times=t, func=fun1, parms=p ,method="ode45") # method="ode45")
  
  return(list("out"=out))
}

# ----------------------------------------------------------------------------
# SimSEIRopt function:
# --------------------
# Function to simulate the spread of infection using the model
# INPUT: input - structure containing all the user entered information
# OUTPUT: named list consisting of df - wide format of the timecourse of each variable, N, Ro, r, and doubling time

SimSEIRopt = function(input){
  
  ParamStruct=GetModelParams(input)
  pModel=ParamStruct$pModel
  Tmax=input$Tmax
  
  mix_matrix<-BBC_matrix # use BBC mixing matrix
  
  pr=(1-pModel$pa)*(1-ph_BBC[,2]) # proportion of those who are symptomatic but will recover without requiring hospitalisation
  ph=(1-pModel$pa)*ph_BBC[,2]*(1-pi_BBC[,2]) # proportion of those who require hospitalisation
  pc=(1-pModel$pa)*ph_BBC[,2]*pi_BBC[,2] # proportion of those who will require critical care
  
  mu=pm_BBC[,2]
  
  pModel=c(list("mix_matrix"=mix_matrix,"pr"=pr,"ph"=ph,"pc"=pc,"mu"=mu),pModel)
  
  # Set initial conditions and time interval
  
  X<- data_regions[which(data_regions$Name==input$Region),4:21] # initial populations in the region
  
  Ntot=sum(X)
  
  IC=rep(0,(pModel$K1+4*pModel$K2+9)*18)
  
  ind1<-18*seq(1:pModel$K1)
  
  for (i in 1:18){
      IC[ind1+i]=pModel$InitInf*as.numeric(X[i])/Ntot; # initial number of exposed in each age group
  }
  
  IC[1:18]=as.numeric(X[1:18])
  for (i in 1:18){
      IC[i]=IC[i]-sum(IC[ind1+i])  # initial number of susceptible in each age group
  }
  
  out=GetSpread_SEIR_opt(pModel,Tmax,IC)
  
  sol=out$out
  
  K1=input$K1
  K2=input$K2
  Tmax=input$Tmax
  
  t=sol[,1]
  
  S=apply(sol[,2:19], 1, "sum")
  E=apply(sol[,20:(18*(K1+1)+1)], 1, "sum")
  A=apply(sol[,(18*(K1+1)+2):(18*(K1+K2+1)+1)], 1, "sum")
  IR=apply(sol[,(18*(K1+K2+1)+2):(18*(K1+2*K2+1)+1)], 1, "sum")
  IH=apply(sol[,(18*(K1+2*K2+1)+2):(18*(K1+3*K2+1)+1)], 1, "sum")
  IC=apply(sol[,(18*(K1+3*K2+1)+2):(18*(K1+4*K2+1)+1)], 1, "sum")
  HH=apply(sol[,(18*(K1+4*K2+1)+2):(18*(K1+4*K2+2)+1)], 1, "sum")
  RH=apply(sol[,(18*(K1+4*K2+2)+2):(18*(K1+4*K2+3)+1)], 1, "sum")
  HC=apply(sol[,(18*(K1+4*K2+3)+2):(18*(K1+4*K2+4)+1)], 1, "sum")
  CC=apply(sol[,(18*(K1+4*K2+4)+2):(18*(K1+4*K2+5)+1)], 1, "sum")
  RC=apply(sol[,(18*(K1+4*K2+5)+2):(18*(K1+4*K2+6)+1)], 1, "sum")
  RA=apply(sol[,(18*(K1+4*K2+6)+2):(18*(K1+4*K2+7)+1)], 1, "sum")
  RR=apply(sol[,(18*(K1+4*K2+7)+2):(18*(K1+4*K2+8)+1)], 1, "sum")
  DD=apply(sol[,(18*(K1+4*K2+8)+2):(18*(K1+4*K2+9)+1)], 1, "sum")
  
  DAges=sol[Tmax+1,(18*(K1+4*K2+8)+2):(18*(K1+4*K2+9)+1)]
  
  df=list("t"=t,"S"=S,"E"=E,"A"=A,"IR"=IR,"IH"=IH,"IC"=IC,"HH"=HH,"RH"=RH,"HC"=HC,"CC"=CC,"RC"=RC,"RA"=RA,"RR"=RR,"DD"=DD,"DAges"=DAges)
  
  return(df)
  
}

# ----------------------------------------------------------------------------
# SimQuarOpt:
# --------------------
# Function to simulate the spread of infection using the model
# INPUT: input - structure containing all the user entered information
# OUTPUT: named list consists of solution

SimQuarOpt = function(input){
  
  ParamStruct=GetModelParams(input)
  pModel=ParamStruct$pModel
  
  mix_matrix<-BBC_matrix # use BBC mixing matrix for the baseline
  
  pr=(1-pModel$pa)*(1-ph_BBC[,2]) # proportion of those who are symptomatic but will recover without requiring hospitalisation
  ph=(1-pModel$pa)*ph_BBC[,2]*(1-pi_BBC[,2]) # proportion of those who require hospitalisation
  pc=(1-pModel$pa)*ph_BBC[,2]*pi_BBC[,2] # proportion of those who will require critical care
  
  mu=pm_BBC[,2]
  
  pModel=c(list("mix_matrix"=mix_matrix,"pr"=pr,"ph"=ph,"pc"=pc,"mu"=mu),pModel)
  Tmax=pModel$Tmax
  
  # Set initial conditions and time interval
  
  X<- data_regions[which(data_regions$Name==pModel$Region),4:21] # initial populations in the region
  
  Ntot=sum(X)
  
  IC=rep(0,(pModel$K1+4*pModel$K2+9)*18)
  
  ind1<-18*seq(1:pModel$K1)
  
  for (i in 1:18){
    IC[ind1+i]=pModel$InitInf*as.numeric(X[i])/Ntot; # initial number of exposed in each age group
  }
  
  IC[1:18]=as.numeric(X[1:18])
  for (i in 1:18){
    IC[i]=IC[i]-sum(IC[ind1+i])  # initial number of susceptible in each age group
  }
  
  #run ODEs
  Trun1=pModel$Tstart
  
  out=GetSpread_SEIR_opt(pModel,Trun1,IC)
  
  sol1=out$out
  
  t1=sol1[,1]
  
  N1=length(t1)
  
  IC2=rep(0,(pModel$K1+4*pModel$K2+9)*8)
  
  for (i in 0:(pModel$K1+4*pModel$K2+8)){
    IC2[8*i+1]=sol1[N1,18*i+2]
    IC2[8*i+2]=sol1[N1,18*i+3]+sol1[N1,18*i+4]+0.6*sol1[N1,18*i+5]
    IC2[8*i+3]=0.4*sol1[N1,18*i+5]+sol1[N1,18*i+6]+sol1[N1,18*i+7]
    IC2[8*i+4]=sol1[N1,18*i+8]+sol1[N1,18*i+9]
    IC2[8*i+5]=sol1[N1,18*i+10]+sol1[N1,18*i+11]
    IC2[8*i+6]=sol1[N1,18*i+12]+sol1[N1,18*i+13]
    IC2[8*i+7]=sol1[N1,18*i+14]+sol1[N1,18*i+15]
    IC2[8*i+8]=sol1[N1,18*i+16]+sol1[N1,18*i+17]+sol1[N1,18*i+18]+sol1[N1,18*i+19]
  }
  
  ParamStruct=GetModelParams(input)
  pModel=ParamStruct$pModel
  
  mix_matrix<-Comix_matrix # use CoMix mixing matrix for the lockdown
  
  pr=(1-pModel$pa)*(1-ph_Com[,2]) # proportion of those who are symptomatic but will recover without requiring hospitalisation
  ph=(1-pModel$pa)*ph_Com[,2]*(1-pi_Com[,2]) # proportion of those who require hospitalisation
  pc=(1-pModel$pa)*ph_Com[,2]*pi_Com[,2] # proportion of those who will require critical care
  
  mu=pm_Com[,2]
  
  pModel=c(list("mix_matrix"=mix_matrix,"pr"=pr,"ph"=ph,"pc"=pc,"mu"=mu),pModel)
  
  Trun2=pModel$Tfinish-pModel$Tstart
  
  out=GetSpread_Comix_opt(pModel,Trun2,IC2)
  
  sol2=out$out
  
  t2=sol2[,1]
  
  N2=length(t2)
  
  IC3=rep(0,(input$K1+4*input$K2+9)*18)
  
  for (i in 0:(pModel$K1+4*pModel$K2+8)){
    IC3[18*i+1]=sol2[N2,8*i+2]
    IC3[18*i+2]=5*sol2[N2,8*i+3]/13
    IC3[18*i+3]=5*sol2[N2,8*i+3]/13
    IC3[18*i+4]=3*sol2[N2,8*i+3]/13+2*sol2[N2,8*i+4]/12
    IC3[18*i+5]=5*sol2[N2,8*i+4]/12
    IC3[18*i+6]=5*sol2[N2,8*i+4]/12
    IC3[18*i+7]=0.5*sol2[N2,8*i+5]
    IC3[18*i+8]=0.5*sol2[N2,8*i+5]
    IC3[18*i+9]=0.5*sol2[N2,8*i+6]
    IC3[18*i+10]=0.5*sol2[N2,8*i+6]
    IC3[18*i+11]=0.5*sol2[N2,8*i+7]
    IC3[18*i+12]=0.5*sol2[N2,8*i+7]
    IC3[18*i+13]=0.5*sol2[N2,8*i+8]
    IC3[18*i+14]=0.5*sol2[N2,8*i+8]
    IC3[18*i+15]=sol2[N2,8*i+9]*as.numeric(X[15])/(as.numeric(X[15])+as.numeric(X[16])+as.numeric(X[17])+as.numeric(X[18]))
    IC3[18*i+16]=sol2[N2,8*i+9]*as.numeric(X[16])/(as.numeric(X[15])+as.numeric(X[16])+as.numeric(X[17])+as.numeric(X[18]))
    IC3[18*i+17]=sol2[N2,8*i+9]*as.numeric(X[17])/(as.numeric(X[15])+as.numeric(X[16])+as.numeric(X[17])+as.numeric(X[18]))
    IC3[18*i+18]=sol2[N2,8*i+9]*as.numeric(X[18])/(as.numeric(X[15])+as.numeric(X[16])+as.numeric(X[17])+as.numeric(X[18]))
  }
  
  ParamStruct=GetModelParams(input)
  pModel=ParamStruct$pModel
  
  mix_matrix<-BBC_matrix # use BBC mixing matrix for the baseline after lockdown
  
  pr=(1-pModel$pa)*(1-ph_BBC[,2]) # proportion of those who are symptomatic but will recover without requiring hospitalisation
  ph=(1-pModel$pa)*ph_BBC[,2]*(1-pi_BBC[,2]) # proportion of those who require hospitalisation
  pc=(1-pModel$pa)*ph_BBC[,2]*pi_BBC[,2] # proportion of those who will require critical care
  
  mu=pm_BBC[,2]
  
  pModel=c(list("mix_matrix"=mix_matrix,"pr"=pr,"ph"=ph,"pc"=pc,"mu"=mu),pModel)
  
  Trun3=pModel$Tmax-pModel$Tfinish
  
  out=GetSpread_SEIR_opt(pModel,Trun3,IC3)
  
  sol3=out$out
  
  K1=pModel$K1
  K2=pModel$K2
  Tmax=pModel$Tmax
  
  t12=sol1[,1]
  N12=length(t12)
  S12=apply(sol1[,2:19], 1, "sum")
  E12=apply(sol1[,20:(18*(K1+1)+1)], 1, "sum")
  A12=apply(sol1[,(18*(K1+1)+2):(18*(K1+K2+1)+1)], 1, "sum")
  IR12=apply(sol1[,(18*(K1+K2+1)+2):(18*(K1+2*K2+1)+1)], 1, "sum")
  IH12=apply(sol1[,(18*(K1+2*K2+1)+2):(18*(K1+3*K2+1)+1)], 1, "sum")
  IC12=apply(sol1[,(18*(K1+3*K2+1)+2):(18*(K1+4*K2+1)+1)], 1, "sum")
  HH12=apply(sol1[,(18*(K1+4*K2+1)+2):(18*(K1+4*K2+2)+1)], 1, "sum")
  RH12=apply(sol1[,(18*(K1+4*K2+2)+2):(18*(K1+4*K2+3)+1)], 1, "sum")
  HC12=apply(sol1[,(18*(K1+4*K2+3)+2):(18*(K1+4*K2+4)+1)], 1, "sum")
  CC12=apply(sol1[,(18*(K1+4*K2+4)+2):(18*(K1+4*K2+5)+1)], 1, "sum")
  RC12=apply(sol1[,(18*(K1+4*K2+5)+2):(18*(K1+4*K2+6)+1)], 1, "sum")
  RA12=apply(sol1[,(18*(K1+4*K2+6)+2):(18*(K1+4*K2+7)+1)], 1, "sum")
  RR12=apply(sol1[,(18*(K1+4*K2+7)+2):(18*(K1+4*K2+8)+1)], 1, "sum")
  DD12=apply(sol1[,(18*(K1+4*K2+8)+2):(18*(K1+4*K2+9)+1)], 1, "sum")
  
  t22=sol2[,1]+7*pModel$Tstart
  N22=length(t22)
  S22=apply(sol2[,2:9], 1, "sum")
  E22=apply(sol2[,10:(8*(K1+1)+1)], 1, "sum")
  A22=apply(sol2[,(8*(K1+1)+2):(8*(K1+K2+1)+1)], 1, "sum")
  IR22=apply(sol2[,(8*(K1+K2+1)+2):(8*(K1+2*K2+1)+1)], 1, "sum")
  IH22=apply(sol2[,(8*(K1+2*K2+1)+2):(8*(K1+3*K2+1)+1)], 1, "sum")
  IC22=apply(sol2[,(8*(K1+3*K2+1)+2):(8*(K1+4*K2+1)+1)], 1, "sum")
  HH22=apply(sol2[,(8*(K1+4*K2+1)+2):(8*(K1+4*K2+2)+1)], 1, "sum")
  RH22=apply(sol2[,(8*(K1+4*K2+2)+2):(8*(K1+4*K2+3)+1)], 1, "sum")
  HC22=apply(sol2[,(8*(K1+4*K2+3)+2):(8*(K1+4*K2+4)+1)], 1, "sum")
  CC22=apply(sol2[,(8*(K1+4*K2+4)+2):(8*(K1+4*K2+5)+1)], 1, "sum")
  RC22=apply(sol2[,(8*(K1+4*K2+5)+2):(8*(K1+4*K2+6)+1)], 1, "sum")
  RA22=apply(sol2[,(8*(K1+4*K2+6)+2):(8*(K1+4*K2+7)+1)], 1, "sum")
  RR22=apply(sol2[,(8*(K1+4*K2+7)+2):(8*(K1+4*K2+8)+1)], 1, "sum")
  DD22=apply(sol2[,(8*(K1+4*K2+8)+2):(8*(K1+4*K2+9)+1)], 1, "sum")
  
  t32=sol3[,1]+7*pModel$Tfinish
  N32=length(t32)
  S32=apply(sol3[,2:19], 1, "sum")
  E32=apply(sol3[,20:(18*(K1+1)+1)], 1, "sum")
  A32=apply(sol3[,(18*(K1+1)+2):(18*(K1+K2+1)+1)], 1, "sum")
  IR32=apply(sol3[,(18*(K1+K2+1)+2):(18*(K1+2*K2+1)+1)], 1, "sum")
  IH32=apply(sol3[,(18*(K1+2*K2+1)+2):(18*(K1+3*K2+1)+1)], 1, "sum")
  IC32=apply(sol3[,(18*(K1+3*K2+1)+2):(18*(K1+4*K2+1)+1)], 1, "sum")
  HH32=apply(sol3[,(18*(K1+4*K2+1)+2):(18*(K1+4*K2+2)+1)], 1, "sum")
  RH32=apply(sol3[,(18*(K1+4*K2+2)+2):(18*(K1+4*K2+3)+1)], 1, "sum")
  HC32=apply(sol3[,(18*(K1+4*K2+3)+2):(18*(K1+4*K2+4)+1)], 1, "sum")
  CC32=apply(sol3[,(18*(K1+4*K2+4)+2):(18*(K1+4*K2+5)+1)], 1, "sum")
  RC32=apply(sol3[,(18*(K1+4*K2+5)+2):(18*(K1+4*K2+6)+1)], 1, "sum")
  RA32=apply(sol3[,(18*(K1+4*K2+6)+2):(18*(K1+4*K2+7)+1)], 1, "sum")
  RR32=apply(sol3[,(18*(K1+4*K2+7)+2):(18*(K1+4*K2+8)+1)], 1, "sum")
  DD32=apply(sol3[,(18*(K1+4*K2+8)+2):(18*(K1+4*K2+9)+1)], 1, "sum")
  DAges=sol3[N32,(18*(K1+4*K2+8)+2):(18*(K1+4*K2+9)+1)]
  
  t=c(t12[1:(N12-1)],t22[1:(N22-1)],t32)
  S=c(S12[1:(N12-1)],S22[1:(N22-1)],S32)
  E=c(E12[1:(N12-1)],E22[1:(N22-1)],E32)
  A=c(A12[1:(N12-1)],A22[1:(N22-1)],A32)
  IR=c(IR12[1:(N12-1)],IR22[1:(N22-1)],IR32)
  IH=c(IH12[1:(N12-1)],IH22[1:(N22-1)],IH32)
  IC=c(IC12[1:(N12-1)],IC22[1:(N22-1)],IC32)
  HH=c(HH12[1:(N12-1)],HH22[1:(N22-1)],HH32)
  RH=c(RH12[1:(N12-1)],RH22[1:(N22-1)],RH32)
  HC=c(HC12[1:(N12-1)],HC22[1:(N22-1)],HC32)
  CC=c(CC12[1:(N12-1)],CC22[1:(N22-1)],CC32)
  RC=c(RC12[1:(N12-1)],RC22[1:(N22-1)],RC32)
  RA=c(RA12[1:(N12-1)],RA22[1:(N22-1)],RA32)
  RR=c(RR12[1:(N12-1)],RR22[1:(N22-1)],RR32)
  DD=c(DD12[1:(N12-1)],DD22[1:(N22-1)],DD32)
  
  df=list("t"=t,"S"=S,"E"=E,"A"=A,"IR"=IR,"IH"=IH,"IC"=IC,"HH"=HH,"RH"=RH,"HC"=HC,"CC"=CC,"RC"=RC,"RA"=RA,"RR"=RR,"DD"=DD,"DAges"=DAges)
  
  return(df)
  
}