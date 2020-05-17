rm(list=ls())
library(deSolve)
library(reshape)
library(googlesheets4)
library(plotly)
sheets_deauth()
setwd("/Applications/COVID/COVID19seir/ICST")
source("code/functions.R")
library(compiler)
library(Rcpp)
library(shiny)
library(rmarkdown)

# ----------------------------------------------------------------------------
# Load data:
# --------------------

function(input, output, session) {
  
  # Reactive function that computes the timecourse of epidemic
  
  sim1<-reactive({
     beta<-input$beta
     IncubPeriod<-input$IncubPeriod
     RecovPeriod<-input$RecovPeriod
     FracAsym<-input$FracAsym
     HospRecPer<-input$HospRecPer
     CritProg<-input$CritProg
     CritRec<-input$CritRec
     K1<-input$K1
     K2<-input$K2
     Tmax<-input$Tmax
     InitInf=input$InitInf
     Region=input$Region
     
     inp1=list("InitInf"=InitInf,"IncubPeriod"=IncubPeriod,"RecovPeriod"=RecovPeriod,"FracAsym"=FracAsym,"HospRecPer"=HospRecPer,"CritProg"=CritProg,"CritRec"=CritRec,"Region"=Region,"K1"=K1,"K2"=K2,"beta"=beta,"Tmax"=Tmax)

     res1<-SimSEIRopt(inp1)
  })
  
  # Reactive function that computes the timecourse of a lockdown
  
  sim2<-reactive({
  beta<-input$beta
  IncubPeriod<-input$IncubPeriod
  RecovPeriod<-input$RecovPeriod
  FracAsym<-input$FracAsym
  HospRecPer<-input$HospRecPer
  CritProg<-input$CritProg
  CritRec<-input$CritRec
  K1<-input$K1
  K2<-input$K2
  Tmax<-input$Tmax
  InitInf=input$InitInf
  Region=input$Region
  Tstart=input$Tstart
  Tfinish=input$Tfinish

  inp2=list("Tstart"=Tstart,"Tfinish"=Tfinish,"InitInf"=InitInf,"IncubPeriod"=IncubPeriod,"RecovPeriod"=RecovPeriod,"FracAsym"=FracAsym,"HospRecPer"=HospRecPer,"CritProg"=CritProg,"CritRec"=CritRec,"Region"=Region,"K1"=K1,"K2"=K2,"beta"=beta,"Tmax"=Tmax)

  res2<-SimQuarOpt(inp2)
  })
  
  output$plot0=renderPlotly({
  
    res<-sim1()
    
    if(input$VarShow1=="S"){
      xvar=res$t
      yvar=res$S
      flag=0
      flag1=1
      ytext="Susceptibles (S)"
    }else if(input$VarShow1=="E"){
      xvar=res$t
      yvar=res$E
      flag=0
      flag1=0
      ytext="Exposed (E)"
    }else if(input$VarShow1=="A"){
      xvar=res$t
      yvar=res$A
      flag=0
      flag1=0
      ytext="Asymptomatic (A)"
    }else if(input$VarShow1=="IR"){
      xvar=res$t
      yvar=res$IR
      flag=0
      flag1=0
      ytext="Infected (mild)"
    }else if(input$VarShow1=="CC"){
      xvar=res$t
      yvar=res$CC
      flag=0
      flag1=0
      ytext="Critical cases"
    }else if(input$VarShow1=="R"){
      xvar=res$t
      yvar=res$RA+res$RR+res$RH+res$RC
      flag=0
      flag1=0
      ytext="Recovered"
    }else if(input$VarShow1=="D"){
      xvar=res$t
      yvar=res$DD
      flag=0
      flag1=0
      ytext="Deaths"
    }else if(input$VarShow1=="Inf"){
      xvar=res$t
      yvar=res$E+res$A+res$IR+res$IH+res$IC
      flag=0
      flag1=0
      ytext="Total infected"
    }else if(input$VarShow1=="Cases"){
      xvar=res$t
      yvar=res$IR+res$IH+res$IC
      flag=0
      flag1=0
      ytext="Symptomatic infected"
    }else if(input$VarShow1=="Hosp"){
      xvar=res$t
      yvar=res$HH+res$HC+res$CC
      flag=0
      flag1=0
      ytext="Total hospitalised"
    }else if(input$VarShow1=="AD"){
      xvar=seq(0,85,5)
      yvar=res$DAges
      flag=1
      flag1=0
      xtext="Age"
      ytext="Deaths"
    }

    data<- data.frame(xvar,yvar)
    
    if (flag==0){
      p=plot_ly(data,x=~xvar, y=~yvar, type = 'scatter', mode = 'lines')
      if (flag1==1){
        p=layout(p,xaxis=list(title="Time since introduction (days)",titlefont = list(size = 20),tickfont = list(size = 15),showline=T),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15)),title=input$Region) #,zeroline=T
      }else{
        p=layout(p,xaxis=list(title="Time since introduction (days)",titlefont = list(size = 20),tickfont = list(size = 15)),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15)),title=input$Region) #,zeroline=T
      }
    }else{
      p=plot_ly(data,x=~xvar, y=~yvar, type = 'bar')
      p=layout(p,xaxis=list(title=xtext,titlefont = list(size = 20),tickfont = list(size = 15)),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15),zeroline=T,showline=T),title=input$Region) # ,showline=T,zeroline=T 
    }
    
    p

  })
  
  #Plot timecourse with a lockdown
  
  output$plot1 = renderPlotly({
    
    res1<-sim1()

    res2<-sim2()

    if(input$VarShow2=="S"){
      xvar=res1$t
      yvar1=res1$S
      yvar2=res2$S
      flag=0
      flag1=1
      ytext="Susceptibles (S)"
    }else if(input$VarShow2=="E"){
      xvar=res1$t
      yvar1=res1$E
      yvar2=res2$E
      flag=0
      flag1=0
      ytext="Exposed (E)"
    }else if(input$VarShow2=="A"){
      xvar=res1$t
      yvar1=res1$A
      yvar2=res2$A
      flag=0
      flag1=0
      ytext="Asymptomatic (A)"
    }else if(input$VarShow2=="IR"){
      xvar=res1$t
      yvar1=res1$IR
      yvar2=res2$IR
      flag=0
      flag1=0
      ytext="Infected (mild)"
    }else if(input$VarShow2=="CC"){
      xvar=res1$t
      yvar1=res1$CC
      yvar2=res2$CC
      flag=0
      flag1=0
      ytext="Critical cases"
    }else if(input$VarShow2=="R"){
      xvar=res1$t
      yvar1=res1$RA+res1$RR+res1$RH+res1$RC
      yvar2=res2$RA+res2$RR+res2$RH+res2$RC
      flag=0
      flag1=0
      ytext="Recovered"
    }else if(input$VarShow2=="D"){
      xvar=res1$t
      yvar1=res1$DD
      yvar2=res2$DD
      flag=0
      flag1=0
      ytext="Deaths"
    }else if(input$VarShow2=="Inf"){
      xvar=res1$t
      yvar1=res1$E+res1$A+res1$IR+res1$IH+res1$IC
      yvar2=res2$E+res2$A+res2$IR+res2$IH+res2$IC
      flag=0
      flag1=0
      ytext="Total infected"
    }else if(input$VarShow2=="Cases"){
      xvar=res1$t
      yvar1=res1$IR+res1$IH+res1$IC
      yvar2=res2$IR+res2$IH+res2$IC
      flag=0
      flag1=0
      ytext="Symptomatic infected"
    }else if(input$VarShow2=="Hosp"){
      xvar=res1$t
      yvar1=res1$HH+res1$HC+res1$CC
      yvar2=res2$HH+res2$HC+res2$CC
      flag=0
      flag1=0
      ytext="Total hospitalised"
    }else if(input$VarShow2=="AD"){
      xvar=seq(0,85,5)
      yvar1=res2$DAges
      yvar2=res2$DAges
      flag=1
      flag1=0
      xtext="Age"
      ytext="Deaths"
    }
    
    data <- cbind.data.frame(xvar,yvar1,yvar2)
    
    if (flag==0){
      p=plot_ly(data,x=~xvar)
      p=p %>% add_trace(y = ~yvar1, name = 'Baseline', type = 'scatter', mode = 'lines') 
      p=p %>% add_trace(y = ~yvar2, name = 'Lockdown', type = 'scatter', mode = 'lines',colors=c("#a50f15","#fc9272"))
      if (flag1==1){
        p=layout(p,xaxis=list(title="Time since introduction (days)",titlefont = list(size = 20),tickfont = list(size = 15),showline=T),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15)),title=input$Region) #,zeroline=T
      }else{
        p=layout(p,xaxis=list(title="Time since introduction (days)",titlefont = list(size = 20),tickfont = list(size = 15)),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15)),title=input$Region) #,zeroline=T
      }
    }else{
      p=plot_ly(data,x=~xvar, y=~yvar1, type = 'bar')
      p=layout(p,xaxis=list(title=xtext,titlefont = list(size = 20),tickfont = list(size = 15)),yaxis=list(title=ytext,titlefont = list(size = 20),tickfont = list(size = 15),zeroline=T,showline=T),title=input$Region)
    }
    
    p
    
  })
  
  # ------------Set the sliders/forms that have dynamic values based on other sliders ----------------------
  
  observeEvent(input$Tstart,{
    if(7*input$Tstart>input$Tmax){
          updateSliderInput(session = session, inputId = "Tstart", value = floor(input$Tmax/7)-1)
        }
  })
  
  observeEvent(input$Tfinish,{
    if(input$Tfinish<=input$Tstart){
      updateSliderInput(session = session, inputId = "Tfinish", value = input$Tstart+1)
    }
    if(7*input$Tfinish>input$Tmax){
      updateSliderInput(session = session, inputId = "Tfinish", value = floor(input$Tmax/7))
    }
  })
  
  observeEvent(input$Tmax,{
    if(7*input$Tstart>input$Tmax){
      updateSliderInput(session = session, inputId = "Tstart", value = floor(input$Tmax/7)-1)
      updateSliderInput(session = session, inputId = "Tfinish", value = floor(input$Tmax/7))
    }
    if(7*input$Tfinish>input$Tmax){
      updateSliderInput(session = session, inputId = "Tfinish", value = floor(input$Tmax/7))
    }
  })
  
  # Reset all parameters if the RESET button is pushed
  observeEvent(input$reset,{
    updateSliderInput(session,'beta',value = 0.16)
    updateSliderInput(session,'IncubPeriod',value = 5)
    updateSliderInput(session,'RecovPeriod',value = 16)
    updateSliderInput(session,'FracAsym',value = 50)
    updateSliderInput(session,'HospRecPer',value = 8)
    updateSliderInput(session,'CritProg',value = 6)
    updateSliderInput(session,'CritRec',value = 10)
    updateSliderInput(session,'K1',value = 6)
    updateSliderInput(session,'K2',value = 6)
    updateSliderInput(session,'Tmax',value = 300)
    updateSliderInput(session,'Tstart',value = 5)
    updateSliderInput(session,'Tfinish',value = 13)
    updateSliderInput(session,'InitInf',value = 1)
  })
}