rm(list=ls())
library(deSolve)
library(plotly)
library(compiler)
library(Rcpp)
source("code/functions.R")

# Set parameters
# (in Shiny app these are input by sliders)

IncubPeriod= 5 #Duration of incubation period
RecovPeriod= 16.6 #Duration of recovery period
FracAsym=50 #Fraction of infections that are asymptomatic
HospRecPer=8 #Duration of hospital recovery for non-critical cases, days
CritProg=6 #Time to progress to critical care in hospital, days
CritRec=10 #Time to recover from critical care, days
Region='West Sussex' #Region being studied
K1=6 # Number of incubation stages
K2=6 # Number of recovery stages
beta=0.1506 # Transmission rate

Tmax= 300 #Initial # infected
InitInf= 100 #Maximum time

data_regions<-read.csv("code/data_new.csv", header = TRUE)
BBC_matrix<-read.table("code/BBC_matrix_all.dat")
ph_BBC<-read.table("code/prop_hosp.dat")
pi_BBC<-read.table("code/prop_icu.dat")
pm_BBC<-read.table("code/mort.dat")

#Put these into an input structure
input=list("InitInf"=InitInf,"IncubPeriod"=IncubPeriod,"RecovPeriod"=RecovPeriod,"FracAsym"=FracAsym,"HospRecPer"=HospRecPer,"CritProg"=CritProg,"CritRec"=CritRec,"Region"=Region,"K1"=K1,"K2"=K2,"beta"=beta,"Tmax"=Tmax)

# Auxiliary C function to compute matrix product

mm_code =
  "NumericVector my_mm(NumericMatrix m, NumericVector v){
   int nRow = m.rows();
   int nCol = m.cols();
   NumericVector ans(nRow);
   double v_j;
   for(int j = 0; j < nCol; j++){
     v_j = v[j];
     for(int i = 0; i < nRow; i++){
       ans[i] += m(i,j) * v_j;
     }
   }
   return(ans);
 }
 "
# Compiling
my_mm = cppFunction(code = mm_code)

# Run simulations

sim=SimSEIRopt(input)

xvar=sim$t # time variable foe a solution
yvar=sim$E+sim$A+sim$IR+sim$IH+sim$IC # total number of infected
data<- data.frame(xvar,yvar)

p=plot_ly(data,x=~xvar, y=~yvar, type = 'scatter', mode = 'lines')
p=layout(p,xaxis=list(title="Time since introduction (days)",titlefont = list(size = 20),tickfont = list(size = 15)),yaxis=list(title="Total infected",titlefont = list(size = 20),tickfont = list(size = 15)),title=input$Region)

p
