library(shiny)
library(shinyWidgets)
library(plotly)
library(reshape)
library(plotly)

data_regions<-read.csv("code/data_new.csv", header = TRUE)

fluidPage(
  titlePanel("Interactive Online Tool for modelling COVID-19 spread and containment"),
  hr(),
  p(div(HTML("<strong>Disclaimer</strong>: This simulation tool is for <strong>research and educational purposes only</strong>. It is not intended to be used for decision-making due to many uncertainties concerning the details of COVID-19 infection and transmission, as well as model limitations. This work is licensed under a <a href=https://creativecommons.org/licenses/by-sa/4.0/> Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License </a>"))),
  
  
  sidebarLayout(
    
    sidebarPanel(
      
      fluidRow(
        chooseSliderSkin(skin="Square",color="Red"),
        
        selectizeInput('Region', 'Select a region', choices = c("choose" = "", levels(data_regions$Name)),selected="Aberdeen City"),
        
        column(width=6,
               
               h4(div(HTML("<em>Choose clinical parameters:</em>"))),
               setSliderColor(c("DeepSkyBlue","#FF4500", "Teal","Red","Green","Magenta","Blue","Purple","DarkGreen","DarkBlue","Crimson","OrangeRed"), c(1, 2, 3,4,5,6,7,8,9,10,11,12)),
               sliderInput("beta", "Transmission rate",0.02, 1, 0.16, step=0.02, post="/day"),
               sliderInput("IncubPeriod", "Duration of incubation period", 1, 20, 5, step=0.5, post = " days"),
               sliderInput("RecovPeriod", "Duration of recovery period", 1, 20, 16, step=0.5, post = " days"),
               sliderInput("FracAsym", "% of asymptomatic infections", 1, 100, 50, step=1, pre="%"),
               sliderInput("HospRecPer", "Duration of hospital recovery for non-critical cases", 1, 10, 8, step=0.5, post = " days"),
               sliderInput("CritProg", "Time to progress to critical care in hospital", 1, 10, 6, step=0.5, post = " days"),
               sliderInput("CritRec", "Time to recover from critical care", 1, 30, 10, step=0.5, post = " days"),
               sliderInput("K1", "Number of incubation stages", 1, 8, 6, step=1),
               sliderInput("K2", "Number of recovery stages", 1, 8, 6, step=1),
               br(),
               
               ),
        column(width=6,
               
               h4(div(HTML("<em>Choose simulation parameters:</em>"))),
               
               numericInput("InitInf","Initial # infected:",value = 100, min = 1, step = 1),
               sliderInput("Tmax", div(HTML("Maximum time")),1, 1000, 300, step=10, post=" days"),
               sliderInput("Tstart", div(HTML("Start of lockdown")),1, 15, 8, step=1, post=" weeks"),
               sliderInput("Tfinish", div(HTML("Finish of lockdown")),1, 30, 16, step=1, post=" weeks"),
               
               br(),
               actionButton("reset", "Reset all"),    
               
               hr(),

        )
      ),
 
    ),
    
    mainPanel(
      
      navbarPage("Output:",
                 
                 tabPanel("Spread",
                          fluidPage(
                            fluidRow(
                              
                              h3("Baseline COVID-19 dynamics"),
                              p(HTML("This is a simulation of a COVID-19 epidemic in a single region <strong>without any interventions</strong>. It is being re-run in real time for each combination of parameters (determined by values on the sliders). For baseline values of parameters, that involves solving a system of 702 differential equqtions, which may take a moment.")),
                              
                              br(),
                              br(),
                              selectInput("VarShow1",
                                          label = "Select variable to show:",
                                          choices = c("Suceptible (S)"="S", "Exposed (E)"="E", "Asymptomatic Infections (A)"="A", "Mild Infections (I_R)"="IR", "Critical Care Cases (CC)"="CC", "Recovered (RA+RR+RH+RC)"="R", "Deaths (D)"="D", "All infected (E + A + all I)"="Inf","All symptomatic (IR+IH+IC)"="Cases","All hospitalized (HH+HC+CC)"="Hosp","Age of deaths"="AD"), 
                                          selected = c("Cases")
                              ),
                              plotlyOutput("plot0"),
                              br(),
                               
                              p(HTML("<b>User instructions:</b> The graph shows the expected numbers of individuals over time who are susceptible, exposed, infected, recovered, hospitalised, in critical care, or dead over time, as well as age distribution of deaths. Infected individuals first pass through an exposed/incubation phase where they are asymptomatic and not infectious, and then move into an infectious stage classified by the clinical status of infection (asymptomatic, mild, severe, or critical). A more detailed description of the model is provided in the Model Description tab. The region, initial condition, and parameter values used to simulate the spread of infection can be specified/adjusted using sliders located in the left-hand panel. Default slider values are equal to estimates taken from the literature (see Model tab). To reset default values, click on the <em>Reset all</em> button. The plot is interactive: hover over it to get values, dragging over a range allows zooming."))
                            )
                          )
                 ),
                 
                 tabPanel("Lockdown",
                          fluidPage(
                            fluidRow(
                              h3("Effect of lockdown on COVID-19 dynamics"),
                              p(HTML("This is a simulation of a COVID-19 epidemic in a single region <strong>with a temporary lockdown</strong>. It is being re-run in real time for each combination of parameters (determined by values on the sliders). For baseline values of parameters, that involves solving a system of 702 differential equqtions, which may take a moment.")),
                              
                              br(),
                              
                              br(),
                              selectInput("VarShow2",
                                          label = "Select variable to show:",
                                          choices = c("Suceptible (S)"="S", "Exposed (E)"="E", "Asymptomatic Infections (A)"="A", "Mild Infections (I_R)"="IR", "Critical Care Cases (CC)"="CC", "Recovered (RA+RR+RH+RC)"="R", "Deaths (D)"="D", "All infected (E + A + all I)"="Inf","All symptomatic (IR+IH+IC)"="Cases","All hospitalized (HH+HC+CC)"="Hosp","Age of deaths"="AD"), 
                                          selected = c("Cases")
                              ),
                              plotlyOutput("plot1"),
                              br(),
                              
                              p(HTML("<b>User instructions:</b> The graph shows the expected numbers of individuals over time who are susceptible, exposed, infected, recovered, hospitalised, in critical care, or dead over time, as well as age distribution of deaths. Infected individuals first pass through an exposed/incubation phase where they are asymptomatic and not infectious, and then move into an infectious stage classified by the clinical status of infection (asymptomatic, mild, severe, or critical). A more detailed description of the model is provided in the Model Description tab. The region, initial condition, and parameter values used to simulate the spread of infection can be specified/adjusted using sliders located in the left-hand panel. Default slider values are equal to estimates taken from the literature (see Model tab). To reset default values, click on the <em>Reset all</em> button. The plot is interactive: hover over it to get values, dragging over a range allows zooming."))
                            )
                          )
                 ),
                 
  
                 tabPanel("Model", br(),
                          fluidRow(column(12,
                                          withMathJax(),
                                          h2("Model Description"),
                                          includeMarkdown("SEIR.Rmd"),
                           ))),
                 
                 tabPanel("About",
                          fluidPage(
                            br(),
                            includeMarkdown("About.Rmd")
                          ))
                 
      ),
      width=7
    )
    
  )
  
)

