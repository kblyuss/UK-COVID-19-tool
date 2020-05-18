# Interactive COVID-19 Simulation Tool

The files in this repository produce the R Shiny application https://kblyuss.shinyapps.io/icst that models the spread of COVID-19 in the UK, and its containment using a lockdown. The app allows one to choose any UK region, the timing of introduction and lifting of lockdown, as well as modify various parameters characterising clinical features of disease transmission.

All functions for simulations are contained in the **code/functions.R** The **code** folder also contains all necessary data files, including the file with an age distribution for all UK regions, as well as all mixing matrices.

To explore the dynamics of disease and lockdown without the need to run the app, you can run the files

+ **runSpread.R** for baseline disease dynamics
+ **runQuar.R** for the dynamics with lockdown

These files can be run independently of the app (described in **server.R**, **ui.R**, and **global.R**), and they rely on the same functions and data files that are contained in the **code** folder.
