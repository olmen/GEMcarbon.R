##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  This file is updated from EGM_fluxfunction_20171205.R which was designed by Cecile Girardin  ++
##  Annotation added by Simone                                                                   ++
##  updated by Huanyuan Zhang 2021 Jan                                                           ++
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##  Codes to estimate soil efflux in umol m-2 sec-1 based on raw EGM-4 data,  ++
##  see GEM protocol for more information                                     ++
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



##----------------------
##  Step 1 Preparation  
##----------------------
rm(list = ls())
# load packages
library (tidyverse)
setwd("F:/Oxford/Chapter_two/soil_respiration/output") #where output file will be

# read in soil respiration auxillary functions from GitHub
source("F:/Oxford/Chapter_two/GEMcarbon.R/soilrespiration_auxfunctions.R")

# Parameters
#  set a DEFAULT collar diameter
collardiameter = 12 # 12 cm in Africa & Andes & TAM, 10.6 cm in Malaysia, 10.143 in JEN 


##--------------------
##  Step 2 Load data  
##--------------------


# If we do have air temperature and collar height, merge the two datasets.

rtot  = read.table("F:/Oxford/Chapter_two/soil_respiration/input/tot_soil_resp_20190116_KOG.csv", sep=",", header=T)
#rpart = read.table("part_soil_resp_20180226.csv", sep=",", header=T)
temp_vwc_ch = read.table("F:/Oxford/Chapter_two/soil_respiration/input/weather_averages_20190130_KOG.csv", sep=",", header=T)

rtot = mutate(rtot, uid = paste(plot_code, sub_plot, year,  sep = '_'))
temp_vwc_ch = mutate(temp_vwc_ch, uid = paste(plot_code, sub_plot, year,  sep = '_'))
datafile = left_join(rtot, temp_vwc_ch[,c('uid','air_temp_c','collar_height_cm')], by = "uid")
#Check key field one by one and make sure they are good


Hmisc::describe(datafile$plot_code)
Hmisc::describe(datafile$sub_plot)

Hmisc::describe(datafile$collar_number)

Hmisc::describe(datafile$replica)

Hmisc::describe(datafile$day)
Hmisc::describe(datafile$month)

Hmisc::describe(datafile$year)

Hmisc::describe(datafile$measurement_code)

Hmisc::describe(datafile$treatment_code_partitioning)

Hmisc::describe(datafile$collar_height_cm)

Hmisc::describe(datafile$air_temp_c)

Hmisc::describe(datafile$co2ref_ppm_sec)

##---------------------------
##  Step 3 Flux calculation  
##---------------------------

# Replace missing temperature and collar height with default values or the median of the whole dataset.
w = which(is.na(datafile$air_temp_c))
Temperature_fill<-ifelse(is.na(median(datafile$air_temp_c)),  25,  median(datafile$air_temp_c,na.rm = T))
datafile$air_temp_c[w] = Temperature_fill

w = which(is.na(datafile$collar_height_cm))
ch_fill<-ifelse(is.na(median(datafile$collar_height_cm)),  5.09,  median(datafile$collar_height_cm,na.rm = T))
datafile$collar_height_cm[w] = ch_fill

## Corrections and conversions
# add a temperature correction from Sotta et al. 2004 Q10=1.8 and k=0.0613
corrsresA = exp(-0.0695*(1))
# Convert units umol m2 s-1 to MgC ha month = 1mo=2592000sec, 10000m2=1ha, 1000000umol = 1 mol, 1mol = 12 g, 1000000g=1Mg
convert = (2592000*10000*12)/(1000000*1000000)
# !! It would be better to Convert to umol m2 s-1 here, from g CO2 m-2 h-1.

data<-datafile
# Set a unique id for each measurement
# Very important that the name of your columns is exactly the same as the GEM example
# to be able to run these functions
data$codew = paste(data$plot_code, data$sub_plot, data$collar_number, data$replica, data$day, data$month, data$year, data$measurement_code, data$treatment_code_partitioning, sep=".") # In most cases you will need to add the data$replica to this unique id.

# Sort by code and time (time)
data = data [order (data $ codew, data $ time),]

# get unique identifier for each measurement


uid = unique(data$codew)
xx  = c()
yy  = c()
zz  = c()
zzz = c()
ddd = c()

for (i in 1:length(uid)) {
  sub=data[data$codew == uid [i],]       # Drop first four measurements: & data$time >= 20 #It was sub=subset(data, subset = (data$codew == uid [i])), but I don't know why this does not work for ANK

  id       = tail(sub$codew, n=1) 
      
  
  P        = tail(sub$atmp_mb, n=1)*100                                           # ambient pressure at t10 (Pa)
  Ta       = tail(sub$air_temp_c, n=1)                                            # air temp at t10 (deg C)
  ch       = tail(sub$collar_height_cm, n=1)                                               # see gap filling function fill.na() in soilrespiration_auxfinctions.r
  codep    = tail(sub$treatment_code_partitioning, n=1)
  plot     = tail(sub$plot_code, n=1)
  A_collar = (pi*((collardiameter/2)/100)^2)                                      # Area of the collar (m2)
  Vd       = 0.0012287                                                            # m3 (constant)
  A        = 0.00950                                                              # m2 (constant)
  Ru       = 8.3144                                                               # J mol-1 K-1 (constant)
  Va       = A*(ch/100)                                                           # additional volume m3
  Vtot     = Vd+Va
  
  # This is the equation from the GEM manual, we have replaced it with the equation below.
  #C10      = tail(ten_co2, n=1)                                                   # last CO2 measurement of last 10 measurements
  #C1       = head(ten_co2, n=1)                                                   # first CO2 measurement of last 10 measurements
  #t10      = tail(ten_time, n=1)                                                  # last time step of 10 last measurements
  #t1       = head(ten_time, n=1)                                                  # first time step of 10 last measurements
  #fl       = ((C10 - C1)/(t10 - t1)) * (P/(Ta + 273.15))*(Vd/A)*((44.01*0.36)/Ru) # CO2 efflux (g CO2 m-2 h-1). This is the equation we have in the GEM manual. We need to update the manal.
  #flux     = (fl*A/Vd*(Va+Vd)/A)*6.312                                           # Convert to umol m-2 s-1, and correct for collar height.
  
  sub <-sub[!is.na(sub$co2ref_ppm_sec), ]
  sub <- sub[!is.na(sub$time), ]
  ten_co2  = tail(sub$co2ref_ppm_sec, n=10)                                                   
  ten_time = tail(sub$time, n=10) 
  
  
  if (sum(is.na(ten_co2)) < 7 & length(ten_co2) >= 7){                   # if at least 7 of the 10 values are different from NA, we apply a linear a regression to get CO2 flux


    fit      = lm(ten_co2~ten_time)
    Co2slope = fit$coefficients[2]                                       # ["ten_time"]
    
    flux    = Co2slope * P * Vtot / (Ru * (Ta + 273.15)) / A_collar                 # output is in umol m-2 s-1. This equation was provided by Terhi Riutta, January 2018.
    NA_note = 'NA'
    if (!is.numeric(flux)) {
      warning(paste0('flux calculation fail in ','id'))
      NA_note = 'flux calculation fail'
    }
    
  }else{
    flux = NA
    NA_note = 'not enough EGM values for fitting curve'
  }
  
  xx        = rbind(xx, id)
  yy        = rbind(yy, flux)
  zz        = rbind(zz, as.character(codep))
  zzz       = rbind(zzz, as.character(plot))
  ddd      =  rbind(ddd, as.character(NA_note))
}

rownames (xx) = NULL
rownames (yy) = NULL
rownames (zz) = NULL
rownames (zzz) = NULL
rownames (ddd) = NULL
Res = data.frame (cbind (xx, yy, zz, zzz,ddd)) #zz, generating a new data set
# from the calculations above
# then you change the name of these new columns
#look at your data set and see what's changing
# note that you now have a new set called Res
colnames (Res) = c ("codew", "Rflux_umolm2sec", "part_code", "plot_code",'NA_note')
# Make sure the flux is numeric
Res $ fluxnum = as.numeric ( (Res $ Rflux_umolm2sec)) 

# Filter outlyers
## Deleting negative values and values above 15 umol m^-2 s^-1 because they don't make sense for soil respiration
Res$NA_note[Res$fluxnum <= 0 | Res$fluxnum >= 15] <- 'flux value exceed limits, should remove'
Res$fluxnum[Res$fluxnum <= 0 | Res$fluxnum >= 15] <- NA


## Here you are just organizing your data with the information that
## need, year, month, plot, converting the efflux into carbon
# split the code into the information we need (plot_code, sub_plot, collar_number, replica, day, month, year).
Res $ codew = as.character (Res $ codew)
Res $ Rflux_MgC_ha_mo = Res $ fluxnum * convert * corrsresA # see the conversion and correction factors above.

temp = (strsplit (Res $ codew, "[.]"))
Res $ plot_code = unlist (lapply (temp, `[[`, 1))
Res $ sub_plot = unlist (lapply (temp, `[[`, 2))
Res $ collar_number = unlist (lapply (temp, `[[`, 3))
Res $ replica = unlist (lapply (temp, `[[`, 4))
Res $ day = unlist (lapply (temp, `[[`, 5))
Res $ month = unlist (lapply (temp, `[[`, 6))
Res $ year = unlist (lapply (temp, `[[`, 7))
Res$collection    = as.Date(paste(Res$year, Res$month, Res$day, sep="."), format="%Y.%m.%d") 
Res $ measurement_code = unlist (lapply (temp, `[[`, 8))
Res $ treatment_code_partitioning = unlist (lapply (temp, `[[`, 9))

# convert character to numeric variables !! CHECK THIS: NAs introduced by coercion

Res[,c('Rflux_MgC_ha_mo','fluxnum')] <- sapply(Res[,c('Rflux_MgC_ha_mo','fluxnum')],as.numeric)
Res<-Res[,!(names(Res) %in% "Rflux_umolm2sec")]
Res<-rename(Res,Rflux_umolm2sec=!!as.name("fluxnum"))
##---------------------------
##  Step 4 Calculate mean  
##---------------------------


# Res - average per collar (average replicas!)
## here you generate a data set with the average of
## efflux already in MgC per ha and month for each pipe
## 
Res1 = Res %>% group_by(plot_code, measurement_code, 
                        treatment_code_partitioning, #for partitioning
                        #disturbance_code_control, #for distubance
                        #litter_code, #for IC
                        #cwd_num, #for cwd
                        year, month, sub_plot, collar_number) %>% 
  dplyr::summarise(collectiondate = min(collection),
                   N = length(Rflux_MgC_ha_mo[!is.na(Rflux_MgC_ha_mo)]),
                   collar_fluxnum_umolm2sec = mean(Rflux_umolm2sec, na.rm = T),
                   collar_Rflux_MgC_ha_mo = mean(Rflux_MgC_ha_mo, na.rm = T)) %>% 
  arrange(plot_code, year, month, sub_plot, collar_number) %>% data.frame(.)

# Average over the whole plot: one value per plot per month
Res2 = Res1 %>% group_by(plot_code, measurement_code, treatment_code_partitioning, year, month) %>% 
  dplyr::summarize(avg = mean(collar_Rflux_MgC_ha_mo, na.rm = T), 
                   sd = sd(collar_Rflux_MgC_ha_mo, na.rm = T),
                   collection_date = max(collectiondate)) %>% data.frame(.)


Res2 %>% group_by(plot_code) %>%
  ggplot(data=., aes(month, avg, colour=year)) + geom_point() +
  facet_wrap(~plot_code)

write.csv (Res1, file = paste0(substr(plot,1,3),"_rsoil_TOTAL_average_per_collar.csv"))
write.csv (Res2, file =paste0(substr(plot,1,3), "_rsoil_TOTAL_one_value_plot_month.csv"))