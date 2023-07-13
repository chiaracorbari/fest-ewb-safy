## Initialization and package loading
options(warn=-1)
library(dplyr)
library(lubridate)
library(data.table)
library(purrr)
library(raster)
library(tiff)
options(warn=0)
source("/home/eoafrica/shared/FINALE/newton.R")
source("/home/eoafrica/shared/FINALE/safy.R")
source("/home/eoafrica/shared/FINALE/waterbalance.R")
# r_source = robjects.r['source']

## input data

main_path = "/home/eoafrica/shared/FINALE/"
air_press = 98000 # Air Pressure [Pa]
water_dens = 1000 # Water Density [kg m-3]
evap_lat_heat = 2450 # Latent Heat of Evaporation [kJ kg-1]
mj_a_j = 1000 # Unit conversion coefficient
dt_safy = 86400 # SAFY works daily
dt_fest = 3600 # FEST works hourly

## meteo data input

filename=paste(main_path,"meteo_data.txt",sep="")
meteo=read.table(filename,header=T)
meteo=meteo%>%mutate(dataorastr=paste(date,hour))
meteo$dataora=as.POSIXct(meteo$dataorastr,format="%d/%m/%Y %H:%M",
                         tz="Africa/Casablanca")
meteo=meteo[,-(1:2)]
meteo=meteo[,-7]  # removal of excess columns
meteo$doy=yday(meteo$dataora) # DOY
meteo$cont=1 # counting column
# temporal continuity maintained
if (year(meteo$dataora[1]) < year(meteo$dataora[dim(meteo)[1]])) {
  # crop season spans two different years
  for (yy in 1:dim(meteo)[1]) {
    if (meteo$doy[yy]<150) {
      meteo$doy[yy]=meteo$doy[yy]+365
    }
  }
}

# average meteorological forcings for daily computations
meteo_sum = meteo %>%
  group_by(doy) %>%
  summarise(across(.cols = where(is.numeric),
                   .fns = ~sum(.x, na.rm = TRUE)), .groups="keep")
doymeteo=meteo_sum$doy
tair_daily = meteo_sum$temp/meteo_sum$cont # ?C
rad_daily = meteo_sum$rad/meteo_sum$cont # W m-2

irrig = meteo$irrig # mm
rain = meteo$prec # mm
water_inputs = (irrig+rain)/1000/dt_fest # m s-1
tair = meteo$temp # ?C
rad = meteo$rad # W m-2
wind_speed = meteo$wind # m s-1
if (meteo$rh[1]<1) {
  rh = meteo$rh*100 # %
} else {
  rh = meteo$rh # %
}

## Soil and Vegetation information
datapath=paste(main_path,sep="")

# input Fv

path_fv = paste(datapath,"FV/",sep="")
listfv = list.files(path=path_fv,pattern="*.tif")
datafv = map(listfv, function(x)
  as.matrix(readTIFF(paste(path_fv, x, sep=""))))

# input albedo

path_alb = paste(datapath,"ALBEDO/",sep="")
listalb = list.files(path=path_alb,pattern="*.tif")
dataalb = map(listalb, function(x)
  as.matrix(readTIFF(paste(path_alb, x, sep=""))))

# soil data

datapath = paste(main_path,"matrices/SOIL_PARAMETERS/",sep="")

soil_type = "loamy soil" # soil definition for capillary rise computation
# possible types: sandy s., loamy s., sandy clay s., silty clay s.
vadose=as.matrix(read.table(paste(datapath,"Groundwater_Level.txt",sep=""),header=F,skip=6))

# shallow layer (0-10 cm)
ks=as.matrix(read.table(paste(datapath,"Conductivity_10cm.txt",sep=""),header=F,skip=6))
# Sat. permeability [m s-1]
bci=as.matrix(read.table(paste(datapath,"Pore_Size_10cm.txt",sep=""),header=F,skip=6))
# Brooks-Corey Index [-]
sm_res=as.matrix(read.table(paste(datapath,"Residual_Water_Content_10cm.txt",sep=""),header=F,skip=6))
# Residual SM [m3 m-3]
sm_res=sm_res+0.1 # correction
sm_sat=as.matrix(read.table(paste(datapath,"Saturated_Water_Content_10cm.txt",sep=""),header=F,skip=6))
sm_sat=sm_sat+0.5
# Saturation SM [m3 m-3]
wp=as.matrix(read.table(paste(datapath,"Wilting_Point_10cm.txt",sep=""),header=F,skip=6))
# Wilting Point [m3 m-3]
fc=as.matrix(read.table(paste(datapath,"Field_Capacity_10cm.txt",sep=""),header=F,skip=6))
# Field Capacity [m3 m-3]
bub_pres=as.matrix(read.table(paste(datapath,"Bubbling_Pressure_10cm.txt",sep=""),header=F,skip=6))
# Bubbling pressure [m]
CN=as.matrix(read.table(paste(datapath,"Curve_Number.txt",sep=""),header=F,skip=6))
# Curve Number [-]
rsmin=as.matrix(read.table(paste(datapath,"Min_Stomatal_Resistance_Area_Morocco.txt",sep=""),header=F,skip=6))
rsmin=rsmin*2.5
# Minimum Stomatal Res. [s m-1]
rsoil=as.matrix(read.table(paste(datapath,"Min_Soil_Resistance_Area_Morocco.txt",sep=""),header=F,skip=6))
rsoil=rsoil+100
# Wet-soil Res. [s m-1]
soil_depth=as.matrix(read.table(paste(datapath,"Soil_Depth.txt",sep=""),header=F,skip=6))
# Soil depth [m]

# deeper layer (10-70 cm)
ks2=as.matrix(read.table(paste(datapath,"Conductivity_70cm.txt",sep=""),header=F,skip=6))
# Sat. permeability [m s-1]
bci2=as.matrix(read.table(paste(datapath,"Pore_Size_70cm.txt",sep=""),header=F,skip=6))
# Brooks-Corey Index [-]
sm_res2=as.matrix(read.table(paste(datapath,"Residual_Water_Content_70cm.txt",sep=""),header=F,skip=6))
# Residual SM [m3 m-3]
sm_sat2=as.matrix(read.table(paste(datapath,"Saturated_Water_Content_70cm.txt",sep=""),header=F,skip=6))
# Saturation SM [m3 m-3]
wp2=as.matrix(read.table(paste(datapath,"Wilting_Point_70cm.txt",sep=""),header=F,skip=6))
# Wilting Point [m3 m-3]
fc2=as.matrix(read.table(paste(datapath,"Field_Capacity_70cm.txt",sep=""),header=F,skip=6))
# Field Capacity [m3 m-3]
soil_depth2=as.matrix(read.table(paste(datapath,"Soil_Depth.txt",sep=""),header=F,skip=6))
# Soil depth [m]


# crop data (SAFY)
datapath = paste(main_path,"matrices//SAFY_PARAMETERS//",sep="")

# Growth parameters
Pgro_R2P=as.matrix(read.table(paste(datapath,"Pgro_R2P.txt",sep=""),header=F,skip=6))
# Global to PAR incident radiation ratio [-]
Pgro_Kex=as.matrix(read.table(paste(datapath,"Pgro_Kex.txt",sep=""),header=F,skip=6))
# Extinction of radiation in canopy (0.3-1)
Pgro_Lue=as.matrix(read.table(paste(datapath,"Pgro_Lue.txt",sep=""),header=F,skip=6))
# Effective light-use efficiency [g.MJ-1]
Pgro_MsZero=as.matrix(read.table(paste(datapath,"Pgro_Ms0.txt",sep=""),header=F,skip=6))
# Emergence Dry Mass Value (g/m2=100 x t/ha)
Pgro_Sla=as.matrix(read.table(paste(datapath,"Pgro_Sla.txt",sep=""),header=F,skip=6))
# Specific Leaf-Area 0.024 (m2 g-1)
Pgro_P2G=as.matrix(read.table(paste(datapath,"Pgro_P2G.txt",sep=""),header=F,skip=6))
# Partition coefficient To Grain

# Phenological parameters
Pfen_PrtA=as.matrix(read.table(paste(datapath,"Pfen_PrtA.txt",sep=""),header=F,skip=6))
# Partitioning to Leaves (after Maas, 1993) / Vary the origin slope of partition
Pfen_PrtB=as.matrix(read.table(paste(datapath,"Pfen_PrtB.txt",sep=""),header=F,skip=6))
# Vary the day of max LAI (partition=0)
Pfen_SenA=as.matrix(read.table(paste(datapath,"Pfen_SenA.txt",sep=""),header=F,skip=6))
# Temperature Threshold to Start Senescence (?C)
Pfen_SenB=as.matrix(read.table(paste(datapath,"Pfen_SenB.txt",sep=""),header=F,skip=6))
# Vary the rate of Senescence (?C)

# Temperature Effect On Development - result in TpS=TempStress[0-1]
Ptfn_Tmin =  5  # Minimum Temperature for Plant Development (?C)
Ptfn_Topt = 22  # Optimum Temperature for Plant Development (?C)
Ptfn_Tmax = 45  # Maximum Temperature for Plant Development (?C)
Ptfn_TpSn = 3   # Make vary the length of plateau around optimum T

p_crop = 0.55 # FAO depletion fraction (for tomatoes)
Pfen_MrgD = 336 # Day of Plant Emergence
fresh = 1 # dry-fresh biomass conversion (1 for maize, 0.055 for tomatoes)
date_root = 1551 # partitioning d-day between ET contributors

## Initial conditions (and variables initialization)

sizet = dim(meteo) # temporal steps
sizet = sizet[1]
sizem = dim(ks) # matrices size
matzero = matrix(0,nrow=sizem[1],ncol=sizem[2])
rain_flag = matzero # Rain flag
inter_time = matzero # non-rain (duration) time
S_ini = matzero # initial saturation
et = matzero # ET
    et_seasonal = et
sm = wp+0.11 # Initial soil moisture
sm2 = wp2+0.11 # Initial deep soil moisture
sm_sat = matzero+0.5
sm_sat2 = matzero+0.5
lai = matzero+0.1 # Initial Leaf Area Index
GLA = matzero
gross_rain = matzero
net_rain = matzero
S = matzero # soil Saturation
tempdegree = matzero # Degree Days (?C days)
SMT = matzero
drymass = matzero # Dry aerial mass [g m-2] (100 ton ha-1)
DAM = matzero
grainmass = matzero # Dry grain mass [g m-2] (100 ton ha-1)
DGM = matzero
partition = matzero # Partition-To-Leaf Index (-)
PRT = matzero
Tsurf = matzero
Rnetta = matzero
G = matzero
H = matzero
LEup = matzero
Tsurfup = matzero
smup = matzero
smup2 = matzero
percolaz = matzero

LEexp = matrix(0,nrow=sizet,ncol=1)
matexport = matrix(0,nrow=sizet,ncol=10)
#safy_data = matrix(0,nrow=sizet,ncol=)

x = 1 # Vegetation Fraction time level
xa = 1 # Albedo time level
xmax = length(datafv) # Total temporal levels available for vegetation
              
lai_data = matrix(NA,nrow=sizet,ncol=3)

xs = 0 # daily meteorological mean counter for SAFY
              
pix = c(20,35) # extraction pixel for the variables of interest

# Progress bar
pb = txtProgressBar(min = 0, max = sizet, style = 3, width = 50, char = "=")

## Start temporal cycle
for (i in 1:sizet) {
  
  # update Fv
  str = listfv[x]
  str1 = substr(str,17,26)
  datefv = as.Date(str1,format="%Y-%m-%d")
  doyfv = yday(datefv)
  if (doyfv<150) {
      doyfv = doyfv + 365
    }
  if (doyfv==meteo$doy[i]) {
      # acquire data only in not many Nans are present within the field
      map = datafv[[x]]
      #count = sum(mask*map/map)
      #count [is.na(count)] = 0
      #if (count==sum(mask)) { # all pixels in the field are storing valid data
          fveg_obs = datafv[[x]]
      #  for (ii in 1:sizem[1]) {
      #      for (jj in 1:sizem[2]) {
      #          lai_obs[ii,jj] = -log(1-min(1,max(0,fveg_obs[ii,jj])))/0.5
      #          lai_obs[ii,jj] = lai_obs[ii,jj]/8*4
      #          fveg_obs[ii,jj] = 1-exp(-0.5*lai_obs[ii,jj])
      #          }
      #      }
      lai_data[i,1] = fveg_obs[pix[1],pix[2]]
          lai_data[i,2] = -log(1-fveg_obs[pix[1],pix[2]])/0.5
      #} else { # otherwise, keep old data
      #    fveg_obs = fveg_obs
      #}
    if (x==xmax) {
      x = xmax
    } else {
      x = x+1
    }
  }
  
  # update albedo
  str = listalb[xa]
  str2 = substr(str,18,27)
  datealb = as.Date(str2,format="%Y-%m-%d")    
  doyalb = yday(datealb)
  if (doyalb<150) {
      doyalb = doyalb + 365
    }
  if (doyalb==meteo$doy[i]) {
      # acquire data only in not many Nans are present within the field
      map = dataalb[[x]]
      #count = sum(mask*map/map)
      #count [is.na(count)] = 0
      #if (count==sum(mask)) { # all pixels in the field are storing valid data
          alb_obs = dataalb[[xa]]
      if (alb_obs[pix[1],pix[2]]<0.10) {
          alb_obs = fveg_obs*0.18 + (1-fveg_obs)*0.23
          }
      lai_data[i,3] = alb_obs[pix[1],pix[2]]
      #} else { # otherwise, keep old data
      #    alb_obs = alb_obs
      #}
    if (xa==xmax) {
      xa = xmax
    } else {
      xa = xa+1
    }
  }
  
  # control of SAFY update
  update_safy=0
  if (i>1) {
    if (meteo$doy[i]>meteo$doy[i-1]) {
      # update meteorological variables for SAFY every midnight
      update_safy=1
    }
  }
   
  for (ii in 1:sizem[1]) {
    for (jj in 1:sizem[2]) {
    # ii=pix[1]
    # jj=pix[2]
      if (ks[ii,jj]>0) {
        
        # SAFY
        if (update_safy==1) {
          
          safy_output = safy(tempdegree[ii,jj],lai[ii,jj],drymass[ii,jj],grainmass[ii,jj],
                             partition[ii,jj],meteo$doy[i],
                             Pgro_Kex[ii,jj],Pgro_Lue[ii,jj],Pgro_MsZero[ii,jj],
                             Pgro_Sla[ii,jj],Pgro_P2G[ii,jj],Pfen_MrgD,
                             Pfen_PrtA[ii,jj],Pfen_PrtB[ii,jj],Pfen_SenA[ii,jj],Pfen_SenB[ii,jj],
                             Ptfn_Tmin,Ptfn_Topt,Ptfn_Tmax,Ptfn_TpSn,Pgro_R2P[ii,jj],
                             mean(meteo[meteo$doy==meteo$doy[i],3]), # mean daily temperature
                             mean(meteo[meteo$doy==meteo$doy[i],2]), # mean daily radiation
                             fc[ii,jj],wp[ii,jj],
                             fc2[ii,jj],wp2[ii,jj],soil_depth[ii,jj],p_crop,
                             sm[ii,jj],sm2[ii,jj],date_root,i)
          # tair_daily[xs],rad_daily[xs]
            #xs=xs+1
            #safy_data[i,1] = mean(meteo[meteo$doy==meteo$doy[i],3])
            #safy_data[i,2] = mean(meteo[meteo$doy==meteo$doy[i],2])
          
          SMT[ii,jj] = safy_output[[1]]
          DAM[ii,jj] = safy_output[[2]]
          DGM[ii,jj] = safy_output[[3]]
          GLA[ii,jj] = safy_output[[4]]
          PRT[ii,jj] = safy_output[[5]]
          
        } else {
          
          SMT[ii,jj] = tempdegree[ii,jj]
          DAM[ii,jj] = drymass[ii,jj]
          DGM[ii,jj] = grainmass[ii,jj]
          GLA[ii,jj] = lai[ii,jj]
          PRT[ii,jj] = partition[ii,jj]
          
        }
        
        Tsurf[ii,jj] = tair[i]
        
        # Energy balance
        newton_output = newton(sm[ii,jj],sm2[ii,jj],tair[i],rad[i],soil_depth[ii,jj],wind_speed[i],air_press,rh[i],
                               rsmin[ii,jj],rsoil[ii,jj],sm_res[ii,jj],sm_sat[ii,jj],bci[ii,jj],
                               wp[ii,jj],bub_pres[ii,jj],fc[ii,jj],
                               min(1,max(0,alb_obs[ii,jj])),
                               min(max(0,fveg_obs[ii,jj]),1),
                               lai[ii,jj],Tsurf[ii,jj],date_root,i)
        
        Rnetta[ii,jj] = newton_output[[1]]
        G[ii,jj] = newton_output[[2]]
        H[ii,jj] = newton_output[[3]]
        LEup[ii,jj] = newton_output[[4]]
        Tsurfup[ii,jj] = newton_output[[5]]
        
        et[ii,jj] = LEup[ii,jj]/(water_dens*evap_lat_heat*mj_a_j) # Evapotranspiration
        
        # Water balance
        watbal_output = waterbalance(sm[ii,jj],sm2[ii,jj],et[ii,jj],
                                     ks[ii,jj],ks2[ii,jj],
                                     sm_res[ii,jj],sm_res2[ii,jj],
                                     sm_sat[ii,jj],sm_sat2[ii,jj],
                                     bci[ii,jj],bci2[ii,jj],
                                     bub_pres[ii,jj],fc[ii,jj],wp[ii,jj],
                                     dt_fest,soil_depth[ii,jj],
                                     soil_depth2[ii,jj],CN[ii,jj],
                                     water_inputs[i],
                                     gross_rain[ii,jj],net_rain[ii,jj],
                                     rain_flag[ii,jj],inter_time[ii,jj],S[ii,jj],
                                     soil_type,date_root,i,vadose[ii,jj])
        
        smup[ii,jj] = watbal_output[[1]]
        smup2[ii,jj] = watbal_output[[2]]
        gross_rain[ii,jj] = watbal_output[[3]]
        net_rain[ii,jj] = watbal_output[[4]]
        rain_flag[ii,jj] = watbal_output[[5]]
        inter_time[ii,jj] = watbal_output[[6]]
        S[ii,jj] = watbal_output[[7]]
        percolaz[ii,jj] = watbal_output[[9]]
        
     } else {
        
        SMT[ii,jj] = -9999
        DAM[ii,jj] = -9999
        DGM[ii,jj] = -9999
        GLA[ii,jj] = -9999
        PRT[ii,jj] = -9999
        
        Tsurf[ii,jj] = -9999
        
        Rnetta[ii,jj] = -9999
        G[ii,jj] = -9999
        H[ii,jj] = -9999
        LEup[ii,jj] = -9999
        Tsurfup[ii,jj] = -9999
        
        et[ii,jj] = -9999
        
        smup[ii,jj] = -9999
        smup2[ii,jj] = -9999
        gross_rain[ii,jj] = -9999
        net_rain[ii,jj] = -9999
        rain_flag[ii,jj] = -9999
        inter_time[ii,jj] = -9999
        S[ii,jj] = -9999
        percolaz[ii,jj] = -9999
        
      }
    }
  }
  
  et_seasonal = et_seasonal + et*dt_fest*1000 # from m s-1 to mm h-1
  # parameters update
  LE     = LEup
  sm     = smup
  sm2    = smup2
  Tsurf  = Tsurfup
  tempdegree = SMT
  partition = PRT
  drymass = DAM
  grainmass = DGM
  lai = GLA
  yield = grainmass/100/fresh
  
  # export
  pix = c(20,35) # extraction pixel for the variables of interest
  matexport[i,1]=LE[pix[1],pix[2]]
  matexport[i,2]=H[pix[1],pix[2]]
  matexport[i,3]=Rnetta[pix[1],pix[2]]
  matexport[i,4]=et[pix[1],pix[2]]
  matexport[i,5]=Tsurf[pix[1],pix[2]]
  matexport[i,6]=sm[pix[1],pix[2]]
  matexport[i,7]=sm2[pix[1],pix[2]]
  matexport[i,8]=lai[pix[1],pix[2]]
  matexport[i,9]=grainmass[pix[1],pix[2]]
  matexport[i,10]=G[pix[1],pix[2]]
    
  # export RET and LAI matrices
  timelst=as.POSIXct("2021-11-22 11:00",format="%Y-%m-%d %H:%M",tz="Africa/Casablanca")
  if (meteo$dataora[i]==timelst) {
      write.table(Tsurf,file="RET_2021-11-22.txt",row.names=FALSE,col.names=FALSE)
    }
  timelai=as.POSIXct("2022-02-15 11:00",format="%Y-%m-%d %H:%M",tz="Africa/Casablanca")
  if (meteo$dataora[i]==timelai) {
      write.table(lai,file="LAI_2022-02-15.txt",row.names=FALSE,col.names=FALSE)
    }
    
  # progress bar update
  setTxtProgressBar(pb,i)
  
}

close(pb)
              
write.table(et_seasonal,file="ET.txt",row.names=FALSE,col.names=FALSE)

varlabels=c("Date","Latent_Heat","Sensible_Heat","Net_Radiation",
            "ET","RET","SM","SM2","LAI","Grainmass","SHF")
              
write.table(data.frame(dates=meteo$dataora,matexport),"/home/eoafrica/shared/FINALE/fest_safy_results.csv",sep=",",
            row.names=FALSE,col.names=varlabels)

#write.table(data.frame(dates=meteo$dataora,lai_data),"/home/eoafrica/shared/FINALE/lai_overpasses.csv",sep=",",
 #           row.names=FALSE,col.names=c("Date","LAI (S2)"))

write.table(data.frame(dates=meteo$dataora,lai_data),"/home/eoafrica/shared/FINALE/veg_data.csv",sep=",",
            row.names=FALSE,col.names=c("Date","fv","lai","alb"))
              
#write.table(data.frame(safy_data),"/home/eoafrica/shared/FINALE/safy_meteo.csv",sep="\t",
 #           row.names=meteo$dataora,col.names=FALSE)

#write.table(fveg_obs,file="fv.txt",row.names=FALSE,col.names=FALSE)

#write.table(alb_obs,file="alb.txt",row.names=FALSE,col.names=FALSE)
            