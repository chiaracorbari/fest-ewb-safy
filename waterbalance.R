waterbalance = function(sm,sm2,et,ks,ks2,sm_res,sm_res2,sm_sat,sm_sat2,
                        Bci,Bci2,bubpress,fc,wp,dt,soil_depth,soil_depth2,CN2,
                        rain,gross_rain,net_rain,rain_flag,
                        inter_time,S,soil,date_root,i,vadose) {
  
  # previous variables initialization
  old_net_rain = net_rain
  old_gross_rain = gross_rain
  old_sm = sm
    old_sm2 = sm2
  saturation = ((old_sm-sm_res)/(sm_sat-sm_res))
    saturation2 = ((old_sm2-sm_res2)/(sm_sat2-sm_res2))
  
  interirriga = 2*dt # ????????????????
  if ((rain_flag==0)&&(rain*1000*3600>0.5)){ # it's starting to rain
    rain_flag = 1
    inter_time = 0
    
    # Curve Number
    CN1=CN2-(20*(100-CN2))/(100-CN2+exp(2.533-0.0636*(100-CN2)))
    if (CN1<0.4*CN2) {
      CN1=0.4*CN2
    }
    CN3=CN2*exp(0.00673*(100-CN2))
    
    # Saturation levels [mm]
    S1 = 254*(100/CN1-1)
    S2 = 254*(100/CN2-1)
    S3 = 254*(100/CN3-1)
    # Weights
    W2 = 2*(log(0.5/(1-(S2/S1))-5)-log(1/(1-(S3/S1))-1))
    W1 = log(1/(1-(S3/S1))-1)+W2
    
    #	!S_cn=S1*(1-(saturation/(saturation+exp(W1-W2*saturation))))
    S = S1-saturation*(S1-S3) # global saturation [mm]
    
  } else if ((rain_flag==1)&&(rain*1000*3600>0.5)) { # it has already started to rain
    rain_flag = 2 # it continues to rain
    
  } else if (rain*1000*3600<0.5) { # it's raining
    inter_time = inter_time+dt
    runoff = 0
    infiltration = rain
  }
  
  if (inter_time>interirriga) {
    rain_flag=0
    gross_rain=0
  }
  if (rain_flag==0) {
    gross_rain=0
    runoff=0
    infiltration=0
  }
  
  if ((rain_flag==1||rain_flag==2)) { # && rain>0
  # if (rain_flag==1 || rain_flag==2) && rain*1000*3600>0.5
    
    # there is rain, soil saturation needs to be updated
    if (S==0) {
      gross_rain=old_gross_rain+rain*1000*dt
      runoff=rain
      infiltration=0
    } else {
      old_gross_rain=gross_rain
      
      if (old_gross_rain >= 0.2*S) { # more than 20% saturation, some runoff (SCS-CN)
        old_net_rain = ((old_gross_rain-0.2*S)^2)/(old_gross_rain+(1-0.2)*S)
      } else { # less than 20% saturation, no net rain
        old_net_rain = 0
      }
      if (rain==0) {
        old_gross_rain = 0
        old_net_rain = 0
      }
      gross_rain = old_gross_rain+rain*1000*dt
      if (gross_rain>=0.2*S) {
        net_rain = ((gross_rain-0.2*S)^2)/(gross_rain+(1-0.2)*S)
        runoff = (net_rain-old_net_rain)/1000/dt
        infiltration = rain-runoff
        if (runoff>rain) {
          runoff = rain
          infiltration = 0
        }
        if (runoff < 0) {
          runoff = 0
          infiltration = rain
        }
      } else {
        runoff = 0
        infiltration = rain
      }
    }
  } else {
    runoff = 0
    infiltration = rain
  }
  
  # capillary
  x_inf=10^(2+0.3*(fc*100-10)/(30-10))/100
  if (infiltration>0) { # no capillary
    ris=0
  } else {
    if (old_sm>=fc) {
      ris=0
    } else {
      
      ks_mm_day=ks*1000*86400; # from m s-1
      
      if (ks>5*10^(-5)) {
        a_cr=-0.3112-ks_mm_day*10^(-5)
        b_cr=-3.4936 + 0.2416*log(ks_mm_day)
      }
      
      if (ks<5*10^(-5) & ks>5*10^(-6)) {
        a_cr=-0.4986 + 9*ks_mm_day*10^(-5)
        b_cr=-4.8320 + 0.4778*log(ks_mm_day)
      }
      
      if (ks<5*10^(-6) & ks>7*10^(-7)) {
        a_cr=-0.5677 - 4*ks_mm_day*10^(-5)
        b_cr=-6.51189 + 0.5922*log(ks_mm_day)
      }
      
      if (ks<7*10^(-7)) {
        a_cr=-0.6366 + 8*ks_mm_day*10^(-4)
        b_cr=-1.9165 + 0.7063*log(ks_mm_day)
      }
      ris_max=exp((log(vadose-soil_depth)-b_cr)/a_cr)
      x=16
      fcr=1-(((old_sm-wp)/(fc-wp))^x)
      ris=ris_max*fcr/1000/dt
    }
  }
  
  if (ris<0) {
    ris=0
  }
  #ris=0
  

  percol = ks*((old_sm-sm_res)/(sm_sat-sm_res))^(2/Bci+3) # percolation [m s-1]
  
  ## Water balance(s)
  if (i<date_root) {
    sm = old_sm+(infiltration-percol-et)*dt/(soil_depth) # sm update
  } else {
    sm = old_sm+(infiltration-percol-1/3*et)*dt/(soil_depth) # sm update
  }
  
  
  # case 1: the soil dries too much
  if (sm<sm_res) {
    ds = (sm-sm_res)*soil_depth/dt
    percol = percol + ds
    sm = sm_res # SM set to minimum
    
    # if percolation is negative, ET becomes too high
    if (percol < 0) {
      ds = -percol
      percol = 0
      et = et-ds
    }
  }
  
  # case 2: the soil saturates too much
  if (sm>sm_sat) {
    ds = (sm-sm_sat)*soil_depth/dt
    infiltration = infiltration - ds
    if (infiltration<0) {
      ds = -infiltration
      infiltration = 0
    }
    sm=sm_sat
  }
  
  
  ## deep layer
  
  if (vadose>soil_depth2) {
    percol2 = ks2*((old_sm2-sm_res2)/(sm_sat2-sm_res2))^(2/Bci2+3) # percolation [m s-1]
  } else {
    percol2=0
  }
  
  if (saturation2>0.01) {
    unsHyCond=ks2*saturation2^(2/Bci2+3)
    matPot=bubpress/saturation2^(2/Bci2+3)
  } else {
    unsHyCond=ks2*0.0001^(2/Bci2+3)
    matPot=bubpress/0.0001^(2/Bci2+3)
  }
  meanHyCond=2*unsHyCond*ks2/(unsHyCond+ks2)
  if (vadose>soil_depth2/2) {
    ris2=meanHyCond*matPot/(vadose-soil_depth2/2)
  } else {
    ris2=ks2
  }
  if (ris2>ks2) {
    ris2=ks2
  }
  #ris2=0
  
  ## Water balance(s)
  if (i<date_root) {
    sm2 = old_sm2+(percol-percol2-ris+ris2)*dt/(soil_depth2) # sm update
  } else {
    sm2 = old_sm2+(percol-percol2-ris+ris2-2/3*et)*dt/(soil_depth2) # sm update
  }
  
  # case 1: the soil dries too much
  if (sm2<sm_res2) {
    ds = (sm2-sm_res2)*soil_depth2/dt
    percol2 = percol2 + ds
    sm2 = sm_res2 # SM set to minimum
    
    # if percolation is negative, ET becomes too high
    if (percol2 < 0) {
      ds = -percol2
      percol2 = 0
      et = et-ds
    }
  }
  
  # case 2: the soil saturates too much
  if (sm2>sm_sat2) {
    ds = (sm2-sm_sat2)*soil_depth2/dt
    sm2=sm_sat2
  }
  
  return(list(sm,sm2,gross_rain,net_rain,rain_flag,
              inter_time,S,percol,percol2,ris,ris2,runoff,infiltration))
}