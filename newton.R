## The crop-energy-water balance FEST-EWB-SAFY model (Corbari et al., 2022) couples 
## the distributed energy-water balance FEST-EWB (Corbari et al, 2011), which allows computing continuously
## in time and distributed in space both soil moisture and evapotranspiration fluxes, 
## and the SAFY (Duchemin et al, 2008), simple model for yield prediction and plant development

## [1] Corbari, C., Ravazzani, G. and Mancini, M. (2011), A distributed thermodynamic model for energy and mass balance 
## computation: FEST-EWB. Hydrol. Process., 25: 1443-1452. https://doi.org/10.1002/hyp.7910
## [2] Corbari, C., Ben Charfi, I., Al Bitar, A., Skokovic, D., Sobrino, J.A., Perelli, C., Branca, G., Mancini, M. (2022),
## A fully coupled crop-water-energy balance model based on satellite data for maize and tomato crops yield estimates: 
##   The FEST-EWB-SAFY model, Agr. Wat. Man., 272: 107850. https://doi.org/10.1016/j.agwat.2022.107850\n
## [3] Duchemin, B., Maisongrande, P., Boulet, G., Benhadj, I., 2008. A simple algorithm for yield estimates: 
##   Evaluation for semi-arid irrigated winter wheat monitored with green leaf area index. Environ. Modell. Softw. 
## 23(7), 876-892. https://doi.org/10.1016/j.envsoft.2007.10.003


##  The code has been written by Nicola Paciolla and Chiara Corbari
## chiara.corbari@polimi.it


newton = function(ums,ums2,Ta,Rsd,dz,wind,P,uma,rsmin,
                  rsoil,smr,sms,bci,wpo,bpr,fca,alb,fveg,
                  laima,Ts,date_root,i) {
  
  iter = 0 # no. of iterations
  TsK = Ts+273.15 # Kelvin surface temperature
  
  #hv = 0.5 # !!!!!!!!!!!!!!!!!!!
  hv = fveg*0.6
  
  elev = 49 # elevation above sea level
  cl_index = 0 # cloud index (0 = no clouds)
  
  ess = 0.98 # soil emissivity [-]
  R = 0.287 # gas-specific constant [kJ kg-1 K-1]
  rhowat = 1000 # water density [kg m-3]
  Bo = 0.0000000567 # Boltzmann constant [W m-2 K-4]
  k = 0.41 # von Karman constant
  # soil-resistance coefficients
    a1 = 3.5
    a2 = 2.3
    a3 = 33.5
  # cloudy-sky atmospheric emissivity coefficients
    b1 = 0.2
    b2 = 2	
  E = 0.622	# fraction of water vapour molecular weight (MW) against dry-air MW
  
  
  hw = 2.5 # height of wind measurement [m above ground level]
  hrh = 2.5 # height of relative humidity measurement [m above ground level]
  
  # Wind minimum setting
  if (wind<1) {
    ws=1
  } else {
    ws=wind
  }
  
  # Aerodynamic variables for bare soil [m]
  z1 = 0.1
  zm1 = hw
  zd1 = 0.6666*z1
  zrw1 = 0.123*z1
  zh1 = hrh
  zrh1 = 0.1*zrw1
  
  # Aerodynamic resistance for bare soil [s m-1]
  rabs = (log((zm1-zd1)/zrw1) * log((zh1-zd1)/zrh1)) / ((k^2)*ws)
  
  # Aerodynamic variables for vegetated soil [m]
  zm = hw-hv
  zd = 0.6666*hv
  zrw = 0.123*hv
  zh = hrh-hv
  zrh = 0.1*zrw
  
  # Aerodynamic resistance allocation
  if (hv==0) {
    ra=rabs
  } else {
    # Aerodynamic resistance [m s-1]
    ra = (log((zm-zd)/zrw) * log((zh-zd)/zrh)) / ((k^2)*ws)
  }
    
  l = 2.501-(0.002361*Ta) # latent heat of vapourization [MJ kg-1]
  P = 101.3*((293-0.0065*elev)/293)^5.26 # air pressure [kPa]
  gpsi = 0.00163*P/l # psychrometric constant [kPa ?C-1]
  cp = gpsi*E*l/P # humid air specific heat [ MJ kg-1 ?C-1]
  roa = 3.486*P/(275+Ta) # air density [kg m-3]
  
  # Soil thermal conductivity [W m-1 K-1]
  if (ums<=smr) {
    gthermc=418.6*0.00041
  } else if (ums>=sms) {
    gthermc=418.6*0.00041
  } else {
    SMrel = (ums-smr)/(sms-smr)
    psi = bpr/(SMrel^(1/bci))
    psicm = psi*100
    pf = log10(psicm)
    if (pf<=5.1) {
      gthermc=418.6*exp(-(pf+2.7))
    } else {
      gthermc=418.6*0.00041
    }
  }
  
  saturation=ums/sms
  if (saturation>=0 && saturation<0.25) {
    gthermc=3
  } else if (saturation>=0.2 && saturation<0.54) {
    gthermc=5.05
  } else if (saturation>=0.54 && saturation<0.6) {
    gthermc=6.2
  } else if (saturation>=0.6 && saturation<0.75) {
    gthermc=7.5
  } else if (saturation>=0.75) {
    gthermc=8.8
  }
  gthermc=1.8*saturation^0.4781
  
  # Tzero day/night computation
  if (Rsd>150) {
    Tzero = 0.6607 * Ta + 5.0789 + 273.15
  } else {
    Tzero = 0.6607 * Ta + 8.5789 + 273.15
  }
  
  # Air emissivity day/night computation
  if (Rsd>50) {
    #eac = 0.000007*((Ta+273.15)^2)
    eac = 0.000010*((Ta+273.15)^2)
  } else {
    # eac=0.0000092*((Ta+273.15)^2) # Zillman
    eac = 0.00001*((Ta+273.15)^2)
  }
    
  # Effective vegetation resistance
    if (i<date_root) {
      ums1=ums
    } else {
      ums1=ums2
    }
  #ums1=ums

  if (fveg>0) {
    if (ums==smr) {
      # rc=(rsmin/laima)
      ums1 = smr+0.01
      rc = (rsmin/laima)*((sms-smr)/(ums1-smr))
    } else {
      rc = (rsmin/laima)*((sms-smr)/(ums1-smr))
    }
  } else {
    rc = 0
  }
  
  # Soil resistance
  rs = (a1*((sms/ums1)^a2))+a3+rsoil
  
  ## Newton-Raphson closure of energy balance
  countiter=0
  while (iter<0.9) {
    
    # Water vapour tension
    es = 0.6108*exp((17.27*Ts)/(Ts+237.3))
    ea = 0.6108*exp((17.27*Ta)/(Ta+237.3))*uma/100
    dftfac=(237*17.27)/((237.3+Ts)^2)
    
    # Energy balance
    ft=Rsd*(1-alb)+
      (eac*Bo*((Ta+273.15)^4))*(1+b1*((1-cl_index)^b2))-
      ess*Bo*(TsK^4)-
      (gthermc/dz)*(TsK-Tzero)-
      fveg*1000000*roa*cp*(es-ea)/(gpsi*(ra+rc))-
      (1-fveg)*1000000*roa*cp*(es-ea)/(gpsi*(rabs+rs))-
      fveg*1000000*roa*cp*(Ts-Ta)/ra-
      (1-fveg)*1000000*roa*cp*(Ts-Ta)/rabs
    
    # Derivative of energy balance
    dftdt=-4*ess*Bo*(TsK^3)-
      gthermc/dz-
      fveg*1000000*roa*cp*es*dftfac/(gpsi*(ra+rc))-
      (1-fveg)*1000000*roa*cp*es*dftfac/(gpsi*(rabs+rs))-
      fveg*1000000*roa*cp/ra-
      (1-fveg)*1000000*roa*cp/rabs
    
    deltnr=ft/dftdt
    TsK=TsK-deltnr
    Ts=TsK-273.15
    
    if ((abs(deltnr)<0.01)&&(abs(ft)<0.01)) {
      iter=1
    } else {
		countiter=countiter+1
	}
	if (countiter>100) {
		break
	}
  }
  
  # Final energy fluxes
  Rldc = eac*Bo*((Ta+273.15)^4) # downwelling longwave radiation [W m-2]
  Rld = Rldc*(1+b1*((1-cl_index)^b2)) # cloudliness correction
  Rlu = ess*Bo*(TsK^4) # upwelling longwave radiation [W m-2]
  Rnet = Rsd*(1-alb)+Rld-Rlu # net radiation [W m-2]
  LE = fveg*1000000*roa*cp*(es-ea)/(gpsi*(ra+rc))+
    (1-fveg)*1000000*roa*cp*(es-ea)/(gpsi*(rabs+rs)) # latent heat [W m-2]
  ETeff = LE/(rhowat*l*1000000) # Evapotranspiration [???]
  
  G = (gthermc/dz)*(TsK-Tzero) # Soil Heat Flux [W m-2]
  Hs = fveg*1000000*roa*cp*(Ts-Ta)/ra+
    (1-fveg)*1000000*roa*cp*(Ts-Ta)/rabs # Sensible heat flux [W m-2]
  
  balance = Rnet-G-LE-Hs
  
  
  if (abs(balance)>0.1) {
    print(balance)
    break
  }
  
  return(list(Rnet,G,Hs,LE,Ts,rc))
}