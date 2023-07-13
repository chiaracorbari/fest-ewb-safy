safy = function(SMT,GLA,DAM,DGM,PRT,jDay,
                Pgro_Kex,Pgro_Lue,Pgro_Ms0,Pgro_Sla,Pgro_P2G,
                Pfen_MrgD,Pfen_PrtA,Pfen_PrtB,Pfen_SenA,Pfen_SenB,
                Ptfn_Tmin,Ptfn_Topt,Ptfn_Tmax,Ptfn_TpSn,
                Pgro_R2P,Tair,Rsin,FC,WP,FC2,WP2,depth,p,SM,SM2,date_root,i)
  {

  LAI_ini=Pgro_Ms0*Pgro_Sla  # Initial Value of LAI  [m2 m-2]
  if (jDay == Pfen_MrgD) { # DAY OF EMERGENCE
    # Initialization of Vegetation Model
    SMT = max(Tair-Ptfn_Tmin,0)
    DAM = Pgro_Ms0
    GLA = LAI_ini
    Day_Of_Emergence = jDay

  } else if (jDay>Pfen_MrgD && GLA>=LAI_ini) { # VEGETATIVE PERIOD

    # Temperature Sum and Stress
    if ((Tair < Ptfn_Tmin)||(Tair > Ptfn_Tmax)) {
      # T outside the functioning range
      TpS=0
    } else if (Tair <= Ptfn_Topt) {
      # T lower than the optimal value
      TpS = 1-((Tair-Ptfn_Topt)/(Ptfn_Tmin-Ptfn_Topt))^Ptfn_TpSn
    } else {
      # T higher than the optimal value
      TpS = 1-((Tair-Ptfn_Topt)/(Ptfn_Tmax-Ptfn_Topt))^Ptfn_TpSn
    }

    SMT = SMT+max(Tair-Ptfn_Tmin,0);

    # Daily Total PAR Absorbed by the Canopy, Daily Dry Mass Production
    if (Rsin>0) {
      Rglb=Rsin*86400/1000000 # [MJ/m2/day]
    } else {
      Rglb=0
    }

    PAR = Pgro_R2P*Rglb*(1-exp(-Pgro_Kex*GLA))

    # Water stress (if not, change to Ks = 1)    
    if (i<date_root) {
      wstress = (SM-WP)/((1-p)*(FC-WP))
    } else {
      smeff=mean(c(SM,SM2))
      smeff=SM2
      wstress = (smeff-WP2)/((1-p)*(FC2-WP2))
    }
    
    wstress = min(wstress,1)
    wstress = max(wstress,0.1)

    ddam=Pgro_Lue*PAR*TpS*wstress

    DAM = DAM + ddam
    # GLA

    # Partitioning, Green LAI Increase (DLP) and Leave Senescence Function (DLM)
    PRT_DAY_minus1 = PRT
    PRT = max(1-Pfen_PrtA*exp(Pfen_PrtB*SMT),0)
    DLP = ddam*PRT*Pgro_Sla

    if (SMT>Pfen_SenA) {
      DLM=GLA*(SMT-Pfen_SenA)/Pfen_SenB
    } else {
      DLM=0
    }

    GLA=GLA+DLP-DLM

    # Yield (Grain Mass increase after the leaf production period)
    if (PRT==0) {
      if (PRT_DAY_minus1>0) { # End of Leaf Growing Period
        SMT_DAY_minus1 = SMT
        Day_Of_Anthesis = jDay
        DGM=0
      } else {
        DGM=DGM+Pgro_P2G*DAM
      }
    }

    # End of Vegetation Modelling if LAI < initial value
    if (GLA<LAI_ini) {
      Day_Of_Senescence=jDay
    }
  }

  return(list(SMT,DAM,DGM,GLA,PRT))
  # SMT = SUM of TEMPERATURE for DEGREE-DAY APPROACH (??C)
  # TpS = TEMPERATURE STRESS INDEX ([0-1],unitless)
  # PAR = ABSORBED PHOTOSYNTHETIC ACTIVE RADIATION (MJ/m2/day)
  # DAM, DGM = DRY AERIAL & GRAIN MASS (g/m2=100 x t/ha)
  # PRT = PARTITION-TO-LEAF INDEX ([0-1],unitless)
  # GLA = GREEN LEAF AREA INDEX (LAI, m2/m2)
  # DLP,DLM = DELTA of GREEN LAI from DAY D to D+1  (P=plus;M=minus;m2/m2)

}
