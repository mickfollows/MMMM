c ......................................................
c DECLARATION OF VARIABLES (INCLUDE FILE) FOR
c macro-molecular model of phototrophic and heterotrophic microbes
c .......................................................
c       implicit none
c declarations
c      integer jmax
c      parameter(jmax=1)
c      integer ntmax
c      parameter(ntmax=500)
c cell size
       real*8 rcell(jmax)
       real*8 cellvol(jmax)
       real*8 init_cellvol
       real*8 Ctotal
c cell state variables
       real*8 QBc_ave(jmax)
       real*8 QBc(jmax)
       real*8 QBn(jmax)
       real*8 QBp(jmax)
       real*8 c(jmax)
       real*8 n(jmax)
       real*8 p(jmax)
       real*8 CH(jmax)
       real*8 AA(jmax)
       real*8 NUC(jmax)
       real*8 PR(jmax)
       real*8 PIG(jmax)
       real*8 Bc(jmax)
       real*8 Bn(jmax)
       real*8 Bp(jmax)
       real*8 Qc(jmax)
       real*8 Qn(jmax)
       real*8 Qp(jmax)
       real*8 QCH(jmax)
       real*8 QAA(jmax)
       real*8 QNUC(jmax)
       real*8 QPR(jmax)
       real*8 QPIG(jmax)
c maximum quotas for inorganic stores
       real*8 Qc_max(jmax)
       real*8 Qn_max(jmax)
       real*8 Qp_max(jmax)
       real*8 QCH_max(jmax)
c trophic strategy
       integer itrophic(jmax)
c cell density
       real*8 X(jmax)

c uptake related
       real*8 uptakec(jmax) 
       real*8 uptaken(jmax) 
       real*8 uptakep(jmax) 
       real*8 uptakeCH(jmax) 
       real*8 uptakeAA(jmax) 
       real*8 uptakeNUC(jmax) 
       real*8 uptakePR(jmax) 
       real*8 uptakePIG(jmax) 
       real*8 sum_uptake_c 
       real*8 sum_uptake_n  
       real*8 sum_uptake_p  
       real*8 sum_uptake_CH 
       real*8 sum_uptake_AA 
       real*8 sum_uptake_NUC
       real*8 sum_uptake_PR 
       real*8 sum_uptake_PIG
       real*8 Vmax_c(jmax) 
       real*8 Vmax_n(jmax) 
       real*8 Vmax_p(jmax) 
       real*8 Vmax_CH(jmax) 
       real*8 Vmax_AA(jmax) 
       real*8 Vmax_NUC(jmax) 
       real*8 Vmax_PR(jmax) 
       real*8 Vmax_PIG(jmax) 
       real*8 Ksat_c_up(jmax) 
       real*8 Ksat_n_up(jmax) 
       real*8 Ksat_p_up(jmax) 
       real*8 Ksat_CH_up(jmax) 
       real*8 Ksat_AA_up(jmax) 
       real*8 Ksat_NUC_up(jmax) 
       real*8 Ksat_PR_up(jmax) 
       real*8 Ksat_PIG_up(jmax) 
       real*8 Ksatc(jmax)
       real*8 Ksatn(jmax)
       real*8 Ksatp(jmax)
       real*8 KsatCH(jmax)
       real*8 KsatAA(jmax)



c synthesis related
       real*8 AAsynth(jmax) 
       real*8 NUCsynth(jmax) 
       real*8 PRsynth(jmax) 
       real*8 PIGsynth(jmax) 

c energy related
       real*8 betaE(jmax) 
       real*8 Kenergy(jmax)
c light/photosynth related
       real*8 photosynth(jmax) 
       real*8 Iharvest(jmax)
       real*8 alphaE(jmax)
       real*8 ephoto(jmax) 

c other parameters
       real*8 D
       real*8 flushtime
       real*8 dummy
       real*8 dummy1 
       real*8 dummy2 
       real*8 dummy3 
 
c timestepping 
       real*8 deltat

c energetic costs
       real*8 epsilon_uptakec(jmax) 
       real*8 epsilon_uptaken(jmax) 
       real*8 epsilon_uptakep(jmax) 
       real*8 epsilon_uptakeCH(jmax) 
       real*8 epsilon_uptakeAA(jmax) 
       real*8 epsilon_AAsynth(jmax) 
       real*8 epsilon_NUCsynth(jmax) 
       real*8 epsilon_PRsynth(jmax) 
       real*8 epsilon_PIGsynth(jmax) 
       real*8 epsilon_maintenance(jmax) 

c "half-saturation quotas"
       real*8 QCH_sat_AAsynth(jmax)     
       real*8 Qn_sat_AAsynth(jmax)     
       real*8 QCH_sat_NUCsynth(jmax)     
       real*8 QAA_sat_NUCsynth(jmax)     
       real*8 Qp_sat_NUCsynth(jmax)     
       real*8 QCH_sat_PRsynth(jmax)     
       real*8 QAA_sat_PRsynth(jmax)     
       real*8 QCH_sat_PIGsynth(jmax)     
       real*8 Qn_sat_PIGsynth(jmax)     
       real*8 Qc_sat_photo(jmax)     
c synthesis rate constants
       real*8 KAAsynth(jmax)
       real*8 KNUCsynth(jmax)
       real*8 KPIGsynth(jmax)
       real*8 GammaPR_synth(jmax)
       real*8 multiplier

c time derivatives
c  ... cell state variables
       real*8 dcdt(jmax)
       real*8 dndt(jmax)
       real*8 dpdt(jmax)
       real*8 dCHdt(jmax)
       real*8 dAAdt(jmax)
       real*8 dNUCdt(jmax)
       real*8 dPRdt(jmax)
       real*8 dPIGdt(jmax)
c ... extra array for predictor-corrector
       real*8 dcdt_o(jmax)
       real*8 dndt_o(jmax)
       real*8 dpdt_o(jmax)
       real*8 dCHdt_o(jmax)
       real*8 dAAdt_o(jmax)
       real*8 dNUCdt_o(jmax)
       real*8 dPRdt_o(jmax)
       real*8 dPIGdt_o(jmax)
c  ... medium variables
       real*8 dcmdt
       real*8 dnmdt
       real*8 dpmdt
       real*8 dCHmdt
       real*8 dAAmdt
       real*8 dNUCmdt
       real*8 dPRmdt
       real*8 dPIGmdt
c ... extra array for predictor-corrector
       real*8 dcmdt_o
       real*8 dnmdt_o
       real*8 dpmdt_o
       real*8 dCHmdt_o
       real*8 dAAmdt_o
       real*8 dNUCmdt_o
       real*8 dPRmdt_o
       real*8 dPIGmdt_o

c medium concentrations
       real*8 c_m
       real*8 n_m
       real*8 p_m
       real*8 CH_m
       real*8 AA_m
       real*8 NUC_m
       real*8 PR_m
       real*8 PIG_m
c incoming medium concentrations
       real*8 c_in_m
       real*8 n_in_m
       real*8 p_in_m
       real*8 CH_in_m
       real*8 AA_in_m
       real*8 NUC_in_m
       real*8 PR_in_m
       real*8 PIG_in_m


c exudation  
       real*8 exCH(jmax)
       real*8 exAA(jmax)
       real*8 exNUC(jmax)
       real*8 exPR(jmax)
       real*8 exPIG(jmax)
       real*8 exude_CH(jmax)
       real*8 exude_AA(jmax)
       real*8 exude_NUC(jmax)
       real*8 exude_PR(jmax)
       real*8 exude_PIG(jmax)

       real*8 sum_exude_CH 
       real*8 sum_exude_AA 
       real*8 sum_exude_NUC
       real*8 sum_exude_PR 
       real*8 sum_exude_PIG



c elemental ratios 
       real*8 Rnc_AA(jmax)
       real*8 Rnc_PR(jmax)
       real*8 Rnc_NUC(jmax)
       real*8 Rnc_PIG(jmax)
       real*8 Rpc_NUC(jmax)


c other stuff
       real*8 approx_growth_timescale
       real*8 cellvol_um


c scaling
      real*8 X_out
      real*8 c_out
      real*8 n_out
      real*8 p_out
      real*8 Bc_out
      real*8 c_m_out
      real*8 n_m_out
      real*8 p_m_out
      real*8 CH_out
      real*8 AA_out
      real*8 NUC_out
      real*8 PR_out
      real*8 PIG_out
      real*8 CH_m_out
      real*8 AA_m_out
      real*8 NUC_m_out
      real*8 PR_m_out
      real*8 PIG_m_out

