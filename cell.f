c.......................................................
c.......................................................
       program cell
c ......................................................
c macro-molecular model of phototrophic and heterotrophic microbes
c density based formulation, simple chemostat environment
c Mick Follows, Jason Bragg, Oliver Jahn 2010-2011
c .......................................................

       implicit none
c declarations
       integer jmax
       parameter(jmax=1)
       real*8 timemax
c one day = 86400.0 seconds
c two days = 172800.0 seconds
c three days 259200.0 seconds
cc      parameter(timemax=1400000.0d0)
       parameter(timemax=518400.0d0)
c      parameter(timemax=350000.0d0)
c      parameter(timemax=3000000.0d0)
c      parameter(timemax=6000000.0d0)
       integer ntmax
c import include files
       include "CELL_ARRAYS.h"
c local variables
       integer j
       integer nt
       integer nout
       real*8 dtout
       real*8 time 
       real*8 time_out
       real*8 ab
       real*8 TotalC
       real*8 TotalC_out
       real*8 TotalN
       real*8 TotalN_out
       real*8 TotalP
       real*8 TotalP_out
       real*8 Bn_out
       real*8 Bp_out
       real*8 total_respiration(jmax)
       real*8 sum_respiration
       real*8 Pphoto 
       real*8 pigfrac
       real*8 coo    
       integer idone
       integer idilute
       integer idilutemax  

c arrays for the Davidson lab data
       integer ic
       integer icmax
       parameter(icmax=19)
       real*8 labt(icmax)
       real*8 labNm(icmax)
       real*8 labX(icmax)
       real*8 labCN(icmax)
c ... and arrays for cost function stuff
       integer icout
       integer ic1
       real*8 modt(icmax)
       real*8 modNm(icmax)
       real*8 modX(icmax)
       real*8 modCN(icmax)
       real*8 timetest
c cost function
       real*8 cdummy
       real*8 cdummy1
       real*8 cost
       real*8 mincost
c factor for max QCH
       real*8 QCH_max_factor
c some initialization variables
       real*8 xinit 
       real*8 GammaPR_synth_init 
       real*8 QCH_max_factor_init 
c some loop variables
       integer ione
       integer itwo
       integer ithree
       integer ifour 
       integer ifive 
       integer isix  
       integer iseven
c photosynth factor
       real*8 pfactor
       real*8 pfactor_init
c other factor
       real*8 KAAfactor
       real*8 VmaxFactor
       real*8 SizeFactor


c READ IN THE DAVIDSON LAB DATA
       write(6,*)'read the Davidson Tp data'
       open(20, file='labdata.dat', status='old')
       read(20,*) 
       read(20,*) 
       write(6,*)' labt    labNm    labX   labCN '
       do ic = 1, icmax 
          read(20,21)labt(ic),labNm(ic),labX(ic),labCN(ic)
          write(6,*)labt(ic),labNm(ic),labX(ic),labCN(ic)
       end do
 21    format(4e14.6)
       close(20)
c set the counter for saving model output/eval cost function
       ic1 = 0

c SET SOME TIMESTEPPING AND OUTPUT VARIABLES
c      timestep (s)
c      deltat =  0.1d0
c      deltat =  0.001d0
       deltat =  0.1d0
       write(6,*)'timestep (s) ',deltat
c max number of timesteps
       ntmax = timemax / deltat

c How often to write out results?
c   ... dtout = time between output (s)
cc       dtout = 300.0d0
c make dtout 1/4 day...
       dtout = 0.25*86400.0
       nout = dtout / deltat 

c set icout - count between eval of cost function
       icout = nout

c initialize time
       time = 0.0d0
c initialize idone - controls writing header to ouput
       idone = 0

c SET TRAIT AND PARAMETER VALUES

c "STANDARD" ......................
c First set control variables - Defaults 
c initial X
c      xinit = 50.0d9
c GammaPRsynth
c      GammaPR_synth_init = 0.9d-4
c QCH_max
c      QCH_max_factor_init = 0.3
c pfactor
c      pfactor_init = 0.01
c KAAfactor
c      KAAfactor = 1.5
c VmaxFactor
c      VmaxFactor = 1.0
c SizeFactor
c      SizeFactor = 1.0
c....................................

c First set control variables - Defaults 
c initial X
       xinit = 20.0d9
c GammaPRsynth
c      GammaPR_synth_init = 0.6d-4
       GammaPR_synth_init = 0.1d-3
c QCH_max
       QCH_max_factor_init = 0.6
c pfactor
       pfactor_init = 0.01
c KAAfactor
       KAAfactor = 1.2
c VmaxFactor
       VmaxFactor = 2.0
c SizeFactor
       SizeFactor = 0.8
      
c set initial value (high) for minimum cost!!!
       mincost = 1.0d12

c open results file
       open(27,file='costs.dat')
c LOOP THROUGH CONTROL VARIABLES
c      do 700 ione = 1,5 
c       xinit = 1.0d10*ione 
c      do 701 itwo = 1,5 
c       GammaPR_synth_init = 2.0d-5*itwo
c      do 702 ithree = 1,5
c       QCH_max_factor_init = 0.2*ithree
c      do 703 ifour = 1,5
c       pfactor_init = 0.002*ifour
c      do 704 ifive = 1,5 
c       KAAfactor = 0.4*ifive 
c      do 705 isix = 1,5 
c       VmaxFactor = 0.4*isix   
c      do 706 iseven = 1,5 
c       SizeFactor = 0.4*iseven 

c initialize state variables
c --- initialize cell radius, volume and number density 
        do j=1,jmax
c specify trophic strategy
c ... 0 = heterotroph,  1 = phototroph,  2 = mixotroph
          itrophic(j) = 1
c specify initial cell number density (cell m-3)
c 1 cell ml-1 = 106 cell m-3
c         X(j) =  5.0d8
c         X(j) =  50.0d9
c          write(6,*)'!!!! setting initial X(j)....' 
          X(j) =  xinit 

c cell radius (m)
cMICK-10july13     rcell(j) = 1.0d-6
          rcell(j) = 1.5d-6
          rcell(j) = rcell(j) * SizeFactor
c         rcell(j) = 3.0d-6
c         rcell(j) = 5.0d-6
c cell volume (m3)
          cellvol(j) = 1.33333d0 * 3.14159 * rcell(j)*rcell(j)*rcell(j)
c average total cellular carbon content, QBc_ave(j)
c based on empirical size relationships 
          init_cellvol = cellvol(j) * 1.0d18 
c needs input in um3 and output in fmol C cell-1
          call carbon_from_vol(init_cellvol,Ctotal)
          QBc_ave(j) = Ctotal * 1.0d-15
c initial Bc(j) 
          Bc(j) = X(j) * QBc_ave(j)
c........................................................
c         write(6,*)'j',j,'   QBc  ',QBc_ave(j),'   Bc   ',
c    &        Bc(j),'   X ',X(j),' init_cellvol ',init_cellvol
c........................................................
        enddo

c --- initial cell state variables: bulk density (mol m-3) 
      do 200 j=1,jmax

c inorg carbon
        if(itrophic(j) .eq. 0)then
           c(j) = Bc(j)*0.001d0 
        else
           c(j) = Bc(j)*0.02d0 
        end if 
c macromolecules (protein as resid to ensure conservation)
        CH(j) = Bc(j) * 0.3d0 
        AA(j) = Bc(j) * 0.3d0 
        NUC(j) = Bc(j) *0.0001d0 
        if(itrophic(j) .eq. 0)then
          PR(j) = Bc(j) - ( c(j) + CH(j) + AA(j) + NUC(j) )
          PIG(j) = 0.0d0 
        else
         pigfrac = 0.05d0
cc         pigfrac = 0.0d0
          PR(j) = (Bc(j) - ( c(j) + CH(j) + AA(j) + NUC(j) )) 
     &             * (1.0d0 - pigfrac)
          PIG(j) = (Bc(j) - ( c(j) + CH(j) + AA(j) + NUC(j) ))*pigfrac
        endif
c inorg nitrogen
        n(j) = c(j)*1.0
c inorg phosphorus
        p(j) = c(j) 
 200  end do

c STOICHIOMETRY OF BIOCHEMICAL COMPONENTS
c        NC_amino(i) = 1.0d0/4.0d0
c        PC_protein(i) = 0.0d0
c        PC_nucleic(i) = 1.0d0/8.0d0
c initialize elemental ratios
       do j=1,jmax
         Rnc_AA(j) = 0.25d0
         Rnc_PR(j) = Rnc_AA(j)
         Rpc_NUC(j) = 0.125d0 
         Rnc_NUC(j) = Rnc_AA(j)
         Rnc_PIG(j) = 0.07d0 
       end do

c set uptake parameters according to cell size
c Vmax ~ E Kh  and  Ksat ~ Kh/Ke then Vmax ~ R^2 and Ksat^R
c where R is cell radius
c base nitrate on Litchman et al (2007) 
        do j=1,jmax
c Litchman et al, volume relationship: consistent with Aksnes and Egge, 1991
c Vmax(N) (umol N cell-1 day-1)
c ... Litchman uses cell vol in micrometres cubed so...
         cellvol_um = cellvol(j)*1.0d18       
         Vmax_n(j) = 9.1d-9 * ( (cellvol_um)**0.67 )
c TEST         Vmax_n(j) = 5.0 *  9.1d-9 * ( (cellvol_um)**0.67 )
         Vmax_n(j) = Vmax_n(j) * VmaxFactor
c ... convert to moles cell-1 s-1
         Vmax_n(j) =  Vmax_n(j) * 1.0d-6 / 86400.0d0
c        Half saturation (umol N ) from Litchman --- NOTE MICROMOLAR!
         Ksat_n_up(j)  = 0.17d0 * (cellvol_um)**0.27 
c convert to moles m-3
         Ksat_n_up(j)  =  Ksat_n_up(j) * 1.0d3 / 1.0d6
c testing ................................................
c        write(6,*)'j, volume (um3) ',j, cellvol_um
c        write(6,*)'j, test Vmax_n (moles cell-1 s-1)',j, vmax_N(j) 
c        write(6,*)'j, test Ksat_n_up (moles m-3)',j, Ksat_n_up(j) 
c testing ................................................

c scale for other substrates according to (Redfield) requirement
c   --- ok for Vmax because scales with density of transporters which
c       could be allocated accordingly(?)
c VERY SKETCHY.....
c         Vmax_c(j) = Vmax_n(j)*(106.0d0/16.0d0)
c heterotroph, no DIC uptake
          if(itrophic(j) .eq. 0)then
              Vmax_c(j) = 0.0d0
          else
              Vmax_c(j) = Vmax_n(j)*(106.0d0/16.0d0)
c             Vmax_c(j) = Vmax_n(j)*2.0d0
c             Vmax_c(j) = Vmax_n(j)*20.0d0
          endif
          Ksat_c_up(j)  = Ksat_n_up(j)*(106.0d0/16.0)
          Vmax_p(j) = Vmax_n(j)*(1.0d0/16.0d0)
          Ksat_p_up(j)  = Ksat_n_up(j)*(1.0d0/16.0d0)
          if(itrophic(j) .ne. 1)then
             Vmax_CH(j) = 1.5 * Vmax_n(j)*(106.0d0/16.0d0)
          else 
             Vmax_CH(j) = 0.0d0
          endif 
          Ksat_CH_up(j)  = Ksat_n_up(j)*(106.0d0/16.0)
c initially no uptake of amino acid, nucleic acid, protein or pigment
          Vmax_AA(j) = 0.0d0
          Ksat_AA_up(j)  = Ksat_n_up(j)
          Vmax_NUC(j) = 0.0d0
          Ksat_NUC_up(j)  = Ksat_n_up(j)
          Vmax_PR(j) = 0.0d0 
          Ksat_PR_up(j)  = Ksat_n_up(j)
          Vmax_PIG(j) = 0.0d0
          Ksat_PIG_up(j)  = Ksat_n_up(j)
c testing ................................................
c        write(6,*)'j, test Vmax_c (moles cell-1 s-1)',j, Vmax_c(j) 
c        write(6,*)'j, test Vmax_p (moles cell-1 s-1)',j, Vmax_p(j) 
c        write(6,*)'j, test Ksat_p_up (moles m-3)',j, Ksat_p_up(j) 
c        write(6,*)'j, test Vmax_CH (moles cell-1 s-1)',j, Vmax_CH(j) 
c        write(6,*)'j, test Ksat_CH_up (moles m-3)',j, Ksat_CH_up(j) 
c testing ................................................
        end do

        do j=1,jmax
c maximum quotas for internal inorganic stores
          Qp_max(j) = QBc_ave(j)*(1.0d0/106.0d0)*0.2
          Qn_max(j) = QBc_ave(j)*(16.0d0/106.0d0)*0.2
c MICK 8 july 2013          Qc_max(j) = QBc_ave(j)*0.3 
          Qc_max(j) = QBc_ave(j)*0.1
c max quota for CH - stop photosynth at some point if CH accumulating
c         QCH_max_factor = 0.3
c         write(6,*)'!!! setting QCH_max_factor_init ... '
          QCH_max_factor = QCH_max_factor_init 

          QCH_max(j) = QBc_ave(j)*QCH_max_factor
c         QCH_max(j) = QBc_ave(j)*0.5
        end do
c testing ...............................................


c synthesis parameters
c half-saturations expressed as cell quotas 
c because cellular internal conc = cell quota / cell vol
        do j=1,jmax
c AA synth
c synthesis rate in moles C cell-1 s-1
c Assume can produce one cells worth of protein per division period.
c Most of cell carbon must pass through AA stage... 
c so estimate 
c              KAAsynth ~ QBc_ave / division period
c where min division period is a few hours
c ...................................................
c         write(6,*)'MAKE SYNTH RATES SIZE DEPENDENT??????'
c ...................................................
          approx_growth_timescale = 24.0d0*3600.0d0 
c TESTING ...................................
          KAAsynth(j) = QBc_ave(j) / approx_growth_timescale  
cMICK-8jul13    KAAsynth(j) = KAAsynth(j)*1.5
c         KAAfactor = 1.5
          KAAsynth(j) = KAAsynth(j)*KAAfactor
c         KAAsynth(j) = KAAsynth(j)*0.3
c         KAAsynth(j) = KAAsynth(j)*0.1
c         write(6,*)'j ',j,'KAAsynth ',KAAsynth(j)
c ...........................................

c "half-saturation quotas" in moles cell-1
c   .... if QBc_ave is average cell carbon content
c        and typical CH carbon fraction of that (25% for phototroph,
c        few % for heterotroph) so make half-saturation CH quota about
c        one percent of average total carbon
cc        QCH_sat_AAsynth(j) = QBc_ave(j)*0.01d0
          QCH_sat_AAsynth(j) = QBc_ave(j)*0.1d0
c    ... Make inorganic nitrogen store small 
c    ... try a small fraction of Redfieldian N:C 
          Qn_sat_AAsynth(j) =  QBc_ave(j)*(16.0d0/106.0d0)*0.05
c NUC synth
c TESTING ................................
          KNUCsynth(j) = KAAsynth(j) * 0.2
c         KNUCsynth(j) = 0.0d0        
c ........................................
c   .... if QBc_ave is average cell carbon content
c        and typical AA carbon just a few percent (?? less) of that,
c        then make min AA carbon quota a fraction of one 
c        percent of average total carbon
          QAA_sat_NUCsynth(j) =  QBc_ave(j)*0.01d0
c    ... Make inorganic nitrogen store small 
c    ... try a small fraction of Redfieldian P:C 
          Qp_sat_NUCsynth(j) =   QBc_ave(j)*(1.0d0/106.0d0)*0.02
c    ... carbo hydrate half sat same as for AA synth
          QCH_sat_NUCsynth(j) =  QCH_sat_AAsynth(j) 
c PR synth
c protein synth rate proportional to NUC
c   ... so GammaPR_synth is s-1
c   ... since protein is half of cell, must be able to reproduce protein
c   ... content on division timescale or greater, accounting for 
c   ... efficiency of NUC (hence "multiplier")
c         multiplier = 5.0d0
c         multiplier = 10.0d0
c         GammaPR_synth(j) = multiplier 
c    &                        * (1.0d0 / approx_growth_timescale)
c TESTING ....................................
c         GammaPR_synth(j) = 4.0d-4
c         GammaPR_synth(j) = 0.4d-4
c         GammaPR_synth(j) = 0.9d-4
c         write(6,*)'!!! setting GammaPR_synth ... '
          GammaPR_synth(j) = GammaPR_synth_init

c ............................................
          QAA_sat_PRsynth(j) =  QAA_sat_NUCsynth(j) 
          QCH_sat_PRsynth(j) =  QCH_sat_AAsynth(j)
c PIG synth
c ... if heterotroph, no pigment synthesis
          KPIGsynth(j) = 0.0d0 
          QCH_sat_PIGsynth(j) =  QCH_sat_AAsynth(j) 
          Qn_sat_PIGsynth(j) =  Qn_sat_AAsynth(j)
c Photosynthesis
c ... carbon limitation
          Qc_sat_photo(j) = QBc_ave(j)*0.01d0
        end do

c energetic costs (mol CH cost per mol C synth etc ...)
        do j=1,jmax
c epsilon_PRsynth approx = 1/yield factor
c     ... ultimately will be function of substrate...
c scale other costs accordingly ...
c         epsilon_PRsynth(j) = 1.0d0
          epsilon_PRsynth(j) = 2.0d0
c         epsilon_PRsynth(j) = 1.5d0
          epsilon_AAsynth(j) = 0.1d0*epsilon_PRsynth(j)
          epsilon_NUCsynth(j) = 0.1d0*epsilon_PRsynth(j)
          epsilon_PIGsynth(j) = 0.5d0*epsilon_PRsynth(j)
c uptake costs
          epsilon_uptakec(j) = 0.01d0*epsilon_PRsynth(j)
          epsilon_uptaken(j) = 0.01d0*epsilon_PRsynth(j)
          epsilon_uptakep(j) = 0.01d0*epsilon_PRsynth(j)
          epsilon_uptakeCH(j) = 0.01d0*epsilon_PRsynth(j)
          epsilon_uptakeAA(j) = 0.01d0*epsilon_PRsynth(j)
c maintenance cost
c         write(6,*)'WARNING: epsilon_maintenance set to 0.0 !!!!!'
c         write(6,*)'WARNING: epsilon_maintenance set to 0.0 !!!!!'
c         write(6,*)'WARNING: epsilon_maintenance set to 0.0 !!!!!'
c         epsilon_maintenance(j) = 0.0
          epsilon_maintenance(j) = 0.5d0 / 86400.0


        end do


c set exudation rate constants
c exCH etc are rate constants (s-1) which determine how rapidly
c macromolecular pools are lost to medium
        do j=1,jmax
c         exCH(j) = 1.0d-7
c         exCH(j) = 5.0d-7
          exCH(j) = 0.0d0 
          exAA(j) = exCH(j) 
          exNUC(j) = exCH(j)
          exPR(j) = exCH(j)
          exPIG(j) = exCH(j)
        end do

c initialize medium concentrations
c 1 micromol kg-1 ~= 10-3 mol m-3 
c       c_m = 30.0d-3 
        c_m = 2000.0d-3 
        n_m = 10.0d-3
cMICK-10jul13        n_m = 8.0d-3
        p_m = 2.0d-3
        CH_m = 30.0d-3 
        AA_m = 0.0d0 
        NUC_m = 0.0d0 
        PR_m = 0.0d0
        PIG_m = 0.0d0
c initialize incoming medium concentrations
c 1 micromol kg-1 ~= 10-3 mol kg-1
        c_in_m = 2000.0d-3 
        n_in_m = 30.0d-3
c       n_in_m = 0.3d-3
c       p_in_m = 0.001d-3
        p_in_m = 10.0d-3
        CH_in_m = 30.0d-3 
c       CH_in_m = 300.0d-3 
        AA_in_m = 0.0d0 
        NUC_in_m = 0.0d0 
        PR_in_m = 0.0d0
        PIG_in_m = 0.0d0

c write out some initial values
c       write(6,*)'INITIAL VALUES........' 
c       do j=1,jmax
c          write(6,900)nt,j,X(j),Bc(j),n_m,CH(j),AA(j),NUC(j),
c    &           PR(j),PIG(j) 
c       end do
c       write(6,*)'......................' 

c initialize timestep arrays
        do j=1,jmax
           dcdt(j) = 0.0d0
           dndt(j) = 0.0d0
           dpdt(j) = 0.0d0
           dCHdt(j) = 0.0d0
           dAAdt(j) = 0.0d0
           dNUCdt(j) = 0.0d0
           dPRdt(j) = 0.0d0
           dPIGdt(j) = 0.0d0
        end do
       dcmdt = 0.0d0
       dnmdt = 0.0d0
       dpmdt = 0.0d0
       dCHmdt = 0.0d0
       dAAmdt = 0.0d0
       dNUCmdt = 0.0d0
       dPRmdt = 0.0d0
       dPIGmdt = 0.0d0
c initialize predictor-corrector timestep arrays
        do j=1,jmax
           dcdt_o(j) = 0.0d0
           dndt_o(j) = 0.0d0
           dpdt_o(j) = 0.0d0
           dCHdt_o(j) = 0.0d0
           dAAdt_o(j) = 0.0d0
           dNUCdt_o(j) = 0.0d0
           dPRdt_o(j) = 0.0d0
           dPIGdt_o(j) = 0.0d0
        end do
        dcmdt_o = 0.0d0
        dnmdt_o = 0.0d0
        dpmdt_o = 0.0d0
        dCHmdt_o = 0.0d0
        dAAmdt_o = 0.0d0
        dNUCmdt_o = 0.0d0
        dPRmdt_o = 0.0d0
        dPIGmdt_o = 0.0d0

c open output file to compile results from several integrations
       open(14,file='chemostat.dat')
 
c how many dilution rates to try?
c      idilutemax = 15 
       idilutemax = 1 
       do 320 idilute = 1, idilutemax

c cycle through flushing time (s) and dilution rate (s-1)
c        flushtime = (0.5d0 + idilute*0.25d0)*86400.0
c        flushtime = 2.5d0*86400.0
c        D=1.0d0/flushtime

c BATCH CULTURE,  set idilutemax=1 and D=0.0
         D=0.0d0
c        D = (0.05 + idilute*0.1d0) / 86400.0

c open output file for single integration
         open(12,file='output.dat')

c TIME STEPPING .................................
         do 300 nt = 1, ntmax

c update time
           time = (nt-1)*deltat
c TEST .......................................
c        write(6,*)'deltat  ',deltat,'  nt ',nt,'   time (s) ',time
c ............................................
c time out is time in days
           time_out = time / 86400.0d0

c BIOLOGICAL PROCESSES................
c step forward macromolecular components
c total cell carbon actually "fixed" but this will give dynamic 
c macromolecular and elemental ratios
          do j = 1,jmax
c total carbon biomass (mol C m-3 in medium) and cell quota (mol C cell-1)
c number density, X(j) (cells m-3)
            Bc(j) = c(j) + CH(j) + AA(j) + NUC(j) + PR(j) + PIG(j)
            X(j) = Bc(j) / QBc_ave(j) 
c total nitrogen biomass (mol N m-3 in medium) and cell quota (mol N cell-1)
            Bn(j) = n(j) + AA(j)*Rnc_AA(j) + NUC(j)*Rnc_NUC(j) 
     &             + PR(j)*Rnc_PR(j) + PIG(j)*Rnc_PIG(j)
            QBn(j) = Bn(j)/X(j)
c total phosphorus biomass (mol P m-3 in medium) and cell quota (mol P cell-1)
            Bp(j) = p(j) + NUC(j)*Rpc_NUC(j)
            QBp(j) = Bp(j)/X(j) 
c cell quotas of macromolecular components
c (mol C cell-1) ....  i.e.   ( mol C m-3 medium ) / ( cell m-3 medium )
            QCH(j)  = CH(j)/X(j) 
            QAA(j)  = AA(j)/X(j) 
            QNUC(j) = NUC(j)/X(j) 
            QPR(j)  = PR(j)/X(j) 
            QPIG(j) = PIG(j)/X(j) 

            Qn(j) = n(j)/X(j)
            Qp(j) = p(j)/X(j)
          end do

c write out some initial values
c       write(6,*)'TEST 1 ...............' 
c       do j=1,jmax
c          write(6,900)nt,j,X(j),Bc(j),n_m,CH(j),AA(j),NUC(j),
c    &           PR(j),PIG(j) 
c       end do
c       write(6,*)'......................' 

c CELLULAR PROCESSES ........................

c NB CAREFUL WITH INORGANIC CARBON
c   ---- phototrophs take it up
c   ---- heterotrophy drives loss 
c   ADD A LOSS IF c_b > c_m OR SOMETHING ????

c uptake terms
          do j = 1,jmax
c simple limitation - if quota reaches max then dont take anything up
c .TESTING..........................................
c        write(6,*)'Qc_max, Qc ',Qc_max(j), Qc(j) 
c        write(6,*)'Qn_max, Qn ',Qn_max(j), Qn(j) 
c        write(6,*)'Qp_max, Qp ',Qp_max(j), Qp(j) 
c .TESTING..........................................
            if(Qc(j) .ge. Qc_max(j))then
              uptakec(j) = 0.0d0
            else
c COULD ASSUME UPTAKE OF C REGULATED BY CO2 WHICH IS BUFFERED AT 30.0 umol l-1
c DODGY, BUT MAY BE MORE APPROPRIATE THAN USING FULL DIC POOL???
cccc           coo = 30.0d-3
            coo = 1000.0d-3
c             uptakec(j) = Vmax_c(j) * coo / (coo + Ksat_c_up(j)) 
              uptakec(j) = Vmax_c(j) * c_m / (c_m + Ksat_c_up(j)) 
            endif
            if(Qn(j) .ge. Qn_max(j))then
               uptaken(j) = 0.0d0
            else
               uptaken(j) = Vmax_n(j) * n_m / (n_m + Ksat_n_up(j)) 
            endif
            if(Qn(j) .ge. Qn_max(j))then
               uptakep(j) = 0.0d0 
            else
               uptakep(j) = Vmax_p(j) * p_m / (p_m + Ksat_p_up(j)) 
            endif
            uptakeCH(j) = Vmax_CH(j) * CH_m / (CH_m + Ksat_CH_up(j)) 
            uptakeAA(j) = Vmax_AA(j) * AA_m / (AA_m + Ksat_AA_up(j)) 
            uptakeNUC(j) = Vmax_NUC(j) 
     &                      * NUC_m / (NUC_m + Ksat_NUC_up(j)) 
            uptakePR(j) = Vmax_PR(j) * PR_m / (PR_m + Ksat_PR_up(j)) 
            uptakePIG(j) = Vmax_PIG(j) 
     &                      * PIG_m / (PIG_m + Ksat_PIG_up(j)) 
          end do

                 
c.......................................................................
  920   format('nt ',i6,' j ',i6,'  upC ',e10.3,'  UpN ',e10.3,
     & '  upP ',e10.3,' upCH ',e10.3,'  upAA ',e10.3, '  upNUC ',e10.3,
     & '  upPR ',e10.3,' upPIG ',e10.3)
c       do j=1,jmax
c          write(6,920)nt,j,uptakec(j),uptaken(j),uptakep(j),
c    &          uptakeCH(j),uptakeAA(j),uptakeNUC(j),
c    &          uptakePR(j),uptakePIG(j) 
c       end do
c.......................................................................

c Photosynthesis
c specify photosynthesis as a function of pigment and light
c C:N chlorophyll ~ 14:1
c Chl-a has 55 carbons per molecule
c Chl accounts for between 1-7% of cell carbon (not carefully researched?)
c Cullen and McKintyre --- culturing article --- PmaxChl = 2 - 8  g C (g Chl)-1 hr-1
c 1 g C (g Chl)-1 hr-1 = 2.1d-2 mol C (mol Chl)-1 s-1
c so, fixing photosynth rate, Pphoto =  8d-2 mol C (mol Chl-1) s-1 
c photosynth = Pphoto * PIG(j)   mol CH m-3 s-1
c         Pphoto = 2.1d-2
          Pphoto = 1.0d-2
          do j = 1,jmax
c if heterotroph (itrophic = 0) then photosynthesis = 0
            if(itrophic(j) .eq. 0)then
               photosynth(j) = 0.0d0
            else
c if QCH exceeds max then photosynth shut down
               if(QCH(j) .ge. QCH_max(j))then 
                 photosynth(j) = 0.0d0
               else
                 Qc(j) = c(j)/X(j)
                 dummy = Qc(j)/(Qc(j) + Qc_sat_photo(j))
c                photosynth(j) = Pphoto*PIG(j)*dummy
c                photosynth(j) = Pphoto*PR(j)*0.1*dummy
                 pfactor = pfactor_init
                 photosynth(j) = Pphoto*PR(j)*pfactor*dummy
cc                photosynth(j) = Pphoto*PR(j)*0.005*dummy
               endif
c.........................................................
c              if(nt .eq. 1)then
c                write(6,*)'!! WARNING: photosyth proportion to PR ...'
c                write(6,*)'!! WARNING: photosyth proportion to PR ...'
c              endif
c.........................................................
c              photosynth(j) = 1.0d-19 * X(j) 
c              photosynth(j) = 1.0d-15 
c TESTING ................
c            write(6,*)'Photosynth 1',photosynth(j)
c            write(6,*)'Uptake_c*X',(uptakec(j)*X(j))
c            write(6,*)'Uptake_n*X',(uptaken(j)*X(j))
c
c            photosynth(j) = 1.0d-18*X(j) 
c            write(6,*)'Photosynth 2',photosynth(j)
c TESTING..................................
c         write(6,*)'j,X,photosynth ',j,X(j),photosynth(j)
c         write(6,*)'Qc,  Qc_sat_photo, dummy ',Qc(j), Qc_sat_photo(j),
c    &                dummy
c         write(6,*)'Qn_sat_AAsynth ',Qn_sat_AAsynth(j)
c............................................
            end if
          end do


c exudation
          do j = 1,jmax
             exude_CH(j)  = exCH(j) * CH(j)
             exude_AA(j)  = exAA(j) * AA(j)
             exude_NUC(j) = exNUC(j) * NUC(j)
             exude_PR(j)  = exPR(j) * PR(j)
             exude_PIG(j) = exPIG(j) * PIG(j)
          end do

c internal cellular reactions depend on cell quotas
c  --- equivalent to using internal concentrations.
c  --- NOT environmental concentrations!
          do j = 1,jmax
c amino acid synthesis
            dummy1 = QCH(j) / (QCH_sat_AAsynth(j) +  QCH(j))
            dummy2 = Qn(j) / (Qn_sat_AAsynth(j) +  Qn(j))
c double check in case slight negative...
            if(Qn(j) .le. 0.0d0)dummy2 = 0.0d0
            if(dummy1 .le. dummy2)then 
              dummy = dummy1
            else 
              dummy = dummy2
            end if
            AAsynth(j) = X(j) * KAAsynth(j) * dummy 
c TEST .................................................
c           write(6,*)'Qn  ',Qn(j),'  AAsynth  ',AAsynth(j)
c TEST .................................................
c nucleic acid synthesis
            dummy1 = QAA(j) / (QAA_sat_NUCsynth(j) +  QAA(j))
            dummy2 = Qp(j) / (Qp_sat_NUCsynth(j) +  Qp(j))
c double check in case slight negative...
            if(Qp(j) .le. 0.0d0)dummy2 = 0.0d0
            dummy3 = QCH(j) / (QCH_sat_NUCsynth(j) +  QCH(j))
            if(dummy1 .le. dummy2)then 
              dummy = dummy1
            else 
              dummy = dummy2
            end if
            if(dummy .gt. dummy3)dummy = dummy3
            NUCsynth(j) = X(j) * KNUCsynth(j) * dummy 
c protein synthesis
            dummy1 = QAA(j) / (QAA_sat_PRsynth(j) +  QAA(j))
            dummy2 = QCH(j) / (QCH_sat_PRsynth(j) +  QCH(j))
            if(dummy1 .le. dummy2)then 
              dummy = dummy1
            else 
              dummy = dummy2
            end if
            PRsynth(j) = X(j) * GammaPR_synth(j) * QNUC(j) * dummy

c pigment synthesis: CH and N limited
c  KPIGsynth(j) should be function of light
            if(itrophic(j) .ne. 0)then
              dummy1 = QCH(j) / (QCH_sat_PIGsynth(j) +  QCH(j))
              dummy2 = Qn(j) / (Qn_sat_PIGsynth(j) +  Qn(j))
              if(Qn(j) .le. 0.0d0)dummy2 = 0.0d0
              if(dummy1 .le. dummy2)then 
                dummy = dummy1
              else 
                dummy = dummy2
              end if
c             PIGsynth(j) = X(j) * KPIGsynth(j) * dummy 
c SIMPLIFIED ALTERNATIVE: Set pigment synthesis as fixed fraction of
c protein synthesis (i.e. approx fixed Chl:C ratio)
c             PIGsynth(j) = X(j) * PRsynth(j) * pigfrac
              PIGsynth(j) = 0.0d0  
            else
              PIGsynth(j) = 0.0d0  
            endif
          end do

c Belt and braces: re-initialize rates of change 
          do j=1,jmax
             dcdt(j) = 0.0d0
             dndt(j) = 0.0d0
             dpdt(j) = 0.0d0
             dCHdt(j) = 0.0d0
             dAAdt(j) = 0.0d0
             dNUCdt(j) = 0.0d0
             dPRdt(j) = 0.0d0
             dPIGdt(j) = 0.0d0
          end do
          dcmdt = 0.0d0
          dnmdt = 0.0d0
          dpmdt = 0.0d0
          dCHmdt = 0.0d0
          dAAmdt = 0.0d0
          dNUCmdt = 0.0d0
          dPRmdt = 0.0d0
          dPIGmdt = 0.0d0

c rates of change from BIOLOGICAL processes only ....
          do j=1,jmax

c cellular organic stores 
c carbon
             dcdt(j) = uptakec(j)*X(j)  - photosynth(j) 
c nitrogen
             dndt(j) = uptaken(j)*X(j)  - AAsynth(j)*Rnc_AA(j) 
     &            - PIGsynth(j)*Rnc_PIG(j) 

c......... TESTING ............................................
c          write(6,*)'dndt ',dndt, 'uptaken ',uptaken(j),
c    &       ' AAsynth ',AAsynth(j),'PIGsynth ',PIGsynth(j)
c......... TESTING ............................................

c phosphorus
             dpdt(j) = uptakep(j)*X(j)  - NUCsynth(j)*Rpc_NUC(j) 

c evaluate total CH consumption by respiration
             total_respiration(j) = 
c          ... and consumption for respiration/energy source
     &                  uptakec(j)*epsilon_uptakec(j) 
     &                + uptaken(j)*epsilon_uptaken(j) 
     &                + uptakep(j)*epsilon_uptakep(j) 
     &                + uptakeCH(j)*epsilon_uptakeCH(j) 
     &                + uptakeAA(j)*epsilon_uptakeAA(j) 
     &                + AAsynth(j)*epsilon_AAsynth(j)
     &                + NUCsynth(j)*epsilon_NUCsynth(j)
c epsilon_PRsynth(j) significant for total yield factor (major contribution)
     &                + PRsynth(j)*epsilon_PRsynth(j)
     &                + PIGsynth(j)*epsilon_PIGsynth(j)
c          ... "maintenance energy" 
     &                + Bc(j)*epsilon_maintenance(j)

c.....................................................
c           write(6,*)'test1 j ',j,'  total resp ', total_respiration(j)
c.....................................................

c cellular organic stores - carbon currency
c carbohydrate / lipid
             dCHdt(j) = 
c sources due to uptake and photosynth
     &                 uptakech(j)*X(j) + photosynth(j)  
c    losses ... direct consumption through synthesis
     &               - AAsynth(j) - PIGsynth(j) 
c          ... total respiration 
     &               - total_respiration(j)
c          ... and exudation
     &               - exude_CH(j)

c amino acids / nucleic acid precursors
             dAAdt(j) = uptakeAA(j)*X(j) + AAsynth(j) - NUCsynth(j) 
     &            - PRsynth(j) - exude_AA(j)

c          write(6,*)'uptakeAA   AAsynth   NUCsynth  PRsynth  exude_AA'
c          write(6,*)uptakeAA(j), AAsynth(j), NUCsynth(j), PRsynth(j),
c    &               exude_AA(j) 

c nucleic acids
             dNUCdt(j) = uptakeNUC(j)*X(j) + NUCsynth(j) - exude_NUC(j) 

c TEST .....................................................
c          write(6,*)'uptakeNUC   NUCsynth  exude_NUC  dNUCdt'
c          write(6,*)uptakeNUC(j),NUCsynth(j),exude_NUC(j),dNUCdt(j) 
c TEST .....................................................

c protein
             dPRdt(j) =  uptakePR(j)*X(j) + PRsynth(j) - exude_PR(j) 
c pigment
             dPIGdt(j) =  uptakePIG(j)*X(j) + PIGsynth(j) - exude_PIG(j) 
          end do

c DILUTION OF CELL PROPERTIES
c --- assumes no cells coming into chamber
          do j=1,jmax
             dcdt(j) = dcdt(j) - D*c(j) 
             dndt(j) = dndt(j) - D*n(j)
             dpdt(j) = dpdt(j) - D*p(j)
             dCHdt(j) = dCHdt(j) - D*CH(j)
             dAAdt(j) = dAAdt(j) - D*AA(j)
             dNUCdt(j) = dNUCdt(j) - D*NUC(j)
             dPRdt(j) =  dPRdt(j)- D*PR(j)
             dPIGdt(j) =  dPIGdt(j)- D*PIG(j)
           end do
c RATES OF CHANGE OF MEDIUM PROPERTIES
c NOTE: uptake SUBTRACTED from medium
c --- first sum all respiration (source of inorg C)
c       ... add respired carbon directly to external inorg pool
c --- and uptake of inorg carbon 

c Currently all respired carbon is immediately lost to external dissolved 
c inorganic pool, regardless of trophic strategy
           sum_respiration = 0.0d0
           sum_uptake_c = 0.0
           do j = 1,jmax
             sum_respiration = sum_respiration + total_respiration(j)
             sum_uptake_c = sum_uptake_c + uptakec(j)*X(j)
           end do
           dcmdt = dcmdt - sum_uptake_c + sum_respiration 
     &                 +  D * (c_in_m - c_m) 
c --- inorg nitrogen
           sum_uptake_n = 0.0
           do j = 1,jmax
             sum_uptake_n = sum_uptake_n + uptaken(j)*X(j) 
           end do
           dnmdt = dnmdt - sum_uptake_n + D * (n_in_m - n_m) 
c --- inorg phosphorus
           sum_uptake_p = 0.0
           do j = 1,jmax
             sum_uptake_p = sum_uptake_p + uptakep(j)*X(j)
           end do
           dpmdt = dpmdt - sum_uptake_p + D * (p_in_m - p_m) 
c --- organic carbon
           sum_uptake_CH = 0.0
           sum_exude_CH = 0.0
           do j = 1,jmax
             sum_uptake_CH = sum_uptake_CH + uptakeCH(j)*X(j)
             sum_exude_CH = sum_exude_CH + exude_CH(j)
           end do
           dCHmdt = dCHmdt - sum_uptake_CH + D * (CH_in_m - CH_m) 
     &           + sum_exude_CH
c --- Amino acids
           sum_uptake_AA = 0.0
           sum_exude_AA = 0.0
           do j = 1,jmax
             sum_uptake_AA = sum_uptake_AA + uptakeAA(j)*X(j)
             sum_exude_AA = sum_exude_AA + exude_AA(j) 
           end do
           dAAmdt = dAAmdt - sum_uptake_AA + D * (AA_in_m - AA_m) 
     &           + sum_exude_AA
c --- Nucleic acids
           sum_uptake_NUC = 0.0
           sum_exude_NUC = 0.0
           do j = 1,jmax
             sum_uptake_NUC = sum_uptake_NUC + uptakeNUC(j)*X(j)
             sum_exude_NUC = sum_exude_NUC + exude_NUC(j) 
           end do
           dNUCmdt = dNUCmdt - sum_uptake_NUC + D * (NUC_in_m - NUC_m) 
     &           + sum_exude_NUC
c --- protein
           sum_uptake_PR = 0.0
           sum_exude_PR = 0.0
           do j = 1,jmax
             sum_uptake_PR = sum_uptake_PR + uptakePR(j)*X(j)
             sum_exude_PR = sum_exude_PR + exude_PR(j) 
           end do
           dPRmdt = dPRmdt - sum_uptake_PR + D * (PR_in_m - PR_m) 
     &           + sum_exude_PR
c --- pigment 
           sum_uptake_PIG = 0.0
           sum_exude_PIG = 0.0
           do j = 1,jmax
             sum_uptake_PIG = sum_uptake_PIG + uptakePIG(j)*X(j)
             sum_exude_PIG = sum_exude_PIG + exude_PIG(j) 
           end do
           dPIGmdt = dPIGmdt - sum_uptake_PIG + D * (PIG_in_m - PIG_m) 
     &           + sum_exude_PIG

c ..................................................................
c check for total mass conservation
c Carbon
          TotalC = 0.0d0
          TotalN = 0.0d0
          TotalP = 0.0d0
          do j=1,jmax
           TotalC = TotalC + c(j) + CH(j) + AA(j) + NUC(j) + PR(j) 
     &                    + PIG(j) 
           Bn(j) = n(j) + AA(j)*Rnc_AA(j) + NUC(j)*Rnc_NUC(j) 
     &             + PR(j)*Rnc_PR(j) + PIG(j)*Rnc_PIG(j)
           Bp(j) = p(j) + NUC(j)*Rpc_NUC(j)
           TotalN = TotalN + Bn(j)
           TotalP = TotalP + Bp(j)
          end do
c ... add medium concentrations
          TotalC = TotalC + c_m + CH_m + AA_m + NUC_m + PR_m + PIG_m
          TotalN = TotalN + n_m + AA_m*Rnc_AA(1) + NUC_m*Rnc_NUC(1) 
     &           + PR_m*Rnc_PR(1) + PIG_m*Rnc_PIG(1)
          TotalP = TotalP + p_m + NUC_m*Rpc_NUC(1)
c ..................................................................

c step forward time
c --- first, cell components
c --- predictor-corrector step
c --- ab = Adams-Bashforth weighting
          ab = 0.7d0
          do j=1,jmax
            c(j) = c(j) + ( dcdt(j)*ab + dcdt_o(j)*(1.0d0-ab) ) * deltat
            n(j) = n(j) + ( dndt(j)*ab + dndt_o(j)*(1.0d0-ab) ) * deltat
            p(j) = p(j) + ( dpdt(j)*ab + dpdt_o(j)*(1.0d0-ab) ) * deltat
            CH(j) = CH(j) + (dCHdt(j)*ab + dCHdt_o(j)*(1.0d0-ab))*deltat
            AA(j) = AA(j) + (dAAdt(j)*ab + dAAdt_o(j)*(1.0d0-ab))*deltat
            NUC(j) = NUC(j) + 
     &              (dNUCdt(j)*ab + dNUCdt_o(j)*(1.0d0-ab))*deltat 
            PR(j) = PR(j) + (dPRdt(j)*ab + dPRdt_o(j)*(1.0d0-ab))*deltat
            PIG(j) = PIG(j) +
     &              (dPIGdt(j)*ab + dPIGdt_o(j)*(1.0d0-ab))*deltat 

c update total carbon biomass (mol C m-3 in medium) and number density
            Bc(j) = c(j) + CH(j) + AA(j) + NUC(j) + PR(j) + PIG(j)
            X(j) = Bc(j) / QBc_ave(j) 
C TESTING ...................
c         if(n(j) .lt. 1.0d-50)n(j) = 1.0d-50
c         if(p(j) .lt. 1.0d-50)p(j) = 1.0d-50
C TESTING ...................
          end do

c --- second, medium concentrations
          c_m = c_m + (dcmdt*ab + dcmdt_o*(1.0d0-ab))*deltat 
          n_m = n_m + (dnmdt*ab + dnmdt_o*(1.0d0-ab))*deltat 
          p_m = p_m + (dpmdt*ab + dpmdt_o*(1.0d0-ab))*deltat 
          CH_m = CH_m + (dCHmdt*ab + dCHmdt_o*(1.0d0-ab))*deltat 

C TESTING ...................
c       if(c_m .lt. 1.0d-80)c_m = 1.0d-80
c       if(n_m .lt. 1.0d-80)n_m = 1.0d-80
c       if(p_m .lt. 1.0d-80)p_m = 1.0d-80
c       if(CH_m .lt. 1.0d-80)CH_m = 1.0d-80
C TESTING ...................

          AA_m = AA_m + (dAAmdt*ab + dAAmdt_o*(1.0d0-ab))*deltat 
          NUC_m = NUC_m + (dNUCmdt*ab + dNUCmdt_o*(1.0d0-ab))*deltat 
          PR_m = PR_m + (dPRmdt*ab + dPRmdt_o*(1.0d0-ab))*deltat 
          PIG_m = PIG_m + (dPIGmdt*ab + dPIGmdt_o*(1.0d0-ab))*deltat 

c ...............................................................
c check for total mass conservation
c Carbon
          TotalC = 0.0d0
          TotalN = 0.0d0
          TotalP = 0.0d0
          do j=1,jmax
           TotalC = TotalC + c(j) + CH(j) + AA(j) + NUC(j) + PR(j) 
     &                    + PIG(j) 
           Bn(j) = n(j) + AA(j)*Rnc_AA(j) + NUC(j)*Rnc_NUC(j) 
     &             + PR(j)*Rnc_PR(j) + PIG(j)*Rnc_PIG(j)
           Bp(j) = p(j) + NUC(j)*Rpc_NUC(j)
           TotalN = TotalN + Bn(j)
           TotalP = TotalP + Bp(j)
          end do
c ... add medium concentrations
          TotalC = TotalC + c_m + CH_m + AA_m + NUC_m + PR_m + PIG_m
          TotalN = TotalN + n_m + AA_m*Rnc_AA(1) + NUC_m*Rnc_NUC(1) 
     &           + PR_m*Rnc_PR(1) + PIG_m*Rnc_PIG(1)
          TotalP = TotalP + p_m + NUC_m*Rpc_NUC(1)
c .................................................................


c --- reset Adams-Bashforth time arrays 
          do j=1,jmax
            dcdt_o(j) =  dcdt(j)
            dndt_o(j) =  dndt(j)
            dpdt_o(j) =  dpdt(j)
            dCHdt_o(j) =  dCHdt(j)
            dAAdt_o(j) =  dAAdt(j)
            dNUCdt_o(j) =  dNUCdt(j)
            dPRdt_o(j) =  dPRdt(j)
            dPIGdt_o(j) =  dPIGdt(j)
          end do
          dcmdt_o = dcmdt 
          dnmdt_o = dnmdt 
          dpmdt_o = dpmdt 
          dCHmdt_o = dCHmdt 
          dAAmdt_o = dAAmdt 
          dNUCmdt_o = dNUCmdt 
          dPRmdt_o = dPRmdt 
          dPIGmdt_o = dPIGmdt 

c write some output
  930     format(f10.2,i6,23e12.4) 
  900     format('nt ',i6,' j ',i6,'  X ',e10.3,'  Bc ',e10.3,
     & '  N_m ',e10.3,' CH ',e10.3,'  AA ',e10.3, '  NUC ',e10.3,
     & '  PR ',e10.3,' PIG ',e10.3)

          if(nt .eq. 1)then
c          write(6,*)'time j  X_out  Bc_out c_out  n_out  p_out  
c    &  CH_out AA_out  NUC_out  PR_out  PIG_out  c_m_out  n_m_out  
c    &  p_m_out  CH_m_out  AA_m_out  NUC_m_out  PR_m_out  PIG_m_out 
c    &  TotalC_out TotalN_out TotalP_out Bn_out Bp_out ' 
c write header to output file
           write(12,*)' time j  X  Bc  c  n  p  CH  AA  NUC  PR  PIG 
     & c_m  n_m  p_m  CH_m  AA_m  NUC_m  PR_m  PIG_m 
     & TotalC TotalN TotalP  Bn_out  Bp_out' 
          end if

          do j=1,jmax
c  X(j) in cells m-3, rewrite as cells ml-1
c  concentrations in moles m-3 medium, rewrite as mol l-1
          X_out = X(j) * 1.0d-6
          Bc_out = Bc(j) * 1.0d3 
          Bn_out = Bn(j) * 1.0d3 
          Bp_out = Bp(j) * 1.0d3 
          TotalC_out = TotalC * 1.0d3 
          TotalN_out = TotalN * 1.0d3 
          TotalP_out = TotalP * 1.0d3 
          c_m_out = c_m * 1.0d3
          if(c_m_out .gt. 0.0d0 .and. c_m_out 
     & .lt. 1.0d-50)c_m_out=1.0d-50
          n_m_out = n_m * 1.0d3
          if(n_m_out .gt. 0.0d0 .and. n_m_out 
     & .lt. 1.0d-50)n_m_out=1.0d-50
          p_m_out = p_m * 1.0d3
          if(p_m_out .gt. 0.0d0 .and. p_m_out 
     & .lt. 1.0d-50)p_m_out=1.0d-50
          CH_m_out = CH_m * 1.0d3
          if(CH_m_out .gt. 0.0d0 .and. CH_m_out 
     & .lt. 1.0d-50)CH_m_out=1.0d-50
          AA_m_out = AA_m * 1.0d3
          PR_m_out = PR_m * 1.0d3
          PIG_m_out = PIG_m * 1.0d3
          NUC_m_out = NUC_m * 1.0d3
          if(NUC_m_out .gt. 0.0d0 .and. NUC_m_out 
     & .lt. 1.0d-50)NUC_m_out=1.0d-50
          c_out = c(j) * 1.0d3
          if(c_out .gt. 0.0d0 .and. c_out 
     & .lt. 1.0d-50)c_out=1.0d-50
          n_out = n(j) * 1.0d3
          if(n_out .gt. 0.0d0 .and. n_out 
     & .lt. 1.0d-50)n_out=1.0d-50
          p_out = p(j) * 1.0d3
          if(p_out .gt. 0.0d0 .and. p_out 
     & .lt. 1.0d-50)p_out=1.0d-50
          CH_out = CH(j) * 1.0d3 
          if(CH_out .gt. 0.0d0 .and. CH_out 
     & .lt. 1.0d-50)CH_out=1.0d-50
          AA_out = AA(j) * 1.0d3 
          if(AA_out .gt. 0.0d0 .and. AA_out 
     & .lt. 1.0d-50)AA_out=1.0d-50
          NUC_out = NUC(j) * 1.0d3 
          PR_out = PR(j) * 1.0d3 
          PIG_out = PIG(j) * 1.0d3 
c          write(6,900)nt,j,X_out,Bc_out,n_m_out,CH_out,AA_out,
c    &            NUC_out, PR_out, PIG_out, TotalC_out

            if(mod(nt, nout) .eq. 0)then
c             write(6,930)time_out,j,X_out,Bc_out,
c    &            c_out, n_out, p_out, CH_out, AA_out, 
c    &            NUC_out, PR_out, PIG_out,
c    &            c_m_out, n_m_out, p_m_out, CH_m_out, AA_m_out, 
c    &            NUC_m_out, PR_m_out, PIG_m_out, 
c    &            TotalC_out,TotalN_out, TotalP_out,
c    &            Bn_out, Bp_out
              write(12,930)time_out,j,X_out,Bc_out,
     &            c_out, n_out, p_out, CH_out, AA_out, 
     &            NUC_out, PR_out, PIG_out,
     &            c_m_out, n_m_out, p_m_out, CH_m_out, AA_m_out, 
     &            NUC_m_out, PR_m_out, PIG_m_out,
     &            TotalC_out,TotalN_out, TotalP_out,
     &            Bn_out, Bp_out
c if appropriate time, then put diagnostics into array 
c for evaluation of cost function
                  do ic = 1,icmax
                    timetest = (time_out - labt(ic))/labt(ic)
                    if(timetest .le. 0.01 .and. timetest .ge. -0.01)then
                      modt(ic) = time_out
                      modNm(ic) = n_m_out
                      modX(ic) = X_out
                      modCN(ic) = Bc_out/Bn_out
                    end if
                   end do

            end if
          end do

 300     end do 
c END OF TIMESTEPPING LOOP
         close(12)
       
c write out some data from the last  timestep assuming equilibrium
c write header only once
         if(idone .eq. 0)then
c        write(6,*)'c_in_m  n_in_m  p_in_m  CH_in_m  AA_in_m  NUC_in_m  
c    & PR_in_m  PIG_in_m  epsilon_PR  D  total_resp  Bc  Bn  Bp'
         write(14,*)'c_in_m  n_in_m  p_in_m  CH_in_m  AA_in_m  NUC_in_m 
     & PR_in_m  PIG_in_m  epsilon_PR  D  total_resp  Bc  Bn  Bp'
c reset idone so header isnt written again
         idone = idone + 1
         end if

c write results
 950     format(14e12.4) 
c        write(6,950)c_in_m, n_in_m, p_in_m, CH_in_m, AA_in_m, NUC_in_m,
c    & PR_in_m, PIG_in_m, epsilon_PRsynth, 
c    & D, total_respiration(1), Bc(1), Bn(1), Bp(1)
         write(14,950)c_in_m, n_in_m, p_in_m, CH_in_m,AA_in_m,NUC_in_m,
     & PR_in_m, PIG_in_m, epsilon_PRsynth, 
     & D, total_respiration(1), Bc(1), Bn(1), Bp(1)

 320   end do
c close compilation results file
c MICK - SHOULD THIS BE close(14)????
c 12 is closed already above????
       close(12)

c EVALUATE COST FUNCTION ..........
c first quick check
c      write(6,*)'ic,     lab time,    model time'
c      do ic = 1,icmax 
c        write(6,*)ic,labt(ic), modt(ic) 
c      end do
c      write(6,*)'ic,     lab Nm,     model Nm'
c      do ic = 1,icmax 
c        write(6,*)ic,labNm(ic), modNm(ic) 
c      end do
c      write(6,*)'ic,     lab X,      model X'
c      do ic = 1,icmax 
c        write(6,*)ic, labX(ic), modX(ic) 
c      end do
c      write(6,*)'ic,    lab CN,     model CN'
c      do ic = 1,icmax 
c        write(6,*)ic, labCN(ic), modCN(ic) 
c      end do

c now evaluate costfunction
c first sum squares of normalized model/data differences
       cost = 0.0d0 
c      write(6,*)'COST = ',cost
       do 399 ic = 1, icmax
          if(labNm(ic) .ne. -999.999 .and. labNm(ic) .gt. 0.15)then
            cdummy = (labNm(ic) - modNm(ic)) / labNm(ic)
c           write(6,*)'ic = ',ic,'labNm = ',
c    &                labNm(ic),'modNm = ',modNm(ic)
            cdummy1 = cdummy*cdummy
            cost = cost + cdummy1
          end if
c         write(6,*)'COST1 = ',cost
          if(modX(ic) .ne. 0.0000)then
            cdummy = (labX(ic) - modX(ic)) / labX(ic)
            cdummy1 = cdummy*cdummy
            cost = cost + cdummy1
          end if
c         write(6,*)'COST2 = ',cost
          if(labCN(ic) .ne. -999.999 .and. modCN(ic) .ne. 0.0000)then
            cdummy = (labCN(ic) - modCN(ic)) / labCN(ic)
            cdummy1 = cdummy*cdummy
            cost = cost + cdummy1
          end if
 399   end do     


C WRITE OUT DETAILS...
       if(ione .eq. 1 .and. itwo .eq. 1 .and. ithree .eq. 1)then
         write(6,*)'xinit  QCH_max_factor   GammaPR_synth  pfac VmaxFac
     & SizeFac  cost'
         write(27,*)'xinit  QCH_max_factor   GammaPR_synth  pfac VmaxFac
     & SizeFac  cost'
       endif
       write(6,750)xinit,QCH_max_factor_init,GammaPR_synth_init,
     &           pfactor_init,KAAfactor,VmaxFactor,SizeFactor,cost
       write(27,750)xinit,QCH_max_factor_init,GammaPR_synth_init,
     &           pfactor_init,KAAfactor,VmaxFactor,SizeFactor,cost
 750   format(4e12.3,4f12.4)
c END VARIABLE LOOPS....

       if(cost .lt. mincost)mincost = cost
       write(6,*)'COST3 = ',cost,'mincost = ',mincost

c706   end do
c705   end do
c704   end do
c703   end do
c702   end do
c701   end do
c700   end do

c close output file
       close(27)

       write(6,*)'MINIMUM COST = ',mincost

c end of program
       end

C ===================================================================
C ===================================================================

