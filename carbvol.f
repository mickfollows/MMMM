C ===================================================================
C ===================================================================

c ......................................................................
c cell volume / carbon relationships ...................................
c Mick Follows, Jason Bragg: June 2009 .................................
c ......................................................................
c cellvol = cell volume (m-3)
c empirical total cell carbon to cell volume from Lovdal et al. 
c Ctot = a*(V)^b     Ctot cell-1 (fmol C)   V (um^3)    
c a=18.7 fmol C  b=0.89 
c ......................................................................
       subroutine carbon_from_vol(cell_vol,total_cell_carb)
c total cell carbon from cell volume
       implicit none
       real*8 total_cell_carb
       real*8 cell_vol
       real*8 a, b
       a=18.7
       b=0.89
c testing ...........
c      write(6,*)'cell_vol',cell_vol
c      write(6,*)'cell_vol**b',cell_vol**b
c testing ...........

       total_cell_carb =  a * (cell_vol**b)
       return
       end
c .....................................................................
       subroutine vol_from_carbon(total_cell_carb,cell_vol)
c cell volume from total cell carbon
       implicit none
       real*8 total_cell_carb
       real*8 cell_vol
       real*8 a, b, bm
       a=18.7
       b=0.89
       bm = 1.0d0/b
       cell_vol =  (total_cell_carb/a)**bm
       return 
       end


C ===================================================================
C ===================================================================

