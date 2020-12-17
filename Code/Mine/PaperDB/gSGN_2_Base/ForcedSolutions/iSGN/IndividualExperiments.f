      program main
   
      implicit none
   
      Integer wdirlen, xbc_len,NormFileI,tlist_len
      PARAMETER(wdirlen= 300,xbc_len=5012,NormFileI = 10,tlist_len=3)
      
      CHARACTER(len =wdirlen) wdir
  
      integer x_len,n_GhstCells
      DOUBLE PRECISION tlist(tlist_len)
     
      DOUBLE PRECISION xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,
     . ga,a0,a1,a2,a3,a4,a5,b1max,b1min,phi0,phi1,
     . b1a6,b1a7,b2max,b2min,b2a6,b2a7,alpha,dtodx
     
      INTEGER effeclenwdir
      
      wdir = "/home/jp/Documents/" // 
     . "Work/PostDoc/Projects/Steve/1DWaves/" //
     . "RegularisedSerre/Data/RAW" //
     . "/Models/gSGNForcedLimAll/ForcedSolution" //
     .  "/exp1/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
            
      n_GhstCells = 6
      x_len = xbc_len - 2 *n_GhstCells
      
      !SWWE equations
      
      ga = 9.81
            
      xstart = -50d0
      xend = 100.0d0
      
      dx = (xend - xstart) / (x_len -1) 
      
      tstart = 0d0
      tend = 10.0d0
            
      ga = 9.81d0
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      
      tlist(1) = 2.5d0
      tlist(2) = 5d0
      tlist(3) = 7.5d0
      
      theta = 1.2
      
      
      ! betamax,betamin,a5,a6,a7 control dispersion properties, size of terms
      b1max = 2d0
      b1min = -2d0/3d0
      a5 = -0.9
      phi0 = xstart - a5*tstart
      phi1 = xend - a5*tend 
      
      b1a6 = (b1max - b1min) / (phi1 - phi0)
      b1a7 = b1min - b1a6*phi0
      
      b2max = 3
      b2min = 0
      
      b2a6 = (b2max - b2min) / (phi1 - phi0)
      b2a7 = b2min - b2a6*phi0
      
      alpha = max(1d0,b2max / (2d0/3d0 + b1max))
      maxwavespeed = a2 + a4 + dsqrt(alpha*ga*(a0 + a1))
      Cr = 0.5
      dtodx = Cr/maxwavespeed
      dt = dx*dtodx 

      open(NormFileI, file = wdir(1:effeclenwdir)//'Norms.dat') 
 
      
      call NumericalSolveForced(xbc_len,n_GhstCells,
     . tstart,tend,xstart,dx,dt,ga,theta, 
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . tlist,tlist_len,wdir(1:effeclenwdir),effeclenwdir,
     . NormFileI)
     
      close(NormFileI) 
      
      end 
