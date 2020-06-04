      program main
   
      implicit none
   
      Integer wdirlen, xbc_len
      PARAMETER(wdirlen= 300,xbc_len=10012)
      
      CHARACTER(len =wdirlen) wdir
  
      integer i,x_len,n_GhstCells
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),Gbc_init(xbc_len),
     . ubc_init(xbc_len),hbc_fin(xbc_len),Gbc_fin(xbc_len),
     . ubc_fin(xbc_len),hbc_fin_a(xbc_len),Gbc_fin_a(xbc_len),
     . ubc_fin_a(xbc_len)
     
      DOUBLE PRECISION hl,hr,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,currenttime,
     . dbalpha
     
      INTEGER effeclenwdir,dtsep
      
      wdir = "/home/jp/Documents/" // 
     . "Work/PostDoc/Projects/Steve/1DWaves/" //
     . "RegularisedSerre/Data/RAW" //
     . "/Models/gSGNForcedLimhG/Beta1Beta2LimitOpt1/"//
     . "SmoothDB/2to1/alpha0p1/timeseries/exp1/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(1, file = wdir(1:effeclenwdir)//'InitVal.dat') 
      open(2, file = wdir(1:effeclenwdir)//'EndVals.dat') 
      open(3, file = wdir(1:effeclenwdir)//'EndAnaVals.dat') 
      open(4, file = wdir(1:effeclenwdir)//'Params.dat') 

      
      n_GhstCells = 6
      x_len = xbc_len - 2 *n_GhstCells
      
      !SWWE equations
      
      ga = 9.81
         
      hl = 2d0
      hr = 1d0
      dbalpha = 0.1d0
      
      xstart = -100d0
      xend = 100d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 20d0
      
      dtsep = 50
      
      
      
      !ensures that x(0) = xstart and x(x_len) = x_end
      dx = (xend - xstart) / (x_len -1) 
      
      Cr = 0.5
      maxwavespeed = dsqrt(ga*(hl))
      dt  = (Cr / maxwavespeed) *dx
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call SmoothDB(xbc,xbc_len,
     . hl,hr,dbalpha,hbc_init,ubc_init,Gbc_init) 
      
      !solve gSGN with beta values until currenttime > tend
      call NumericalSolveTSPrint(tstart,tend,
     . ga,theta,dx,dt,n_GhstCells,xbc,xbc_len,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . dtsep,wdir(1:effeclenwdir),effeclenwdir)

     
      ! get analytic values of h,u,G
      call Dambreak(xbc,xbc_len,currenttime,ga,hl,hr,
     . hbc_fin_a,ubc_fin_a,Gbc_fin_a) 


      !write out initial values and end values
      do i = 1,xbc_len
         write(1,*) xbc(i),hbc_init(i),Gbc_init(i),ubc_init(i)
         write(2,*) xbc(i),hbc_fin(i),Gbc_fin(i),ubc_fin(i)
         write(3,*) xbc(i),hbc_fin_a(i),Gbc_fin_a(i),ubc_fin_a(i)
      end do
      
      !write out parameters      write(5,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(4,*) 'xstart :',xstart
      write(4,*) 'xend :',xend
      write(4,*) 'x_len :',x_len
      write(4,*) 'dx = (x_end - x_start) / (x_len -1)  :' , dx
      write(4,*) 'n_GhstCells :',n_GhstCells
      write(4,*) 'tstart :', tstart
      write(4,*) 'tend :',tend 
      write(4,*) 'actual end time :', currenttime
      write(4,*) 'dt/dx :' , (Cr / maxwavespeed)
      write(4,*) 'dt = dx*(dt/dx)  :' , dt
      write(4,*) 'gravity :' , ga
      write(4,*) 'hl :' , hl
      write(4,*) 'hr :' , hr
  
      
      
      close(1)
      close(2)
      close(3)
      close(4)

      
      end 
