c =========================
c DB single experiment
c
c ====================   
      subroutine SingleGaussian(a0,a1,a2,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,EnergFile,ExpWdir,ExpWdir_len,dtsep)
     
      DOUBLE PRECISION a0,a1,a2,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
     
      integer xbc_len,n_GhstCells,EnergFile,ExpWdir_len,dtsep
      
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len), ubc_init(xbc_len),  hbc_fin(xbc_len),
     . Gbc_fin(xbc_len), ubc_fin(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      integer i
      DOUBLE PRECISION currenttime
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call Gaussian(xbc,xbc_len,
     . a0,a1,a2,hbc_init,ubc_init,Gbc_init) 
      
      !solve gSGN with beta values until currenttime > tend    
      call NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc,xbc_len,
     . hbc_init,Gbc_init,ubc_init,
     . Energs_init,currenttime,hbc_fin,Gbc_fin,ubc_fin,Energs_fin,
     . dtsep,ExpWdir,ExpWdir_len)
     
c Conservation Norm Tests  

      do i = 1,4
         !if denominator small just return absolute error
         if (dabs(Energs_init(i)) .lt. 10d0**(-10)) then
            Energ_Err(i) = dabs(Energs_fin(i) - Energs_init(i))
         else
            Energ_Err(i) = dabs(Energs_fin(i) - Energs_init(i))/
     .       dabs(Energs_init(i))
         end if
      end do   

           
      ! write out individual experiment details
      open(1, file = ExpWdir//'InitVal.dat') 
      open(2, file = ExpWdir//'EndVals.dat') 
      open(4, file = ExpWdir//'Params.dat') 
      
      !write out initial,end and analytic values 
      do i = 1,xbc_len
         write(1,*) xbc(i),hbc_init(i),Gbc_init(i),ubc_init(i)
         write(2,*) xbc(i),hbc_fin(i),Gbc_fin(i),ubc_fin(i)
      end do
      
      !write out parameters      write(5,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(4,*) 'xstart :',xstart
      write(4,*) 'xend :',xbc(xbc_len)
      write(4,*) 'x_len :',xbc_len - 2*n_GhstCells
      write(4,*) 'n_GhstCells :',n_GhstCells
      write(4,*) 'xbc_len :',xbc_len
      write(4,*) 'dx = (x_end - x_start) / (x_len -1)  :' , dx
      write(4,*) 'tstart :', tstart
      write(4,*) 'tend :',tend 
      write(4,*) 'actual end time :', currenttime
      write(4,*) 'dt:' , dt
      write(4,*) 'gravity :' , ga
      write(4,*) 'hl :' , hl
      write(4,*) 'hr :' , hr
      write(4,*) 'beta1 :' , beta1
      write(4,*) 'beta2 :' , beta2
      
      
      write(EnergFile,*) dx,dt,beta1,beta2,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3),Energ_Err(4)
      
      end 


 
      program main
         
      implicit none
   
      Integer wdirlen,EnergFile
      PARAMETER(wdirlen= 300,EnergFile = 99)
      
      CHARACTER(len =wdirlen) wdir
     
      CHARACTER(len=2) strdiri
     
  
      integer expi,x_len,xbc_len,n_GhstCells,dtsep
      DOUBLE PRECISION a0,a1,a2,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha
     
      INTEGER effeclenwdir
      
      wdir = "/home/jp/Documents/" // 
     . "Work/PostDoc/Projects/Steve/1DWaves/" //
     . "RegularisedSerre/Data/RAW" //
     . "/Models/gSGNCA/VaryBeta/ImpDisp/Gaussian/timeseries/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(EnergFile, file = wdir(1:effeclenwdir)//'Energy.dat') 

      
      n_GhstCells = 6
      
      !SWWE equations
      ga = 9.81d0
      
      
      a0 = 1.0d0
      a1 = 0.5
      a2 = 10.0d0
      
      xstart = -250d0
      xend = 250d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 50d0
      
      dtsep = 50
      
      x_len = 20000
      xbc_len = x_len + 2 *n_GhstCells

      dx = (xend - xstart) / (x_len -1)
           
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 0,40
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
      
         !beta1 = -2d0/3d0 + (expi - 5)*0.1
         !beta2 = beta1 + 2d0/3d0
         
         beta1 = expi/100d0 - 0.1d0
         beta2 = expi/100d0 - 0.1d0
         
         !alpha is a factor on g*h, that determines wavespeed
         !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
         ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
         if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
            alpha = 1d0
         else
            alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
         end if
         
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
           
         !have to trim charachter string
         call SingleGaussian(a0,a1,a2,ga,beta1,beta2,xstart,
     .      xbc_len,n_GhstCells,dx,tstart,tend,dt,theta,
     .      EnergFile,wdir(1:effeclenwdir)//strdiri//'/',
     .      effeclenwdir+3,dtsep)
      
      end do
      
      close(EnergFile)
      
      end
