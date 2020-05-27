c =========================
c DB single experiment
c
c ====================   
      subroutine SingleSWWEDB(hl,hr,ga,beta1,beta2
     . ,xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,ExpWdir,ExpWdir_len)
     
      DOUBLE PRECISION hl,hr,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
     
      integer xbc_len,n_GhstCells,NormFile,EnergFile,ExpWdir_len
      
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len), ubc_init(xbc_len),  hbc_fin(xbc_len),
     . Gbc_fin(xbc_len), ubc_fin(xbc_len), hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      DOUBLE PRECISION Norms(3),sumerr(3),suma(3)
      integer i
      DOUBLE PRECISION currenttime
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call Dambreak(xbc,xbc_len,tstart,
     . ga,hl,hr,hbc_init,ubc_init,Gbc_init) 
      
      !solve gSGN with beta values until currenttime > tend
      call NumericalSolve(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc_len,
     . hbc_init,Gbc_init,ubc_init,
     . Energs_init,currenttime,hbc_fin,Gbc_fin,ubc_fin,Energs_fin)
     
      ! get analytic values of h,u,G
      call Dambreak(xbc,xbc_len,currenttime,ga,hl,hr,
     . hbc_fin_a,ubc_fin_a,Gbc_fin_a) 
     
     
c Convergence Norm Tests     
      !calculate norm values for h,u,G
      !Convergence Norms
      do i = 1,3
         sumerr(i) = 0
         suma(i) = 0
      end do
      
      !sum L2 norm of q - q*, and q*
      do i = 1, xbc_len
         sumerr(1) = sumerr(1) + (hbc_fin_a(i) - hbc_fin(i))**2
         suma(1) = suma(1) + (hbc_fin_a(i))**2
         
         sumerr(2) = sumerr(2) + (ubc_fin_a(i) - ubc_fin(i))**2
         suma(2) = suma(2) + (ubc_fin_a(i))**2
         
         sumerr(3) = sumerr(3) + (Gbc_fin_a(i) - Gbc_fin(i))**2
         suma(3) = suma(3) + (Gbc_fin_a(i))**2 
      end do
      
      do i = 1,3
         !if sum is very small, just return absolute error
         if (suma(i) .lt. 10d0**(-10)) then
            Norms(i) = sqrt(sumerr(i))
         else
            Norms(i) = sqrt(sumerr(i)) / sqrt(suma(i))
         end if
      end do

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
      open(3, file = ExpWdir//'EndAnaVals.dat') 
      open(4, file = ExpWdir//'Params.dat') 
      
      !write out initial,end and analytic values 
      do i = 1,xbc_len
         write(1,*) xbc(i),hbc_init(i),Gbc_init(i),ubc_init(i)
         write(2,*) xbc(i),hbc_fin(i),Gbc_fin(i),ubc_fin(i)
         write(3,*) xbc(i),hbc_fin_a(i),Gbc_fin_a(i),ubc_fin_a(i)
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
      
      !write out information for group of experiments (Norms, Energy)
      write(NormFile,*) dx,dt,beta1,beta2,Norms(1),Norms(2),Norms(3)
      
      write(EnergFile,*) dx,dt,beta1,beta2,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3),Energ_Err(4)
      
      end 


 
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,EnergFile
      PARAMETER(wdirlen= 300,NormFile = 98, EnergFile = 99)
      
      CHARACTER(len =wdirlen) wdir
     
      CHARACTER(len=2) strdiri
     
  
      integer expi,x_len,xbc_len,n_GhstCells
      DOUBLE PRECISION hl,hr,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha
     
      INTEGER effeclenwdir
      
      wdir = "/home/jp/Documents/" // 
     . "Work/PostDoc/Projects/Steve/1DWaves/" //
     . "RegularisedSerre/Data/RAW" //
     . "/Models/gSGN/VaryBeta/ImpDisp/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(EnergFile, file = wdir(1:effeclenwdir)//'Energy.dat') 
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 

      
      n_GhstCells = 6
      
      !SWWE equations
      ga = 9.81d0
      
      
      hl = 2.0d0
      hr = 1.0d0
      
      xstart = -50d0
      xend = 50d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 7d0
      
      x_len = 10000
      xbc_len = x_len + 2 *n_GhstCells

      dx = (xend - xstart) / (x_len -1)
           
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 0,30
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
      
         !beta1 = -2d0/3d0 + expi*0.1
         !beta2 = beta1 + 2d0/3d0
         beta1 = expi*1d0/90d0
         beta2 = expi*1d0/90d0
         
         !alpha is a factor on g*h, that determines wavespeed
         !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
         ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
         if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
            alpha = 1d0
         else
            alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
         end if
         
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(hl))
         dt  = (Cr / maxwavespeed) *dx
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
           
         !have to trim charachter string
         call SingleSWWEDB(hl,hr,ga,beta1,beta2,xstart,
     .      xbc_len,n_GhstCells,dx,tstart,tend,dt,theta,
     .      NormFile,EnergFile,
     .      wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)
      
      end do
      
      close(EnergFile)
      close(NormFile)
      
      end
