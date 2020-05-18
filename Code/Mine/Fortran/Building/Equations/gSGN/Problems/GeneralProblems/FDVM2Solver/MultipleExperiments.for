
      module IndividualExperiments
      use FDVM2Procedures
      use InitialConditionProblems
      implicit none
      
      contains
      
      subroutine SingleSerreSoliton(a0,a1,ga,beta1,beta2
     . ,xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,ExpWdir)
     
      real*8, intent(in) :: a0,a1,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
      integer, intent(in) :: xbc_len,n_GhstCells
      CHARACTER(len=*), intent(in) :: ExpWdir
      integer, intent(in) :: NormFile,EnergFile
      
      real*8, dimension(xbc_len) ::xbc,hbc_init,Gbc_init,ubc_init,
     . hbc_fin,Gbc_fin,ubc_fin,hbc_fin_a,Gbc_fin_a,ubc_fin_a
     
      
      real*8, dimension(4) :: Energs_init, Energs_fin,Energ_Err
      real*8, dimension(3) :: Norms,sumerr,suma
      integer i
      real*8 currenttime
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call SerreSoliton(xbc,tstart,ga,a0,a1,hbc_init,ubc_init,Gbc_init)
      
      !solve gSGN with beta values until currenttime > tend
      call NumericalSolve(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . Energs_init,currenttime,hbc_fin,Gbc_fin,ubc_fin,Energs_fin)
     
      ! get analytic values of h,u,G
      call SerreSoliton(xbc,currenttime,ga,a0,a1,
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
      write(4,*) 'a0 :' , a0
      write(4,*) 'a1 :' , a1
      write(4,*) 'beta1 :' , beta1
      write(4,*) 'beta2 :' , beta2
      
      !write out information for group of experiments (Norms, Energy)
      write(NormFile,*) dx,Norms(1),Norms(2),Norms(3)
      
      write(EnergFile,*) dx,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3),Energ_Err(4)
      
      end subroutine SingleSerreSoliton
      
      end module IndividualExperiments
 
      program main
      use IndividualExperiments
         
      implicit none
   
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/Models/gSGN/ConstantBetas/Serre/LoopSoliton/"
     
      CHARACTER(len=2) :: strdiri
     
      integer, parameter:: NormFile = 98, EnergFile = 99
  
      integer expi,x_len,xbc_len,n_GhstCells, lowestresx
      real*8 a0,a1,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(EnergFile, file = wdir//'Energy.dat') 
      open(NormFile, file = wdir//'Norms.dat') 

      
      n_GhstCells = 6
      
      !Serre equations
      beta1 = 0d0
      beta2 = 0d0
      
      a0 = 1.0
      a1 = 0.7
      ga = 9.81
      
      xstart = -50d0
      xend = 100d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 10d0
      
      !alpha is a factor on g*h, that determines wavespeed
      !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
      ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
      if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
         alpha = 1
      else
         alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
      end if
      
      lowestresx = 100
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 0,10
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir //'/'//strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         call SingleSerreSoliton(a0,a1,ga,beta1,beta2,xstart,
     .      xbc_len,n_GhstCells,dx,tstart,tend,dt,theta,
     .      NormFile,EnergFile,wdir//'/'//strdiri//'/')
      
      end do
      
      close(EnergFile)
      close(NormFile)
      

      
      end program main
