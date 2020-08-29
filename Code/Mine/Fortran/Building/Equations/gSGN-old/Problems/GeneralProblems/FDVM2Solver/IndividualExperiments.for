      program main
      use FDVM2Procedures
      use InitialConditionProblems
   
      implicit none
   
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/Models/gSGN/ConstantBetas/Serre/Soliton/"
  
      integer i,x_len,xbc_len,n_GhstCells
      real*8, dimension(5012) :: xbc,hbc_init,Gbc_init,ubc_init,
     . hbc_fin,Gbc_fin,ubc_fin,hbc_fin_a,Gbc_fin_a,ubc_fin_a
      real*8, dimension(4) :: Energs_init, Energs_fin
      real*8 a0,a1,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha,currenttime
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(1, file = wdir//'InitVal.dat') 
      open(2, file = wdir//'EndVals.dat') 
      open(3, file = wdir//'EndAnaVals.dat') 
      open(4, file = wdir//'Params.dat') 
      open(5, file = wdir//'Energy.dat')
      
      n_GhstCells = 6
      xbc_len = size(xbc)
      x_len = xbc_len - 2 *n_GhstCells
      
      !Serre equations
      beta1 = 0d0
      beta2 = 0d0
      
      a0 = 1.0
      a1 = 0.7
      ga = 9.81
      
      xstart = -50d0
      xend = 100d0
      
      theta = 1.2d0
      
      !ensures that x(0) = xstart and x(x_len) = x_end
      dx = (xend - xstart) / (x_len -1) 
      
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
      
      Cr = 0.5
      maxwavespeed = dsqrt(alpha*ga*(a0 + a1))
      dt  = (Cr / maxwavespeed) *dx
      
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
      write(4,*) 'a0 :' , a0
      write(4,*) 'a1 :' , a1
      write(4,*) 'beta1 :' , beta1
      write(4,*) 'beta2 :' , beta2
  
      
      !write out energies
      write(5,*) 'initial h  G  uh  H'
      write(5,*) Energs_init(1), Energs_init(2), 
     . Energs_init(3), Energs_init(4)
      write(5,*) ''
      write(5,*) 'End h  G  uh  H'
      write(5,*) Energs_fin(1), Energs_fin(2), 
     . Energs_fin(3), Energs_fin(4)   
      write(5,*) ''  
      write(5,*) 'Relative h  G  uh  H'
      write(5,*) 
     . dabs(Energs_fin(1) - Energs_init(1))/ dabs(Energs_init(1)), 
     . dabs(Energs_fin(2) - Energs_init(2))/ dabs(Energs_init(2)),
     . dabs(Energs_fin(3) - Energs_init(3))/ dabs(Energs_init(3)),
     . dabs(Energs_fin(4) - Energs_init(4))/ dabs(Energs_init(4))
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)

      
      end program main
