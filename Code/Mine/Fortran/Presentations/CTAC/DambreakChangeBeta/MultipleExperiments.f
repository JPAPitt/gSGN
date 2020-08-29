      subroutine SWWEDBBetaChange(h0,h1,ga,beta1list,beta2list,
     . beta1_len,beta2_len, xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)
     
      implicit none
      DOUBLE PRECISION h0,h1,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
      integer xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile,
     . tlist_len,beta1_len,beta2_len
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      DOUBLE PRECISION Norms(3),sumerr(3),suma(3)
      DOUBLE PRECISION tlist(tlist_len),beta1list(beta1_len),
     . beta2list(beta2_len)
      integer i
      DOUBLE PRECISION currenttime
      
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes     
      call Dambreak(xbc,xbc_len,tstart,ga,h0,h1,
     . hbc_init,ubc_init,Gbc_init)
     
      
      !solve gSGN with beta values until currenttime > tend     
      call NumericalSolveTSPrint(tstart,tend,
     . ga,beta1list,beta2list,beta1_len,beta2_len,
     . theta,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energs_init, Energs_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len)
     
      ! get analytic values of h,u,G     
      call Dambreak(xbc,xbc_len,currenttime,ga,h0,h1,
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
            Norms(i) = dsqrt(sumerr(i))
         else
            Norms(i) = dsqrt(sumerr(i)) / dsqrt(suma(i))
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
      open(1, file = ExpWdir//'Init.dat') 
      open(2, file = ExpWdir//'End.dat') 
      open(3, file = ExpWdir//'EndAna.dat') 
      open(4, file = ExpWdir//'Params.dat') 
      open(5, file = ExpWdir//'Energy.dat') 
      
      !write out initial,end and analytic values 
      do i = 1,xbc_len
         write(1,*) tstart,xbc(i),hbc_init(i),Gbc_init(i),ubc_init(i)
         write(2,*) currenttime,xbc(i),hbc_fin(i),Gbc_fin(i),ubc_fin(i)
         write(3,*) currenttime,xbc(i),hbc_fin_a(i),Gbc_fin_a(i),
     .               ubc_fin_a(i)
      end do
      
      !write out parameters      write(5,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(4,*) 'xstart',xstart
      write(4,*) 'xend',xbc(xbc_len)
      write(4,*) 'x_len',xbc_len - 2*n_GhstCells
      write(4,*) 'n_GhstCells',n_GhstCells
      write(4,*) 'xbc_len',xbc_len
      write(4,*) 'dx' , dx
      write(4,*) 'tstart', tstart
      write(4,*) 'tend',tend 
      write(4,*) 'actual_end_time', currenttime
      write(4,*) 'dt' , dt
      write(4,*) 'theta' , theta
      write(4,*) 'gravity' , ga
      write(4,*) 'h0' , h0
      write(4,*) 'h1' , h1
      
      !write out enery
      write (5,*) 'When , h , G ,  uh , Energy'
      write(5,*) 'Initial ',Energs_init(1),Energs_init(2),
     .   Energs_init(3),Energs_init(4)
      write(5,*) 'End ',Energs_fin(1),Energs_fin(2),
     .   Energs_fin(3),Energs_fin(4)
      
      !write out information for group of experiments (Norms, Energy)
      write(NormFile,*) dx,Norms(1),Norms(2),Norms(3)
      
      write(EnergFile,*) dx,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3),Energ_Err(4)
      
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      end
      

 
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,EnergFile,tlist_len,beta1_len,beta2_len
      PARAMETER(wdirlen= 300,NormFile = 98, EnergFile = 99,
     . tlist_len=2001,beta1_len=4, beta2_len = 4)
      
      CHARACTER(len =wdirlen) wdir
     
      CHARACTER(len=2) strdiri
     
  
      integer expi,x_len,xbc_len,n_GhstCells, lowestresx
      DOUBLE PRECISION h0,h1,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha
     
      DOUBLE PRECISION tlist(tlist_len),beta1list(beta1_len),
     . beta2list(beta2_len)
     
      INTEGER effeclenwdir
      
      wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/"//
     . "RegularisedSerre/Data/RAW/Presentations/CTAC/"//
     . "Animation/DB6/" 
     
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
      
      beta1list(1) = 0d0
      beta1list(2) = 2d0/15d0
      beta1list(3) = 2d0/15d0 -2d0/3d0
      beta1list(4) = -2d0/3d0
      
      beta2list(1) = 0d0
      beta2list(2) = 2d0/15d0
      beta2list(3) = 2d0/15d0
      beta2list(4) = 0d0
      
      h0 = 2.0
      h1 = 1.0
      
      xstart = -1000d0
      xend = 1000d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 200d0
      
      beta1 = beta1list(1)
      beta2 = beta2list(1)
      
      call EqualSpaced(tstart,tend,tlist_len,0.1, tlist)
      
      !alpha is a factor on g*h, that determines wavespeed
      !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
      ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
      if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
         alpha = 1d0
      else
         alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
      end if
      
      lowestresx = 100
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 0,0
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(7)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(h0 + h1))
         dt  = (Cr / maxwavespeed) *dx
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         
         call SWWEDBBetaChange(h0,h1,ga,beta1list,beta2list,beta1_len,
     .   beta2_len,xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,tlist,tlist_len,
     .      wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

      
      end do
      
      close(EnergFile)
      close(NormFile)
      
      end
