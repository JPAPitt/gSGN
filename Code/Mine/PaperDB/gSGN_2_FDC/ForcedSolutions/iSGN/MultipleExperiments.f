      
      subroutine ForceSol(ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,NormFile,EnergFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)
     
      implicit none
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,dx,tstart,tend,dt
      integer xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile,
     . tlist_len
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . ubc_init(xbc_len),hbc_fin(xbc_len),
     . ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . ubc_fin_a(xbc_len),plchld(xbc_len),phbc_init(xbc_len),
     . pubc_init(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(3), Energs_fin(3),Energ_Err(3)
      DOUBLE PRECISION Norms(2),sumerr(2),suma(2)
      DOUBLE PRECISION tlist(tlist_len)
      integer i
      DOUBLE PRECISION currenttime
      
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      print *, 1
      
      !get initial conditions at all cell nodes     
      call ForcedSol(xbc,xbc_len,tstart,a0,a1,a2,a3,a4,
     . hbc_init,ubc_init) 
      call ForcedSol(xbc,xbc_len,tstart-dt,a0,a1,a2,a3,a4,
     . phbc_init,pubc_init) 
      
      print *, 2
      !solve gSGN with beta values until currenttime > tend     
      call NumericalSolveTSPrint(tstart,tend,
     . ga,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,ubc_init,
     . phbc_init,pubc_init,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . currenttime,hbc_fin,ubc_fin,
     . Energs_init, Energs_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len)
     
     
      ! get analytic values of h,u,G
      call ForcedSol(xbc,xbc_len,currenttime,a0,a1,a2,a3,a4,
     . hbc_fin_a,ubc_fin_a) 
     
     
c Convergence Norm Tests     
      !calculate norm values for h,u,G
      !Convergence Norms
      do i = 1,2
         sumerr(i) = 0
         suma(i) = 0
      end do
      
      !sum L2 norm of q - q*, and q*
      do i = 1, xbc_len
         sumerr(1) = sumerr(1) + (hbc_fin_a(i) - hbc_fin(i))**2
         suma(1) = suma(1) + (hbc_fin_a(i))**2
         
         sumerr(2) = sumerr(2) + (ubc_fin_a(i) - ubc_fin(i))**2
         suma(2) = suma(2) + (ubc_fin_a(i))**2
      end do
      
      do i = 1,2
         !if sum is very small, just return absolute error
         if (suma(i) .lt. 10d0**(-10)) then
            Norms(i) = dsqrt(sumerr(i))
         else
            Norms(i) = dsqrt(sumerr(i)) / dsqrt(suma(i))
         end if
      end do

c Conservation Norm Tests  

      do i = 1,3
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
         write(1,*) tstart,xbc(i),hbc_init(i),ubc_init(i)
         write(2,*) currenttime,xbc(i),hbc_fin(i),ubc_fin(i)
         write(3,*) currenttime,xbc(i),hbc_fin_a(i),ubc_fin_a(i)
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
      write(4,*) 'gravity' , ga
      write(4,*) 'a0' , a0
      write(4,*) 'a1' , a1
      write(4,*) 'a2' , a2
      write(4,*) 'a3' , a3
      write(4,*) 'a4' , a4
      write(4,*) 'a5' , a5
      write(4,*) 'b1a6' , b1a6
      write(4,*) 'b1a7' , b1a7   
      write(4,*) 'b2a6' , b2a6
      write(4,*) 'b2a7' , b2a7  
      
      !write out enery
      write (5,*) 'When , h , G ,  uh , Energy'
      write(5,*) 'Initial ',Energs_init(1),Energs_init(2),
     .   Energs_init(3)
      write(5,*) 'End ',Energs_fin(1),Energs_fin(2),
     .   Energs_fin(3)
      
      !write out information for group of experiments (Norms, Energy)
      write(NormFile,*) dx,Norms(1),Norms(2)
      
      write(EnergFile,*) dx,Energ_Err(1),Energ_Err(2),
     .   Energ_Err(3)
      
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      end
      

 
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,EnergFile,tlist_len
      PARAMETER(wdirlen= 300,NormFile = 98, EnergFile = 99,
     . tlist_len=31)
      
      CHARACTER(len =wdirlen) wdir
     
      CHARACTER(len=2) strdiri
     
  
      integer expi,x_len,xbc_len,n_GhstCells, lowestresx
      DOUBLE PRECISION ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,dtodx
     
      DOUBLE PRECISION tlist(tlist_len)
     
      INTEGER effeclenwdir
      
      wdir = "../../Results/FSol/iSGN/"  
     
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
      
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      
      a5 = 0.0

      b1a6 = 0.0
      b1a7 = 2.0/3.0 + 2.0 / 15.0
      
      
      b2a6 = 0.0
      b2a7 = 2.0/ 15.0
      
      xstart = -100d0
      xend = 100d0
      
      
      tstart = 0d0
      tend = 10.0d0
      
      call EqualSpaced(tstart,tend,tlist_len,1.0, tlist)
      
      
      lowestresx = 100
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 0,11
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)
         
         alpha = max(1d0,b2a7 / (b1a7))
         maxwavespeed = a2 + a4 + dsqrt(alpha*ga*(a0 + a1))
         Cr = 0.5
         dtodx = Cr/maxwavespeed
         dt = dx*dtodx 

         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
                  
         call ForceSol(ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,NormFile,EnergFile,tlist,tlist_len,
     .      wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

      
      end do
      
      close(EnergFile)
      close(NormFile)
      
      end
