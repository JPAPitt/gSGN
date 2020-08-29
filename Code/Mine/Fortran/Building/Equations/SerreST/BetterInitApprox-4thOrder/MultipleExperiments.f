      
      subroutine SingleSerrePeakon(a0,a1,ga,tau
     . ,xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,ExpWdir,ExpWdir_len)
     
      implicit none
      DOUBLE PRECISION a0,a1,ga,tau,xstart,dx,
     . tstart,tend,dt,theta
      integer xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      DOUBLE PRECISION Norms(3),sumerr(3),suma(3)
      DOUBLE PRECISION tlist(1)
      integer i
      DOUBLE PRECISION currenttime
      
      tlist(1) = (tend - tstart) / 2d0
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call Peakon(xbc,xbc_len,tstart,
     . dx,ga,a0,a1,hbc_init,ubc_init,Gbc_init)
     
      
      !solve gSGN with beta values until currenttime > tend     
      call NumericalSolveTSPrint(tstart,tend,
     . ga,tau,theta,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energs_init, Energs_fin,
     . tlist,1,ExpWdir,ExpWdir_len)
     
      ! get analytic values of h,u,G
      call Peakon(xbc,xbc_len,currenttime,dx,ga,a0,a1,
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
         sumerr(1) = sumerr(1) + dx*dabs(hbc_fin_a(i) - hbc_fin(i))
         suma(1) = suma(1) + dx*dabs(hbc_fin_a(i))
         
         sumerr(2) = sumerr(2) + dx*dabs(ubc_fin_a(i) - ubc_fin(i))
         suma(2) = suma(2) + dx*dabs(ubc_fin_a(i))
         
         sumerr(3) = sumerr(3) + dx*dabs(Gbc_fin_a(i) - Gbc_fin(i))
         suma(3) = suma(3) + dx*dabs(Gbc_fin_a(i))
      end do
      
      do i = 1,3
         !if sum is very small, just return absolute error
         if (suma(i) .lt. 10d0**(-10)) then
            Norms(i) = (sumerr(i))
         else
            Norms(i) = (sumerr(i)) / (suma(i))
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
      write(4,*) 'tau :' , tau
      
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
   
      Integer wdirlen,NormFile,EnergFile
      PARAMETER(wdirlen= 300,NormFile = 98, EnergFile = 99)
      
      CHARACTER(len =wdirlen) wdir
     
      CHARACTER(len=2) strdiri
     
  
      integer expi,x_len,xbc_len,n_GhstCells, lowestresx
      DOUBLE PRECISION a0,a1,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,tau
     
      INTEGER effeclenwdir
      
      wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/" //
     . "1DWaves/RegularisedSerre/Data/RAW/Models/SerreST" //
     . "/Peakon4-L1/";
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open output files
      open(EnergFile, file = wdir(1:effeclenwdir)//'Energy.dat') 
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 

      
      n_GhstCells = 6
      
      !SWWE equations
      ga = 1d0
      
      tau = 1d0/3d0
      
      a0 = 1.0
      a1 = 0.5
      
      xstart = -50d0
      xend = 50d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 0d0
      
      
      lowestresx = 100
      
      !perform the soliton experiment a number of times, decreasing \Delta x each time
      do expi = 1,12
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         
         call SingleSerrePeakon(a0,a1,ga,tau,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,
     .      wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

      
      end do
      
      close(EnergFile)
      close(NormFile)
      
      end
