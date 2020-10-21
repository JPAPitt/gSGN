      subroutine DBSolve(hl,hr,grav,xstart,xbc_len,n_GhstCells,dx,
     . dt,startt,endt,currentt,
     . tolcalc,NormFile,ConsFile,ExpWdir,ExpWdir_len)
      
     
      implicit none
      integer xbc_len,n_GhstCells,i,ExpWdir_len,NormFile,ConsFile
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hAN(xbc_len),
     . uAN(xbc_len),GAN(xbc_len),hAI(xbc_len),
     . uAI(xbc_len),GAI(xbc_len),
     . hAE(xbc_len),uAE(xbc_len),GAE(xbc_len),
     . hl,hr,xstart,tolcalc,L2h,L2G,dx,
     . grav, dt, startt,endt,currentt,
     . TotalhN, TotalGN,TotalhI, TotalGI,TotalGFlux
     
     
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !Get initial conditions
      call Dambreak(xbc,xbc_len,startt,grav,hl,hr,hAI,uAI,GAI)
      
      do i = 1,xbc_len
         hAN(i) = hAI(i)
         uAN(i) = uAI(i)
         GAN(i) = GAI(i)
      end do
     
      ! solver
      call SWWESolver(hAN,GAN,xbc_len,n_GhstCells,grav,dx,dt,
     . tolcalc,startt,endt,currentt)
     
      call Dambreak(xbc,xbc_len,currentt,grav,hl,hr,hAE,uAE,GAE)
      
      
      call L2Norm(xbc_len,hAN,hAE,L2h) 
      call L2Norm(xbc_len,GAN,GAE,L2G) 
      
      TotalhN = sum(hAN)*dx 
      TotalGN = sum(GAN)*dx 
      TotalhI = sum(hAI)*dx 
      TotalGI = sum(GAI)*dx 
      TotalGFlux = grav*currentt/2.0*(hl**2 - hr**2)   

      write(ConsFile,*)  dx, (TotalhN -TotalhI)/TotalhI,
     .  (TotalGN -TotalGFlux)/TotalGFlux
      write(NormFile,*)  dx, L2h, L2G
      
      open(11, file = ExpWdir//'xhhGG.dat') 
      do i = 1,xbc_len
         write(11,*) xbc(i), hAN(i), hAE(i),GAN(i), GAE(i)
      end do
      close(11)

      open(12, file = ExpWdir//'xhiGi.dat') 
      do i = 1,xbc_len
         write(12,*) xbc(i), hAI(i),GAI(i)
      end do
      close(12)
      
      end
      
      
      
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,ConsFile
      PARAMETER(wdirlen= 300,NormFile=1,ConsFile=2)
      
      INTEGER effeclenwdir,n_GhstCells,lowestresx,x_len,xbc_len,
     . expi
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) strdiri
      
      DOUBLE PRECISION hl,hr,grav,xstart,xend,dx,
     . tolcalc,startt,endt,dt,currentt,Cr,maxwavespeed
     
      
      wdir = "../Results/Validation/Solver/DB21/15sCT/"
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      CALL SYSTEM('rm -rf '//wdir(1:effeclenwdir))
      CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir)) 
      
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 
      open(ConsFile, file = wdir(1:effeclenwdir)//'Cons.dat') 
      
      tolcalc = 10d0**(-10)
      
      hl = 2
      hr = 1
      grav = 9.81
      
      xstart = -100.0
      xend = 100.0
      
      n_GhstCells = 6
      lowestresx = 100
            
      startt = 0.0
      endt = 15.0
      
      
      do expi = 1,8
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(grav*max(hl,hr))
         dt  = (Cr / maxwavespeed) *dx

         print *, 'Experiment', expi
         call DBSolve(hl,hr,grav,xstart,xbc_len,n_GhstCells,dx,
     . dt,startt,endt,currentt,
     . tolcalc,NormFile,ConsFile,wdir(1:effeclenwdir)//strdiri//'/',
     . effeclenwdir+3)
     


      end do

      close(NormFile)
      close(ConsFile)

      
     
      end
 
