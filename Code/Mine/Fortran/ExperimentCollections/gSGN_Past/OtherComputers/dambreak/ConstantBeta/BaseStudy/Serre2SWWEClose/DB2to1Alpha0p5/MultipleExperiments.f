      
      subroutine SingleSWWESDB(hl,hr,dbalpha,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,NormFile,EnergFile,tlist,tlist_len,ExpWdir,ExpWdir_len)       
      implicit none
      DOUBLE PRECISION hl,hr,dbalpha,ga,beta1,beta2,xstart,dx,
     . tstart,tend,dt,theta
      integer xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile,
     . tlist_len
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len)
     
      
      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      DOUBLE PRECISION Norms(3),sumerr(3),suma(3)
      DOUBLE PRECISION tlist(tlist_len)
      integer i
      DOUBLE PRECISION currenttime
      

      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !get initial conditions at all cell nodes
      call SmoothDB(xbc,xbc_len,hl,hr,dbalpha,
     . hbc_init,ubc_init,Gbc_init)

      
      !solve gSGN with beta values until currenttime > tend     
      call NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energs_init, Energs_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len)
     
      ! get DB analytic values of h,u,G
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
      write(4,*) 'dbalpha :' , dbalpha
      write(4,*) 'beta1 :' , beta1
      write(4,*) 'beta2 :' , beta2
      
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
   
      Integer wdirlen,NormFile,EnergFile,tlist_len
      PARAMETER(wdirlen= 300,NormFile = 98, 
     .   EnergFile = 99,tlist_len=35)
      
      CHARACTER(len =wdirlen) wdir, betawdir,dxbetawdir
     
      CHARACTER(len=2) strdirdxi,strdirbetai
     
  
      integer x_len,xbc_len,n_GhstCells, lowestresx 
      DOUBLE PRECISION hl,hr,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,beta1,beta2,alpha,dbalpha
     
      INTEGER effeclenwdir,effeclenbetawdir,effeclendxbetawdir,
     . betaexpi,dxexpi,dxres_num,betapair_num
     
      DOUBLE PRECISION tlist(tlist_len)
      
      wdir = "../../../ConstantBeta/SmoothDBInvestigation/" //
     . "Serre2SWWEClose/AspRat2to1/DBalpha0p5/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
     

      n_GhstCells = 6
      dxres_num = 6
      betapair_num = 20
      
      !SWWE equations
      ga = 9.81d0
      hl = 2.0
      hr = 1.0
      dbalpha = 0.5
      
      xstart = -250d0
      xend = 250d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 35d0
      
      call EqualSpaced(tstart,tend,tlist_len, tlist)
      
      
      lowestresx = 1000
      
      
      do betaexpi = 0,betapair_num
      
         write (strdirbetai,'(I2.2)') betaexpi
         
         betawdir = wdir(1:effeclenwdir) // 'Beta' //strdirbetai//'/'
         call LenTrim(betawdir,wdirlen,effeclenbetawdir)
            
         CALL SYSTEM('mkdir -p '// betawdir(1:effeclenbetawdir))
         

         beta1 = -2d0/3d0 + 0.1*( (1.0*betaexpi) / (betapair_num))
         beta2 = 0d0             
      
         
         !alpha is a factor on g*h, that determines wavespeed
         !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
         ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
         if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
            alpha = 1d0
         else
            alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
         end if
         
         !open output files
         open(EnergFile, file = betawdir(1:effeclenbetawdir) //
     .      'Energy.dat') 
         open(NormFile, file = betawdir(1:effeclenbetawdir) //
     .      'Norms.dat') 
      
         !perform the soliton experiment a number of times, decreasing \Delta x each time
         do dxexpi = 0,dxres_num

            write (strdirdxi,'(I2.2)') dxexpi
         
            dxbetawdir = betawdir(1:effeclenbetawdir) // 'dx' //
     .            strdirdxi //'/'
     
            call LenTrim(dxbetawdir,wdirlen,effeclendxbetawdir)
               
            CALL SYSTEM('mkdir -p '// dxbetawdir(1:effeclendxbetawdir))
                        
            
            x_len = lowestresx *2**(dxexpi)
            xbc_len = x_len + 2 *n_GhstCells

            dx = (xend - xstart) / (x_len -1)
            Cr = 0.5
            maxwavespeed = dsqrt(alpha*ga*(hl + hr))
            dt  = (Cr / maxwavespeed) *dx
                        
            print *, '++++++++++++ Beta Num : ', 
     .         betaexpi ,'|| dx num : ' , dxexpi, '++++++++++++'
            
            
            call SingleSWWESDB(hl,hr,dbalpha,ga,beta1,beta2,
     .         xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     .         dt,theta,NormFile,EnergFile,tlist,tlist_len,
     .         dxbetawdir(1:effeclendxbetawdir),effeclendxbetawdir)  

         
         end do
         
         close(EnergFile)
         close(NormFile)
      
      end do
      
      
      end
