      
      subroutine SingleSWWESDB(hl,hr,dbalpha,ga,b10,b11,b20,b21,t0,
     . xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     . dt,theta,tlist,tlist_len,ExpWdir,ExpWdir_len)     
       
      implicit none
      DOUBLE PRECISION hl,hr,dbalpha,ga,b10,b11,b20,b21,t0,xstart,dx,
     . tstart,tend,dt,theta
      integer xbc_len,n_GhstCells,ExpWdir_len,
     . tlist_len
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . Gbc_init(xbc_len),ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len),ubc_fin(xbc_len),hbc_fin_a(xbc_len),
     . Gbc_fin_a(xbc_len),ubc_fin_a(xbc_len)
     
      
      DOUBLE PRECISION Energs0_init(4), Energs0_fin(4)
      DOUBLE PRECISION Energs1_init(4), Energs1_fin(4)
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
     . ga,b10,b11,b20,b21,t0,theta,dx,dt,n_GhstCells,xbc_len,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . Energs0_init, Energs0_fin,
     . Energs1_init, Energs1_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len,hl,hr)
     
      ! get DB analytic values of h,u,G
      call Dambreak(xbc,xbc_len,currenttime,ga,hl,hr,
     . hbc_fin_a,ubc_fin_a,Gbc_fin_a)
     
           
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
      write(4,*) 'b10 :' , b10
      write(4,*) 'b11 :' , b11
      write(4,*) 'b20 :' , b20
      write(4,*) 'b21 :' , b21
      write(4,*) 't0 :' , t0
      
      !write out enery
      write (5,*) 'When , h , G ,  uh , Energy'
      write(5,*) 'Initial 0 : ',Energs0_init(1),Energs0_init(2),
     .   Energs0_init(3),Energs0_init(4)
      write(5,*) 'End 0 : ',Energs0_fin(1),Energs0_fin(2),
     .   Energs0_fin(3),Energs0_fin(4)
      write(5,*) 'Initial 1 : ',Energs1_init(1),Energs1_init(2),
     .   Energs1_init(3),Energs1_init(4)
      write(5,*) 'End 1 : ',Energs1_fin(1),Energs1_fin(2),
     .   Energs1_fin(3),Energs1_fin(4)      
      
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      
      end
      

 
      program main
         
      implicit none
   
      Integer wdirlen,tlist_len
      PARAMETER(wdirlen= 300,tlist_len=200)
      
      CHARACTER(len =wdirlen) wdir
     
      integer x_len,xbc_len,n_GhstCells
      DOUBLE PRECISION hl,hr,ga,xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,dbalpha,
     . b10,b11,b20,b21,t0,alpha0,alpha1,alpha
     
      INTEGER effeclenwdir
     
      DOUBLE PRECISION tlist(tlist_len)
      
      wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/" //
     . "1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll" //
     . "/BetaVary/SmoothDBInvestigation/" //
     . "OneOffTests/10to1/dbalpha0p5/SWWE2Serre/"
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
     

      n_GhstCells = 6
      
      !SWWE equations
      ga = 9.81d0
      hl = 10.0
      hr = 1.0
      dbalpha = 0.5
      
      xstart = -1000d0
      xend = 1000d0
      
      theta = 1.2d0
      
      tstart = 0d0
      tend = 100d0
      
      call EqualSpaced(tstart,tend,tlist_len, tlist)

      x_len = 1000 *2**(4)
      xbc_len = x_len + 2 *n_GhstCells
      
      b20 = 0.0
      b10 = b20 -2d0/3d0
      
      b11 = 0d0
      b21 = 0d0
      
      t0 = 20d0
      
      !alpha is a factor on g*h, that determines wavespeed
      !when beta1 ~ -2/3, then this ratio would go to infinity unless beta1 = 0
      ! thus we limit ourselves to beta1 ~ -2/3 only when beta1 = 0
      if  (dabs(2d0/3d0 + b10) < 10d0**(-10))  then
         alpha0 = 1d0
      else
         alpha0 = max(1d0,b20 / (2d0/3d0 + b10))
      end if
      
      if  (dabs(2d0/3d0 + b11) < 10d0**(-10))  then
         alpha1 = 1d0
      else
         alpha1 = max(1d0,b21 / (2d0/3d0 + b11))
      end if
      
      alpha = max(alpha0,alpha1)

      dx = (xend - xstart) / (x_len -1)
      Cr = 0.5
      maxwavespeed = dsqrt(alpha*ga*(hl + hr))
      dt  = (Cr / maxwavespeed) *dx
      
      
      call SingleSWWESDB(hl,hr,dbalpha,ga,b10,b11,b20,b21,t0,
     .         xstart,xbc_len,n_GhstCells,dx,tstart,tend,
     .         dt,theta,tlist,tlist_len,
     .         wdir(1:effeclenwdir),effeclenwdir)  

      end
