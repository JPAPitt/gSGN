      program main
   
      implicit none
   
      Integer wdirlen, xbc_len,NormFileI,tlist_len
      PARAMETER(wdirlen= 300,NormFileI = 10,tlist_len=3)
      
      CHARACTER(len =wdirlen) wdir, expdir
      CHARACTER(len=2) expnumstr
  
      integer x_len,n_GhstCells
      DOUBLE PRECISION tlist(tlist_len)
     
      DOUBLE PRECISION xstart,xend,tstart,tend,
     . dx,dt,theta,Cr,maxwavespeed,
     . ga,a0,a1,a2,a3,a4,a5,a6,a7,dtodx
     
      INTEGER effeclenwdir,effeclenexpdir,n,lowi,i
      
      wdir = "/home/jp/Documents/" // 
     . "Work/PostDoc/Projects/Steve/1DWaves/" //
     . "RegularisedSerre/Data/RAW" //
     . "/Models/SerreST/ForcedSolutionShort/"
     
      expdir = wdir
     
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      !Remove previous runs, and make directory to dump data
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      n = 12
      lowi = 100
      n_GhstCells = 6
            
      xstart = -50d0
      xend = 100.0d0
      
      tstart = 0d0
      tend = 1d0 !10.0d0
            
      ga = 1d0
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      a5 = 0d0
      a6 = 0d0
      a7 = 1d0/3d0
      
      tlist(1) = 2.5d0
      tlist(2) = 5d0
      tlist(3) = 7.5d0
      
      theta = 1.2
   
      open(NormFileI, file = wdir(1:effeclenwdir)//'Norms.dat') 
      
      do i = 0,n
         write (expnumstr,'(I2.2)') i
         
         expdir = wdir(1:effeclenwdir)//expnumstr//'/'
         call LenTrim(expdir,wdirlen,effeclenexpdir)
         
         CALL SYSTEM('mkdir -p '// expdir(1:effeclenexpdir))
      
      
         x_len = lowi*2**(i)
         xbc_len= x_len + 2 *n_GhstCells
      
         dx = (xend - xstart) / (x_len -1) 

         maxwavespeed = a2 + a4 + dsqrt(ga*(a0 + a1))
         Cr = 0.5
         dtodx = Cr/maxwavespeed
         dt = dx*dtodx 

    
         print *,'++++++++++++ Experiment : ',i ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         call NumericalSolveForced(xbc_len,n_GhstCells,
     .      tstart,tend,xstart,dx,dt,ga,theta, 
     .      a0,a1,a2,a3,a4,a5,a6,a7,
     .      tlist,tlist_len,expdir(1:effeclenexpdir),effeclenexpdir,
     .      NormFileI)
     
      end do
      close(NormFileI) 
      
      end 
