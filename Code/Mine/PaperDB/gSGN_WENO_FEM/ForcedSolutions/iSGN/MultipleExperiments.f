      subroutine SingleForced(a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,xbc_len,n_GhstCells,cubic_len,dx,tstart,tend,
     . dt,tolcalc,NormFile,AvgNormFile,EnergFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)
     
      implicit none
      DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,dx,
     . tstart,tend,dt,tolcalc
      integer xbc_len,n_GhstCells,ExpWdir_len,NormFile,EnergFile,
     . tlist_len,cubic_len,AvgNormFile
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hbccub_init(cubic_len),
     . Gbccub_init(cubic_len),ubccub_init(cubic_len),
     . hbccub_fin(cubic_len),Gbccub_fin(cubic_len),
     . ubccub_fin(cubic_len),hbccub_fina(cubic_len),
     . Gbccub_fina(cubic_len),ubccub_fina(cubic_len),
     . habc(xbc_len),Gabc(xbc_len),uabc(xbc_len),
     . habca(xbc_len),Gabca(xbc_len),uabca(xbc_len)

      DOUBLE PRECISION Energs_init(4), Energs_fin(4),Energ_Err(4)
      DOUBLE PRECISION Norms(3),AvgNorms(3)
      DOUBLE PRECISION tlist(tlist_len)
      integer i
      DOUBLE PRECISION currenttime
      
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      

      !get initial conditions at all reconstruction points    
      call ForcedSolution(xbc,xbc_len,dx,tstart,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,hbccub_init,Gbccub_init,ubccub_init,cubic_len)
     


      do i = 1,cubic_len
        hbccub_fin(i)= hbccub_init(i)
        ubccub_fin(i)= ubccub_init(i)
        Gbccub_fin(i)= Gbccub_init(i)
      end do
      
     
      call NumericalSolveTSPrint(tstart,tend,
     . ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,tolcalc,dx,dt,
     . n_GhstCells,xbc_len,cubic_len,
     . xbc,hbccub_init,Gbccub_init,ubccub_init,
     . hbccub_fin,Gbccub_fin,ubccub_fin,
     . currenttime,Energs_init, Energs_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len)


      call ForcedSolution(xbc,xbc_len,dx,currenttime,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,hbccub_fina,Gbccub_fina,ubccub_fina,cubic_len)  

     
     
      call CubicToCellAverage(cubic_len,xbc_len,hbccub_fin,
     .   dx,habc)
      call CubicToCellAverage(cubic_len,xbc_len,hbccub_fina,
     .   dx,habca)
      call CubicToCellAverage(cubic_len,xbc_len,Gbccub_fin,
     .   dx,Gabc)
      call CubicToCellAverage(cubic_len,xbc_len,Gbccub_fina,
     .   dx,Gabca)      
      call CubicToCellAverage(cubic_len,xbc_len,ubccub_fin,
     .   dx,uabc)
      call CubicToCellAverage(cubic_len,xbc_len,ubccub_fina,
     .   dx,uabca)

      ! calculate norms
      call L2Norm(xbc_len,habca,habc,AvgNorms(1)) 
      call L2Norm(xbc_len,Gabca,Gabc,AvgNorms(2)) 
      call L2Norm(xbc_len,uabca,uabc,AvgNorms(3)) 
 
     
      ! calculate norms
      call L2Norm(cubic_len,hbccub_fina,hbccub_fin,Norms(1)) 
      call L2Norm(cubic_len,Gbccub_fina,Gbccub_fin,Norms(2)) 
      call L2Norm(cubic_len,ubccub_fina,ubccub_fin,Norms(3)) 
      
      
      
      
      open(11, file = ExpWdir//'xhuGInit.dat') 
      do i = 1,xbc_len
         write(11,*) xbc(i) -0.5*dx, hbccub_init(4*i -3), 
     .    ubccub_init(4*i -3) , Gbccub_init(4*i -3)
         write(11,*) xbc(i) - dx/6.0, hbccub_init(4*i -2), 
     .    ubccub_init(4*i -2) , Gbccub_init(4*i -2)
         write(11,*) xbc(i) + dx/6.0, hbccub_init(4*i -1), 
     .    ubccub_init(4*i -1) , Gbccub_init(4*i -1)
         write(11,*) xbc(i) + 0.5*dx, hbccub_init(4*i), 
     .    ubccub_init(4*i) , Gbccub_init(4*i)
      end do
      
      close(11)
      
      open(12, file = ExpWdir//'xhuGFin.dat') 
      do i = 1,xbc_len
         write(12,*) xbc(i) -0.5*dx, hbccub_fin(4*i -3), 
     .    ubccub_fin(4*i -3) , Gbccub_fin(4*i -3)
         write(12,*) xbc(i) - dx/6.0, hbccub_fin(4*i -2), 
     .    ubccub_fin(4*i -2) , Gbccub_fin(4*i -2)
         write(12,*) xbc(i) + dx/6.0, hbccub_fin(4*i -1), 
     .    ubccub_fin(4*i -1) , Gbccub_fin(4*i -1)
         write(12,*) xbc(i) + 0.5*dx, hbccub_fin(4*i), 
     .    ubccub_fin(4*i) , Gbccub_fin(4*i)
      end do
     
      close(12)
      
      open(12, file = ExpWdir//'xhuGFinA.dat') 
      do i = 1,xbc_len
         write(12,*) xbc(i) -0.5*dx, hbccub_fina(4*i -3), 
     .    ubccub_fina(4*i -3) , Gbccub_fina(4*i -3)
         write(12,*) xbc(i) - dx/6.0, hbccub_fina(4*i -2), 
     .    ubccub_fina(4*i -2) , Gbccub_fina(4*i -2)
         write(12,*) xbc(i) + dx/6.0, hbccub_fina(4*i -1), 
     .    ubccub_fina(4*i -1) , Gbccub_fina(4*i -1)
        write(12,*) xbc(i) + 0.5*dx, hbccub_fina(4*i), 
     .    ubccub_fina(4*i) , Gbccub_fina(4*i)
      end do
      
      close(12)

     
      open(14, file = ExpWdir//'Params.dat') 
      write(14,*) 'xstart',xstart
      write(14,*) 'xend',xbc(xbc_len)
      write(14,*) 'x_len',xbc_len - 2*n_GhstCells
      write(14,*) 'n_GhstCells',n_GhstCells
      write(14,*) 'xbc_len',xbc_len
      write(14,*) 'dx' , dx
      write(14,*) 'tstart', tstart
      write(14,*) 'tend',tend 
      write(14,*) 'actual_end_time', currenttime
      write(14,*) 'dt' , dt
      write(14,*) 'tolcalc' , tolcalc
      write(14,*) 'gravity' , ga
      write(14,*) 'a0' , a0
      write(14,*) 'a1' , a1
      write(14,*) 'a2' , a2
      write(14,*) 'a3' , a3      
      write(14,*) 'a4' , a4
      write(14,*) 'a5' , a5
      write(14,*) 'b1a6' , b1a6
      write(14,*) 'b1a7' , b1a7   
      write(14,*) 'b2a6' , b2a6
      write(14,*) 'b2a7' , b2a7
      close(14)
      
      write(AvgNormFile,*)  dx, AvgNorms(1), AvgNorms(2), AvgNorms(3)
      write(NormFile,*)  dx, Norms(1), Norms(2), Norms(3)
        
      
      end
      
      
      
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,AvgNormFile,EnergFile,tlist_len
      PARAMETER(wdirlen= 300,NormFile=1,AvgNormFile=2,
     . EnergFile=3,tlist_len=30)
      
      INTEGER effeclenwdir,n_GhstCells,lowestresx,x_len,xbc_len,
     . cubic_len, expi
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) strdiri
      
      DOUBLE PRECISION ga,beta2,beta1,xstart,xend,dx,
     . tolcalc,alpha,Cr,maxwavespeed,dt,tend,tstart,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7
     
      DOUBLE PRECISION tlist(tlist_len)
      
      wdir = "../Results/Validation/Forced/Run/1sNoLim/"
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      CALL SYSTEM('rm -rf '//wdir(1:effeclenwdir))
      CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir)) 
      
      open(EnergFile, file = wdir(1:effeclenwdir)//'Energy.dat') 
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 
      open(AvgNormFile, file = wdir(1:effeclenwdir)//'AvgNorms.dat') 
      
      tolcalc = 10.0**(-15)
      
      a0= 1
      a1= 0.5
      a2= 5.0
      a3= 20
      a4= 0.3
      a5= 0.0
      b1a6= 0.0
      b1a7= 2.0/15.0
      b2a6= 0.0
      b2a7= 2.0/15.0
      
      
      ga = 9.81

      
      xstart = -100d0
      xend = 100d0
      
      tstart = 0d0
      tend = 1.0
      
      n_GhstCells = 6
      lowestresx = 100
      
      if  (dabs(b1a7) < 10d0**(-10))  then
         alpha = 1d0
      else
         alpha = max(1d0,b2a7 / (b1a7))
      end if
      
      call EqualSpaced(tstart,tend,tlist_len,0.5, tlist)
      
      do expi = 0,10
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells
         cubic_len = 4*xbc_len 

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*(a0 + a1))
         dt  = (Cr / maxwavespeed) *dx

         call SingleForced(a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . ga,xstart,xbc_len,n_GhstCells,cubic_len,dx,tstart,tend,
     . dt,tolcalc,NormFile,AvgNormFile,EnergFile,tlist,tlist_len,
     . wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)

         


      end do
      
      close(EnergFile)
      close(NormFile)
      
     
      end
 
