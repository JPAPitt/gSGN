      subroutine SingleSWWEDB(h0,h1,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,cubic_len,dx,tstart,tend,
     . dt,tolcalc,NormFile,AvgNormFile,EnergFile,tlist,tlist_len,
     . ExpWdir,ExpWdir_len)
     
      implicit none
      DOUBLE PRECISION h0,h1,ga,beta1,beta2,xstart,dx,
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
     . , Energs_fina(4) , Energ_ErrA(4) 
      DOUBLE PRECISION Norms(3),AvgNorms(3)
      DOUBLE PRECISION tlist(tlist_len)
      integer i
      DOUBLE PRECISION currenttime
      
      
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      

      !get initial conditions at all reconstruction points     
      call Dambreak(xbc,xbc_len,tstart,ga,h0,h1,
     . dx,hbccub_init,ubccub_init,Gbccub_init,cubic_len)


      do i = 1,cubic_len
        hbccub_fin(i)= hbccub_init(i)
        ubccub_fin(i)= ubccub_init(i)
        Gbccub_fin(i)= Gbccub_init(i)
      end do
      
      
      call NumericalSolveTSPrint(tstart,tend,
     . ga,beta1,beta2,tolcalc,dx,dt,n_GhstCells,xbc_len,cubic_len,
     . xbc,hbccub_init,Gbccub_init,ubccub_init,
     . hbccub_fin,Gbccub_fin,ubccub_fin,
     . currenttime,Energs_init, Energs_fin,
     . tlist,tlist_len,ExpWdir,ExpWdir_len)
  
      !get initial conditions at all reconstruction points     
      call Dambreak(xbc,xbc_len,currenttime,ga,h0,h1,
     . dx,hbccub_fina,ubccub_fina,Gbccub_fina,cubic_len)  
     
      ! get energy in analytic solution
      call TotalEnergy(xbc_len,cubic_len,hbccub_fina,ubccub_fina,
     . Gbccub_fina,ga, beta1,beta2,n_GhstCells,dx,Energs_fina)   
     
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
      
      ! Conservation Norm Tests  
      do i = 1,4
         !if denominator small just return absolute error
         if (dabs(Energs_init(i)) .lt. 10d0**(-10)) then
            Energ_Err(i) = dabs(Energs_fin(i) - Energs_init(i))
         else
            Energ_Err(i) = dabs(Energs_fin(i) - Energs_init(i))/
     .       dabs(Energs_init(i))
         end if
         
         if (dabs(Energs_fina(i)) .lt. 10d0**(-10)) then
            Energ_ErrA(i) = dabs(Energs_fin(i) - Energs_fina(i))
         else
            Energ_ErrA(i) = dabs(Energs_fin(i) - Energs_fina(i))/
     .       dabs(Energs_fina(i))
         end if
      end do 
      
      
      
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
      
      open(12, file = ExpWdir//'CAxhuGFin.dat') 
      do i = 1,xbc_len
         write(12,*) xbc(i), habc(i), uabc(i), Gabc(i)
      end do
      close(12)

      open(12, file = ExpWdir//'CAxhuGFinA.dat') 
      do i = 1,xbc_len
         write(12,*) xbc(i), habca(i), uabca(i), Gabca(i)
      end do
      close(12)
      
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


      
      open(13, file = ExpWdir//'Energy.dat')       
      write (13,*) 'When , h , G ,  uh , Energy'
      write(13,*) 'Initial ',Energs_init(1),Energs_init(2),
     .   Energs_init(3),Energs_init(4)
      write(13,*) 'End ',Energs_fin(1),Energs_fin(2),
     .   Energs_fin(3),Energs_fin(4)
      close(13)
     
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
      write(14,*) 'h0' , h0
      write(14,*) 'h1' , h1
      write(14,*) 'beta1' , beta1
      write(14,*) 'beta2' , beta2
      close(14)
      
      write(AvgNormFile,*)  dx, AvgNorms(1), AvgNorms(2), AvgNorms(3)
      write(NormFile,*)  dx, Norms(1), Norms(2), Norms(3)
      write(EnergFile,*) dx, Energ_Err(1), Energ_Err(2),
     .  Energ_Err(3), Energ_Err(4),Energ_ErrA(1), Energ_ErrA(2),
     .  Energ_ErrA(3), Energ_ErrA(4)
        
      
      end subroutine
      
      
      
      program main
         
      implicit none
   
      Integer wdirlen,NormFile,AvgNormFile,EnergFile,tlist_len
      PARAMETER(wdirlen= 300,NormFile=1,AvgNormFile=2,
     . EnergFile=3,tlist_len=30)
      
      INTEGER effeclenwdir,n_GhstCells,lowestresx,x_len,xbc_len,
     . cubic_len, expi
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) strdiri
      
      DOUBLE PRECISION h0,h1,ga,beta2,beta1,xstart,xend,dx,
     . tolcalc,alpha,Cr,maxwavespeed,dt,tend,tstart
     
      DOUBLE PRECISION tlist(tlist_len)
      
      wdir = "../Results/Validation/Run/iSGN/10s/"
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      CALL SYSTEM('rm -rf '//wdir(1:effeclenwdir))
      CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir)) 
      
      open(EnergFile, file = wdir(1:effeclenwdir)//'Energy.dat') 
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 
      open(AvgNormFile, file = wdir(1:effeclenwdir)//'AvgNorms.dat') 
      
      tolcalc = 10.0**(-12)
      
      h0 = 2.0
      h1 = 1.0
      ga = 9.81

      beta1 = 0.8
      beta2 = 0.133333333
      
      xstart = -100d0
      xend = 100d0
      
      tstart = 0d0
      tend = 10.0
      
      n_GhstCells = 6
      lowestresx = 100
      
      if  (dabs(beta1) < 10d0**(-10))  then
         alpha = 1d0
      else
         alpha = max(1d0,beta2 / (beta1))
      end if
      
      call EqualSpaced(tstart,tend,tlist_len,0.5, tlist)
      
      do expi = 6,6
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells
         cubic_len = 4*xbc_len 

         dx = (xend - xstart) / (x_len -1)
         Cr = 0.5
         maxwavespeed = dsqrt(alpha*ga*max(h0,h1))
         dt  = (Cr / maxwavespeed) *dx




         call  SingleSWWEDB(h0,h1,ga,beta1,beta2,
     . xstart,xbc_len,n_GhstCells,cubic_len,dx,tstart,tend,
     . dt,tolcalc,NormFile,AvgNormFile,EnergFile,tlist,tlist_len,
     . wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)
         


      end do
      
      close(EnergFile)
      close(NormFile)
      close(AvgNormFile)
      
     
      end program
 
