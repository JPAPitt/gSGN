c ===
c Function to make time series of Energies
c ===
      subroutine EnergyTimeSeries(wdir,wdir_len,paramfilename,
     . paramfilename_len,ts_len,ts,Energys,dx)
      INTEGER ts_len,wdir_len,paramfilename_len
      CHARACTER(len=wdir_len) wdir
      CHARACTER(len=paramfilename_len) paramfilename
      CHARACTER(len=2) strct
      
      DOUBLE PRECISION ts(ts_len), Energys(ts_len,4)
      
      DOUBLE PRECISION dx,ga,beta1,beta2
      INTEGER nGhstcell,xbc_len
      
      integer j
      
      call  ReadParamFile(wdir//paramfilename,
     .   wdir_len + paramfilename_len,nGhstcell,
     .   xbc_len,dx,ga,beta1,beta2)
      do j = 1,ts_len -1
         write (strct,'(I2)') j
         
         call ReadResFile(wdir//strct //'.dat',wdir_len +6,
     .   xbc_len,nGhstcell,ga, beta1,beta2,dx,ts(j),Energys(j,:))
     
      end do   
      
      call ReadResFile(wdir//'End.dat',wdir_len +7,
     .   xbc_len,nGhstcell,ga, beta1,beta2,dx,ts(j),Energys(j,:))
      end

c ===
c Function to
c
c ===
      subroutine ReadParamFile(FileLoc,FileLoc_len,nGhstcell,
     . xbc_len,dx,ga,beta1,beta2)
      integer FileLoc_len,xbc_len,nGhstcell
      CHARACTER(len=FileLoc_len) FileLoc
      DOUBLE PRECISION ga,beta1,beta2,dx
      
      CHARACTER(len=1) str1
      
      open(1, file =FileLoc, action = 'read')
      do j = 1,16
      if (j .eq. 4) then
         READ(1,*) str1,nGhstcell
      else if (j .eq. 5) then
         READ(1,*) str1, xbc_len
      else if (j .eq. 6) then
         READ(1,*) str1,dx
      else if (j .eq. 12) then
         READ(1,*) str1, ga
      else if (j .eq. 15) then
         READ(1,*) str1, beta1
      else if (j .eq. 16) then
         READ(1,*) str1,beta2
      else
         READ(1,*) str1
      end if
     
      end do
      close(1)
            
      end



c===
c function to read output files
c
c===
      subroutine ReadResFile(FileLoc,FileLoc_len,xbc_len,nGhst,ga,
     .   beta1,beta2,dx,t,Energy)
      integer FileLoc_len,xbc_len,nGhst
      CHARACTER(len=FileLoc_len) FileLoc
      DOUBLE PRECISION t,ga,beta1,beta2,dx
      DOUBLE PRECISION x(xbc_len),h(xbc_len),G(xbc_len),u(xbc_len)
      DOUBLE PRECISION Energy(4)
      
      open(1, file =FileLoc, action = 'read')
      do j = 1,  xbc_len
         READ(1,*) t,x(j),h(j),G(j),u(j)
      end do
      close(1)
   
      call TotalEnergy(xbc_len,h,u,G,ga,beta1,beta2,
     . nGhst,dx,Energy)
      end


c =====
c Function to generate quartic coefficients using
c f(x_{j-2}),f(x_{j-1}),f(x_{j}),f(x_{j+1}),f(x_{j+2})
c =====      
      subroutine QuarticInterp(qjm2,qjm1,qj,qjp1,qjp2,dx,QuartCoeff)
     
      DOUBLE PRECISION qjm2,qjm1,qj,qjp1,qjp2,dx
      DOUBLE PRECISION QuartCoeff(5)
      
      QuartCoeff(1) = (qjp2 - 4*qjp1 + 6*qj - 4*qjm1 + qjm2) /
     . (24* (dx**4))
      QuartCoeff(2) = (qjp2 - 2*qjp1 + 2*qjm1 - qjm2) /
     . (12* (dx**3))  
      QuartCoeff(3) = (-qjp2 + 16*qjp1 - 30*qj + 16*qjm1 - qjm2) /
     . (24* (dx**2))  
      QuartCoeff(4) = (-qjp2 + 8*qjp1 - 8*qjm1 + qjm2) /
     . (12*dx)
      QuartCoeff(5) = qj
      end


c ========================
c Functions to evaluate Quartic at x, with quartic centered around x_j
c xmxj = x - x_j 
c ====================      
      subroutine QuarticCoeffEvalxj(QuarticCoeff,xmxj,qatxj)
      
      DOUBLE PRECISION QuarticCoeff(5)
      DOUBLE PRECISION xmxj,qatxj
      
      qatxj = QuarticCoeff(1)*xmxj**4 + QuarticCoeff(2)*xmxj**3 +
     . QuarticCoeff(3)*xmxj**2 + QuarticCoeff(4)*xmxj +
     . QuarticCoeff(5)
     
      end
      
      subroutine QuarticCoeffEvalGradxj(QuarticCoeff,xmxj,dqatxj)
      
      DOUBLE PRECISION QuarticCoeff(5)
      DOUBLE PRECISION xmxj,dqatxj
      
      dqatxj = 4*QuarticCoeff(1)*xmxj**3 + 3*QuarticCoeff(2)*xmxj**2 +
     . 2*QuarticCoeff(3)*xmxj + QuarticCoeff(4)
     
      end

c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      subroutine AllEnergiesIntegralCell(xbc_len,h,u,G,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      integer j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len),G(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2,beta1p2o3
      DOUBLE PRECISION CellEnergies(4)
      
      integer i
      DOUBLE PRECISION fGPe(4),sGPe(4),tGPe(4)
      
      DOUBLE PRECISION GPmxj,hGP,GGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5), GCoeff(5)
      
      if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
         beta1p2o3 = 0
      else
         beta1p2o3 = beta1 + 2d0/3d0
      end if
      
      call QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      call QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      call QuarticInterp(G(j-2),G(j-1),G(j),G(j+1),G(j+2),dx,GCoeff)
      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      fGPe(1) = hGP
      fGPe(2) = GGP
      fGPe(3) = hGP*uGP
      fGPe(4) = (hGP*uGP**2 + (1d0/2d0)*beta1p2o3*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !second gauss point
      GPmxj = 0.0 
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      sGPe(1) = hGP
      sGPe(2) = GGP
      sGPe(3) = hGP*uGP
      sGPe(4) = (hGP*uGP**2 + (1d0/2d0)*beta1p2o3*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !third gauss point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + (1d0/2d0)*beta1p2o3*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !weight the values at gauss points to get approximate integral over cell
      do i = 1,4
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .+ (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      end do
      
      end
      
      !Function to sum all energies
      subroutine TotalEnergy(xbc_len,hbc,ubc,Gbc,ga,beta1,beta2,
     . n_GhstCells,dx,TotEnergVals)
      
      integer xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),Gbc(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION TotEnergVals(4)
      
      DOUBLE PRECISION CellEnergVals(4)
      integer i,j
            
      !running totals for energy values, start at 0
      do i = 1,4
         TotEnergVals(i) = 0.0
      end do
      
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      do j= n_GhstCells + 1, xbc_len - n_GhstCells
         call  AllEnergiesIntegralCell(xbc_len,hbc,ubc,Gbc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,4
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end
      
c get length of string without trailing whitespace      
      subroutine LenTrim(String,n,resultval)
      CHARACTER(len=n) String

      integer n, resultval

      resultval = 1

      do while ((String(resultval:resultval) .NE. ' ' ) .AND.
     . (resultval .LE. n )  )
         resultval = resultval + 1
      end do

      resultval = resultval - 1
      end 
      
      
      program main
         
      implicit none
      Integer wdirlen,ts_len
      PARAMETER(wdirlen= 200,ts_len=35)
      
      INTEGER j,effeclenwdir,i
      DOUBLE PRECISION ts(ts_len),Energys(ts_len,4),dx
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) stri
      
      wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/" //
     . "1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/" //
     . "Validation/AnalyticSolutions/DBSWWE/"     
      
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      do i = 0,12
         write (stri,'(I2.2)') i
         call EnergyTimeSeries(wdir(1:effeclenwdir) // stri // '/',
     .      effeclenwdir +3, "Params.dat", 10,ts_len,ts,Energys,dx)
     
     
         open(2, file = 'DB_DX'//stri //'.dat') 
         do j = 1,ts_len
            write(2,*) ts(j),dx,Energys(j,1),Energys(j,2),Energys(j,3),
     .         Energys(j,4) 
         end do
      
      end do
      
      

      
      end

      


