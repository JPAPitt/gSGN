c ====================================================================================
c Module of fortran subroutines use FDVM2 to solve generalised Serre - Green -Naghdi equations (gSGN)
c it takes initial conditions and returns solutions Q(x,t) where t > tend (Q = [h,G])
c the subroutine also calculates total conserved quantities initially and in the final solution
c ====================================================================================

c=====================================
c Numerical solver program, evolving system through time and producing Q(x,t) where t > tend
c The solver uses arrays xbc,hbc_init,Gbc_init and ubc_init
c to produce: 
c initial energy - Energs_init
c solution for h at final time hbc_fin
c solution for G at final time Gbc_fin
c solution for u at final time ubc_fin
c currenttime -time of produced solution currenttime = tstart + dt*number of timesteps performed
c Notes:
c 1. ubc_init only needs to have u defined at ghost cells
c function will update interior using up to date values of h and G.
c
c 2. Solver assumes beta's are constant and ghost cell values are constant- constant dirichlet boundary conditions
c=====================================
      subroutine NumericalSolveTSPrint(tstart,tend,
     . ga,theta,dx,dt,n_GhstCells,xbc_len,
     . xbc,hbc_init,uhbc_init,ubc_init,
     . currenttime,hbc_fin,uhbc_fin,ubc_fin,
     . Energ_Init, Energ_Fin,
     . tlist,tlist_len,ExpWdir,expwdirlen)
     
     
      implicit none
      
      integer n_GhstCells,xbc_len,expwdirlen,tlist_len
      CHARACTER(len=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,
     . theta,dx,dt,currenttime
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . uhbc_init(xbc_len),
     . ubc_init(xbc_len),hbc_fin(xbc_len),
     . uhbc_fin(xbc_len), ubc_fin(xbc_len)
     
      DOUBLE PRECISION tlist(tlist_len)
      
      DOUBLE PRECISION Energ_Init(3), Energ_Fin(3)
      
      integer i,ileft,iright,dtcount,filecount
      CHARACTER(len=2) strct
      
      !initial time
      currenttime  = tstart
      dtcount = 0
      filecount = 1
      
      !loop over and set hbc_fin,Gbc_fin to initial conditions
      ileft = 1
      iright = xbc_len
      do i = ileft,iright
         hbc_fin(i) = hbc_init(i) 
         uhbc_fin(i) = uhbc_init(i) 
         ubc_fin(i) = ubc_init(i)
      end do
      
      !calculate initial Energies
      call Getufromhuh(xbc_len,hbc_fin,uhbc_fin,ubc_fin,
     . n_GhstCells)
     
      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,uhbc_fin,ga,
     . n_GhstCells,dx,Energ_Init)

      !evolve the system through time
      do while (currenttime  .LT. tend )   
      
         if ( dabs(currenttime - tlist(filecount)) .LT. 0.9*dt  )  then
              
            write (strct,'(I2)') filecount
            open(9, file = ExpWdir// strct //'.dat') 
            do i = 1,xbc_len
               write(9,*) currenttime,xbc(i),hbc_fin(i),
     .            uhbc_fin(i),ubc_fin(i)
            end do
            close(9)
            
            filecount = filecount + 1
         end if

         
         call EvolveStepWrap(xbc_len,hbc_fin,uhbc_fin,ubc_fin,ga,
     .    theta,n_GhstCells,dx,dt)
     
         currenttime  = currenttime  + dt
         dtcount = dtcount + 1
         print *, 'Current Time : ', currenttime 
      end do
       
      !calculate end energies  
      call Getufromhuh(xbc_len,hbc_fin,uhbc_fin,ubc_fin,n_GhstCells)
      
      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,uhbc_fin,ga,
     . n_GhstCells,dx,Energ_Fin)
                 
      end


c  ********************************************************************************
c  Part 2 : Finite Difference Solver
c  Functions thet solve finite difference matrix Au = G
c  for u. A is tridiagonal, so we solve using Thomas algorithm
c
c
c
c ********************************************************************************

c =====================
c Function that builds the three diagonals of the matrix produced by the finite difference approximation to elliptic equation
c It then uses a matrix solve on the diagonals, and the value of G to solve for u (for gSGN equations)
c ==================      
      subroutine Getufromhuh(xbc_len,hbc,uhbc,ubc,n_GhstCells)
      
      implicit none
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),uhbc(xbc_len),ubc(xbc_len)

      integer i
      
      do i=n_GhstCells+1,xbc_len - n_GhstCells 
         ubc(i) = uhbc(i) / hbc(i)
      end do 

      
      end


c  ********************************************************************************
c  Part 3 : Evolution
c  Functions that evolve h^n, G^n  -> h^{n+1}, G^{n+1}
c  Using u solve, Kurganov method and Rk2 step
c
c
c ********************************************************************************


c ====
c Evolve Step Wrap, is the wrapper function that takes h^n, G^n and produces
c our second order approximation to h^n+1 and G^n+1, using RK time stepping
c ====
      subroutine EvolveStepWrap(xbc_len,hbc,uhbc,ubc,ga,theta
     . ,n_GhstCells,dx,dt)
          
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),uhbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION ga,theta,dx,dt
     
      !local variables
      DOUBLE PRECISION hpbc(xbc_len),uhpbc(xbc_len),
     . hppbc(xbc_len),uhppbc(xbc_len)
     
      integer i
         
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call Getufromhuh(xbc_len,hbc,uhbc,ubc,n_GhstCells)
            
      !Update cell averages (first order approximation to h^{n+1}, G^{n+1})
      call SingleEulerStep(xbc_len,hbc,uhbc,ubc,ga,
     . theta,n_GhstCells,dt,dx,hpbc,uhpbc)
     
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call Getufromhuh(xbc_len,hpbc,uhpbc,ubc,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+2}, G^{n+2})
      call SingleEulerStep(xbc_len,hpbc,uhpbc,ubc,ga,
     . theta,n_GhstCells,dt,dx,hppbc,uhppbc)
      
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      do i= 1,xbc_len
         hbc(i) = ( hbc(i ) + hppbc(i))/2d0
         uhbc(i) = ( uhbc(i ) + uhppbc(i))/2d0
      end do
          
      
      end 
    

c ====
c SingleEulerStep
c produces new cell averages of h,G using forward Euler step,
c with flux approximated using Kurganov's method
c ====
      subroutine SingleEulerStep(xbc_len,hbc,uhbc,ubc,ga,
     . theta,n_GhstCells,dt,dx,hpbc,uhpbc)
     
     
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION dt,dx,ga,theta
      DOUBLE PRECISION hbc(xbc_len),uhbc(xbc_len),ubc(xbc_len),
     . hpbc(xbc_len),uhpbc(xbc_len)
     
     
      DOUBLE PRECISION cdhi,cduhi, fih,fiuh,foh,fouh
      
      integer i,ileft,iright
      
      !we assume boundary conditions are constant
      ileft = 1
      iright = n_GhstCells
      do i = ileft,iright
         hpbc(i) = hbc(i)
         uhpbc(i) = uhbc(i)
         
         hpbc(xbc_len - n_GhstCells + i) = 
     .      hbc(xbc_len - n_GhstCells + i)
         uhpbc(xbc_len - n_GhstCells + i) = 
     .      uhbc(xbc_len - n_GhstCells + i)
      end do
      
      !Now we update interior, first calculating flux across left interior boundary 
      !then loop over cells to get flux in/out and thus updating cell average values
      
      !Do left boundary
      i = n_GhstCells
      
      !initial gradient of h,G across cell i
      call ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta,cdhi)
      call ReconLinLimGrad(uhbc(i-1),uhbc(i),uhbc(i+1),theta,cduhi)
      call ReconLinLimGrad(ubc(i-1),ubc(i),ubc(i+1),theta,cdui)
      
      !calculates foh,foG which is flux across x_{i +1/2}
      ! it also updates cdhi,cdGi to be gradient across cell i + 1
      call Fluxxiph(xbc_len,hbc,uhbc,ubc,ga,theta,i,
     . foh,fouh,cdhi,cduhi,cdui)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiuh = fouh
      !loop over interior cells (do not update ghost cells)
      ileft = n_GhstCells + 1
      iright = xbc_len - n_GhstCells
      do i = ileft, iright
      
         !calculates foh,foG which is flux across x_{i +1/2}
         ! it also updates cdhi,cdGi to be gradient across cell i + 1
         call Fluxxiph(xbc_len,hbc,uhbc,ubc,ga,theta,i,
     . foh,fouh,cdhi,cduhi,cdui)
             
         hpbc(i) = hbc(i) - dt*(foh - fih)/dx 
     
         uhpbc(i) = uhbc(i) - dt*(fouh - fiuh)/dx 

         !flux out becomes flux in on next cell
         fih = foh
         fiuh = fouh
         
      end do 
           
      end
      
c subroutine that given arrays, and i calculates flux across x_{i+1/2} for gSGN equations  
c note that - cdhi,cdGi is inout, on the way in cdhi,cdGi are our approximation to the gradient
c of h and G across cell i, on the way out cdhi,cdGi are our approximation to the gradient of h and G
c across cell i +1 (which will be cell i as the method moves to the next cell)    
      subroutine Fluxxiph(xbc_len,hbc,uhbc,ubc,ga,theta,
     . i,foh,fouh,cdhi,cduhi,cdui)
     
       integer i,xbc_len
       DOUBLE PRECISION hbc(xbc_len),uhbc(xbc_len),ubc(xbc_len)
       DOUBLE PRECISION ga,theta,
     . cdhi,cduhi,cdui,foh,fouh
     
       DOUBLE PRECISION feluh,felh,feruh,ferh,sr,sl,isrmsl,
     . hir,uhir,uir,uhip1l,uip1l,hip1l
                
c     centered approximation
      !use slopes from last calculation         
      hir = hbc(i) + cdhi/2d0
      uhir = uhbc(i) + cduhi/2d0
      uir = ubc(i) + cdui/2d0  

      !reconstruct values on right side of edge x_{i+1/2}
      call ReconLinLimGrad(hbc(i),hbc(i+1),hbc(i+2),theta,cdhip1)
      call ReconLinLimGrad(uhbc(i),uhbc(i+1),uhbc(i+2),theta,cduhip1)
      call ReconLinLimGrad(ubc(i),ubc(i+1),ubc(i+2),theta,cduip1)
      
      hip1l = hbc(i + 1) - cdhip1/2d0
      uhip1l = uhbc(i + 1) - cduhip1/2d0
      uip1l = ubc(i + 1) - cduip1/2d0
      
      sl  = min(0d0, uir - dsqrt(ga*hir) ,
     . uip1l - dsqrt(ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(ga*hir) ,
     . uip1l + dsqrt(ga*hip1l)  )
      
      !left and right flux
      felh = uhir
      feluh = uir*uhir + ga*(hir**2)/2d0 
      
      ferh = uhip1l
      feruh = uip1l*uhip1l + ga*(hip1l**2)/2d0 
     
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if 
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      fouh = isrmsl*(sr*feluh - sl*feruh + sl*sr*(uhip1l - uhir))
      
      !return gradient of h,G across cell i + 1 (for next iteration)
      cduhi = cduhip1
      cdhi = cdhip1   
      cdui = cduip1 
      end

c ====
c ReconLinLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconLinLimGrad(qim1,qi,qip1,theta,cdq)
      DOUBLE PRECISION qim1,qi,qip1,theta,cdq,
     . fdq,mdq,bdq
      
      fdq = qip1 - qi
      mdq = 0.5*(qip1 - qim1)
      bdq = qi - qim1
      
      call minmod(theta*fdq,mdq,theta*bdq,cdq)
      
      end

c ===
c minmod is used for limiting and returns smallest size argument if all arguments have same sign,
c otherwise it returns 0.
c ===      
      subroutine minmod(a,b,c,d)
      implicit none
      DOUBLE PRECISION a,b,c, d
      
      if ((a .gt. 0d0) .and. (b .gt. 0d0) .and. (c .gt. 0d0)) then
         d = min(a,b,c)
      else if ((a .lt. 0d0) .and. (b .lt. 0d0) .and. (c .lt. 0d0)) then
         d = max(a,b,c)
      else
         d = 0d0
      end if
      
      end 

c  ********************************************************************************
c  Part 4 : Analysis Functions
c  Functions that do analyses on solutions
c  1. For Conservation -h, uh, G and Energy  (was at least second order accurate (FD solve, holding it back to second order for G))
c  
c
c
c ********************************************************************************
      
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


c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      subroutine AllEnergiesIntegralCell(xbc_len,h,u,uh,ga,j
     . ,dx,CellEnergies)
      
      integer j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len),uh(xbc_len)
      DOUBLE PRECISION dx,ga
      DOUBLE PRECISION CellEnergies(3)
      
      integer i
      DOUBLE PRECISION fGPe(3),sGPe(3),tGPe(3)
      
      DOUBLE PRECISION GPmxj,hGP,uhGP,uGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5), uhCoeff(5)
      
      call QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      call QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      call QuarticInterp(uh(j-2),uh(j-1),uh(j),uh(j+1),
     . uh(j+2),dx,uhCoeff)
      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalxj(uhCoeff,GPmxj,uhGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      
      fGPe(1) = hGP
      fGPe(2) = uhGP
      fGPe(3) = (uhGP*uGP + ga*hGP**2 )/2d0

      !second gauss point
      GPmxj = 0.0 
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalxj(uhCoeff,GPmxj,uhGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      
      sGPe(1) = hGP
      sGPe(2) = uhGP
      sGPe(3) = (uhGP*uGP + ga*hGP**2 )/2d0
      
      !third gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalxj(uhCoeff,GPmxj,uhGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      
      tGPe(1) = hGP
      tGPe(2) = uhGP
      tGPe(3) = (uhGP*uGP + ga*hGP**2 )/2d0
      
      !weight the values at gauss points to get approximate integral over cell
      do i = 1,3
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .+ (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      end do
      
      end
      
      !Function to sum all energies
      subroutine TotalEnergy(xbc_len,hbc,ubc,uhbc,ga,
     . n_GhstCells,dx,TotEnergVals)
      
      integer xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),uhbc(xbc_len)
      DOUBLE PRECISION dx,ga
      DOUBLE PRECISION TotEnergVals(3)
      
      DOUBLE PRECISION CellEnergVals(3)
      integer i,j
            
      !running totals for energy values, start at 0
      do i = 1,3
         TotEnergVals(i) = 0.0
      end do
      
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      do j= n_GhstCells + 1, xbc_len - n_GhstCells
         call  AllEnergiesIntegralCell(xbc_len,hbc,ubc,uhbc,ga,
     .         j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,3
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end

      


