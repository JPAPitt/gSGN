c ====================================================================================
c Module of fortran subroutines use FDVM2 to solve generalised Serre - Green -Naghdi equations (gSGN)
c it takes initial conditions and returns solutions Q(x,t) where t > tend (Q = [h,G])
c the subroutine also calculates total conserved quantities initially and in the final solution
c ====================================================================================

      
c **********************************************
c Numerical solve that prints time series
c
c
c *******************************************
      subroutine NumericalSolveTSPrint(tstart,tend,
     . ga,beta2,theta,dx,dt,n_GhstCells,xbc,xbc_len,
     . hbc_init,Gbc_init,ubc_init,
     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
     . dtsep,ExpWdir,expwdirlen)
     
      implicit none
      
      integer n_GhstCells,xbc_len,expwdirlen,dtsep
      CHARACTER(len=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,beta2,theta,dx,dt,
     . currenttime
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),Gbc_init(xbc_len),
     . ubc_init(xbc_len),hbc_fin(xbc_len),
     . Gbc_fin(xbc_len), ubc_fin(xbc_len)
     
      
      integer i,ileft,iright,dtcount,filecount
      CHARACTER(len=10) strct
      
      !initial time
      currenttime  = tstart
      dtcount = 0
      filecount = 0
      
      !loop over and set hbc_fin,Gbc_fin to initial conditions
      ileft = 1
      iright = xbc_len
      do i = ileft,iright
         hbc_fin(i) = hbc_init(i) 
         Gbc_fin(i) = Gbc_init(i) 
         ubc_fin(i) = ubc_init(i)
      end do

      
      !calculate initial Energies
      call GetufromhG(xbc_len,hbc_fin,Gbc_fin,ubc_fin,
     . dx,n_GhstCells)

     
      
      !evolve the system through time
      do while (currenttime  .LT. tend )     

         if ((MOD(dtcount,dtsep) .EQ. 0)
     .      .OR. (currenttime + dt  .GE. tend ))  then
            write (strct,'(I10)') filecount
            open(5, file = ExpWdir// strct //'.dat') 
            
            
            do i = 1,xbc_len
               write(5,*) currenttime,xbc(i),hbc_fin(i),
     .            Gbc_fin(i),ubc_fin(i)
            end do
            close(5)
            filecount = filecount + 1
         end if
         
         call EvolveStepWrap(xbc_len,hbc_fin,Gbc_fin,ubc_fin,ga,
     .   beta2,theta,n_GhstCells,dx,dt)
     
         currenttime  = currenttime  + dt
         dtcount = dtcount + 1
         print *, 'Current Time : ', currenttime 
      end do
      
      !calculate end energies  
      call GetufromhG(xbc_len,hbc_fin,Gbc_fin,ubc_fin,
     . dx,n_GhstCells) 
                 
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
      subroutine GetufromhG(xbc_len,hbc,Gbc,ubc,dx,n_GhstCells)
      
      implicit none
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION dx,beta1
                   
      DOUBLE PRECISION subdiag1(xbc_len),
     . diag(xbc_len),
     . supdiag1(xbc_len),
     . RHS(xbc_len)
          
      DOUBLE PRECISION ht1,ht2,dGc,a0
     
      integer i
                  
      !calculate diagonals in interior
      !set RHS B
            
      do i=n_GhstCells+1,xbc_len - n_GhstCells 
      
         a0 = 0.5d0
         
         dGc =dabs(Gbc(i+1) - Gbc(i-1))/(2*dx)
         
         if (dGc .LE. a0) then
            beta1 = 0d0
         else 
            beta1 = -2d0/3d0
         end if 
         
         
         ht1 = (beta1/2d0 + 1d0/3d0)*(hbc(i)**3/(dx*dx))
         ht2 = (3*beta1/2d0 + 1d0)*
     .      (hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1)))
         
         
         subdiag1(i)  = -ht1 + ht2
         diag(i) = hbc(i) + (2d0*ht1)
         supdiag1(i)  = -ht1 - ht2 
                  
         RHS(i) = Gbc(i)
      end do 
      
      !boundary conditions
      !first and last n_GhstCells  x n_GhstCells in tridiagmatrix
      !Should be identity
      do i = 1, n_GhstCells
         ! left
         subdiag1(i) = 0d0
         diag(i) = 1d0
         supdiag1(i) = 0d0
         RHS(i) = ubc(i)
         
         !right
         subdiag1(xbc_len - n_GhstCells + i) = 0d0
         diag(xbc_len - n_GhstCells + i) = 1d0
         supdiag1(xbc_len - n_GhstCells + i)= 0d0

         RHS(xbc_len - n_GhstCells + i) = 
     .      ubc(xbc_len - n_GhstCells + i)
      end do

      call ThomasTriSolve(xbc_len,subdiag1,diag,supdiag1,RHS,ubc)
      
      end


c ================
c ThomasTriSolve solves a tridiagonal matrix equation Ax = b, where A is tridiagonal
c =================   
      subroutine ThomasTriSolve(n,a,b,c,d,x)
      implicit none
      
      integer n
      DOUBLE PRECISION a(n), c(n), b(n), d(n),x(n)
      
      !  --- Local variables ---
      integer i
      DOUBLE PRECISION q
      !  --- Elimination ---
      do i = 2,n
         q = a(i)/b(i - 1)
         b(i) = b(i) - c(i - 1)*q
         d(i) = d(i) - d(i - 1)*q
      end do
      ! --- Backsubstitution ---
      q = d(n)/b(n)
      x(n) = q
      do i = n - 1,1,-1
         q = (d(i) - c(i)*q)/b(i)
         x(i) = q
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
      subroutine EvolveStepWrap(xbc_len,hbc,Gbc,ubc,ga,beta2,theta
     . ,n_GhstCells,dx,dt)
          
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION ga,beta2,theta,dx,dt
     
      !local variables
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len),
     . hppbc(xbc_len),Gppbc(xbc_len)
     
      integer i
         
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(xbc_len,hbc,Gbc,ubc,dx,n_GhstCells)
            
      !Update cell averages (first order approximation to h^{n+1}, G^{n+1})
      call SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,beta2
     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(xbc_len,hpbc,Gpbc,ubc,dx,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+2}, G^{n+2})
      call SingleEulerStep(xbc_len,hpbc,Gpbc,ubc,ga,beta2
     . ,theta,n_GhstCells,dt,dx,hppbc,Gppbc)
      
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      do i= 1,xbc_len
         hbc(i) = ( hbc(i ) + hppbc(i))/2d0
         Gbc(i) = ( Gbc(i ) + Gppbc(i))/2d0
      end do
          
      
      end 


c ====
c SingleEulerStep
c produces new cell averages of h,G using forward Euler step,
c with flux approximated using Kurganov's method
c ====
      subroutine SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,beta2
     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
     
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION dt,dx,ga,beta2,theta
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len),
     . hpbc(xbc_len),Gpbc(xbc_len)
     
     
      DOUBLE PRECISION cdhi,cdGi, fih,fiG,foh,foG
      
      integer i,ileft,iright
      
      !we assume boundary conditions are constant
      ileft = 1
      iright = n_GhstCells
      do i = ileft,iright
         hpbc(i) = hbc(i)
         Gpbc(i) = Gbc(i)
         
         hpbc(xbc_len - n_GhstCells + i) = 
     .      hbc(xbc_len - n_GhstCells + i)
         Gpbc(xbc_len - n_GhstCells + i) = 
     .      Gbc(xbc_len - n_GhstCells + i)
      end do
      
      !Now we update interior, first calculating flux across left interior boundary 
      !then loop over cells to get flux in/out and thus updating cell average values
      
      !Do left boundary
      i = n_GhstCells
      
      !initial gradient of h,G across cell i
      call ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta,cdhi)
      call ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta,cdGi)
      
      !calculates foh,foG which is flux across x_{i +1/2}
      ! it also updates cdhi,cdGi to be gradient across cell i + 1
      call Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta2,theta,dx,i,
     . foh,foG,cdhi,cdGi)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      
      !loop over interior cells (do not update ghost cells)
      ileft = n_GhstCells + 1
      iright = xbc_len - n_GhstCells
      do i = ileft, iright
      
         !calculates foh,foG which is flux across x_{i +1/2}
         ! it also updates cdhi,cdGi to be gradient across cell i + 1
         call Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta2,theta,dx,i,
     . foh,foG,cdhi,cdGi)
             
         hpbc(i) = hbc(i) - dt*(foh - fih)/dx 
     
         Gpbc(i) = Gbc(i) - dt*(foG - fiG)/dx 

         !flux out becomes flux in on next cell
         fih = foh
         fiG = foG
         
      end do 
           
      end
      
c subroutine that given arrays, and i calculates flux across x_{i+1/2} for gSGN equations  
c note that - cdhi,cdGi is inout, on the way in cdhi,cdGi are our approximation to the gradient
c of h and G across cell i, on the way out cdhi,cdGi are our approximation to the gradient of h and G
c across cell i +1 (which will be cell i as the method moves to the next cell)    
      subroutine Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,beta2,theta,dx,i,
     . foh,foG,cdhi,cdGi)
     
       integer i,xbc_len
       DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
       DOUBLE PRECISION ga,beta1,beta2,theta,dx,
     . cdhi,cdGi,foh,foG
     
       DOUBLE PRECISION cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l,
     . dhir,ddhir,dhip1l,ddhip1l, alpha,dGc,a0
           
c     centered approximation
      !use slopes from last calculation         
      hir = hbc(i) + cdhi/2d0
      Gir = Gbc(i) + cdGi/2d0
      uir = (ubc(i+1)+ubc(i))/2d0
      duir = (ubc(i+1) - ubc(i)) /dx
      dhir = (hbc(i+1) - hbc(i)) /dx
      ddhir = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )/(2*dx**2)

      !reconstruct values on right side of edge x_{i+1/2}
      call ReconLinLimGrad(hbc(i),hbc(i+1),hbc(i+2),theta,cdhip1)
      call ReconLinLimGrad(Gbc(i),Gbc(i+1),Gbc(i+2),theta,cdGip1)
      
      hip1l = hbc(i + 1) - cdhip1/2d0
      Gip1l = Gbc(i + 1) - cdGip1/2d0
      uip1l = (ubc(i+1)+ubc(i))/2
      duip1l = (ubc(i+1) - ubc(i)) /dx
      dhip1l = (hbc(i+1) - hbc(i)) /dx
      ddhip1l = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )
     .     /(2*dx**2)
     
      a0 = 0.5d0
      
      dGc =dabs(Gbc(i+1) - Gbc(i-1))/(2*dx)
      
      if (dGc .LE. a0) then
         beta1 = 0d0
      else 
         beta1 = -2d0/3d0
      end if 
         
      ! speed bounds
      !alpha = max(1,beta2/ (2/3 + beta1) )
      ! only works if beta1 != 2/3
      ! if beta1 == 2/3, then must have beta2 - enforced at user level
      if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))  then
         alpha = 1
      else
         alpha = max(1d0,beta2 / (2d0/3d0 + beta1))
      end if
      
      sl  = min(0d0, uir - dsqrt(alpha*ga*hir) ,
     . uip1l - dsqrt(alpha*ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(alpha*ga*hir) ,
     . uip1l + dsqrt(alpha*ga*hip1l)  )
      
      !left and right flux
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 
     .      - 2d0/3d0*(1d0 + 3d0*beta1/2d0)*hir**3*duir**2
     .      - beta2/2*ga*hir**2*(hir*ddhir + dhir**2/2)
           
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - 2d0/3d0*(1d0 + 3d0*beta1/2d0)*hip1l**3*duip1l**2
     .      - beta2/2*ga*hip1l**2*(hip1l*ddhip1l + dhip1l**2/2)
     
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if 
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
      
      !return gradient of h,G across cell i + 1 (for next iteration)
      cdGi = cdGip1
      cdhi = cdhip1   
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
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION CellEnergies(4)
      
      integer i
      DOUBLE PRECISION fGPe(4),sGPe(4),tGPe(4)
      
      DOUBLE PRECISION GPmxj,hGP,GGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5), GCoeff(5)
      
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
      fGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
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
      sGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !third gauss point
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
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

      


