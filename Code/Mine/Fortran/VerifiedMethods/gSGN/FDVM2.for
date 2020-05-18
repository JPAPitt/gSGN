c ====================================================================================
c Module of fortran subroutines use FDVM2 to solve generalised Serre - Green -Naghdi equations (gSGN)
c it takes initial conditions and returns solutions Q(x,t) where t > tend (Q = [h,G])
c the subroutine also calculates total conserved quantities initially and in the final solution
c ====================================================================================



      module FDVM2Procedures
      implicit none

      contains

c=====================================
c Numerical solver program, evolving system through time and producing Q(x,t) where t > tend
c The solver uses arrays xbc,hbc_init,Gbc_init and ubc_init
c to produce: 
c initial energy - Energs_init
c solution for h at final time hbc_fin
c solution for G at final time Gbc_fin
c solution for u at final time ubc_fin
c Notes:
c 1. ubc_init only needs to have u defined at ghost cells
c function will update interior using up to date values of h and G.
c
c 2. Solver assumes beta's are constant and ghost cell values are constant- constant dirichlet boundary conditions
c=====================================
      subroutine NumericalSolve(tstart,tend,
     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc,
     . hbc_init,Gbc_init,ubc_init,
     . Energs_init,currenttime,hbc_fin,Gbc_fin,ubc_fin,Energs_fin)
     
      integer, intent(in) :: n_GhstCells
      real*8,intent(in) :: tstart,tend,ga,beta1,beta2,theta,dx,dt
      real*8,dimension(:),intent(in) :: xbc,hbc_init,Gbc_init,ubc_init
      real*8,dimension(:),intent(out) :: Energs_init,hbc_fin,Gbc_fin,
     . ubc_fin,Energs_fin
      real*8,intent(out) :: currenttime 
      
      integer i,ileft,iright
      
      !initial time
      currenttime  = tstart
      
      !loop over and set hbc_fin,Gbc_fin to initial conditions
      ileft = 1
      iright = size(xbc)
      do i = ileft,iright
         hbc_fin(i) = hbc_init(i) 
         Gbc_fin(i) = Gbc_init(i) 
         ubc_fin(i) = ubc_init(i)
      end do
      
      !calculate initial Energies
      call GetufromhG(hbc_fin,Gbc_fin,ubc_fin,beta1,dx,n_GhstCells)
      call TotalEnergy(hbc_fin,ubc_fin,Gbc_fin,ga,beta1,beta2,
     . n_GhstCells,dx,Energs_init)
     
      
      !evolve the system through time
      do while (currenttime  < tend ) 
      
         call EvolveStepWrap(hbc_fin,Gbc_fin,ubc_fin,ga,beta1,beta2
     . ,theta,n_GhstCells,dx,dt)
     
         currenttime  = currenttime  + dt
         print *, 'Current Time : ', currenttime 
      end do
      
      !calculate end energies  
      call GetufromhG(hbc_fin,Gbc_fin,ubc_fin,beta1,dx,n_GhstCells)
      call TotalEnergy(hbc_fin,ubc_fin,Gbc_fin,ga,beta1,beta2,
     . n_GhstCells,dx,Energs_fin)
      
                 
    
      end subroutine NumericalSolve


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
      subroutine GetufromhG(hbc,Gbc,ubc,beta1,dx,n_GhstCells)
      
      real*8,dimension(:),intent(in) :: hbc,Gbc
      real*8,dimension(:), intent(inout)   :: ubc
      real*8,intent(in) :: dx,beta1
      integer,intent(in) :: n_GhstCells
                   
      real*8 subdiag1(size(hbc)),
     . diag(size(hbc)),
     . supdiag1(size(hbc)),
     . RHS(size(hbc))
          
      real*8 ht1,ht2
     
      integer i,n
      
      n = size(hbc)
                  
      !calculate diagonals in interior
      !set RHS B
            
      do i=n_GhstCells+1,n - n_GhstCells 
         ht1 = (beta1/2d0 + 1d0/3d0)*(hbc(i)**3/(dx*dx))
         ht2 = (3*beta1/2d0 + 1d0)*
     .      (hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1)))
         
         
         subdiag1(i)  = -ht1 + ht2
         diag(i) = hbc(i) + (2d0*ht1)
         supdiag1(i)  = -ht1 - ht2 
                  
         RHS(i) = Gbc(i)
      end do 
      
      !first and last n_GhstCells  x n_GhstCells in tridiagmatrix
      !Should be identity
      subdiag1(:n_GhstCells) = 0d0
      diag(:n_GhstCells) = 1d0
      supdiag1(:n_GhstCells) = 0d0
      
      RHS(:n_GhstCells) = ubc(:n_GhstCells)
      
      subdiag1(n - (n_GhstCells-1 ) : n) = 0d0
      diag(n - (n_GhstCells-1 ) : n) = 1d0
      supdiag1(n - (n_GhstCells-1 ) : n) = 0d0
      
      RHS(n - (n_GhstCells-1 ) : n) = ubc(n - (n_GhstCells-1 ) : n)
      
   
      call ThomasTriSolve(n,subdiag1,diag,supdiag1,RHS,ubc)
      
      end subroutine GetufromhG


c ================
c ThomasTriSolve solves a tridiagonal matrix equation Ax = b, where A is tridiagonal
c =================   
      subroutine ThomasTriSolve(n,a,b,c,d,x)
      integer, intent(in) :: n
      real*8, intent(in) :: a(n), c(n)
      real*8, intent(inout), dimension(n) :: b, d
      real*8, intent(out) :: x(n)
      !  --- Local variables ---
      integer :: i
      real*8 :: q
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
      
      end subroutine ThomasTriSolve


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
      subroutine EvolveStepWrap(hbc,Gbc,ubc,ga,beta1,beta2,theta
     . ,n_GhstCells,dx,dt)
     
      real*8,dimension(:),intent(inout) :: hbc,Gbc,ubc
      real*8, intent(in) :: ga,beta1,beta2,theta,dx,dt
      integer, intent(in) :: n_GhstCells
     
      real*8, dimension(size(hbc)) :: hpbc,Gpbc,hppbc,Gppbc
      integer :: i,xbc_len
      
      xbc_len = size(hbc)      
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(hbc,Gbc,ubc,beta1,dx,n_GhstCells)
            
      !Update cell averages (first order approximation to h^{n+1}, G^{n+1})
      call SingleEulerStep(hbc,Gbc,ubc,ga,beta1,beta2
     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(hpbc,Gpbc,ubc,beta1,dx,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+2}, G^{n+2})
      call SingleEulerStep(hpbc,Gpbc,ubc,ga,beta1,beta2
     . ,theta,n_GhstCells,dt,dx,hppbc,Gppbc)
      
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      do i= 1,xbc_len
         hbc(i) = ( hbc(i ) + hppbc(i))/2d0
         Gbc(i) = ( Gbc(i ) + Gppbc(i))/2d0
      end do
          
      
      end subroutine EvolveStepWrap


c ====
c SingleEulerStep
c produces new cell averages of h,G using forward Euler step,
c with flux approximated using Kurganov's method
c ====
      subroutine SingleEulerStep(hbc,Gbc,ubc,ga,beta1,beta2
     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
     
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: dt,dx,ga,beta1,beta2,theta
      real*8,dimension(:),intent(in) :: hbc,Gbc,ubc
      
      real*8, dimension(size(hbc)),intent(out) :: hpbc,Gpbc
     
     
      real*8 cdhi,cdGi, fih,fiG,foh,foG
      
      integer i,xbc_len,ileft,iright
      
      xbc_len = size(hbc)
      
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
      call Fluxxiph(hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
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
         call Fluxxiph(hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
     . foh,foG,cdhi,cdGi)
             
         hpbc(i) = hbc(i) - dt*(foh - fih)/dx 
     
         Gpbc(i) = Gbc(i) - dt*(foG - fiG)/dx 

         !flux out becomes flux in on next cell
         fih = foh
         fiG = foG
         
      end do 
           
      end subroutine SingleEulerStep
      
c subroutine that given arrays, and i calculates flux across x_{i+1/2} for gSGN equations  
c note that - cdhi,cdGi is inout, on the way in cdhi,cdGi are our approximation to the gradient
c of h and G across cell i, on the way out cdhi,cdGi are our approximation to the gradient of h and G
c across cell i +1 (which will be cell i as the method moves to the next cell)    
      subroutine Fluxxiph(hbc,Gbc,ubc,ga,beta1,beta2,theta,dx,i,
     . foh,foG,cdhi,cdGi)
     
       real*8, dimension(:), intent(in) :: hbc,Gbc,ubc
       real*8, intent(in) :: ga,beta1,beta2,theta,dx
       integer, intent(in) :: i
       real*8, intent(inout) :: cdhi,cdGi
       real*8, intent(out) :: foh,foG
     
       real*8 cdhip1,cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l,
     . dhir,ddhir,dhip1l,ddhip1l, alpha
     
      !reconstruct values on left side of edge x_{i+1/2}
      hir = hbc(i) + cdhi/2
      Gir = Gbc(i) + cdGi/2
      uir = (ubc(i+1)+ubc(i))/2
      duir = (2*ubc(i) - 3*ubc(i-1) + ubc(i-2)) /dx
      dhir = (2*hbc(i) - 3*hbc(i-1) + hbc(i-2))  /dx
      ddhir = (5*hbc(i) - 13*hbc(i-1) + 11*hbc(i-2) - 3*hbc(i-3))
     .      /(2*dx**2)

      !reconstruct values on right side of edge x_{i+1/2}
      call ReconLinLimGrad(hbc(i),hbc(i+1),hbc(i+2),theta,cdhip1)
      call ReconLinLimGrad(Gbc(i),Gbc(i+1),Gbc(i+2),theta,cdGip1)
      
      hip1l = hbc(i + 1) - cdhip1/2
      Gip1l = Gbc(i + 1) - cdGip1/2
      uip1l = (ubc(i+1)+ubc(i))/2
      duip1l = (-2*ubc(i+1) + 3*ubc(i+2) - ubc(i+3)) /dx
      dhip1l = (-2*hbc(i+1) + 3*hbc(i+2) - hbc(i+3)) /dx
      ddhip1l = (5*hbc(i+1) - 13*hbc(i+2) + 11*hbc(i+3) 
     .    - 3*hbc(i+4)) /(2*dx**2)
  
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
      cdhi = cdhip1 
      cdGi = cdGip1
         
      end subroutine  Fluxxiph
      
c ====
c ReconLinLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconLinLimGrad(qim1,qi,qip1,theta,cdq)
      real*8, intent(in) :: qim1,qi,qip1,theta
      real*8, intent(out) :: cdq
      real*8 :: fdq,mdq,bdq
      
      fdq = qip1 - qi
      mdq = 0.5*(qip1 - qim1)
      bdq = qi - qim1
      
      call minmod(theta*fdq,mdq,theta*bdq,cdq)
      
      end subroutine ReconLinLimGrad

c ===
c minmod is used for limiting and returns smallest size argument if all arguments have same sign,
c otherwise it returns 0.
c ===      
      subroutine minmod(a,b,c,d)
      implicit none
      real*8, intent(in) :: a,b,c
      real*8,intent(out) :: d
      
      if ((a .gt. 0) .and. (b .gt. 0) .and. (c .gt. 0)) then
         d = min(a,b,c)
      else if ((a .lt. 0) .and. (b .lt. 0) .and. (c .lt. 0)) then
         d = max(a,b,c)
      else
         d = 0.0
      end if
      
      end subroutine minmod

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
     
      real*8, intent(in) :: qjm2,qjm1,qj,qjp1,qjp2,dx
      real*8,dimension(5), intent(out) :: QuartCoeff
      
      QuartCoeff(1) = (qjp2 - 4*qjp1 + 6*qj - 4*qjm1 + qjm2) /
     . (24* (dx**4))
      QuartCoeff(2) = (qjp2 - 2*qjp1 + 2*qjm1 - qjm2) /
     . (12* (dx**3))  
      QuartCoeff(3) = (-qjp2 + 16*qjp1 - 30*qj + 16*qjm1 - qjm2) /
     . (24* (dx**2))  
      QuartCoeff(4) = (-qjp2 + 8*qjp1 - 8*qjm1 + qjm2) /
     . (12*dx)
      QuartCoeff(5) = qj
      end subroutine QuarticInterp


c ========================
c Functions to evaluate Quartic at x, with quartic centered around x_j
c xmxj = x - x_j 
c ====================      
      subroutine QuarticCoeffEvalxj(QuarticCoeff,xmxj,qatxj)
      
      real*8, dimension(5),intent(in) :: QuarticCoeff
      real*8, intent(in) :: xmxj
      real*8, intent(out) :: qatxj
      
      qatxj = QuarticCoeff(1)*xmxj**4 + QuarticCoeff(2)*xmxj**3 +
     . QuarticCoeff(3)*xmxj**2 + QuarticCoeff(4)*xmxj +
     . QuarticCoeff(5)
     
      end subroutine QuarticCoeffEvalxj
      
      subroutine QuarticCoeffEvalGradxj(QuarticCoeff,xmxj,dqatxj)
      
      real*8, dimension(5),intent(in) :: QuarticCoeff
      real*8, intent(in) :: xmxj
      real*8, intent(out) :: dqatxj
      
      dqatxj = 4*QuarticCoeff(1)*xmxj**3 + 3*QuarticCoeff(2)*xmxj**2 +
     . 2*QuarticCoeff(3)*xmxj + QuarticCoeff(4)
     
      end subroutine QuarticCoeffEvalGradxj

c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      subroutine AllEnergiesIntegralCell(h,u,G,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      real*8, dimension(:),intent(in) :: h,u,G
      real*8, intent(in) :: dx,ga,beta1,beta2
      integer, intent(in) :: j
      real*8, dimension(4),intent(out) :: CellEnergies
      
      integer ::  i
      real*8 :: fGPe(4),sGPe(4),tGPe(4)
      
      real*8 :: GPmxj,hGP,GGP,uGP,uxGP,hxGP
      real*8 :: hCoeff(5), uCoeff(5), GCoeff(5)
      
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
      
      end subroutine AllEnergiesIntegralCell
      
      !Function to sum all energies
      subroutine TotalEnergy(hbc,ubc,Gbc,ga,beta1,beta2,
     . n_GhstCells,dx,TotEnergVals)
      
      real*8, dimension(:),intent(in) :: hbc,ubc,Gbc
      real*8, intent(in) :: dx,ga,beta1,beta2
      integer, intent(in) :: n_GhstCells
      real*8, dimension(4),intent(out) :: TotEnergVals
      
      real*8 CellEnergVals(4)
      integer i,j,xbc_len
      
      xbc_len = size(hbc)
      
      !running totals for energy values, start at 0
      do i = 1,4
         TotEnergVals(i) = 0.0
      end do
      
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      do j= n_GhstCells + 1, xbc_len - n_GhstCells
         call  AllEnergiesIntegralCell(hbc,ubc,Gbc,ga,beta1,beta2
     .      ,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,4
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end subroutine TotalEnergy
      

      
      end module FDVM2Procedures
      


