c ====================================================================================
c Fortran program that uses FDVM2 to solve Serre equations
c Convergence and Conservation properites for soliton agree with thesis results
c ready to try and solve Surface Tension and Generalised SWWE
c
c ====================================================================================



      module FDVM2Procedures
      implicit none

      contains

c=====================================
c Program that given array of h calculates h* the reconstructed array
c Analysis shows second order until 10^-10 by which point we get round-off error effects
c=====================================
      subroutine Solve(x_start,x_end,x_len,n_GhstCells,tstart,tend,
     . dtodx,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
     
      
      real*8,intent(in) :: x_start,x_end,
     . ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,tstart,tend,dtodx
      integer, intent(in) :: x_len, n_GhstCells
      
                         
      ! locations for x,u,h,G (cell nodes in interior and bc)
      real*8 xbc(x_len + 2*n_GhstCells),
     . hbc(x_len + 2*n_GhstCells),
     . Gbc(x_len + 2*n_GhstCells),
     . ubc(x_len + 2*n_GhstCells),
     . hAEbc(x_len + 2*n_GhstCells),
     . GAEbc(x_len + 2*n_GhstCells),
     . uAEbc(x_len + 2*n_GhstCells),
     . beta1bc(x_len + 2*n_GhstCells),
     . beta2bc(x_len + 2*n_GhstCells)
     
      real*8 ct,dt,norm(3),dx
    
      integer xbc_len,i,cti
      xbc_len = x_len + 2*n_GhstCells

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
      dt = dtodx*dx
            
      !Get locations of reconstruction plus Ghost Cells
      call Generatexbc(x_start,dx,x_len,n_GhstCells,xbc)
      
            
      !analytic values for h,u,G at t = tstart
      call Initloc(xbc,tstart,a0,a1,a2,a3,a4,a5,b1a6,b1a7,hbc,Gbc,ubc)
      
      
      beta1bc(:) = b1a6*(xbc(:) - a5*tstart) + b1a7
      beta2bc(:) = b2a6*(xbc(:) - a5*tstart) + b2a7
      
      do i=1,xbc_len
         write(2,*) xbc(i),hbc(i),Gbc(i),ubc(i),beta1bc(i),beta2bc(i)
      end do  
      
      ct = tstart
      cti = 0
      do while (ct < tend ) 
      !do while (cti < 100) 
         call EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,
     .      dx,dt,xbc,ct,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
         ct = ct + dt
         cti = cti + 1
         print *, 'Current Time : ', ct
      end do
      
      !get u at end
      call GetufromhG(hbc,Gbc,ubc,xbc,ct,a5,b1a6,b1a7,dx,n_GhstCells)
      
      !analytic values for h,u,G at t = tend
      call Initloc(xbc,ct,a0,a1,a2,a3,a4,a5,b1a6,b1a7,
     . hAEbc,GAEbc,uAEbc)
      
      beta1bc(:) = b1a6*(xbc(:) - a5*ct) + b1a7
      beta2bc(:) = b2a6*(xbc(:) - a5*ct) + b2a7

      do i=1,xbc_len
         write(3,*) xbc(i),hAEbc(i),hbc(i),GAEbc(i),Gbc(i),
     .    uAEbc(i),ubc(i),beta1bc(i),beta2bc(i) 
      end do  
      
      !write out parameters of experiment
      write(4,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(4,*) 'x_start :',x_start
      write(4,*) 'x_end :',x_end
      write(4,*) 'x_len :',x_len
      write(4,*) 'dx = (x_end - x_start) / (x_len -1)  :' , dx
      write(4,*) 'n_GhstCells :',n_GhstCells
      write(4,*) 'tstart :', tstart
      write(4,*) 'tend :',tend 
      write(4,*) 'dt/dx :' , dtodx
      write(4,*) 'dt = dx*(dt/dx)  :' , dt
      write(4,*) 'gravity :' , ga
      write(4,*) 'a0 :' , a0
      write(4,*) 'a1 :' , a1
      write(4,*) 'a2 :' , a2
      write(4,*) 'a3 :' , a3
      write(4,*) 'a4 :' , a4
      write(4,*) 'a5 :' , a5
      write(4,*) 'b1a6 :' , b1a6 
      write(4,*) 'b1a7 :' , b1a7
      write(4,*) 'b2a6 :' , b2a6 
      write(4,*) 'b2a7 :' , b2a7    

      !Convergence Norms
      norm(1) = sqrt(sum((hAEbc(:) - hbc(:)) **2)) / 
     . sqrt(sum(hAEbc(:)**2))       
      norm(2) = sqrt(sum((GAEbc(:) - Gbc(:)) **2)) / 
     . sqrt(sum(GAEbc(:)**2)) 
      norm(3) = sqrt(sum((uAEbc(:) - ubc(:)) **2)) / 
     . sqrt(sum(uAEbc(:)**2)) 
     
      !write out norms
      write(1,*) dx,norm(1),norm(2),norm(3)
                 
    
      end subroutine Solve


c  ********************************************************************************
c  Part 1 : Initial Conditions
c  Functions that generate Initial Conditions
c
c
c
c
c ********************************************************************************


c==============================
c Function to generate x values at all points: All points in cell 
c for reconstruction + ghost cells
c==============================
      subroutine Generatexbc(x_start,dx,x_len,n_GhstCells,xbc)
          
      implicit none
      real*8,intent(in) :: x_start,dx
      integer,intent(in) :: x_len,n_GhstCells
      real*8, intent(out) :: xbc(x_len + 2*n_GhstCells)
      
      integer i,xbc_len
      
      xbc_len = x_len + 2*n_GhstCells
      
      !Left boundary
      do i=1,x_len + 2*n_GhstCells  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end subroutine Generatexbc


c ====================
c Function to generate analytic h,u,G at all locations x
c ====================      
      subroutine Initloc(x,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,h,G,u)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,intent(in) :: t,a0,a1,a2,a3,a4,a5,b1a6,b1a7
      real*8,dimension(size(x)), intent(out)   :: h,G,u
                   
      real*8 PHI,EXPPHI1,dhdx,dudx,d2udx2,beta1
      
      integer i,n
      
           
      n= size(x)      
      do i=1,n 
         beta1 = b1a6*(x(i) - a5*t) + b1a7
         PHI  = x(i) - a2*t
         EXPPHI1 = dexp(-PHI**2 / (2*a3))
         h(i) = EXPPHI1*a1 + a0
         u(i) = EXPPHI1*a4
         dhdx = -EXPPHI1*PHI*a1/a3
         dudx = -EXPPHI1*PHI*a4/a3
         d2udx2 = -a4*(-PHI**2 + a3)*dexp(-PHI**2/(2*a3))/a3**2
         
         G(i) = h(i)*u(i) -(beta1/2d0 + 1d0/3d0)*
     .      (d2udx2*h(i)**3 + 3*dhdx*dudx*h(i)**2)
      end do 
      
            
      end subroutine Initloc



c  ********************************************************************************
c  Part 2 : Finite Difference Solver
c  Functions thet solve finite difference matrix Au = G
c  for u. A is tridiagonal, so we solve using Thomas algorithm
c
c
c
c ********************************************************************************

c =====================
c Function to perform FD for u
c ==================      
      subroutine GetufromhG(hbc,Gbc,ubc,xbc,t,a5,b1a6,b1a7
     . ,dx,n_GhstCells)
      implicit none
      
      real*8,dimension(:),intent(in) :: hbc,Gbc,xbc
      real*8,dimension(:), intent(inout)   :: ubc
      real*8,intent(in) :: dx,t,a5,b1a6,b1a7
      integer,intent(in) :: n_GhstCells
                   
      real*8 subdiag1(size(hbc)),
     . diag(size(hbc)),
     . supdiag1(size(hbc)),
     . RHS(size(hbc))
          
      real*8 ht1,ht2,beta1
     
      integer i,n
      
      n = size(hbc)
                  
      !calculate diagonals in interior
      !set RHS B
            
      do i=n_GhstCells+1,n - n_GhstCells 
         beta1 = b1a6*(xbc(i) - a5*t) + b1a7
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

      subroutine ThomasTriSolve(n,a,b,c,d,x)
      implicit none
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

      subroutine FluxStep(hbc,Gbc,ubc,n_GhstCells,dt,dx,xbc,t,ga,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,hp,Gp)
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: dt,dx,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7
      real*8,dimension(:),intent(in) :: xbc,hbc,Gbc,ubc
      
      real*8, dimension(size(xbc) - 2*n_GhstCells),
     .   intent(out) :: hp,Gp
     
     
      real*8 cdhi,cdGi,cdhip1,cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l, fih,fiG,foh,foG,
     . dhir,ddhir,dhip1l,ddhip1l,
     . dhdt, dGdt, dfluxhdx, dfluxGdx, alpha,beta1,beta2
      
      integer n,nbc,i
      
      nbc = size(xbc)
      n = nbc - 2*n_GhstCells
      
      !Do left boundary
      i = n_GhstCells
      
      cdhi = (hbc(i+1) - hbc(i-1))/2
      cdGi = (Gbc(i+1) - Gbc(i-1))/2
      
c     centered approximation
      !use slopes from last calculation         
      hir = hbc(i) + cdhi/2
      Gir = Gbc(i) + cdGi/2
      uir = (ubc(i+1)+ubc(i))/2
      duir = (ubc(i+1) - ubc(i)) /dx
      dhir = (hbc(i+1) - hbc(i)) /dx
      ddhir = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )/(2*dx**2)

      !reconstruct values on right side of edge x_{i-1/2}
      cdhip1 =  (hbc(i+2) - hbc(i))/2
      cdGip1 = (Gbc(i+2) - Gbc(i))/2
      
      hip1l = hbc(i + 1) - cdhip1/2
      Gip1l = Gbc(i + 1) - cdGip1/2
      uip1l = (ubc(i+1)+ubc(i))/2
      duip1l = (ubc(i+1) - ubc(i)) /dx
      dhip1l = (hbc(i+1) - hbc(i)) /dx
      ddhip1l = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )
     .     /(2*dx**2)
          
      ! speed bounds
      !alpha = max(1,beta2/ (2/3 + beta1) )
      ! only works if beta1 != 2/3
      ! if beta1 == 2/3, then must have beta2 - enforced at user level
      beta1 = b1a6*(xbc(i) + 0.5*dx - a5*t) + b1a7
      beta2 = b2a6*(xbc(i) + 0.5*dx - a5*t) + b2a7
      if  (dabs(2d0/3d0 + beta1) < 10d0**(-10))   then
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
      
      !The flux from cell i -1 to cell i, in the next iteration where i -> i +1
      ! is the same as the flux we just calculated
      fih = foh
      fiG = foG
      
      !similarly, cdhip1 and cdGip1 become cdhi, cdgi when i -> i +1
      cdhi = cdhip1
      cdGi = cdGip1
      
      ! i is numbering by xbc (cells and ghost cells)
      do i = n_GhstCells + 1, nbc - n_GhstCells
      
         hir = hbc(i) + cdhi/2
         Gir = Gbc(i) + cdGi/2
         uir = (ubc(i+1)+ubc(i))/2
         duir = (ubc(i+1) - ubc(i)) /dx
         dhir = (hbc(i+1) - hbc(i)) /dx
         ddhir = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )/(2*dx**2)

         !reconstruct values on right side of edge x_{i-1/2}
         cdhip1 =  (hbc(i+2) - hbc(i))/2
         cdGip1 = (Gbc(i+2) - Gbc(i))/2
         
         hip1l = hbc(i + 1) - cdhip1/2
         Gip1l = Gbc(i + 1) - cdGip1/2
         uip1l = (ubc(i+1)+ubc(i))/2
         duip1l = (ubc(i+1) - ubc(i)) /dx
         dhip1l = (hbc(i+1) - hbc(i)) /dx
         ddhip1l = (hbc(i+2)  - hbc(i+1) - hbc(i) + hbc(i-1) )
     .         /(2*dx**2)
     
         ! speed bounds
         !alpha = max(1,beta2/ (2/3 + beta1) )
         ! only works if beta1 != 2/3
         ! if beta1 == 2/3, then must have beta2 - enforced at user level
         beta1 = b1a6*(xbc(i) + 0.5*dx - a5*t) + b1a7
         beta2 = b2a6*(xbc(i) + 0.5*dx - a5*t) + b2a7
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
         
         call Forcedht(xbc(i),t,a1,a2,a3,dhdt)
         
         call ForcedGt(xbc(i),t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,dGdt)
         
         call Forcedfluxhx(xbc(i),t,a0,a1,a2,a3,a4,dfluxhdx) 
         call ForcedfluxGx(xbc(i),t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxGdx)      
         
         hp(i - n_GhstCells ) = hbc(i) + dt*dhdt 
     .    + dt*( dfluxhdx  - (foh - fih)/dx ) 
     
         Gp(i - n_GhstCells ) =  Gbc(i)+ dt*dGdt 
     .    + dt*( dfluxGdx  - (foG - fiG)/dx ) 

         
         !The flux from cell i -1 to cell i, in the next iteration where i -> i +1
         ! is the same as the flux we just calculated
         fih = foh
         fiG = foG
         
         cdhi = cdhip1
         cdGi = cdGip1
      end do      
      end subroutine FluxStep 
      
      subroutine Forcedht(x,t,a1,a2,a3,dhdt) 
      real*8, intent(out) :: dhdt
      real*8, intent(in) :: x,t,a1,a2,a3
      real*8 :: EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      dhdt = EXPPHI1*PHI*a1*a2/a3
      
      end subroutine Forcedht
      
      subroutine ForcedGt(x,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,dGdt)
      real*8, intent(out) :: dGdt 
      real*8, intent(in) :: x,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7
      real*8 :: EXPPHI1,PHI,hi,ui
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1

      dGdt = 3.0*EXPPHI1**2*PHI**3*a1**2*a2*a5*b1a6*hi*t*ui/a3**3 
     . - 3.0*EXPPHI1**2*PHI**3*a1**2*a2*b1a6*hi*ui*x/a3**3 
     . - 3.0*EXPPHI1**2*PHI**3*a1**2*a2*b1a7*hi*ui/a3**3 
     . - 2*EXPPHI1**2*PHI**3*a1**2*a2*hi*ui/a3**3 
     . + 1.5*EXPPHI1*PHI**2*a1*a5*b1a6*hi**2*ui/a3**2 
     . + EXPPHI1*PHI*a1*a2*ui/a3 
     . + 4.5*PHI**3*a1*a2*a5*b1a6*hi**2*t*ui*dexp(-PHI**2/(2*a3))/a3**3 
     . - 4.5*PHI**3*a1*a2*b1a6*hi**2*ui*x*dexp(-PHI**2/(2*a3))/a3**3 
     . - 4.5*PHI**3*a1*a2*b1a7*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**3 
     . - 3*PHI**3*a1*a2*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**3 
     . + 0.5*PHI**3*a2*a4*a5*b1a6*hi**3*t*dexp(-PHI**2/(2*a3))/a3**3 
     . - 0.5*PHI**3*a2*a4*b1a6*hi**3*x*dexp(-PHI**2/(2*a3))/a3**3 
     . - 0.5*PHI**3*a2*a4*b1a7*hi**3*dexp(-PHI**2/(2*a3))/a3**3 
     . - PHI**3*a2*a4*hi**3*dexp(-PHI**2/(2*a3))/(3*a3**3) 
     . + 0.5*PHI**2*a4*a5*b1a6*hi**3*dexp(-PHI**2/(2*a3))/a3**2 
     . - 4.5*PHI*a1*a2*a5*b1a6*hi**2*t*ui*dexp(-PHI**2/(2*a3))/a3**2 
     . + 4.5*PHI*a1*a2*b1a6*hi**2*ui*x*dexp(-PHI**2/(2*a3))/a3**2 
     . + 4.5*PHI*a1*a2*b1a7*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**2 
     . + 3*PHI*a1*a2*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**2 
     . + PHI*a2*hi*ui/a3 
     . - 1.5*PHI*a2*a4*a5*b1a6*hi**3*t*dexp(-PHI**2/(2*a3))/a3**2 
     . + 1.5*PHI*a2*a4*b1a6*hi**3*x*dexp(-PHI**2/(2*a3))/a3**2 
     . + 1.5*PHI*a2*a4*b1a7*hi**3*dexp(-PHI**2/(2*a3))/a3**2 
     . + PHI*a2*a4*hi**3*dexp(-PHI**2/(2*a3))/a3**2 
     . - 0.5*a4*a5*b1a6*hi**3*dexp(-PHI**2/(2*a3))/a3
      end subroutine ForcedGt
      
      subroutine Forcedfluxhx(x,t,a0,a1,a2,a3,a4,dfluxhdx) 
      real*8, intent(out) :: dfluxhdx
      real*8, intent(in) :: x,t,a0,a1,a2,a3,a4
      real*8 :: EXPPHI1,PHI,hi,ui,dhdx,dudx
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      
      dfluxhdx = ui*dhdx + hi*dudx
      
      end subroutine Forcedfluxhx
      
      subroutine ForcedfluxGx(x,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,dfluxGdx)
      real*8, intent(out) :: dfluxGdx
      real*8, intent(in) :: x,t,ga,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7
      real*8 :: EXPPHI1,PHI,hi,ui
      real*8 :: dhdx,d2hdx2,d3hdx3,dudx,d2udx2,d3udx3,beta1,beta2,
     . db1x, db2x
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      
      dhdx = -EXPPHI1*PHI*a1/a3
      d2hdx2 = -a1*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3hdx3 = PHI*a1*(-PHI**2 + 3*a3)*EXPPHI1/a3**3
      dudx = -EXPPHI1*PHI*a4/a3
      d2udx2 = -a4*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3udx3 = PHI*a4*(-PHI**2 + 3*a3)*EXPPHI1/a3**3
      
      beta1 = b1a6*(x - a5*t) + b1a7
      db1x = b1a6
      beta2 = b2a6*(x - a5*t) + b2a7 
      db2x = b2a6
     
      dfluxGdx = -2*d2udx2*dudx*hi**3*(beta1 + 2d0/3d0)
     .  - 3*dhdx*dudx**2*hi**2*(beta1 + 2d0/3d0) 
     . - dhdx*ga*hi*(d2hdx2*hi + dhdx**2/2)*beta2 
     . + dhdx*ga*hi - dudx**2*hi**3*db1x
     . + dudx*(hi*ui - (d2udx2*hi**3 + 3*dhdx*dudx*hi**2)*
     .   (0.5*beta1 + 1d0/3d0)) 
     . - ga*hi**2*(2*d2hdx2*dhdx 
     . + d3hdx3*hi)*beta2/2d0 
     . - ga*hi**2*(d2hdx2*hi + dhdx**2/2d0)*db2x/2 
     . + ui*(dhdx*ui + dudx*hi - 0.5*(d2udx2*hi**3 + 3*dhdx*dudx*hi**2)
     .  *db1x
     . + (-0.5*beta1 - 1d0/3d0)*
     . (3*d2hdx2*dudx*hi**2 + 6*d2udx2*dhdx*hi**2 
     . + d3udx3*hi**3 + 6*dhdx**2*dudx*hi))
     
      end subroutine ForcedfluxGx
      
      subroutine EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,
     . dx,dt,xbc,t,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
     
           
      real*8,dimension(:),intent(inout) :: hbc,Gbc,ubc
      real*8, dimension(:),intent(in) :: xbc
      real*8, intent(in) :: dx,dt,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,t,ga
     
      real*8, dimension(size(hbc)) :: hpbc,Gpbc
      real*8, dimension(size(hbc) - 2*n_GhstCells) :: hp,Gp
      integer :: n_GhstCells,nbc,n
      
      nbc = size(hbc)
      n = nbc - 2*n_GhstCells
      
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(hbc,Gbc,ubc,xbc,t,a5,b1a6,b1a7
     . ,dx,n_GhstCells)
            
      !Update cell averages
      call FluxStep(hbc,Gbc,ubc,n_GhstCells,dt,dx,xbc,t,ga,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,hp,Gp)
     
      !BC's
      call Initloc(xbc(1:n_GhstCells),
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,hpbc(1:n_GhstCells),
     .   Gpbc(1:n_GhstCells),ubc(1:n_GhstCells))
      
      ! fill interior
      hpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = hp(:)
      Gpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = Gp(:)
      
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,
     .   hpbc(nbc - n_GhstCells + 1:nbc),
     .   Gpbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc))
          
      
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(hpbc,Gpbc,ubc,xbc,t + dt,a5,b1a6,b1a7
     . ,dx,n_GhstCells)
      
      !Update cell averages    
      call FluxStep(hpbc,Gpbc,ubc,n_GhstCells,dt,dx,xbc,t+dt,ga,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,hp,Gp)
      
      !boundary conditions of hbc, Gbc stay same, just have to update interior
      hbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     . ( hbc(n_GhstCells + 1 : nbc - n_GhstCells ) + hp(:))/2d0
     
      Gbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     .  ( Gbc(n_GhstCells + 1 : nbc - n_GhstCells) + Gp(:) )/2d0
     
      !BC's
      call Initloc(xbc(1:n_GhstCells),
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,hpbc(1:n_GhstCells),
     .   Gpbc(1:n_GhstCells),ubc(1:n_GhstCells))
      
      
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,a0,a1,a2,a3,a4,a5,b1a6,b1a7,
     .   hpbc(nbc - n_GhstCells + 1:nbc),
     .   Gpbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc))
          
      
      end subroutine EvolveStepWrap

      
      end module FDVM2Procedures
      

c  ********************************************************************************
c  Main Program
c  Functions that loops over many dx (x_len), to investigate convergence
c  
c  The wdir is the working directory - all data printed from solver will end up in there
c
c
c ********************************************************************************
      program main
      
      use FDVM2Procedures
            
      implicit none
      
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/Models/gSGN/Forced/NV-beta1beta2Vary/"
  

      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 startx,endx,ga,tstart,tend,Cr,
     . maxwavespeed,dtodx,a0,a1,a2,a3,a4,a5,alpha,
     . b1max,b1min,b2max,b2min,phi0,phi1,b1a6,b1a7,b2a6,b2a7
      
      
      n=15
      n_gc = 6
      
      lowxlen = 100
      
      startx = -50d0
      endx = 100.0d0
      
      tstart = 0d0
      tend = 10.0d0
            
      ga = 9.81d0
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      
      ! betamax,betamin,a5,a6,a7 control dispersion properties, size of terms
      b1max = 2d0
      b1min = -2d0/3d0
      a5 = -0.9
      phi0 = startx - a5*tstart
      phi1 = endx - a5*tend 
      
      b1a6 = (b1max - b1min) / (phi1 - phi0)
      b1a7 = b1min - b1a6*phi0
      
      b2max = 3
      b2min = 0
      
      b2a6 = (b2max - b2min) / (phi1 - phi0)
      b2a7 = b2min - b2a6*phi0
      
      !alpha is the maximum ratio between beta2 / (2/3 +beta1 )
      alpha = 1.2
      maxwavespeed = a2 + a4 + dsqrt(alpha*ga*(a0 + a1))
      Cr = 0.5
      dtodx = Cr/maxwavespeed
            
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      !open files comparing multiple dx values
      open(1, file = wdir//'Norms.dat') 
      
      
      do i=1,n
         
         !create directory to stored output data
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         
         !open files for single experiment
         open(2, file = trim(fileloci)//'InitVars.dat')  
         open(3, file = trim(fileloci)//'OutVars.dat')  
         open(4, file = trim(fileloci)//'Params.dat')  
         
         
         print *,'Experiment : ',i ,' || ', '# Cells :',
     .      lowxlen*(2**(i-1))

         call Solve(startx,endx,lowxlen*(2**(i-1)),n_gc,tstart,tend,
     .   dtodx,ga,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
         
      end do
      
      close(1)  
      close(2)
      close(3)
      close(4)

      end program main
      

      
