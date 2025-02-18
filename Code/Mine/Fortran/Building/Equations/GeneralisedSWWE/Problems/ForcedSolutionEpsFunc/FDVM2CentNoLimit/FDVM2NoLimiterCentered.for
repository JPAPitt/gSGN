c ====================================================================================
c Fortran program that uses FDVM2 to solve Serre equations
c Convergence and Conservation properites for soliton agree with thesis results
c ready to try and solve Surface Tension and Generalised SWWE
c
c ====================================================================================



c=====================================
c Program that given array of h calculates h* the reconstructed array
c Analysis shows second order until 10^-10 by which point we get round-off error effects
c=====================================
      subroutine Solve(x_start,x_end,x_len, n_GhstCells,tstart,tend,
     . dtodx,ga,a0,a1,a2,a3,a4,a5,a6,a7,FilePre,norm,dx)
      implicit none
      
      
      real*8,intent(in) :: x_start,x_end,
     . ga,a0,a1,a2,a3,a4,a5,a6,a7,tstart,tend,dtodx
      integer, intent(in) :: x_len, n_GhstCells
      character (LEN=*), intent(in)   :: FilePre
      character (len = 5) :: strti
      real*8, intent(out)   :: norm(3),dx
      
                         
      ! locations for x,u,h,G (cell nodes in interior and bc)
      real*8 xbc(x_len + 2*n_GhstCells),
     . hbc(x_len + 2*n_GhstCells),
     . Gbc(x_len + 2*n_GhstCells),
     . ubc(x_len + 2*n_GhstCells),
     . hAEbc(x_len + 2*n_GhstCells),
     . GAEbc(x_len + 2*n_GhstCells),
     . uAEbc(x_len + 2*n_GhstCells),
     . Epsbc(x_len + 2*n_GhstCells)
     
      real*8 ct,dt
    
      integer xbc_len,i,cti
      xbc_len = x_len + 2*n_GhstCells
      

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
      dt = dtodx*dx
            
      !Get locations of reconstruction plus Ghost Cells
      xbc = Generatexbc(x_start,dx,x_len,n_GhstCells)
      
            
      !analytic values for h,u,G at t = tstart
      call Initloc(xbc,tstart,a0,a1,a2,a3,a4,a5,a6,a7,
     .   hbc,Gbc,ubc,Epsbc)
      
      open(1, file = FilePre//'InitVars.dat', status = 'unknown')  
      do i=1,xbc_len
         write(1,*) xbc(i),hbc(i),Gbc(i),ubc(i),Epsbc(i)
      end do  
      close(1)
      
      ct = tstart
      cti = 0
      do while (ct < tend ) 
      !do while (cti < 10)
         call EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,dx,dt,
     .     xbc,ct,a0,a1,a2,a3,a4,a5,a6,a7)
         ct = ct + dt
         cti = cti + 1
         print *, 'Current Time : ', ct
      end do

      
      !get u at end
      call GetufromhG(xbc,hbc,Gbc,ubc,ct,a5,a6,a7,dx,n_GhstCells)
      
      
      !analytic values for h,u,G at t = tend
      call Initloc(xbc,ct,a0,a1,a2,a3,a4,a5,a6,a7,
     . hAEbc,GAEbc,uAEbc,Epsbc)
      
   
      open(2, file = FilePre//'AllVarsEnd.dat', status = 'unknown')  
      do i=1,xbc_len
         write(2,*) xbc(i),hAEbc(i),hbc(i),GAEbc(i),Gbc(i),
     .    uAEbc(i),ubc(i),Epsbc(i)
      end do  
      close(2)
      
      !write out parameters of experiment
      open(5, file = FilePre//'Parameters.dat', status = 'unknown')  
      write(5,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(5,*) 'x_start :',x_start
      write(5,*) 'x_end :',x_end
      write(5,*) 'x_len :',x_len
      write(5,*) 'dx = (x_end - x_start) / (x_len -1)  :' , dx
      write(5,*) 'n_GhstCells :',n_GhstCells
      write(5,*) 'tstart :', tstart
      write(5,*) 'tend :',tend 
      write(5,*) 'dt/dx :' , dtodx
      write(5,*) 'dt = dx*(dt/dx)  :' , dt
      write(5,*) 'gravity :' , ga
      write(5,*) 'a0 :' , a0
      write(5,*) 'a1 :' , a1
      write(5,*) 'a2 :' , a2
      write(5,*) 'a3 :' , a3
      write(5,*) 'a4 :' , a4
      write(5,*) 'a5 :' , a5
      write(5,*) 'a6 :' , a6
      write(5,*) 'a7 :' , a7
      close(5)

      !Convergence Norms
      norm(1) = sqrt(sum((hAEbc(:) - hbc(:)) **2)) / 
     . sqrt(sum(hAEbc(:)**2))       
      norm(2) = sqrt(sum((GAEbc(:) - Gbc(:)) **2)) / 
     . sqrt(sum(GAEbc(:)**2)) 
      norm(3) = sqrt(sum((uAEbc(:) - ubc(:)) **2)) / 
     . sqrt(sum(uAEbc(:)**2)) 
                 
    
      contains

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
      function Generatexbc(x_start,dx,x_len,n_GhstCells) result (xbc)
          
      implicit none
      real*8,intent(in) :: x_start,dx
      integer,intent(in) :: x_len,n_GhstCells
      real*8 xbc(x_len + 2*n_GhstCells)
      
      integer i,xbc_len
      
      xbc_len = x_len + 2*n_GhstCells
      
      !Left boundary
      do i=1,x_len + 2*n_GhstCells  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end function Generatexbc


c ====================
c Function to generate analytic h,u,G at all locations x
c ====================      
      subroutine Initloc(x,t,a0,a1,a2,a3,a4,a5,a6,a7,h,G,u,e)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,intent(in) :: t,a0,a1,a2,a3,a4,a5,a6,a7
      real*8,dimension(size(x)), intent(out)   :: h,G,u,e
                   
      real*8 PHI,EXPPHI1,
     . eps,dhdx,dudx,d2udx2
      
      integer i,n
           
      n= size(x)      
      do i=1,n 
         
         PHI  = x(i) - a2*t
         EXPPHI1 = dexp(-PHI**2 / (2*a3))
         h(i) = EXPPHI1*a1 + a0
         u(i) = EXPPHI1*a4
         e(i) = a5*DCOS(a7*(x(i) - a6*t))
         dhdx = -EXPPHI1*PHI*a1/a3
         dudx = -EXPPHI1*PHI*a4/a3
         d2udx2 = -a4*(-PHI**2 + a3)*exp(-PHI**2/(2*a3))/a3**2
         G(i) = h(i)*u(i) -e(i)*(d2udx2*h(i)**3 + 3*dhdx*dudx*h(i)**2)
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
      subroutine GetufromhG(xbc,hbc,Gbc,ubc,t,a5,a6,a7,dx,n_GhstCells)
      implicit none
      
      real*8,dimension(:),intent(in) :: xbc,hbc,Gbc
      real*8,dimension(:), intent(inout)   :: ubc
      real*8,intent(in) :: dx,t,a5,a6,a7
      integer,intent(in) :: n_GhstCells
                   
      real*8 subdiag1(size(hbc)),
     . diag(size(hbc)),
     . supdiag1(size(hbc)),
     . RHS(size(hbc))
          
      real*8 ht1,ht2,eps
     
      integer n
      
      n = size(hbc)
                  
      !calculate diagonals in interior
      !set RHS B
            
      do i=n_GhstCells+1,n - n_GhstCells 
         eps = a5*DCOS(a7*(xbc(i) - a6*t))
         ht1 = eps*(hbc(i)**3/(dx*dx))
         ht2 = 3d0*eps*(hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1)))
         
         
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
      subroutine FluxStep(hbc,Gbc,ubc,n_GhstCells,dt,dx,
     . xbc,t,a0,a1,a2,a3,a4,a5,a6,a7,hp,Gp)
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: dt,dx,t,a0,a1,a2,a3,a4,a5,a6,a7
      real*8,dimension(:),intent(in) :: xbc,hbc,Gbc,ubc
      
      real*8, dimension(size(xbc) - 2*n_GhstCells),
     .   intent(out) :: hp,Gp
     
     
      real*8 cdhi,cdGi,cdhip1,cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l, fih,fiG,foh,foG,
     . dhir,ddhir,dhip1l,ddhip1l,
     . dhdt, dGdt, dfluxhdx, dfluxGdx,SourceG,
     . PHI, EXPPHI1,eps
      
      integer n,nbc,i
      
      nbc = size(xbc)
      n = nbc - 2*n_GhstCells
      
      !Do left boundary
      i = n_GhstCells
      
      !reconstruct values on left side of edge x_{j-1/2}
      cdhi =  (hbc(i+1) - hbc(i-1))/2
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
      sl  = min(0d0, uir - dsqrt(ga*hir) , uip1l - dsqrt(ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(ga*hir) , uip1l + dsqrt(ga*hip1l)  )
      
      eps = a5*DCOS(a7*(xbc(i) + 0.5*dx - a6*t))
      
      !left and right flux
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 
     .      - 2*eps*hir**3*duir**2
     .      - ga*eps*hir**3*ddhir 
     .      - (ga/2d0)*eps*(hir**2)*dhir**2
           
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - 2*eps*hip1l**3*duip1l**2
     .      - ga*eps*hip1l**3*ddhip1l 
     .    - (ga/2d0)*eps*(hip1l**2)*dhip1l**2
     
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
      
c     centered approximation - forced solution seems to require it (solution is smooth but limiting  makes it discontinous, particularly in gradient)
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
     .         /(2*dx**2)
                   
         ! speed bounds
         sl  = min(0d0, uir - dsqrt(ga*hir) , uip1l - dsqrt(ga*hip1l))
         sr  = max(0d0, uir + dsqrt(ga*hir) , uip1l + dsqrt(ga*hip1l))
         
         eps = a5*DCOS(a7*(xbc(i) + 0.5*dx - a6*t))
         
         !left and right flux
         felh = uir*hir
         felG = uir*Gir + ga*(hir**2)/2d0 
     .      - 2*eps*hir**3*duir**2
     .      - ga*eps*hir**3*ddhir 
     .      - (ga/2d0)*eps*(hir**2)*dhir**2
              
         ferh = uip1l*hip1l
         ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - 2d0*eps*hip1l**3*duip1l**2
     .      - ga*eps*hip1l**3*ddhip1l 
     .      - (ga/2d0)*eps*(hip1l**2)*dhip1l**2
        
         if (sr == sl) then
            isrmsl = 0.0
         else
            isrmsl = 1.0 / (sr - sl)
         end if 
         
         !calculate flux from cell i to cell i + 1 (Kurganov)
         foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
         foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
         
         dhdt = Forcedht(xbc(i),t,a0,a1,a2,a3) 
         dGdt = ForcedGt(xbc(i),t,a0,a1,a2,a3,a4,a5,a6,a7) 
         
         dfluxhdx = Forcedfluxhx(xbc(i),t,a0,a1,a2,a3,a4) 
         dfluxGdx = ForcedfluxGx(xbc(i),t,ga,eps,a0,a1,a2,a3
     .      ,a4,a5,a6,a7) 
         
         SourceG = ForcedSourceG(xbc(i),t,a0,a1,a2,a3,a4,a5,a6,a7)  
               
         
         hp(i - n_GhstCells ) = hbc(i) + dt*dhdt 
     .    + dt*( dfluxhdx  - (foh - fih)/dx ) 
     
         Gp(i - n_GhstCells ) = Gbc(i) + dt*dGdt !+ dt*SourceG 
     .    + dt*( dfluxGdx  - (foG - fiG)/dx ) 

         
         !The flux from cell i -1 to cell i, in the next iteration where i -> i +1
         ! is the same as the flux we just calculated
         fih = foh
         fiG = foG
         
         cdhi = cdhip1
         cdGi = cdGip1
      end do
      
      
      end subroutine FluxStep 
      
      function Forcedht(x,t,a0,a1,a2,a3) result (dhdt)
      
      real*8, intent(in) :: x,t,a0,a1,a2,a3
      real*8 :: dhdt,EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      dhdt = EXPPHI1*PHI*a1*a2/a3
      
      end function Forcedht
      
      function ForcedGt(x,t,a0,a1,a2,a3,a4,a5,a6,a7) result (dGdt)
      
      real*8, intent(in) :: x,t,a0,a1,a2,a3,a4,a5,a6,a7
      real*8 :: dGdt,EXPPHI1,PHI,hi,ui,ei
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      ei = a5*DCOS(a7*(x - a6*t))
 
      dGdt = - 6*EXPPHI1**2*PHI**3*a1**2*a2*ei*hi*ui/a3**3 
     . + 3*EXPPHI1*PHI**2*a1*a5*a6*a7*hi**2*ui
     .   *dsin(a6*a7*t - a7*x)/a3**2
     . + EXPPHI1*PHI*a1*a2*ui/a3 
     . - 9*PHI**3*a1*a2*ei*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**3 
     . - PHI**3*a2*a4*ei*hi**3*dexp(-PHI**2/(2*a3))/a3**3 
     . + PHI**2*a4*a5*a6*a7*hi**3*dexp(-PHI**2/(2*a3))
     .   *dsin(a6*a7*t - a7*x)/a3**2 
     . + 9*PHI*a1*a2*ei*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**2 
     . + PHI*a2*hi*ui/a3 
     . + 3*PHI*a2*a4*ei*hi**3*dexp(-PHI**2/(2*a3))/a3**2 
     . - a4*a5*a6*a7*hi**3*dexp(-PHI**2/(2*a3))
     .   *dsin(a6*a7*t - a7*x)/a3 
     
      end function ForcedGt
      
      function Forcedfluxhx(x,t,a0,a1,a2,a3,a4) 
     . result (dfluxhdx)
      
      real*8, intent(in) :: x,t,a0,a1,a2,a3,a4
      real*8 :: dfluxhdx,EXPPHI1,PHI,hi,ui,dhdx,dudx
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      
      dfluxhdx = ui*dhdx + hi*dudx
      
      end function Forcedfluxhx
      
      function ForcedfluxGx(x,t,ga,eps,a0,a1,a2,a3,a4,a5,a6,a7) 
     . result (dfluxGdx)
      
      real*8, intent(in) :: x,t,ga,eps,a0,a1,a2,a3,a4,a5,a6,a7
      real*8 :: dfluxGdx,EXPPHI1,PHI,hi,ui,ei,dedx
      real*8 :: dhdx,d2hdx2,d3hdx3,dudx,d2udx2,d3udx3
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      ei = a5*DCOS(a7*(x - a6*t))
      dedx = -a5*a7*DSIN(a7*(x - a6*t))
      
      dhdx = -EXPPHI1*PHI*a1/a3
      d2hdx2 = -a1*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3hdx3 = PHI*a1*(-PHI**2 + 3*a3)*EXPPHI1/a3**3
      dudx = -EXPPHI1*PHI*a4/a3
      d2udx2 = -a4*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3udx3 = PHI*a4*(-PHI**2 + 3*a3)*EXPPHI1/a3**3

      dfluxGdx = -dedx*hi**2*(d2hdx2*ga*hi 
     . + dhdx**2*ga/2 + 2*dudx**2*hi) 
     . - 2*dhdx*ei*hi*(d2hdx2*ga*hi 
     . + dhdx**2*ga/2 + 2*dudx**2*hi) 
     . + dhdx*ga*hi 
     . + dudx*(-ei*(d2udx2*hi**3 + 3*dhdx*dudx*hi**2) + hi*ui) 
     . - ei*hi**2*(2*d2hdx2*dhdx*ga + 4*d2udx2*dudx*hi 
     .      + d3hdx3*ga*hi + 2*dhdx*dudx**2) 
     . + ui*(dedx*(-d2udx2*hi**3 - 3*dhdx*dudx*hi**2) 
     . + dhdx*ui + dudx*hi 
     . + ei*(-3*d2hdx2*dudx*hi**2 - 6*d2udx2*dhdx*hi**2 
     .      - d3udx3*hi**3 - 6*dhdx**2*dudx*hi))
     
      end function ForcedfluxGx
      
      function ForcedSourceG(x,t,a0,a1,a2,a3,a4,a5,a6,a7) result (Sa)
      
      real*8, intent(in) :: x,t,a0,a1,a2,a3,a4,a5,a6,a7
      real*8 :: dGdt,EXPPHI1,PHI,hi,ui,ei,Sa
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      ei = a5*DCOS(a7*(x - a6*t))
      
      Sa = -3*EXPPHI1*PHI**2*a1*a5*a6*a7*hi**2*ui
     .      *dsin(a6*a7*t - a7*x)/a3**2
     . + 3*EXPPHI1*PHI**2*a1*a5*a7*hi**2*ui**2
     .      *dsin(a6*a7*t - a7*x)/a3**2 
     . - PHI**2*a4*a5*a6*a7*hi**3*dexp(-PHI**2/(2*a3))
     .      *dsin(a6*a7*t - a7*x)/a3**2 
     . + PHI**2*a4*a5*a7*hi**3*ui*dexp(-PHI**2/(2*a3))
     .      *dsin(a6*a7*t - a7*x)/a3**2 
     . + 2*PHI**2*a5*a7*hi**3*ui**2*dsin(a6*a7*t - a7*x)/a3**2 
     . + a4*a5*a6*a7*hi**3*dexp(-PHI**2/(2*a3))*dsin(a6*a7*t - a7*x)/a3 
     . - a4*a5*a7*hi**3*ui*dexp(-PHI**2/(2*a3))*dsin(a6*a7*t - a7*x)/a3

     
      end function ForcedSourceG
      
      subroutine EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,dx,dt,
     . xbc,t,a0,a1,a2,a3,a4,a5,a6,a7)
           
      real*8,dimension(:),intent(inout) :: hbc,Gbc,ubc
      real*8, dimension(:),intent(in) :: xbc
      real*8, intent(in) :: dx,dt,a0,a1,a2,a3,a4,a5,a6,a7,t
      real*8, dimension(size(hbc)) :: hpbc,Gpbc,epbc
      real*8, dimension(size(hbc) - 2*n_GhstCells) :: hp,Gp
      integer :: n_GhstCells,nbc,n
      
      nbc = size(hbc)
      n = nbc - 2*n_GhstCells
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(xbc,hbc,Gbc,ubc,t,a5,a6,a7,dx,n_GhstCells)
      
      !Update cell averages
      call FluxStep(hbc,Gbc,ubc,n_GhstCells,dt,dx,xbc
     . ,t,a0,a1,a2,a3,a4,a5,a6,a7,hp,Gp)
    
      !BC's
      call Initloc(xbc(:n_GhstCells),
     .   t + dt,a0,a1,a2,a3,a4,a5,a6,a7,hpbc(:n_GhstCells),
     .   Gpbc(:n_GhstCells),ubc(:n_GhstCells),epbc(:n_GhstCells))
      
      hpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = hp(:)
      Gpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = Gp(:)
      
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,a0,a1,a2,a3,a4,a5,a6,a7,
     .   hpbc(nbc - n_GhstCells + 1:nbc),
     .   Gpbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc),
     .   epbc(nbc - n_GhstCells + 1:nbc))
     
      
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(xbc,hpbc,Gpbc,ubc,t + dt,a5,a6,a7,dx,n_GhstCells)
      
      !Update cell averages
      call FluxStep(hpbc,Gpbc,ubc,n_GhstCells,dt,dx,xbc
     . ,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,hp,Gp)
      
      !boundary conditions of hbc, Gbc stay same, just have to update interior
      hbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     . ( hbc(n_GhstCells + 1 : nbc - n_GhstCells ) + hp(:))/2
     
      Gbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     . ( Gbc(n_GhstCells + 1 : nbc - n_GhstCells ) + Gp(:))/2
     
      call Initloc(xbc(:n_GhstCells),
     .   t + dt,a0,a1,a2,a3,a4,a5,a6,a7,hbc(:n_GhstCells),
     .   Gbc(:n_GhstCells),ubc(:n_GhstCells),epbc(:n_GhstCells))
            
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,a0,a1,a2,a3,a4,a5,a6,a7,
     .   hbc(nbc - n_GhstCells + 1:nbc),
     .   Gbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc),
     .   epbc(nbc - n_GhstCells + 1:nbc))
      
      
      end subroutine EvolveStepWrap
       
      
      function ReconLinLimGrad(qim1,qi,qip1,theta) result (cdq)
      real*8, intent(in) :: qim1,qi,qip1,theta
      real*8 :: fdq,mdq,bdq,cdq
      
      fdq = qip1 - qi
      mdq = 0.5*(qip1 - qim1)
      bdq = qi - qim1
      
      cdq = minmod(theta*fdq,mdq,theta*bdq)
      
      end function ReconLinLimGrad
      
      function minmod(a,b,c) result (d)
      implicit none
      real*8, intent(in) :: a,b,c
      real*8 d
      
      if ((a > 0) .and. (b>0) .and. (c>0)) then
         d = min(a,b,c)
      else if ((a < 0) .and. (b<0) .and. (c<0)) then
         d = max(a,b,c)
      else
         d = 0.0
      end if
      
      end function


c  ********************************************************************************
c  Part 4 : Analysis Functions
c  Functions that do analyses on solutions
c  1. For Conservation -h, uh, G and Energy  (was at least second order accurate (FD solve, holding it back to second order for G))
c  2. TVD
c
c
c ********************************************************************************
      
c =====
c Function to generate quartic coefficients using
c f(x_{j-2}),f(x_{j-1}),f(x_{j}),f(x_{j+1}),f(x_{j+2})
c =====      
      function QuarticInterp(qjm2,qjm1,qj,qjp1,qjp2,dx)
     .  result (QuartCoeff)
     
      real*8, intent(in) :: qjm2,qjm1,qj,qjp1,qjp2,dx
      real*8,dimension(5) :: QuartCoeff
      
      QuartCoeff(1) = (qjp2 - 4*qjp1 + 6*qj - 4*qjm1 + qjm2) /
     . (24* (dx**4))
      QuartCoeff(2) = (qjp2 - 2*qjp1 + 2*qjm1 - qjm2) /
     . (12* (dx**3))  
      QuartCoeff(3) = (-qjp2 + 16*qjp1 - 30*qj + 16*qjm1 - qjm2) /
     . (24* (dx**2))  
      QuartCoeff(4) = (-qjp2 + 8*qjp1 - 8*qjm1 + qjm2) /
     . (12*dx)
      QuartCoeff(5) = qj
      end function QuarticInterp


c ========================
c Functions to evaluate Quartic at x, with quartic centered around x_j
c xmxj = x - x_j 
c ====================      
      function QuarticCoeffEvalxj(QuarticCoeff,xmxj) result (Val)
      
      real*8, dimension(5),intent(in) :: QuarticCoeff
      real*8, intent(in) :: xmxj
      real*8 :: Val
      
      Val = QuarticCoeff(1)*xmxj**4 + QuarticCoeff(2)*xmxj**3 +
     . QuarticCoeff(3)*xmxj**2 + QuarticCoeff(4)*xmxj +
     . QuarticCoeff(5)
     
      end function QuarticCoeffEvalxj
      
      function QuarticCoeffEvalGradxj(QuarticCoeff,xmxj) result (Val)
      
      real*8, dimension(5),intent(in) :: QuarticCoeff
      real*8, intent(in) :: xmxj
      real*8 :: Val
      
      Val = 4*QuarticCoeff(1)*xmxj**3 + 3*QuarticCoeff(2)*xmxj**2 +
     . 2*QuarticCoeff(3)*xmxj + QuarticCoeff(4)
     
      end function QuarticCoeffEvalGradxj

c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      function AllEnergiesIntegralCell(h,u,G,ga,j,dx) result (Val)
      
      real*8, dimension(:),intent(in) :: h,u,G
      real*8, intent(in) :: dx,ga
      integer, intent(in) :: j
      real*8 :: Val(4),fGPe(4),sGPe(4),tGPe(4)
      
      real*8 :: GPmxj,hGP,GGP,uGP,uxGP
      real*8 :: hCoeff(5), uCoeff(5), GCoeff(5)
      
      hCoeff = QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx)
      uCoeff = QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx)
      GCoeff = QuarticInterp(G(j-2),G(j-1),G(j),G(j+1),G(j+2),dx)
      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      hGP = QuarticCoeffEvalxj(hCoeff,GPmxj)
      GGP = QuarticCoeffEvalxj(GCoeff,GPmxj)
      uGP = QuarticCoeffEvalxj(uCoeff,GPmxj)
      uxGP = QuarticCoeffEvalGradxj(uCoeff,GPmxj)
      
      fGPe(1) = hGP
      fGPe(2) = GGP
      fGPe(3) = hGP*uGP
      fGPe(4) = (hGP*uGP**2 + ga*hGP**2 + (uxGP**2)*(hGP**3)/3d0)/2d0

      !second gauss point
      GPmxj = 0.0 
      hGP = QuarticCoeffEvalxj(hCoeff,GPmxj)
      GGP = QuarticCoeffEvalxj(GCoeff,GPmxj)
      uGP = QuarticCoeffEvalxj(uCoeff,GPmxj)
      uxGP = QuarticCoeffEvalGradxj(uCoeff,GPmxj)
      
      sGPe(1) = hGP
      sGPe(2) = GGP
      sGPe(3) = hGP*uGP
      sGPe(4) = (hGP*uGP**2 + ga*hGP**2 + (uxGP**2)*(hGP**3)/3d0)/2d0
      
      !first gauss point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2d0
      hGP = QuarticCoeffEvalxj(hCoeff,GPmxj)
      GGP = QuarticCoeffEvalxj(GCoeff,GPmxj)
      uGP = QuarticCoeffEvalxj(uCoeff,GPmxj)
      uxGP = QuarticCoeffEvalGradxj(uCoeff,GPmxj)
      
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + ga*hGP**2 + (uxGP**2)*(hGP**3)/3d0)/2d0
      
      Val(:) = (dx /2d0)*( (5.0/9.0)*fgpe(:) + (8.0/9.0)*sgpe(:) +
     . (5.0/9.0)*tgpe(:))
      
      end function AllEnergiesIntegralCell
      
      !Function to sum all energies
      function SumEnergy(h,u,G,ga,n_GhstCells,dx) result (SumEnergVals)
      
      real*8, dimension(:),intent(in) :: h,u,G
      real*8, intent(in) :: dx,ga
      integer, intent(in) :: n_GhstCells
      
      real*8 SumEnergVals(4),EnergVals(4)
      integer i,n
      
      n = size(h)
      SumEnergVals(:) = 0.0
      
      do i= n_GhstCells + 1, n - n_GhstCells
         EnergVals = AllEnergiesIntegralCell(h,u,G,ga,i,dx)
         SumEnergVals(:)  = SumEnergVals(:) + EnergVals(:)
      end do
      
      
      end function SumEnergy
      
      
 
      end subroutine Solve
      

c  ********************************************************************************
c  Main Program
c  Functions that loops over many dx (x_len), to investigate convergence
c  
c  The wdir is the working directory - all data printed from solver will end up in there
c
c
c ********************************************************************************
      program main
            
      implicit none
      
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/Models/GenSWWE/ForcedEpsFunc/expTimeDeriv1/"
  

      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 norm(3),C1EngMeas(4) 
      real*8 dx,startx,endx,ga,tstart,tend,Cr,
     . maxwavespeed,dtodx,a0,a1,a2,a3,a4,a5,a6,a7,
     . PI
      
      
      n=12
      n_gc = 6
      
      lowxlen = 100
      
      startx = -50d0
      endx = 100.0d0
      
      tstart = 0d0
      tend = 10.0d0
      
      PI=DACOS(-1.D0)
     
      ga = 9.81d0
      a0 = 1d0
      a1 = 0.5d0
      a2 = 5d0
      a3 = 20d0
      a4 = 0.3d0
      a5 = 2.0d0
      a6 = -0.5d0
      a7 = PI / 125.0d0
      
      
      maxwavespeed = a2 + a4 + dsqrt(ga*(a0 + a1))
      Cr = 0.5
      dtodx = Cr/maxwavespeed
            
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      open(3, file = wdir//'Norms.dat', status = 'unknown') 
      do i=1,n
         print *,'Experiment : ',i ,' || ', '# Cells :',
     .      lowxlen*(2**(i-1))
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         call Solve(startx,endx,lowxlen*(2**(i-1)),n_gc,tstart,tend,
     .   dtodx,ga,a0,a1,a2,a3,a4,a5,a6,a7,
     .   trim(fileloci),norm,dx)
         !write out
         write(3,*) dx,norm(1), norm(2), norm(3)
         
      end do
      
      close(3)  

      end program main
      
