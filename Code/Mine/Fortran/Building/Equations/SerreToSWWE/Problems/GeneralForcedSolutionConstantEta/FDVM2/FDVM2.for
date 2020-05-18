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
      subroutine Solve(x_start,x_end,x_len, n_GhstCells,tstart,tend,eta,
     . dtodx,theta,ga,a0,a1,a2,a3,a4,FilePre,norm,C1EngMeas,dx)
      implicit none
      
      
      real*8,intent(in) :: x_start,x_end,theta,
     . ga,a0,a1,a2,a3,a4,tstart,tend,dtodx,eta
      integer, intent(in) :: x_len, n_GhstCells
      character (LEN=*), intent(in)   :: FilePre
      character (len = 5) :: strti
      real*8, intent(out)   :: norm(3),C1EngMeas(4),dx
      
                         
      ! locations for x,u,h,G (cell nodes in interior and bc)
      real*8 xbc(x_len + 2*n_GhstCells),
     . hbc(x_len + 2*n_GhstCells),
     . Gbc(x_len + 2*n_GhstCells),
     . ubc(x_len + 2*n_GhstCells),
     . hAEbc(x_len + 2*n_GhstCells),
     . GAEbc(x_len + 2*n_GhstCells),
     . uAEbc(x_len + 2*n_GhstCells)
     
     
      real*8 AnaEnergy(4),CalcEnergy(4)
      real*8 ct,dt
    
      integer xbc_len,i,cti
      xbc_len = x_len + 2*n_GhstCells
      

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
      dt = dtodx*dx
            
      !Get locations of reconstruction plus Ghost Cells
      xbc = Generatexbc(x_start,dx,x_len,n_GhstCells)
      
            
      !analytic values for h,u,G at t = tstart
      call Initloc(xbc,tstart,eta,a0,a1,a2,a3,a4,hbc,Gbc,ubc)
      
      open(1, file = FilePre//'InitVars.dat', status = 'unknown')  
      do i=1,xbc_len
         write(1,*) xbc(i),hbc(i),Gbc(i),ubc(i)
      end do  
      close(1)
      
      ct = tstart
      cti = 0
      do while (ct < tend ) 
         call EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,eta,dx,dt,theta,
     .     xbc,ct,a0,a1,a2,a3,a4)
         ct = ct + dt
         cti = cti + 1
         print *, 'Current Time : ', ct
      end do
      
      !get u at end
      call GetufromhG(hbc,Gbc,ubc,eta,dx,n_GhstCells)
      
      !analytic values for h,u,G at t = tend
      call Initloc(xbc,ct,eta,a0,a1,a2,a3,a4,hAEbc,GAEbc,uAEbc)
      
   
      open(2, file = FilePre//'AllVarsEnd.dat', status = 'unknown')  
      do i=1,xbc_len
         write(2,*) xbc(i),hAEbc(i),hbc(i),GAEbc(i),Gbc(i),
     .    uAEbc(i),ubc(i)
      end do  
      close(2)
      
      !write out parameters of experiment
      open(4, file = FilePre//'Parameters.dat', status = 'unknown')  
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
      write(4,*) 'eta :' , eta
      write(4,*) 'gravity :' , ga
      write(4,*) 'a0 :' , a0
      write(4,*) 'a1 :' , a1
      write(4,*) 'a2 :' , a2
      write(4,*) 'a3 :' , a3
      write(4,*) 'a4 :' , a4
      close(4)

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
      subroutine Initloc(x,t,eta,a0,a1,a2,a3,a4,h,G,u)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,intent(in) :: t,eta,a0,a1,a2,a3,a4
      real*8,dimension(size(x)), intent(out)   :: h,G,u
                   
      real*8 PHI,EXPPHI1,dh3uxo3dx,dhdx,dudx,d2udx2
      
      integer i,n
           
      n= size(x)      
      do i=1,n 
         
         PHI  = x(i) - a2*t
         EXPPHI1 = dexp(-PHI**2 / (2*a3))
         h(i)  = a0 + a1*EXPPHI1
         u(i)  = a4*EXPPHI1
         
c         dhdx = -a1*PHI*EXPPHI1/a3
c         dudx = -a4*PHI*EXPPHI1/a3
c         d2udx2 = a4*( (PHI**2)*EXPPHI1/(a3**2) - EXPPHI1/(a3) )
c         dh3uxo3dx = h(i)**2*( dhdx*dudx + h(i)*d2udx2/3)
c                             
c         G(i) = u(i)*h(i) - eta*dh3uxo3dx

          G(i) = (3*a3**2 - eta*(PHI**2*(3*EXPPHI1*a1 + h(i)) -
     .    a3*h(i))*h(i))*h(i)*u(i)/(3*a3**2)

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
      subroutine GetufromhG(hbc,Gbc,ubc,eta,dx,n_GhstCells)
      implicit none
      
      real*8,dimension(:),intent(in) :: hbc,Gbc
      real*8,dimension(:), intent(inout)   :: ubc
      real*8,intent(in) :: dx,eta
      integer,intent(in) :: n_GhstCells
                   
      real*8 subdiag1(size(hbc)),
     . diag(size(hbc)),
     . supdiag1(size(hbc)),
     . RHS(size(hbc))
          
      real*8 ht1,ht2
     
      integer n
      
      n = size(hbc)
                  
      !calculate diagonals in interior
      !set RHS B
            
      do i=n_GhstCells+1,n - n_GhstCells 
         ht1 = eta*(hbc(i)**3/(3d0*dx*dx))
         ht2 = eta*(hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1)))
         
         
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

      subroutine FluxStep(hbc,Gbc,ubc,n_GhstCells,eta,dt,dx,theta,
     . xbc,t,a0,a1,a2,a3,a4,hp,Gp)
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: theta,dt,dx,eta,t,a0,a1,a2,a3,a4
      real*8,dimension(:),intent(in) :: xbc,hbc,Gbc,ubc
      
      real*8, dimension(size(xbc) - 2*n_GhstCells),
     .   intent(out) :: hp,Gp
     
     
      real*8 cdhi,cdGi,cdhip1,cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l, fih,fiG,foh,foG,
     . dhdt, dGdt, dfluxhdx, dfluxGdx
      
      integer n,nbc,i
      
      nbc = size(xbc)
      n = nbc - 2*n_GhstCells
      
      !Do left boundary
      i = n_GhstCells
      
      !reconstruct values on left side of edge x_{j-1/2}
      cdhi = ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta)
      cdGi = ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta)
      
      hir = hbc(i) + cdhi/2
      Gir = Gbc(i) + cdGi/2
      uir = (ubc(i+1)+ubc(i))/2
      duir = (2*ubc(i) - 3*ubc(i-1) + ubc(i-2)) /dx
      
      !reconstruct values on right side of edge x_{j-1/2}
      cdhip1 = ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta)
      cdGip1 = ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta)
      
      hip1l = hbc(i + 1) - cdhip1/2
      Gip1l = Gbc(i + 1) - cdGi /2
      uip1l = (ubc(i+1)+ubc(i))/2
      duip1l = (-2*ubc(i+1) + 3*ubc(i+2) - ubc(i+3)) /dx
      
      ! speed bounds
      sl  = min(0d0, uir - dsqrt(ga*hir) , uip1l - dsqrt(ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(ga*hir) , uip1l + dsqrt(ga*hip1l)  )
      
      !left and right flux
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 - eta*(2d0/3d0)*hir**3*duir**2
      
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 - eta*
     . (2d0/3d0)*hip1l**3*duip1l**2
     
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
      
         !use slopes from last calculation         
         hir = hbc(i) + cdhi/2
         Gir = Gbc(i) + cdGi/2
         uir = (ubc(i+1)+ubc(i))/2
         duir = (2*ubc(i) - 3*ubc(i-1) + ubc(i-2)) /dx
         
         !reconstruct values on right side of edge x_{i-1/2}
         cdhip1 = ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta)
         cdGip1 = ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta)
         
         hip1l = hbc(i + 1) - cdhip1/2
         Gip1l = Gbc(i + 1) - cdGi /2
         uip1l = (ubc(i+1)+ubc(i))/2
         duip1l = (-2*ubc(i+1) + 3*ubc(i+2) - ubc(i+3)) /dx
         
         ! speed bounds
         sl  = min(0d0, uir - dsqrt(ga*hir) , uip1l - dsqrt(ga*hip1l)  )
         sr  = max(0d0, uir + dsqrt(ga*hir) , uip1l + dsqrt(ga*hip1l)  )
         
         !left and right flux
         felh = uir*hir
         felG = uir*Gir + ga*(hir**2)/2d0 - eta*(2d0/3d0)*hir**3*duir**2
         
         ferh = uip1l*hip1l
         ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 -
     .      eta*(2d0/3d0)*hip1l**3*duip1l**2
        
         if (sr == sl) then
            isrmsl = 0.0
         else
            isrmsl = 1.0 / (sr - sl)
         end if 
         
         !calculate flux from cell i to cell i + 1 (Kurganov)
         foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
         foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
         
         dhdt = Forcedht(xbc(i),t,eta,a0,a1,a2,a3,a4) 
         dGdt = ForcedGt(xbc(i),t,eta,a0,a1,a2,a3,a4) 
         
         dfluxhdx = Forcedfluxhx(xbc(i),t,eta,a0,a1,a2,a3,a4) 
         dfluxGdx = ForcedfluxGx(xbc(i),t,ga,eta,a0,a1,a2,a3,a4)         
         
         !Update cell averages
c         hp(i - n_GhstCells ) = hbc(i) - dt/dx*(foh - fih) 
c     .         - dt*(dhdt + dfluxhdx  )
c         Gp(i - n_GhstCells ) = Gbc(i) - dt/dx*(foG - fiG)
c     .         - dt*(dGdt + dfluxGdx  )

         hp(i - n_GhstCells ) = hbc(i) - dt/dx*(foh - fih)  
     .      + dt*(dhdt + dfluxhdx ) 
         Gp(i - n_GhstCells ) = Gbc(i) - dt/dx*(foG - fiG) 
     .      + dt*(dGdt + dfluxGdx )       
         
         !The flux from cell i -1 to cell i, in the next iteration where i -> i +1
         ! is the same as the flux we just calculated
         fih = foh
         fiG = foG
         
         cdhi = cdhip1
         cdGi = cdGip1
      end do

      end subroutine FluxStep 
      
      function Forcedht(x,t,eta,a0,a1,a2,a3,a4) result (dhdt)
      
      real*8, intent(in) :: x,t,eta,a0,a1,a2,a3,a4
      real*8 :: dhdt,EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      
      dhdt = EXPPHI1*PHI*a1*a2/a3
      
      end function Forcedht
      
      function ForcedGt(x,t,eta,a0,a1,a2,a3,a4) result (dGdt)
      
      real*8, intent(in) :: x,t,eta,a0,a1,a2,a3,a4
      real*8 :: dGdt,EXPPHI1,PHI,hi,ui
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1

      
      dGdt = a2*ui*(3*PHI*a3**2*(EXPPHI1*a1 + hi) - 
     .   eta*hi*(PHI**3*(6*EXPPHI1**2*a1**2 + 9*EXPPHI1*a1*hi + hi**2) 
     .   + a3*hi*(-9*EXPPHI1*PHI*a1 - 3*PHI*hi)))/(3*a3**3)

      end function ForcedGt
      
      function Forcedfluxhx(x,t,eta,a0,a1,a2,a3,a4) 
     . result (dfluxhdx)
      
      real*8, intent(in) :: x,t,eta,a0,a1,a2,a3,a4
      real*8 :: dfluxhdx,EXPPHI1,PHI,hi,ui
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      
      dfluxhdx = -PHI*ui*(EXPPHI1*a1 + hi)/a3
      
      end function Forcedfluxhx
      
      function ForcedfluxGx(x,t,ga,eta,a0,a1,a2,a3,a4) 
     . result (dfluxGdx)
      
      real*8, intent(in) :: x,t,ga,eta,a0,a1,a2,a3,a4
      real*8 :: dfluxGdx,EXPPHI1,PHI,hi,ui
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      
      dfluxGdx = -PHI*(3*EXPPHI1*a1*a3**2*ga*hi
     .   - 2*PHI**2*eta*hi**2*ui**2*(3*EXPPHI1*a1 + 2*hi) 
     .   + 4*a3*eta*hi**3*ui**2 - hi*ui**2*(-3*a3**2 + 
     .   eta*hi*(PHI**2*(3*EXPPHI1*a1 + hi) - a3*hi)) +
     .   ui**2*(3*a3**2*(EXPPHI1*a1 + hi) - eta*hi*
     .   (PHI**2*(6*EXPPHI1**2*a1**2 + 9*EXPPHI1*a1*hi + hi**2)
     .    - 3*a3*hi*(3*EXPPHI1*a1 + hi))))/(3*a3**3)
     
      end function ForcedfluxGx
      
      subroutine EvolveStepWrap(hbc,Gbc,ubc,n_GhstCells,eta,dx,dt,
     . theta,xbc,t,a0,a1,a2,a3,a4)
           
      real*8,dimension(:),intent(inout) :: hbc,Gbc,ubc
      real*8, dimension(:),intent(in) :: xbc
      real*8, intent(in) :: dx,dt,theta,eta,a0,a1,a2,a3,a4,t
      real*8, dimension(size(hbc)) :: hpbc,Gpbc
      real*8, dimension(size(hbc) - 2*n_GhstCells) :: hp,Gp
      integer :: n_GhstCells,nbc,n
      
      nbc = size(hbc)
      n = nbc - 2*n_GhstCells
      
      
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(hbc,Gbc,ubc,eta,dx,n_GhstCells)
      
      !Update cell averages
      call FluxStep(hbc,Gbc,ubc,n_GhstCells,eta,dt,dx,theta,xbc
     . ,t,a0,a1,a2,a3,a4,hp,Gp)
    
      
      !get boundary conditions for hpbc = h' with ghost cells      
c      hpbc(:n_GhstCells) = hbc(:n_GhstCells)
c      Gpbc(:n_GhstCells) = Gbc(:n_GhstCells)

      !BC's
      call Initloc(xbc(:n_GhstCells),
     .   t + dt,eta,a0,a1,a2,a3,a4,hpbc(:n_GhstCells),
     .   Gpbc(:n_GhstCells),ubc(:n_GhstCells))
      
      hpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = hp(:)
      Gpbc(n_GhstCells + 1 : nbc - n_GhstCells ) = Gp(:)
      
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,eta,a0,a1,a2,a3,a4,
     .   hpbc(nbc - n_GhstCells + 1:nbc),
     .   Gpbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc))
     
      !get boundary conditions for hpbc = h' with ghost cells
c      hpbc(nbc - n_GhstCells + 1:nbc) = hbc(nbc - n_GhstCells + 1:nbc)
c      Gpbc(nbc - n_GhstCells + 1:nbc) = Gbc(nbc - n_GhstCells + 1:nbc)
      
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(hpbc,Gpbc,ubc,eta,dx,n_GhstCells)
      
      !Update cell averages
      call FluxStep(hpbc,Gpbc,ubc,n_GhstCells,eta,dt,dx,theta,xbc
     . ,t + dt,a0,a1,a2,a3,a4,hp,Gp)
      
      !boundary conditions of hbc, Gbc stay same, just have to update interior
      hbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     . ( hbc(n_GhstCells + 1 : nbc - n_GhstCells ) + hp(:))/2
     
      Gbc(n_GhstCells + 1 : nbc - n_GhstCells ) = 
     . ( Gbc(n_GhstCells + 1 : nbc - n_GhstCells ) + Gp(:))/2
     
      call Initloc(xbc(:n_GhstCells),
     .   t + dt,eta,a0,a1,a2,a3,a4,hbc(:n_GhstCells),
     .   Gbc(:n_GhstCells),ubc(:n_GhstCells))
            
      !BC's
      call Initloc(xbc(nbc - n_GhstCells + 1:nbc),
     .   t + dt,eta,a0,a1,a2,a3,a4,
     .   hbc(nbc - n_GhstCells + 1:nbc),
     .   Gbc(nbc - n_GhstCells + 1:nbc),
     .   ubc(nbc - n_GhstCells + 1:nbc))
      
      
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
      
      CHARACTER(*), PARAMETER :: wdir = "~/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/Models/Serre2SWWE/Forced/eta0/"
  
      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 norm(3),C1EngMeas(4) 
      real*8 dx,theta,startx,endx,ga,tstart,tend,Cr,
     . maxwavespeed,dtodx,eta,a0,a1,a2,a3,a4
      
      n=12
      n_gc = 3
      
      lowxlen = 100
      
      startx = -50d0
      endx = 100.0d0
      
      tstart = 0.0d0
      tend = 10d0
      
      !eta controls dispersion (1 = Serre, 0 is SWWE)
      eta = 0.0
      
      theta = 1.2d0
      
      ga = 9.81d0
      a0 = 1d0
      a1 = 1d0
      a2 = 5d0
      a3 = 10d0
      a4 = 1d0
      
      maxwavespeed = a2 + a4 + dsqrt(ga*(a0 + a1))
      Cr = 0.5
      dtodx = Cr/maxwavespeed
            
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      open(3, file = wdir//'Norms.dat', status = 'unknown') 
      open(4, file = wdir//'Energy.dat', status = 'unknown') 
      do i=1,n
         print *,'Experiment : ',i ,' || ', '# Cells :',
     .      lowxlen*(2**(i-1))
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         call Solve(startx,endx,lowxlen*(2**(i-1)),n_gc,tstart,tend,
     .   eta,dtodx,theta,ga,a0,a1,a2,a3,a4,
     .   trim(fileloci),norm,C1EngMeas,dx)
         !write out
         write(3,*) dx,norm(1), norm(2), norm(3)
         write(4,*) dx,C1EngMeas(1), C1EngMeas(2), C1EngMeas(3),
     .    C1EngMeas(4)
         
      end do
      
      close(3)  
      close(4) 

      end program main
      
