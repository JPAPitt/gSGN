c=====================================
c Program that given array of h calculates h* the reconstructed array
c Analysis shows second order until 10^-10 by which point we get round-off error effects
c=====================================
      subroutine ReconAll(x_start,x_end,x_len, n_GhstCells,
     . theta,ga,a0,a1,a2,a3,a4,FilePre,norm,C1EngMeas,dx)
      implicit none
      
      real*8,intent(in) :: x_start,x_end,theta,ga,a0,a1,a2,a3,a4
      integer, intent(in) :: x_len, n_GhstCells
      character (LEN=*), intent(in)   :: FilePre
      real*8, intent(out)   :: norm(4),C1EngMeas(4),dx
                         
      ! locations for x,u,h,G (cell nodes in interior and bc)
      real*8 xbc(x_len + 2*n_GhstCells),
     . hbc(x_len + 2*n_GhstCells),
     . Gbc(x_len + 2*n_GhstCells),
     . uabc(x_len + 2*n_GhstCells),
     . ubc(x_len + 2*n_GhstCells)
     
      !reconstructions
      real*8 xR(3*x_len),
     . hR(3*x_len),
     . GR(3*x_len),
     . uR(3*x_len),
     . duR(3*x_len),
     . hRA(3*x_len),
     . GRA(3*x_len),
     . uRA(3*x_len),
     . duRA(3*x_len)
     
      real*8 AnaEnergy(4),CalcEnergy(4)
    
      integer xr_len,xbc_len,i
      xr_len = 3*x_len
      xbc_len = x_len + 2*n_GhstCells

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
            
      !Get locations of reconstruction plus Ghost Cells
      xbc = Generatexbc(x_start,dx,x_len,n_GhstCells)
      
      !Get h,u and G at all these cell nodes
      call GenInitloc(xbc,dx,a0,a1,a2,a3,a4,hbc,Gbc,uabc)
      
      !use analytic values
      !call Initloc(xbc,a0,a1,a2,a3,a4,hbc,Gbc,ubc)
      
      !get initial u (mainly to get BC)
      ubc(:) = uabc(:)
      
      !Get ubc, from h and G from FD method
      call GetufromhG(hbc,Gbc,ubc,dx,n_GhstCells)
      
      !DoReconstruction at cell edges
      call ReconsAll(xbc,hbc,Gbc,ubc,n_GhstCells,dx,theta,
     . xR,hR,GR,uR,duR)
     
      !Get analytic values for recons
      call Initloc(xR,a0,a1,a2,a3,a4,hRA,GRA,uRA)
      
      !get du, at all locations
      duRA(:) = a3*a4*DCOS(a4*xR(:))
      
      open(1, file = FilePre//'All.dat', status = 'unknown')  
      do i=1,xr_len
         write(1,*) xR(i),hRA(i),hR(i),GRA(i),GR(i),uRA(i),uR(i),
     . duRA(i),duR(i)    
      end do  
      close(1)

      !Convergence Norms
      norm(1) = sqrt(sum((hRA(:) - hR(:)) **2)) / 
     . sqrt(sum(hRA(:)**2))       
      norm(2) = sqrt(sum((GRA(:) - GR(:)) **2)) / 
     . sqrt(sum(GRA(:)**2)) 
      norm(3) = sqrt(sum((uRA(:) - uR(:)) **2)) / 
     . sqrt(sum(uRA(:)**2)) 
      norm(4) = sqrt(sum((duRA(:) - duR(:)) **2)) / 
     . sqrt(sum(duRA(:)**2)) 
     
      !Analytic Energy Values
      AnaEnergy(:) = IntegralCons(x_start - dx/2,x_end + dx/2
     .,a0,a1,a2,a3,a4) 
            
      !Conservation Measures
      CalcEnergy(:) = SumEnergy(hbc,ubc,Gbc,ga,n_GhstCells,dx)
      
      C1EngMeas(:) = abs( AnaEnergy(:) - CalcEnergy(:) ) /
     . abs(AnaEnergy(:))
     
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
c Function to generate h,u,G at all locations x
c ====================      
      subroutine Initloc(x,a0,a1,a2,a3,a4,h,G,u)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,dimension(size(x)), intent(out)   :: h,G,u
      real*8 :: hx,ux,uxx
                   
      real*8 a0,a1,a2,a3,a4
      
      integer i,n
      
      n= size(x)
      
      do i=1,n 
         h(i)  = a0 + a1*DSIN(a2*x(i) ) 
         u(i)  = a3*DSIN(a4*x(i) ) 
         
         hx = a1*a2*DCOS(a2*x(i))
         ux = a3*a4*DCOS(a4*x(i))

         uxx = -a3*a4**2*DSIN(a4*x(i))
               
         G(i) = u(i)*h(i) - h(i)**2*hx*ux - h(i)**3*uxx/3.0
         
      end do 
            
      end subroutine Initloc

c ====================
c Function to generate h,u,G at all locations x
c ====================      
      subroutine GenInitloc(x,dx,a0,a1,a2,a3,a4,h,G,u)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,dimension(size(x)), intent(out)   :: h,G,u
                   
      real*8 dx,a0,a1,a2,a3,a4,hbeg,ubeg,hend,uend
      
      integer i,n
      
      n= size(x)
      
      do i=1,n 
         h(i)  = a0 + a1*DSIN(a2*x(i) ) 
         u(i)  = a3*DSIN(a4*x(i) )         
      end do 
      
      !left boundary
      hbeg = a0 + a1*DSIN(a2*(x(1) - dx) )
      ubeg  = a3*DSIN(a4*(x(1) - dx) ) 
      
      !right boundary
      hend = a0 + a1*DSIN(a2*(x(n) + dx) )
      uend  = a3*DSIN(a4*(x(n) + dx) ) 
      
      G(:) = GetGfromuh(h,u,hbeg,ubeg,hend,uend)
            
      end subroutine GenInitloc


c ====================
c Function to generate h,u at all locations x and then calculate G
c ====================      
      function GetGfromuh(h,u,hbeg,ubeg,hend,uend) result (G)
      implicit none
      
      real*8,dimension(:),intent(in) :: h,u
      real*8,intent(in)  :: hbeg,ubeg,hend,uend
      real*8,dimension(size(h))  :: G
      
      !i is for loop
      !n is length of h,u
      integer i,n
      
      real*8 ai,bi,ci,ht1,ht2
      
      n= size(h)
      
      !Left boundary
      i = 1
      ht1 = hbc(i)**3/(3d0*dx*dx)
      ht2 = hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbeg)
         
         
      ai  = -ht1 + ht2
      bi = hbc(i) + 2d0*ht1
      ci  = -ht1 - ht2   
                  
      G(i) = ai*ubeg + bi*u(i) + ci*u(i+1)
      
      !interior
      do i=2,n-1    
         ht1 = hbc(i)**3/(3d0*dx*dx)
         ht2 = hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1))
         
         
         ai  = -ht1 + ht2
         bi = hbc(i) + 2d0*ht1
         ci  = -ht1 - ht2   
                     
         G(i) = ai*u(i-1) + bi*u(i) + ci*u(i+1)
         
      end do 
      
      !Right boundary
      i = n
      ht1 = hbc(i)**3/(3d0*dx*dx)
      ht2 = hbc(i)**2/(4.d0*dx*dx)*(hend - hbc(i-1))
         
         
      ai  = -ht1 + ht2
      bi = hbc(i) + 2d0*ht1
      ci  = -ht1 - ht2   
                  
      G(i) = ai*u(i-1) + bi*u(i) + ci*uend
            
      end function GetGfromuh

c ============
c Analytic values for Integrals of initial conditions
c ==========
      function IntegralCons(xbeg,xend,a0,a1,a2,a3,a4) 
     . result (Integrals)
      implicit none
      
      real*8,intent(in) :: xbeg,xend,a0,a1,a2,a3,a4
      real*8 Integrals(4),h3ux,uh,Hp1,Hp2,Hp3

            
      ! h integral      
      Integrals(1)  =  a0*(xend) -a1/a2*DCOS(a2* (xend)) 
     . - a0*(xbeg) + a1/a2*DCOS(a2* (xbeg))
     
      ! uh integral
      uh = -a0*a3/a4*DCOS(a4*xend)
     . + a1*a3/(a4**2 - a2**2)*(a2*DSIN(a4*xend)
     . *DCOS(a2*xend) - a4*DCOS(a4*xend)
     . *DSIN(a2*xend) )
     . + a0*a3/a4*DCOS(a4*xbeg)
     . - a1*a3/(a4**2 - a2**2)*(a2*DSIN(a4*xbeg)
     . *DCOS(a2*xbeg) - a4*DCOS(a4*xbeg)
     . *DSIN(a2*xbeg))
      !G integral
     
      ! integrating (h^3ux/3)x  gives h^3ux/3
      h3ux = (a0 + a1*DSIN(a2*xend))**3 *
     . a3*a4*DCOS(a4*xend) -
     . (a0 + a1*DSIN(a2*xbeg ))**3 *
     . a3*a4*DCOS(a4*xbeg)
     
      Integrals(2) = uh - h3ux/3 
      
      Integrals(3) = uh
      
      Hp1 = a3**2*(a0*a2**3*a4*xend*dsin(a4*xend)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) + a0*a2**3*a4*xend*
     . dcos(a4*xend)**2/(2*a2**3*a4 - 8*a2*a4**3) -
     . a0*a2**3*dsin(a4*xend)*dcos(a4*xend)/
     . (2*a2**3*a4 - 8*a2*a4**3) - 4*a0*a2*a4**3*xend*dsin(a4*xend)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) - 4*a0*a2*a4**3*xend*cos(a4*xend)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) + 4*a0*a2*a4**2*dsin(a4*xend)*
     . dcos(a4*xend)/(2*a2**3*a4 - 8*a2*a4**3) - 2*a1*a2**2*a4*
     . dsin(a4*xend)**2*dcos(a2*xend)/(2*a2**3*a4 - 8*a2*a4**3) +
     . 4*a1*a2*a4**2*dsin(a2*xend)*dsin(a4*xend)*dcos(a4*xend)
     . /(2*a2**3*a4 - 8*a2*a4**3) + 4*a1*a4**3*dsin(a4*xend)**2*
     . dcos(a2*xend)/(2*a2**3*a4 - 8*a2*a4**3) +
     . 4*a1*a4**3*dcos(a2*xend)*dcos(a4*xend)**2/
     . (2*a2**3*a4 - 8*a2*a4**3)) -
     . a3**2*(a0*a2**3*a4*xbeg*dsin(a4*xbeg)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) + a0*a2**3*a4*xbeg*dcos(a4*xbeg)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) - a0*a2**3*dsin(a4*xbeg)*dcos(a4*xbeg)/
     . (2*a2**3*a4 - 8*a2*a4**3) - 4*a0*a2*a4**3*xbeg*dsin(a4*xbeg)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) - 4*a0*a2*a4**3*xbeg*dcos(a4*xbeg)**2/
     . (2*a2**3*a4 - 8*a2*a4**3) + 4*a0*a2*a4**2*dsin(a4*xbeg)*
     . dcos(a4*xbeg)/(2*a2**3*a4 - 8*a2*a4**3) -
     . 2*a1*a2**2*a4*dsin(a4*xbeg)**2*dcos(a2*xbeg)/
     . (2*a2**3*a4 - 8*a2*a4**3) + 4*a1*a2*a4**2*dsin(a2*xbeg)*
     . dsin(a4*xbeg)*dcos(a4*xbeg)/(2*a2**3*a4 - 8*a2*a4**3) +
     . 4*a1*a4**3*dsin(a4*xbeg)**2*dcos(a2*xbeg)/
     . (2*a2**3*a4 - 8*a2*a4**3) + 4*a1*a4**3*dcos(a2*xbeg)*
     . dcos(a4*xbeg)**2/(2*a2**3*a4 - 8*a2*a4**3))
     
      Hp2 = -ga*(a0**2*xbeg - 2*a0*a1*dcos(a2*xbeg)/a2 +
     . a1**2*xbeg*dsin(a2*xbeg)**2/2 + a1**2*xbeg*cos(a2*xbeg)**2/2 -
     . a1**2*dsin(a2*xbeg)*dcos(a2*xbeg)/(2*a2)) + 
     . ga*(a0**2*xend - 2*a0*a1*dcos(a2*xend)/a2 +
     . a1**2*xend*dsin(a2*xend)**2/2 + a1**2*xend*dcos(a2*xend)**2/2 -
     . a1**2*dsin(a2*xend)*dcos(a2*xend)/(2*a2))
     

      
      Integrals(4) = (Hp1 + Hp2)/2.0

      
      end function IntegralCons   


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
      subroutine GetufromhG(hbc,Gbc,ubc,dx,n_GhstCells)
      implicit none
      
      real*8,dimension(:),intent(in) :: hbc,Gbc
      real*8,dimension(:), intent(inout)   :: ubc
      real*8,intent(in) :: dx
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
         ht1 = hbc(i)**3/(3d0*dx*dx)
         ht2 = hbc(i)**2/(4.d0*dx*dx)*(hbc(i+1) - hbc(i-1))
         
         
         subdiag1(i)  = -ht1 + ht2
         diag(i) = hbc(i) + 2d0*ht1
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
c  Part 3 : Reconstruction Functions
c  Functions that Reconstruct h,u,G,du/dx at all points inside cell
c
c
c
c ********************************************************************************
      
      subroutine ReconsAll(xbc,hbc,Gbc,ubc,n_GhstCells,dx,theta
     . ,xR,hR,GR,uR,duR)
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: dx,theta
      real*8,dimension(:),intent(in) :: xbc,hbc,Gbc,ubc
      
      real*8, dimension(3*(size(xbc) - 2*n_GhstCells)),
     .   intent(out) :: xR,hR,GR,uR ,duR
     
     
      real*8 cdh,cdG
      integer n,nbc,nR,i,j
      
      nbc = size(xbc)
      n = nbc - 2*n_GhstCells
      nR = 3*n
      
      ! i is numbering by xbc (cells and ghost cells)
      do i = n_GhstCells + 1, nbc - n_GhstCells
      
         !j is numbering by equivlent x (cells)
         j = (i -(n_GhstCells + 1) )
         
         !x recon
         xR(3*j + 1) = xbc(i) - 0.5*dx
         xR(3*j + 2) = xbc(i)
         xR(3*j + 3) = xbc(i) + 0.5*dx
         
         !h recon
         cdh = ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta)
         hR(3*j + 1) = hbc(i) - 0.5*cdh
         hR(3*j + 2) = hbc(i)
         hR(3*j + 3) = hbc(i) + 0.5*cdh
         
         !G recon
         cdG = ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta)
         GR(3*j + 1) = Gbc(i) - 0.5*cdG
         GR(3*j + 2) = Gbc(i)
         GR(3*j + 3) = Gbc(i) + 0.5*cdG
         
         !u recon
         uR(3*j + 1) = 0.5*(ubc(i-1) + ubc(i))
         uR(3*j + 2) = ubc(i)
         uR(3*j + 3) = 0.5*(ubc(i) + ubc(i+1))
         
         !du recon
         !forwards 2nd order 1st derivative
         duR(3*j + 1) = (-2*ubc(i) + 3*ubc(i+1) - ubc(i+2)) /dx
         
         duR(3*j + 2) = (ubc(i+1) - ubc(i-1)) / (2d0*dx)
         
         !backwards 2nd order 1st derivative
         duR(3*j + 3) = (2*ubc(i) - 3*ubc(i-1) + ubc(i-2)) /dx
      end do
     
      
      end subroutine ReconsAll
      
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
      
      
 
      end subroutine ReconAll
      

c  ********************************************************************************
c  Main Program
c  Functions that loops over many dx (x_len), to investigate convergence
c
c
c
c
c ********************************************************************************
      program main
            
      implicit none
      
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/FortranTests/EnerguhtoGtou/"
  
      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 norm(4),C1EngMeas(4) 
      real*8 PI,a0,a1,a2,a3,a4,dx,theta,startx,endx,ga
      
      n= 18
      n_gc = 3
      
      lowxlen = 100
      
      startx = 0.0
      endx = 100.0
      
      theta = 1.2
      
      ga = 9.81
      
      a0 = 1d0
      a1 = 0.5d0
      a3 = 0.1d0
      
      PI = 2*DASIN(1.0d0)
      a2 = 2*PI/ 20.0
      a4 =  0!2*PI/ 10.0
      
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      open(3, file = wdir//'Norms.dat', status = 'unknown') 
      open(4, file = wdir//'Energy.dat', status = 'unknown') 
      do i=1,n
         print *,'Reconstruction',i,lowxlen*(2**(i-1))
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         call ReconAll(startx,endx,lowxlen*(2**(i-1)),n_gc,theta,ga,a0,
     .   a1,a2,a3,a4,trim(fileloci),norm,C1EngMeas,dx)


         !write out
         write(3,*) dx,norm(1), norm(2), norm(3), norm(4)
         write(4,*) dx,C1EngMeas(1), C1EngMeas(2), C1EngMeas(3),
     .    C1EngMeas(4)
         
      end do
      
      close(3)  
      close(4) 

      end program main
      
