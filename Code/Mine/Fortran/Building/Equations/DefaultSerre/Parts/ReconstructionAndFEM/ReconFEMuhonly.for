c=====================================
c Program that given array of h calculates h* the reconstructed array
c Analysis shows second order until 10^-10 by which point we get round-off error effects
c=====================================
      subroutine ReconAll(x_start,x_end,x_len,n_PtsCell, n_GhstCells,
     . theta,a0,a1,a2,a3,a4,FilePre,norm,dx)
      implicit none
      
      real*8,intent(in) :: x_start,x_end,theta,a0,a1,a2,a3,a4
      integer, intent(in) :: x_len,n_PtsCell, n_GhstCells
      character (LEN=*), intent(in)   :: FilePre
      real*8, intent(out)   :: norm(3),dx
                   
      real*8 x(x_len),hA(x_len),GA(x_len)
      
      ! locations for x,h,G
      real*8 xhGR(n_PtsCell*(x_len + 2*n_GhstCells)),
     . haR(n_PtsCell*(x_len + 2*n_GhstCells)),
     . GaR(n_PtsCell*(x_len + 2*n_GhstCells)),
     . hR(n_PtsCell*(x_len + 2*n_GhstCells)),
     . GR(n_PtsCell*(x_len + 2*n_GhstCells))
     
      !location arrays for x,u
      !length of u, is  different, because edges of cells are equal
      real*8 xuR((n_PtsCell -1)*x_len -1 + 2*n_PtsCell*n_GhstCells),
     . uaR((n_PtsCell -1)*x_len -1 + 2*n_PtsCell*n_GhstCells),
     . uR((n_PtsCell -1)*x_len -1 + 2*n_PtsCell*n_GhstCells) 

     
      real*8 hbeg(n_PtsCell * n_GhstCells),
     . hend(n_PtsCell * n_GhstCells ),
     . ubeg(n_PtsCell * n_GhstCells),
     . uend(n_PtsCell * n_GhstCells ),
     . Gbeg(n_PtsCell * n_GhstCells),
     . Gend(n_PtsCell * n_GhstCells )
      
      integer xhGR_len,xuR_len,nGP, i
      
            
      xhGR_len = n_PtsCell*(x_len + 2*n_GhstCells)
      xuR_len = (n_PtsCell -1)*x_len -1 + 2*n_PtsCell*n_GhstCells
      nGP = n_PtsCell * n_GhstCells 

      !Initial data
      dx = (x_end - x_start) / (x_len) 
            
      !Get locations of reconstruction plus Ghost Cells
      xhGR = GenerateAllLoc(x_start,dx,x_len,xhGR_len,n_PtsCell,
     . n_GhstCells)
     
      !get cell nodes
      x = xhGR(nGP + 2 : xhGR_len-nGP : 3)
      
      !get locations of u from locations of h,G
      !Get left boundary
      xuR(:nGP) = xhGR(:nGP)
      !Get interior
      !cell nodes
      xuR(nGP + 1: xuR_len - nGP:2) =  xhGR(nGP + 2 : xhGR_len-nGP : 3)
      
      !cell edges
      xuR(nGP + 2: xuR_len - (nGP + 1):2) =
     .   xhGR(nGP + 3 : xhGR_len-(nGP+1) : 3)
     
      ! right boundary

      xuR(xuR_len-(nGP-1):xuR_len) = xhGR(xhGR_len- (nGP-1):xhGR_len)
      
      
      !Calculate Cell Averages for initial conditions      
      call InitCellAvg(x,dx,a0,a1,a2,a3,a4,hA,GA)
      
      !Get h,u and G at all locations for initial conditions
      call Initloc(xhGR,haR,GaR,xuR,uaR,a0,a1,a2,a3,a4)
      
      !use midpoints instead of cell avg
      hA = haR(nGP + 2 : xhGR_len-nGP : 3)
      GA = GaR(nGP + 2 : xhGR_len-nGP : 3)
      
      !Write Out input data
      open(1, file = FilePre//'CellAvg.dat', status = 'unknown')  
      do i=1,x_len  
         write(1,*) x(i),hA(i), GA(i)
      end do  
      close(1)  
      
      !Do Recontruction
      !Reconstruct each cell to giver x^+_j-1/2, x_j, x^-_j+1/2
      
      !get hbeg and hend
      hbeg(:) = haR(1:nGP)
      hend(:) = haR(xhGR_len- (nGP-1):xhGR_len)
      
      !get ubeg and uend
      ubeg(:) = uaR(1:nGP)
      uend(:) = uaR(xuR_len- (nGP-1):xuR_len)
      
      !get ubeg and uend
      Gbeg(:) = GaR(1:nGP)
      Gend(:) = GaR(xhGR_len- (nGP-1):xhGR_len)
      
      !Reconstruct h,G
      hR = PerformReconstruction(hA,theta,hbeg,hend,
     . n_PtsCell, n_GhstCells)
     
      GR = PerformReconstruction(GA,theta,Gbeg,Gend,
     . n_PtsCell, n_GhstCells)
     
      !Solve for u
      uR = FEM(hR,GR,ubeg,uend,x_len,n_PtsCell,n_GhstCells)
      

      !Write Out At Reconstructed Points
      open(2, file = FilePre//'ReconHG.dat', status = 'unknown')  
      do i=1,xhGR_len
         write(2,*) xhGR(i),haR(i), hR(i),GaR(i), GR(i)
      end do  
      close(2)
      
      open(2, file = FilePre//'Reconu.dat', status = 'unknown')  
      do i=1,xuR_len
         write(2,*) xuR(i),uaR(i),uR(i)
      end do  
      close(2)     
      
      !Calculate error L2
      norm(1) = sqrt(sum((haR(:) - hR(:)) **2)) / 
     . sqrt(sum(haR(:)**2))      
      norm(2) = sqrt(sum((GaR(:) - GR(:)) **2)) / 
     . sqrt(sum(GaR(:)**2))   
      norm(3) = sqrt(sum((uaR(:) - uR(:)) **2)) / 
     . sqrt(sum(uaR(:)**2))  
       
      !Calculate error L1
c     norm = sum(abs(hR(:) - hRec(:))) / sum(abs(hR(:)))    
 
      contains


c==============================
c Function to generate x values at all points: All points in cell 
c for reconstruction + ghost cells
c==============================
      function GenerateAllLoc(x_start,dx,x_len,xR_len,n_PtsCell,
     . n_GhstCells) result (xR)
          
      implicit none
      real*8,intent(in) :: x_start,dx
      integer,intent(in) :: x_len,xR_len,n_PtsCell,n_GhstCells
      real*8 xR(xR_len)
      
      integer i
      
      !Left boundary
      do i=1,n_GhstCells 
         xR(n_PtsCell*(i) -2)  =  x_start
     .   - (n_GhstCells- i +1)*dx
     
         xR(n_PtsCell*(i) -1)  =  x_start
     .   - (n_GhstCells- i +1)*dx + 0.5*dx
     
         xR(n_PtsCell*(i))  = x_start 
     .   - (n_GhstCells- i + 1)*dx + dx
     
      end do 

      !Interior
      do i=1,x_len  
         xR(n_PtsCell*(i +n_GhstCells )- 2)  = x_start
     .      + (i-1)*dx
         xR(n_PtsCell*(i +n_GhstCells )- 1)  = x_start
     .      + (i-1)*dx + 0.5*dx
         xR(n_PtsCell*(i +n_GhstCells ))  = x_start
     .      + (i)*dx
      end do  
      
      !Right boundary
      do i=1,n_GhstCells 
         xR(n_PtsCell*(x_len +n_GhstCells + i) -2)  = x_start
     .      + (x_len + i -1)*dx
         xR(n_PtsCell*(x_len +n_GhstCells + i) -1)  = x_start
     .      + (x_len + i -1)*dx + 0.5*dx
         xR(n_PtsCell*(x_len +n_GhstCells + i))  = x_start
     .      + (x_len + i)*dx
      end do 
               
      end function
      
      
c ====================
c Function to generate h,u,G at all locations x
c ====================
      subroutine Initloc(xhG,h,G,xu,u,a0,a1,a2,a3,a4)
      implicit none
      
      real*8,dimension(:),intent(in) :: xhG,xu
      real*8,dimension(size(xu)), intent(out)   :: u
      real*8,dimension(size(xhG)), intent(out)   :: h,G
      real*8,dimension(size(xhG)) :: hx,ux,uxx
                   
      real*8 a0,a1,a2,a3,a4
      
            
      h(:)  = a0 + a1*DSIN(a2*xhG(:) ) 
      u(:)  = a3*DSIN(a4*xu(:) )
      
      hx(:) = a1*a2*DCOS(a2*xhG(:))
      ux(:) = a3*a4*DCOS(a4*xhG(:))

      uxx(:) = -a3*a4*a4*DSIN(a4*xhG(:))
      
      G(:) = a3*DSIN(a4*xhG(:) )*h(:) !- h(:)*h(:)*hx(:)*ux(:) -
c     . h(:)*h(:)*h(:)/3.0*uxx(:)
      
      end subroutine Initloc

c ====================
c Function to generate cell averages of h,G at all cell nodes x
c cells are dx wide
c ====================
      subroutine InitCellAvg(x,dx,a0,a1,a2,a3,a4,hA,GA)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,dimension(size(x)), intent(out)   :: hA,GA
       real*8,dimension(size(x)) :: uh,h3ux
                   
      real*8 a0,a1,a2,a3,a4,dx
            
            
      hA(:)  =  (a0*(x(:) + 0.5*dx) -a1/a2*DCOS(a2* (x(:) + 0.5*dx)) 
     . - a0*(x(:) - 0.5*dx) + a1/a2*DCOS(a2* (x(:) - 0.5*dx)) ) /dx
     
      ! integrating (h^3ux/3)x  gives h^3ux/3
      h3ux(:) = (a0 + a1*DSIN(a2*(x(:) + 0.5*dx) ))**3 *
     . a3*a4*DCOS(a4*(x(:) + 0.5*dx)) -
     . (a0 + a1*DSIN(a2*(x(:) - 0.5*dx) ))**3 *
     . a3*a4*DCOS(a4*(x(:) - 0.5*dx))
     
      uh(:) = -a0*a3/a4*DCOS(a4*(x(:) + 0.5*dx))
     . + a1*a3/(a4**2 - a2**2)*(a2*DSIN(a4*(x(:) + 0.5*dx))
     . *DCOS(a2*(x(:) + 0.5*dx)) - a4*DCOS(a4*(x(:) + 0.5*dx))
     . *DSIN(a2*(x(:) + 0.5*dx))) 
     . +a0*a3/a4*DCOS(a4*(x(:) - 0.5*dx))
     . - a1*a3/(a4**2 - a2**2)*(a2*DSIN(a4*(x(:) - 0.5*dx))
     . *DCOS(a2*(x(:) - 0.5*dx)) - a4*DCOS(a4*(x(:) - 0.5*dx))
     . *DSIN(a2*(x(:) - 0.5*dx)))
      
      
      GA(:) =  uh(:)/ dx!(uh(:) - h3ux(:)/3 )/dx
      
      end subroutine InitCellAvg      


c==============================
c Function to do Reconstruction - get h at x^+_j-1/2, x_j, x^-_j+1/2
c==============================
      function PerformReconstruction(qa,theta,qbeg,qend,
     . n_PtsCell, n_GhstCells) result (qR)
     
      implicit none
      
      integer, intent(in) :: n_PtsCell, n_GhstCells
      
      real*8,intent(in) :: theta
      real*8, dimension(:), intent(in)   :: qa
      real*8,intent(in) :: qbeg(n_PtsCell*n_GhstCells),
     . qend(n_PtsCell*n_GhstCells)
      
      real*8, dimension( n_PtsCell*(size(qa)+ 2*n_GhstCells) ) :: qR
      
      integer i,qr_len,qa_len
      real*8 fdq,mdq,bdq,cdq
      
      qa_len = size(qa)
      qr_len = n_PtsCell*(qa_len+ 2*n_GhstCells)
      
      !Left boundary
      qR(1: n_PtsCell*n_GhstCells) = qbeg(:)
      
      !Left most cell - needs hbeg to calculate derivative
      !cell 1 by i count, left cell point is n_PtsCell*n_GhstCells
      i = 1
      fdq = qa(2) - qa(1)
      mdq = 0.5*(qa(2) - qbeg(2))
      bdq = qa(1) - qbeg(2)
      
      cdq = minmod(theta*fdq,mdq,theta*bdq)
      
      qR(n_PtsCell*(i+n_GhstCells)-2) = qa(1) - 0.5*cdq
      qR(n_PtsCell*(i+n_GhstCells)-1) = qa(1)
      qR(n_PtsCell*(i+n_GhstCells)) = qa(1) + 0.5*cdq
      
      !Interior
      do i=2,qa_len-1
         fdq = qa(i+1) - qa(i)
         mdq = 0.5*(qa(i+1) - qa(i-1))
         bdq = qa(i) - qa(i-1)
         
         cdq = minmod(theta*fdq,mdq,theta*bdq)
         
         qR(n_PtsCell*(i+n_GhstCells)-2) = qa(i) - 0.5*cdq
         qR(n_PtsCell*(i+n_GhstCells)-1) = qa(i)
         qR(n_PtsCell*(i+n_GhstCells)) = qa(i) + 0.5*cdq
      end do
      
      !Right most cell - needs hend to calculate derivative
      !cell x_len by i count, left cell point is n_PtsCell*n_GhstCells
      i=qa_len
      fdq = qend(2) - qa(qa_len)
      mdq = 0.5*(qend(2) - qa(qa_len-1))
      bdq = qa(qa_len) - qa(qa_len-1)
      
      cdq = minmod(theta*fdq,mdq,theta*bdq)
      
      qR(n_PtsCell*(i+n_GhstCells)-2) = qa(i) - 0.5*cdq
      qR(n_PtsCell*(i+n_GhstCells)-1) = qa(i)
      qR(n_PtsCell*(i+n_GhstCells)) = qa(i) + 0.5*cdq
         
      !Right Boundary
      qR(qr_len- (n_PtsCell*n_GhstCells-1):qr_len) = qend(:)
      
      end function
      
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
      
c =========================
c FEM solve
c ========================
      function FEM(hR,GR,ubeg,uend,xc_len,n_PtsCell,n_GhstCells)
     .  result (uR)
      
      implicit none
   
      integer, intent(in) :: xc_len,n_PtsCell,n_GhstCells
      real*8, dimension(:), intent(in)   :: hR,GR
      real*8, dimension(:), intent(in)  :: ubeg, uend
      
      !Number of points for ur, is
      !The number of interior points between (x_-1/2 , x _m+1/2)
      ! which is 2*m -1 
      !Plus the number of boundary points, which includes the boundaries
      ! x_-1/2, x_m+1/2
      real*8, dimension(2*xc_len -1 + 2*size(ubeg))   :: uR,Gg
      
      !Pivot array, used by solver
      integer, dimension(2*xc_len -1 + 2*size(ubeg)) ::IPIV
      
      real*8, dimension(3,3) :: Ae
      real*8, dimension(3) :: Ge
      
      real*8, dimension(2*2 + 2 + 1,
     . 2*xc_len -1 + 2*size(ubeg)) :: Ag
     
      
      integer uR_len,i,nGP,INFO
      
      !values of G and h at cell edges
      !qjmhp -> q^{+}_{j - 1/2},qjphm -> q^{-}_{j + 1/2}
      real*8 Gjmhp,Gjphm,hjmhp,hjphm
      
      
      uR_len = 2*xc_len -1 + 2*size(ubeg)
      nGP = n_PtsCell*n_GhstCells
      !Make global arrays 0, so we can just add elementwise contributions
      AG(:,:) = 0d0
      Gg(:) = 0d0
      
      !Construct Global FEM matrices and arrays
      do i = 1,xc_len + 2*n_GhstCells
         Gjphm= GR(n_PtsCell*i )
         Gjmhp= GR(n_PtsCell*i -2 )
         
         hjphm= hR(n_PtsCell*i )
         hjmhp= hR(n_PtsCell*i -2 )
         
         !Generate element wise matrices Ae,Ge
         call EWFEM(Gjmhp,Gjphm,hjmhp,hjphm,dx,Ae,Ge)
      
         !Assemble the elementwise matrices into global matrices Ag, Gg
         call EContrG(Ae,Ge,i-1,Ag,Gg)

      end do
            
      !fix boundaries, to satisfy the boundary conditions
      !left boundary, first 3x3 matrix in Ag should be identity
      !first 3 elements of Gg should be ubeg
      Gg(:nGP) = ubeg(:)
      
      !0 elements, and then make diagonal 1
      Ag(:,:nGP) = 0d0
      Ag(5,:nGP) = 1d0
      
      !right boundary, last 3x3 matrix in Ag should be identity
      !last 3 elements of Gg should be uend
      Gg(uR_len-(nGP-1):uR_len) = uend(:)
      Ag(:,uR_len-(nGP-1):uR_len) = 0d0
      Ag(5,uR_len-(nGP-1):uR_len) = 1d0
    
c      !http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7
c      
      call dgbsv(uR_len, 2,2, 1,Ag,7,IPIV,Gg, uR_len,INFO)
      uR(:) = Gg(:)
      end function FEM

c ====================
c Function to generate elementwise matrices
c given Gjmhp,Gjphm,hjmhp,hjphm and dx
c ====================
      subroutine EWFEM(Gjmhp,Gjphm,hjmhp,hjphm,dx,Ae,Ge)
      implicit none
      
      real*8,intent(in) :: Gjmhp, Gjphm, hjmhp, hjphm, dx
      real*8,dimension(3,3), intent(out)   :: Ae
      real*8,dimension(3), intent(out)   :: Ge
      real*8, dimension(3,3) :: uhe, h3uxe
                   
      !Calculate G element contribution
      Ge(1) = Gjmhp
      Ge(2) = 2*(Gjmhp + Gjphm)
      Ge(3) = Gjphm
      Ge(:) = dx/6.0*Ge(:)
      
      !calculate uh element contribution
      uhe(1,1) = (7*hjmhp + hjphm)
      uhe(1,2) = (4*hjmhp )
      uhe(1,3) = (-hjmhp - hjphm)
      
      uhe(2,1) = (4*hjmhp)
      uhe(2,2) = (16*hjmhp + 16*hjphm)
      uhe(2,3) = (4*hjphm)
       
      uhe(3,1) = (-hjmhp - hjphm)
      uhe(3,2) = (4*hjphm)
      uhe(3,3) = (hjmhp + 7*hjphm)
      
      
      uhe(:,:) = dx/60.0*uhe(:,:)
      
      !calculate uh element contribution
c      h3uxe(1,1) = 79.0/120*hjmhp*hjmhp*hjmhp +  
c     . 39.0/120*hjmhp*hjmhp*hjphm + 3.0/24*hjmhp*hjphm*hjphm + 
c     . 7.0/120*hjphm*hjphm*hjphm
c      h3uxe(1,2) = -23.0/30*hjmhp*hjmhp*hjmhp + 
c     . -3.0/10*hjmhp*hjmhp*hjphm + -3.0/30*hjmhp*hjphm*hjphm + 
c     . -1.0/6*hjphm*hjphm*hjphm)
c      h3uxe(1,3) = 13.0/120*hjmhp*hjmhp*hjmhp +
c     . -3.0/120*hjmhp*hjmhp*hjphm + -3.0/120*hjmhp*hjphm*hjphm + 
c     . 13.0/120*hjphm*hjphm*hjphm
c      
c      h3uxe(2,1) = -23.0/30*hjmhp*hjmhp*hjmhp + 
c     . -3.0/10*hjmhp*hjmhp*hjphm + -3.0/30*hjmhp*hjphm*hjphm + 
c     . -1.0/6*hjphm*hjphm*hjphm    
c      h3uxe(2,2) = 14.0/15*hjmhp*hjmhp*hjmhp +  
c     . 6.0/15*hjmhp*hjmhp*hjphm + 6.0/15*hjmhp*hjphm*hjphm + 
c     . 14.0/15*hjphm*hjphm*hjphm  
c      h3uxe(2,3) = -1.0/6*hjmhp*hjmhp*hjmhp + 
c     . -3.0/30*hjmhp*hjmhp*hjphm + -3.0/10*hjmhp*hjphm*hjphm +
c     . -23.0/30*hjphm*hjphm*hjphm
c       
c      h3uxe(3,1) = 13.0/120*hjmhp*hjmhp*hjmhp +
c     .  -3.0/120*hjmhp*hjmhp*hjphm + -3.0/120*hjmhp*hjphm*hjphm +  
c     .  13.0/120*hjphm*hjphm*hjphm      
c      h3uxe(3,2) = -1.0/6*hjmhp*hjmhp*hjmhp +  
c     . -3.0/30*hjmhp*hjmhp*hjphm + -3.0/10*hjmhp*hjphm*hjphm + 
c     . -23.0/30*hjphm*hjphm*hjphm)      
c      h3uxe(3,3) = 7.0/120*hjmhp*hjmhp*hjmhp +  
c     . 3.0/24*hjmhp*hjmhp*hjphm + 39.0/120*hjmhp*hjphm*hjphm +
c     . 79.0/120*hjphm*hjphm*hjphm)

c      h3uxe(:,:) = 2.0/(3*dx)*h3uxe(:,:)
      
      Ae(:,:) = uhe(:,:) !+ h3uxe(:,:)
      
      end subroutine EWFEM       
      
c ====
c Add elementwise contribution to global matrix
c ===
      subroutine EContrG(Ae,Ge,i,Ag,Gg)
      implicit none
      !for interior - i ~ (i + n_GhstCells - 1)
      
      integer,intent(in) :: i
      
      real*8,dimension(:,:),intent(in) :: Ae
      real*8,dimension(:,:), intent(out)   :: Ag

      real*8,dimension(:),intent(in) :: Ge
      real*8,dimension(:), intent(out)   :: Gg          
      !add element contribution to global matrix
      ! Ag is in banded form, it is LDAB * N -https://www.netlib.org/lapack/lug/node124.html
      !first super diagonal
      Ag(3, 2*i + 3 ) = 
     . Ag(3, 2*i + 3 ) + Ae(1,3)
     
      !second super diagonal
      Ag(4, 2*i + 2 ) = 
     . Ag(4, 2*i +2 ) + Ae(1,2)
      Ag(4, 2*i + 3 ) = 
     . Ag(4, 2*i +3 ) + Ae(2,3)
     
      !diagonal elements
      Ag(5, 2*i + 1 ) = Ag(5, 2*i +1 )
     . + Ae(1,1)
      Ag(5, 2*i + 2 ) = 
     . Ag(5, 2*i + 2 ) + Ae(2,2)
      Ag(5, 2*i + 3 ) = 
     . Ag(5, 2*i + 3 ) + Ae(3,3)
     
      !first subdiagonal
      Ag(6, 2*i +1 ) = Ag(6, 2*i +1 )
     . + Ae(2,1)
      Ag(6, 2*i + 2 ) = 
     . Ag(6, 2*i + 2 ) + Ae(3,2)
     
      !second subdiagonal
      Ag(7, 2*i  + 1) = Ag(7, 2*i +1 )
     . + Ae(3,1)
     
      !G array
      Gg(2*i + 1 : 2*i + 3) = Gg(2*i + 1 : 2*i + 3) + Ge
      end subroutine EContrG       
      
            
      end subroutine ReconAll
      
c =========================
c Main program to actually run the LinearRecon stuff for various dx/n values
c ========================
      program main
            
      implicit none
      
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/FortranTests/ReconFEM/"
  
      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n
      real*8 norm(3)    
      real*8 PI,a0,a1,a2,a3,a4,dx
      
      n= 16
      
      a0 = 1
      a1 = 0.5
      a3 = 0.1
      
      PI = 2*DASIN(1.0)
      a2 = 2*PI/ 20.0
      a4 = 2*PI/ 10.0
      
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      open(3, file = wdir//'Norms.dat', status = 'unknown') 
      
      do i=1,n
         print *,'Reconstruction',i,100*(2**(i-1))
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         call ReconAll(0.0d0,100.0d0,100*(2**(i-1)),3,1,1.2d0,a0,
     .   a1,a2,a3,a4,trim(fileloci),norm,dx)


         !write out
         write(3,*) dx,norm(1),norm(2),norm(3)
         
      end do
      
      close(3)  
      

      end program main
      
