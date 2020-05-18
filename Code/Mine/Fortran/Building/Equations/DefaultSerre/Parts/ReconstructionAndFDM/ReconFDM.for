c=====================================
c Program that given array of h calculates h* the reconstructed array
c Analysis shows second order until 10^-10 by which point we get round-off error effects
c=====================================
      subroutine ReconAll(x_start,x_end,x_len, n_GhstCells,
     . theta,a0,a1,a2,a3,a4,FilePre,norm,dx)
      implicit none
      
      real*8,intent(in) :: x_start,x_end,theta,a0,a1,a2,a3,a4
      integer, intent(in) :: x_len, n_GhstCells
      character (LEN=*), intent(in)   :: FilePre
      real*8, intent(out)   :: norm(4),dx
                         
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
    
      integer xr_len,xbc_len,i
      xr_len = 3*x_len
      xbc_len = x_len + 2*n_GhstCells

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
            
      !Get locations of reconstruction plus Ghost Cells
      xbc = Generatexbc(x_start,dx,x_len,n_GhstCells)
      
      !Get h,u and G at all these cell nodes
      call Initloc(xbc,a0,a1,a2,a3,a4,hbc,Gbc,uabc)
      
      
      !get initial u (mainly to get BC)
      ubc(:) = uabc(:)
      
      !Get ubc, from h and G from FD method
      call FDforu(hbc,Gbc,ubc,dx,n_GhstCells)
      
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

      norm(1) = sqrt(sum((hRA(:) - hR(:)) **2)) / 
     . sqrt(sum(hRA(:)**2))       
      norm(2) = sqrt(sum((GRA(:) - GR(:)) **2)) / 
     . sqrt(sum(GRA(:)**2)) 
      norm(3) = sqrt(sum((uRA(:) - uR(:)) **2)) / 
     . sqrt(sum(uRA(:)**2)) 
      norm(4) = sqrt(sum((duRA(:) - duR(:)) **2)) / 
     . sqrt(sum(duRA(:)**2)) 
     
      contains


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
      
      subroutine FDforu(hbc,Gbc,ubc,dx,n_GhstCells)
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
      
      end subroutine FDforu

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
      
           
      
      end subroutine ReconAll
      

c =========================
c Main program to actually run the LinearRecon stuff for various dx/n values
c ========================
      program main
            
      implicit none
      
      CHARACTER(*), PARAMETER :: wdir = "/home/jp/Documents/"//
     .   "Work/PostDoc/Projects/Steve/1DWaves/"//
     .   "RegularisedSerre/CodeAndData/Data/RAW"//
     .   "/FortranTests/ReconFDM/"
  
      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 norm(4)    
      real*8 PI,a0,a1,a2,a3,a4,dx,theta,startx,endx
      
      n= 116
      n_gc = 3
      
      lowxlen = 100
      
      startx = 0.0
      endx = 100.0
      
      theta = 1.2
      
      a0 = 1d0
      a1 = 0.5d0
      a3 = 0.1d0
      
      PI = 2*DASIN(1.0d0)
      a2 = 2*PI/ 20.0
      a4 = 2*PI/ 10.0
      
      !Remove previous runs
      CALL SYSTEM('rm -rf '//wdir)
      CALL SYSTEM('mkdir -p '// wdir)
      
      open(3, file = wdir//'Norms.dat', status = 'unknown') 
      
      do i=1,n
         print *,'Reconstruction',i,lowxlen*(2**(i-1))
         write (stri,'(I2.2)') i
         fileloci = wdir // stri // '/'
         CALL SYSTEM('mkdir -p '// fileloci)
         call ReconAll(startx,endx,lowxlen*(2**(i-1)),n_gc,theta,a0,
     .   a1,a2,a3,a4,trim(fileloci),norm,dx)


         !write out
         write(3,*) dx,norm(1), norm(2), norm(3), norm(4)
         
      end do
      
      close(3)  
      

      end program main
      
