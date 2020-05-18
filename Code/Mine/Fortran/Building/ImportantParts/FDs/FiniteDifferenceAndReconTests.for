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
      real*8, intent(out)   :: norm(3),dx
                         
      ! locations for x,u,h,G (cell nodes in interior and bc)
      real*8 xbc(x_len + 2*n_GhstCells),
     . qbc(x_len + 2*n_GhstCells),
     . dqbc(x_len + 2*n_GhstCells),
     . ddqbc(x_len + 2*n_GhstCells)
     
      !reconstructions
      real*8 xR(3*x_len),
     . qR(3*x_len),
     . dqR(3*x_len),
     . ddqR(3*x_len),
     . qRA(3*x_len),
     . dqRA(3*x_len),
     . ddqRA(3*x_len)
    
      integer xr_len,xbc_len,i
      xr_len = 2*x_len
      xbc_len = x_len + 2*n_GhstCells

      !Initial data
      dx = (x_end - x_start) / (x_len -1) 
            
      !Get locations of reconstruction plus Ghost Cells
      xbc = Generatexbc(x_start,dx,x_len,n_GhstCells)
      
      !Get h,u and G at all these cell nodes
      call Initloc(xbc,a0,a1,qbc,dqbc,ddqbc)
      
      !DoReconstruction at cell edges
      call ReconsAll(xbc,qbc,n_GhstCells,dx,theta,
     . xR,qR,dqR,ddqR)
     
      !Get analytic values for recons
      call Initloc(xR,a0,a1,qRA,dqRA,ddqRA)
      
      open(1, file = FilePre//'All.dat', status = 'unknown')  
      do i=1,xr_len
         write(1,*) xR(i),qRA(i),qR(i),dqRA(i),dqR(i),ddqRA(i),ddqR(i)  
      end do  
      close(1)

      norm(1) = sqrt(sum((qRA(:) - qR(:)) **2)) / 
     . sqrt(sum(qRA(:)**2))       
      norm(2) = sqrt(sum((dqRA(:) - dqR(:)) **2)) / 
     . sqrt(sum(dqRA(:)**2)) 
      norm(3) = sqrt(sum((ddqRA(:) - ddqR(:)) **2)) / 
     . sqrt(sum(ddqRA(:)**2)) 
     
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
      subroutine Initloc(x,a0,a1,q,dq,ddq)
      implicit none
      
      real*8,dimension(:),intent(in) :: x
      real*8,intent(in) :: a0,a1
      real*8,dimension(size(x)), intent(out)   :: q,dq,ddq
      
      integer i,n
      
      n= size(x)
      
      q(:)  = a0 + a1*DSIN(a2*x(:) ) 
      dq(:) =  a2*a1*DCOS(a2*x(:) ) 
      ddq(:) =  -(a2**2)*a1*DSIN(a2*x(:) ) 
            
      end subroutine Initloc
      
      
      subroutine ReconsAll(xbc,qbc,n_GhstCells,dx,theta
     . ,xR,qR,dqR,ddqR)
      integer, intent(in) :: n_GhstCells
      real*8, intent(in) :: dx,theta
      real*8,dimension(:),intent(in) :: xbc,qbc
      
      real*8, dimension(3*(size(xbc) - 2*n_GhstCells)),
     .   intent(out) :: xR,qR,dqR,ddqR
     
     
      real*8 cdq
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
         
         !q,dq,ddq recon
         cdq = ReconLinLimGrad(qbc(i-1),qbc(i),qbc(i+1),theta)
         qR(3*j + 1) = qbc(i) - 0.5*cdq
         qR(3*j + 2) = qbc(i)
         qR(3*j + 3) = qbc(i) + 0.5*cdq
         
         dqR(3*j + 1) = (-2*qbc(i) + 3*qbc(i+1) - qbc(i+2)) /dx
         dqR(3*j + 2) = (qbc(i + 1) - qbc(i-1))/(2*dx)
         dqR(3*j + 3) = (2*qbc(i) - 3*qbc(i-1) + qbc(i-2)) /dx

         ddqR(3*j + 1) = (5*qbc(i) - 13*qbc(i+1) + 
     .   11*qbc(i+2) - 3*qbc(i+3) )/(2*dx**2)
     
         ddqR(3*j + 2) = (qbc(i + 1) - 2*qbc(i) + qbc(i-1))/(dx*dx)
         
         ddqR(3*j + 3) = (5*qbc(i) - 13*qbc(i-1) + 
     .   11*qbc(i-2) - 3*qbc(i-3) )/(2*dx**2)       
         
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
     .   "/MethodParts/Reconstructions/"
  
      CHARACTER(200) :: fileloci
      CHARACTER(2) :: stri
      integer i,n,n_gc,lowxlen
      real*8 norm(4)    
      real*8 PI,a0,a1,a2,a3,a4,dx,theta,startx,endx
      
      n= 16
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
      
