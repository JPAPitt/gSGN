c ====================================================================================
c Module of fortran subroutines to generate the initial condition problems
c to be handed to the FDVM2 module to produce the desired numerical solutions
c ====================================================================================        
      module InitialConditionProblems
      implicit none

      contains

c ====
c Subroutine that generates the cell centres, which are our locations for the initial conditions h,G and u
c we generate both the interior and ghost cells. When dx = (xend - xstart) / (x_len -1) , the first cell centre
c in the interior will be xstart and the last cell centre will be xend.
c ====      
      subroutine Generatexbc(x_start,dx,xbc_len,n_GhstCells,xbc)
          
      implicit none
      real*8,intent(in) :: x_start,dx
      integer,intent(in) :: xbc_len,n_GhstCells
      real*8,intent(out) :: xbc(xbc_len)
      
      integer i
      
      !Left boundary
      do i=1,xbc_len  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end subroutine Generatexbc

c ====
c Subroutine that generates the Serre soliton solution - this solution only valid when beta1 = beta2 = 0
c =====      
      subroutine SerreSoliton(x,t,ga,a0,a1,h,u,G)
      
      real*8,dimension(:),intent(in) :: x
      real*8,intent(in) :: t,a0,a1,ga
      real*8,dimension(size(x)), intent(out)   :: h,G,u
      
      real*8 k,c,phi,sechkphi
      integer i,x_len
      
      x_len = size(x)
      k = dsqrt(3*a1) / (2*a0*sqrt(a0 + a1))
      c = dsqrt(ga*(a0 + a1))
      
      do i=1,x_len 
         
         phi  = x(i) - c*t
         sechkphi = 1.0 / dcosh(k*phi)
         h(i)  = a0 + a1*sechkphi**2
         u(i)  = c*(1 - a0 / h(i))
                        
         G(i) = u(i)*h(i) + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*h(i)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*dtanh(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*h(i)*dtanh(k*phi)**2 
   
      end do 
      
      end subroutine SerreSoliton
      
      end module InitialConditionProblems

