c ====================================================================================
c Module of fortran subroutines to generate the initial condition problems
c to be handed to the FDVM2 module to produce the desired numerical solutions
c ====================================================================================        

c ====
c Subroutine that generates the cell centres, which are our locations for the initial conditions h,G and u
c we generate both the interior and ghost cells. When dx = (xend - xstart) / (x_len -1) , the first cell centre
c in the interior will be xstart and the last cell centre will be xend.
c ====      
      subroutine Generatexbc(x_start,dx,xbc_len,n_GhstCells,xbc)
          
      implicit none
      real*8 x_start,dx
      integer xbc_len,n_GhstCells
      real*8 xbc(xbc_len)
      
      integer i
      
      !Left boundary
      do i=1,xbc_len  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end 

c ====
c Subroutine that generates the Serre soliton solution - this solution only valid when beta1 = beta2 = 0
c =====      
      subroutine SerreSoliton(x,x_len,t,ga,a0,a1,h,u,G)
      implicit none
      
      integer x_len
      real*8 x(x_len)
      real*8 t,a0,a1,ga
      real*8 h(x_len),G(x_len),u(x_len)
      
      real*8 k,c,phi,sechkphi
      integer i
      
      k = dsqrt(3*a1) / (2*a0*dsqrt(a0 + a1))
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
      
      end

      
c ====
c Subroutine that generates the hyperbolic tanh
c =====      
      subroutine SmoothDB(x,x_len,hl,hr,alpha,h,u,G)
      implicit none
      
      integer x_len
      real*8 x(x_len),h(x_len),G(x_len),u(x_len)
      real*8 hl,hr,alpha
      
      integer i
      
      do i=1,x_len 
         h(i) = hl + (hr - hl)/2d0*(1d0 + dtanh(x(i)/alpha))
         u(i) = 0d0
         G(i) = 0d0
      end do 
      
      end
      
c ====
c Subroutine that generates gaussian (want really smooth initial conditions to see what might be causing failure)
c =====      
      subroutine Gaussian(x,x_len,a0,a1,a2,h,u,G)
      implicit none
      
      integer x_len
      real*8 x(x_len),h(x_len),G(x_len),u(x_len)
      real*8 a0,a1,a2
      
      integer i
      
      do i=1,x_len 
         h(i) = a0 + a1*dexp( - x(i)**2 / (2*a2**2) )
         u(i) = 0d0
         G(i) = 0d0
      end do 
      
      end


