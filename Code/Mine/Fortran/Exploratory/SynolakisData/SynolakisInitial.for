      program main
      
      integer n
      real*8 x_1, x_xlen,dx
      real*8 x(1000),h(1000),w(1000),u(1000), b(1000)
      real*8 t,h0,amp,beta,xpeak,g
      
      n = size(x)
      
      !location vars
      x_1 = -50
      x_xlen = 150
      dx = (x_xlen - x_1) / (n -1)
      
      !experiment vars
      t = 0
      h0 = 1
      amp = 0.0185
      beta = datan(1.0/19.85)
      xpeak = 38.5
      g =1
      
      
      
      x = GenerateCellNodes(x_1,dx,n)
      
      call Synolakis(x,t,h0,amp,beta,xpeak,g,h,w,u,b) 
      
      open(1, file = 'CellNodesAll.dat', status = 'unknown')  
      do i=1,n  
         write(1,*) x(i),h(i),w(i),u(i),b(i)
      end do  
      close(1)  
      
      contains

c ============
c Function to generate cell nodes, x1,x2,...,x_(x_len)
c ============      
      function GenerateCellNodes(x_1,dx,x_len) result (x)
          
      implicit none
      real*8,intent(in) :: x_1,dx
      integer,intent(in) :: x_len
      real*8 x(x_len)
      
      integer i
      
      !Interior
      do i=1,x_len  
         x(i)  = x_1   + (i-1)*dx
      end do  
      
      end function GenerateCellNodes


c ============
c Subroutine to generate h,u,b at cell nodes
c h - height, w - stage, u - velocity, b - bed
c h0 - height of water before slope
c amp - amplitude of synolakis wave
c beta - slope of beach
c xpeak is initial position of synolakis wave at t = 0
c g - acceleartion due to gravity
c ============      
      subroutine Synolakis(x,t,h0,amp,beta,xpeak,g,h,w,u,b) 
          
      implicit none
      real*8,intent(in) :: t,h0,amp,beta,xpeak,g
      real*8,dimension(:),intent(in) :: x
      real*8,dimension(size(x)),intent(out) :: h,w,u,b
      
      integer i,n
      real*8 c,k
      
      n = size(x)
      
      ! c- speed of synolakis wave (assuming soliton)
      ! k - parameter for sech^2, in soliton
      c = sqrt(g*(h0 + amp))
      k = sqrt(3.0*amp) / (2.0*h0 *sqrt(h0 + amp))
      !Interior
      do i=1,n  
         ! 1.0/dtan(beta) = cot(beta)
         if (x(i) <= 1.0/dtan(beta)) then
            b(i) = - x(i)*dtan(beta)
            
        else
            b(i) = -1
        
        end if
        
        if (b(i) >= 0) then
            h(i) = 0
            u(i) = 0
            w(i) = b(i)
        else
            w(i) = amp*sechsq(k*( (x(i) - xpeak ) - c*t))
            h(i) = w(i) - b(i)
            !wave travels to left, so -c
            u(i)= -c* (1 - h0/ (w(i) + h0))
        end if
           
      end do  
      
      end subroutine Synolakis
      
      !return sech^2(a)
      function sechsq(a) result (b)
      implicit none
      real*8, intent(in) :: a
      real*8 b
      
      b = (2.0/(dexp(a) + dexp(-a)))**2
      
      end function sechsq
      
       
      end program
   
   
