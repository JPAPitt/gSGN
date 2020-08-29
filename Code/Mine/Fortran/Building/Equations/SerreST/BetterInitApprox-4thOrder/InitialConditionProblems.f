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
c Subroutine that generates the Peakon solution for Bo = 1/3 (g =1),
c Left and right approximations
c =====      
      subroutine Peakon(x,x_len,t,dx,ga,a0,a1,h,u,G)
      implicit none
      
      integer x_len
      real*8 x(x_len)
      real*8 t,a0,a1,ga
      real*8 h(x_len),G(x_len),u(x_len)
      
      real*8 c,phi,ni,nip1,nim1,nx,nxx,ux,uxx
      real*8 nip2,nim2,nip3,nim3
      real*8 phim1, phip1,phim2, phip2,phim3, phip3,dx
      integer i

      c = dsqrt(1 + a1)
      
      h(1) = 1
      u(1) = 0
      G(1) = 0
      h(2) = 1
      u(2) = 0
      G(2) = 0
      do i=3,x_len-2
         phi = x(i) - c*t
         phim2 = x(i-2) - c*t  
         phim1 = x(i-1) - c*t        
         phip2 = x(i+2) - c*t  
         phip1 = x(i+1) - c*t    
             
         ni= a1*dexp(-dsqrt(3d0)*dabs(phi))
         nip1= a1*dexp(-dsqrt(3d0)*dabs(phip1))
         nip2= a1*dexp(-dsqrt(3d0)*dabs(phip2))
         nim1= a1*dexp(-dsqrt(3d0)*dabs(phim1))
         nim2= a1*dexp(-dsqrt(3d0)*dabs(phim2))
         nx = (1*nim2 - 8*nim1 + 8*nip1 - nip2)/ (12*dx)
         nxx = (-nim2 + 16*nim1 - 30*ni + 16*nip1 - nip2)/ (12*dx*dx)
         
         ux = c*nx/ (1 + ni)**2
         uxx = c*((1 + ni)*nxx - 2*nx**2)/(ni + 1)**3
         
         h(i) = 1 + ni  
         u(i) = c*(ni / (1+ni))
         G(i) = u(i)*h(i) - (h(i)**3)*uxx/3d0 - h(i)**2*nx*ux
   
   
      end do 
      
      h(x_len) = 1
      u(x_len) = 0
      G(x_len) = 0
      
      h(x_len-1) = 1
      u(x_len-1) = 0
      G(x_len-1) = 0
      
      end

c ====
c Subroutine that generates the SWWE Dambreak Solution
c =====      
      subroutine Dambreak(x,x_len,t,ga,hl,hr,h,u,G)
      implicit none
      
      integer x_len
      real*8 x(x_len),h(x_len),G(x_len),u(x_len)
      real*8 t,hl,hr,ga
      
      real*8 h2,u2,S2,zmin,fzmin,z,fz,fzmax,zmax
      integer i
      
      if (dabs(t) .LE. 10d0**(-10)) then
         do i=1,x_len 

            if (x(i) .LE. 0) then
               h(i) = hl
               u(i) = 0
            else
               h(i) = hr
               u(i) = 0
            end if
            
            G(i) = u(i)*h(i)


         end do 
      
      
      else
      
         !calculate h2
         !h2 = hr/2*dsqrt( dsqrt(1 + 8*((2*h2/(h2 - hr))* ((dsqrt(ga*hl) - dsqrt(ga*h2))/ sqrt(ga*hr)) )**2 ) -1)
         
         !bisection method to calculate h2
         zmin = hl
         fzmin = zmin - hr/2*( dsqrt(1 + 8*((2*zmin/(zmin - hr))*
     .      ((dsqrt(ga*hl) - dsqrt(ga*zmin))/ sqrt(ga*hr)) )**2 ) -1)
         zmax = hr
         fzmax = zmax - hr/2*( dsqrt(1 + 8*((2*zmax/(zmax - hr))*
     .      ((dsqrt(ga*hl) - dsqrt(ga*zmax))/ sqrt(ga*hr)) )**2 ) -1)
         do while (abs(zmax - zmin) .gt. 10d0**(-14))
           z = (zmin + zmax)/2.d0
           fz = z - hr/2*( dsqrt(1 + 8*((2*z/(z - hr))*
     .      ((dsqrt(ga*hl) - dsqrt(ga*z))/ sqrt(ga*hr)) )**2 ) -1)
        
           if(fz .gt. 0.d0)then
             zmin = z
             fzmin = fz
           else
             zmax = z
             fzmax = fz
           end if
           
         end do
         
         h2 = z
         
         S2 = 2d0*h2 / (h2 - hr) *(sqrt(ga*hl) -  dsqrt(ga*h2))
         u2 = 2*(dsqrt(ga*hl) - dsqrt(ga*h2))
               
         do i=1,x_len 
         
            if (x(i) .LE. -t*dsqrt(ga*hl)) then
               h(i) = hl
               u(i) = 0
            else if ((x(i) .GT. -t*dsqrt(ga*hl)) 
     .         .AND. (x(i) .LE. t*(u2 - dsqrt(ga*h2))))  then
               h(i) = 4d0/(9d0*ga) *(dsqrt(ga*hl) - x(i)/(2d0*t))**2
               u(i) = 2d0/3d0*(dsqrt(ga*hl) + x(i)/t)
            else if ((x(i) .GT. t*(u2 - dsqrt(ga*h2)) ) 
     .         .AND. (x(i) .LT. t*S2))  then
               h(i) = h2    
               u(i) = u2  
            else
               h(i) = hr
               u(i) = 0
            end if
            
            G(i) = u(i)*h(i)
         
      
         end do 
      
      end if
     
      
      
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


