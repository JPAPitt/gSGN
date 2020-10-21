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
      DOUBLE PRECISION x_start,dx
      integer xbc_len,n_GhstCells
      DOUBLE PRECISION xbc(xbc_len)
      
      integer i
      
      !Left boundary
      do i=1,xbc_len  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end 

      subroutine ReconTests(x,xbc_len,dx,hA,hjph,a0,a1,a2)
      integer xbc_len
      DOUBLE PRECISION x(xbc_len),hA(xbc_len),hjph(xbc_len),dx,
     . a0, a1,a2,xjmh,xjph
     
      do i=1,xbc_len 
         xjmh = x(i) - dx/2
         xjph = x(i) + dx/2
         hA(i) = (a0*(xjph - xjmh) - 
     .     a1/a2*(cos(a2*xjph) - cos(a2*xjmh))) / dx
         hjph(i) = a0 + a1*sin(a2*xjph)
      
      end do
     
     
      end


c ====
c Subroutine that generates the SWWE Dambreak Solution
c =====      
      subroutine Dambreak(x,xbc_len,t,grav,hl,hr,
     . hA,uA,GA)
      implicit none
      
      
      integer xbc_len
      DOUBLE PRECISION x(xbc_len),hA(xbc_len),uA(xbc_len),
     . GA(xbc_len)
      DOUBLE PRECISION t,hl,hr,grav,cxi
      
      DOUBLE PRECISION h2,u2,S2,zmin,fzmin,z,fz,fzmax,zmax
      integer i
      
      if (dabs(t) .LE. 10d0**(-10)) then
         do i=1,xbc_len 

            ! x_j-1/2
            !print *, i, x(i)
            if (x(i)  .LT. 0) then
               hA(i) = hl
            else
               hA(i) = hr
            end if
            
            uA(i) = 0
            GA(i) = 0            
         end do 
      
      
      else
      
         !calculate h2
         !h2 = hr/2*dsqrt( dsqrt(1 + 8*((2*h2/(h2 - hr))* ((dsqrt(ga*hl) - dsqrt(ga*h2))/ sqrt(ga*hr)) )**2 ) -1)
         
         !bisection method to calculate h2
         zmin = hl
         fzmin = zmin - hr/2*( dsqrt(1 + 8*((2*zmin/(zmin - hr))*
     .      ((dsqrt(grav*hl) - dsqrt(grav*zmin))/ sqrt(grav*hr)) )**2 )
     . -1)
         zmax = hr
         fzmax = zmax - hr/2*( dsqrt(1 + 8*((2*zmax/(zmax - hr))*
     .      ((dsqrt(grav*hl) - dsqrt(grav*zmax))/ sqrt(grav*hr)) )**2 ) 
     . -1)
         do while (abs(zmax - zmin) .gt. 10d0**(-14))
           z = (zmin + zmax)/2.d0
           fz = z - hr/2*( dsqrt(1 + 8*((2*z/(z - hr))*
     .      ((dsqrt(grav*hl) - dsqrt(grav*z))/ sqrt(grav*hr)) )**2 ) 
     . -1)
        
           if(fz .gt. 0.d0)then
             zmin = z
             fzmin = fz
           else
             zmax = z
             fzmax = fz
           end if
           
         end do
         
         h2 = z
         
         S2 = 2d0*h2 / (h2 - hr) *(sqrt(grav*hl) -  dsqrt(grav*h2))
         u2 = 2*(dsqrt(grav*hl) - dsqrt(grav*h2))
               
         do i=1,xbc_len
         
            ! x_i-1/2
            cxi = x(i)
            if (cxi .LE. -t*dsqrt(grav*hl)) then
               hA(i) = hl
               uA(i) = 0
            else if ((cxi .GT. -t*dsqrt(grav*hl)) 
     .         .AND. (cxi .LE. t*(u2 - dsqrt(grav*h2))))  then
               hA(i) = 4d0/(9d0*grav)*(dsqrt(grav*hl) - cxi/(2d0*t))**2
               uA(i) = 2d0/3d0*(dsqrt(grav*hl) + cxi/t)
            else if ((cxi .GT. t*(u2 - dsqrt(grav*h2)) ) 
     .         .AND. (cxi .LT. t*S2))  then
               hA(i) = h2    
               uA(i) = u2  
            else
               hA(i) = hr
               uA(i) = 0
            end if
            
            GA(i) = uA(i)*hA(i)

         end do 
      
      end if
     
      
      
      end
      


