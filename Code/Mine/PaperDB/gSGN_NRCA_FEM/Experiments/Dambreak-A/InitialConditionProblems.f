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

c ====
c Subroutine that generates the Serre soliton solution - this solution only valid when beta1 = beta2 = 0
c =====      
      subroutine SerreSoliton(x_bc,xbc_len,t,ga,a0,a1,dx,hcub_bc,
     .   ucub_bc,Gcub_bc,cubbc_len)
     
      implicit none
      
      
      integer xbc_len,cubbc_len
      DOUBLE PRECISION x_bc(xbc_len)
      DOUBLE PRECISION t,a0,a1,ga,dx
      DOUBLE PRECISION hcub_bc(cubbc_len),Gcub_bc(cubbc_len),
     . ucub_bc(cubbc_len)
      
      DOUBLE PRECISION k,c,phi,sechkphi
      integer i
      
      k = dsqrt(3*a1) / (2*a0*dsqrt(a0 + a1))
      c = dsqrt(ga*(a0 + a1))
      
      do i=1,xbc_len 
         
         ! x_j-1/2
         phi  = x_bc(i) - 0.5*dx - c*t
         sechkphi = 1.0 / dcosh(k*phi)
         hcub_bc(4*i - 3)  = a0 + a1*sechkphi**2
         ucub_bc(4*i - 3)  = c*(1 - a0 / hcub_bc(4*i - 3))
                        
         Gcub_bc(4*i - 3) = ucub_bc(4*i - 3)*hcub_bc(4*i - 3)
     . + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*hcub_bc(4*i - 3)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*dtanh(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*hcub_bc(4*i - 3)
     . * dtanh(k*phi)**2 


         ! x_j-1/6
         phi  = x_bc(i) - dx/6.0 - c*t
         sechkphi = 1.0 / dcosh(k*phi)
         hcub_bc(4*i - 2)  = a0 + a1*sechkphi**2
         ucub_bc(4*i - 2)  = c*(1 - a0 / hcub_bc(4*i - 2))
                        
         Gcub_bc(4*i - 2) = ucub_bc(4*i - 2)*hcub_bc(4*i - 2)
     . + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*hcub_bc(4*i - 2)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*dtanh(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*hcub_bc(4*i - 2)
     . * dtanh(k*phi)**2 

         ! x_j+1/6
         phi  = x_bc(i) + dx/6.0 - c*t
         sechkphi = 1.0 / dcosh(k*phi)
         hcub_bc(4*i - 1)  = a0 + a1*sechkphi**2
         ucub_bc(4*i - 1)  = c*(1 - a0 / hcub_bc(4*i - 1))
                        
         Gcub_bc(4*i - 1) = ucub_bc(4*i - 1)*hcub_bc(4*i - 1)
     . + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*hcub_bc(4*i - 1)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*dtanh(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*hcub_bc(4*i - 1)
     . *dtanh(k*phi)**2 

         ! x_j+1/2
         phi  = x_bc(i) + 0.5*dx - c*t
         sechkphi = 1.0 / dcosh(k*phi)
         hcub_bc(4*i)  = a0 + a1*sechkphi**2
         ucub_bc(4*i)  = c*(1 - a0 / hcub_bc(4*i))
                        
         Gcub_bc(4*i ) = ucub_bc(4*i)*hcub_bc(4*i)
     . + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*hcub_bc(4*i)
     . - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*dtanh(k*phi)**2 
     . - 4.0/3*a0*a1*c*k**2*sechkphi**2*hcub_bc(4*i)
     . * dtanh(k*phi)**2 
        
      end do 
      
      end

c ====
c Subroutine that generates the SWWE Dambreak Solution
c =====      
      subroutine Dambreak(x,xbc_len,t,ga,hl,hr,
     . dx,h,u,G,cubbc_len)
      implicit none
      
      
      integer xbc_len,cubbc_len
      DOUBLE PRECISION x(xbc_len),h(cubbc_len),u(cubbc_len),
     . G(cubbc_len)
      DOUBLE PRECISION t,hl,hr,ga,dx,cxi
      
      DOUBLE PRECISION h2,u2,S2,zmin,fzmin,z,fz,fzmax,zmax
      integer i
      
      if (dabs(t) .LE. 10d0**(-10)) then
         do i=1,xbc_len 

            ! x_j-1/2
            if (x(i) + 0.5*dx .LE. 0) then
               h(4*i-3) = hl
               h(4*i-2) = hl
               h(4*i-1) = hl
               h(4*i) = hl
               
            else
               h(4*i-3) = hr
               h(4*i-2) = hr
               h(4*i-1) = hr
               h(4*i) = hr
            end if
            
            u(4*i-3) = 0
            G(4*i-3) = 0

            
            u(4*i-2) = 0
            G(4*i-2) = 0            
            
            u(4*i-1) = 0
            G(4*i-1) = 0
            
            u(4*i) = 0
            G(4*i) = 0

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
               
         do i=1,xbc_len
         
            ! x_i-1/2
            cxi = x(i) - 0.5*dx
            if (cxi .LE. -t*dsqrt(ga*hl)) then
               h(4*i-3) = hl
               u(4*i-3) = 0
            else if ((cxi .GT. -t*dsqrt(ga*hl)) 
     .         .AND. (cxi .LE. t*(u2 - dsqrt(ga*h2))))  then
               h(4*i-3) = 4d0/(9d0*ga) *(dsqrt(ga*hl) - cxi/(2d0*t))**2
               u(4*i-3) = 2d0/3d0*(dsqrt(ga*hl) + cxi/t)
            else if ((cxi .GT. t*(u2 - dsqrt(ga*h2)) ) 
     .         .AND. (cxi .LT. t*S2))  then
               h(4*i-3) = h2    
               u(4*i-3) = u2  
            else
               h(4*i-3) = hr
               u(4*i-3) = 0
            end if
            
            G(4*i-3) = u(4*i-3)*h(4*i-3)
            
            ! x_i-1/6
            cxi = x(i) - dx/6.0
            if (cxi .LE. -t*dsqrt(ga*hl)) then
               h(4*i-2) = hl
               u(4*i-2) = 0
            else if ((cxi .GT. -t*dsqrt(ga*hl)) 
     .         .AND. (cxi .LE. t*(u2 - dsqrt(ga*h2))))  then
               h(4*i-2) = 4d0/(9d0*ga) *(dsqrt(ga*hl) - cxi/(2d0*t))**2
               u(4*i-2) = 2d0/3d0*(dsqrt(ga*hl) + cxi/t)
            else if ((cxi .GT. t*(u2 - dsqrt(ga*h2)) ) 
     .         .AND. (cxi .LT. t*S2))  then
               h(4*i-2) = h2    
               u(4*i-2) = u2  
            else
               h(4*i-2) = hr
               u(4*i-2) = 0
            end if
            
            G(4*i-2) = u(4*i-2)*h(4*i-2)         


            ! x_i+1/6
            cxi = x(i) + dx/6.0
            if (cxi .LE. -t*dsqrt(ga*hl)) then
               h(4*i-1) = hl
               u(4*i-1) = 0
            else if ((cxi .GT. -t*dsqrt(ga*hl)) 
     .         .AND. (cxi .LE. t*(u2 - dsqrt(ga*h2))))  then
               h(4*i-1) = 4d0/(9d0*ga) *(dsqrt(ga*hl) - cxi/(2d0*t))**2
               u(4*i-1) = 2d0/3d0*(dsqrt(ga*hl) + cxi/t)
            else if ((cxi .GT. t*(u2 - dsqrt(ga*h2)) ) 
     .         .AND. (cxi .LT. t*S2))  then
               h(4*i-1) = h2    
               u(4*i-1) = u2  
            else
               h(4*i-1) = hr
               u(4*i-1) = 0
            end if
            
            G(4*i-1) = u(4*i-1)*h(4*i-1)           

            ! x_i+1/2
            cxi = x(i) + 0.5*dx
            if (cxi .LE. -t*dsqrt(ga*hl)) then
               h(4*i) = hl
               u(4*i) = 0
            else if ((cxi .GT. -t*dsqrt(ga*hl)) 
     .         .AND. (cxi .LE. t*(u2 - dsqrt(ga*h2))))  then
               h(4*i) = 4d0/(9d0*ga) *(dsqrt(ga*hl) - cxi/(2d0*t))**2
               u(4*i) = 2d0/3d0*(dsqrt(ga*hl) + cxi/t)
            else if ((cxi .GT. t*(u2 - dsqrt(ga*h2)) ) 
     .         .AND. (cxi .LT. t*S2))  then
               h(4*i) = h2    
               u(4*i) = u2  
            else
               h(4*i) = hr
               u(4*i) = 0
            end if
            
            G(4*i) = u(4*i)*h(4*i)      

         end do 
      
      end if
     
      
      
      end
      


