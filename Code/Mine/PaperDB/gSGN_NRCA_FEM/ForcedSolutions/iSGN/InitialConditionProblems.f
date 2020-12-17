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
      subroutine SerreSoliton(x_bc,xbc_len,t,ga,a0,a1,dx,hcub_bc,
     .   ucub_bc,Gcub_bc,cubbc_len)
     
      implicit none
      
      
      integer xbc_len,cubbc_len
      real*8 x_bc(xbc_len)
      real*8 t,a0,a1,ga,dx
      real*8 hcub_bc(cubbc_len),Gcub_bc(cubbc_len),ucub_bc(cubbc_len)
      
      real*8 k,c,phi,sechkphi
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
      
c ====================
c Function to generate analytic h,u,G at all locations x
c ====================      



