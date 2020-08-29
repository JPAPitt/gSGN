c ====================================================================================
c Module of fortran subroutines use FDVM2 to solve generalised Serre - Green -Naghdi equations (gSGN)
c it takes initial conditions and returns solutions Q(x,t) where t > tend (Q = [h,G])
c the subroutine also calculates total conserved quantities initially and in the final solution
c ====================================================================================

      subroutine NumericalSolveForced(xbc_len,n_GhstCells,
     . tstart,tend,xstart,dx,dt,ga,theta, 
     . a0,a1,a2,a3,a4,a5,a6,a7,
     . tlist,tlist_len,IndivExpWdir,Indivexpwdirlen,
     . NormFileI)

      implicit none
      
      integer n_GhstCells,xbc_len,Indivexpwdirlen,
     . tlist_len,NormFileI
      
      DOUBLE PRECISION tstart,tend,xstart,dx,dt,ga,theta,
     . a0,a1,a2,a3,a4,a5,a6,a7
     
      CHARACTER(len=Indivexpwdirlen) IndivExpWdir
      
      DOUBLE PRECISION tlist(tlist_len)
     
      DOUBLE PRECISION xbc(xbc_len),
     . hbc(xbc_len),Gbc(xbc_len), ubc(xbc_len),
     . hbc_ana(xbc_len), Gbc_ana(xbc_len), ubc_ana(xbc_len)
     
      DOUBLE PRECISION currenttime,normG,normh,normu
      
      integer i,filecount
      CHARACTER(len=2) strct
      
      print *, xbc_len,n_GhstCells,tstart,tend,xstart,dx,dt,ga,theta,
     . tlist_len,Indivexpwdirlen
      
      !initial time
      currenttime  = tstart
      filecount = 1
      
      !Get Initial Conditions
      call Generatexbc(xstart,dx,xbc_len - 2*n_GhstCells,
     . n_GhstCells,xbc)
      call Initloc(xbc,xbc_len,currenttime,a0,a1,a2,a3,a4,
     . hbc,Gbc,ubc)
     
      
      ! Initial Conditions Out
      open(1, file = IndivExpWdir//'Init.dat') 
      do i = 1,xbc_len
         write(1,*) currenttime,xbc(i),hbc(i),
     .            Gbc(i),ubc(i)
      end do
      close(1)
            

      !evolve the system through time
      do while (currenttime  .LT. tend ) 
      
         if ( dabs(currenttime - tlist(filecount)) .LT. 0.51*dt  )  then
         
            call Initloc(xbc,xbc_len,currenttime,
     .         a0,a1,a2,a3,a4,
     .         hbc_ana,Gbc_ana,ubc_ana)
     
            write (strct,'(I2)') filecount
            open(3, file = IndivExpWdir// strct //'.dat') 
            do i = 1,xbc_len
               write(3,*) currenttime,xbc(i),hbc(i),hbc_ana(i),
     .            Gbc(i),Gbc_ana(i),ubc(i),ubc_ana(i)
            end do
            close(3)
            
            filecount = filecount + 1
         end if

         call EvolveStepWrap(xbc_len,hbc,Gbc,ubc,ga,theta,
     .      n_GhstCells,dx,dt,
     .      xbc,currenttime,a0,a1,a2,a3,a4,a5,a6,a7)
     
     
         currenttime  = currenttime  + dt
         print *, 'Current Time : ', currenttime 
      end do
      
      !calculate  
      call GetufromhG(xbc_len,hbc,Gbc,ubc,dx,n_GhstCells)

      call Initloc(xbc,xbc_len,currenttime,
     .         a0,a1,a2,a3,a4,
     .         hbc_ana,Gbc_ana,ubc_ana)
     
      open(4, file = IndivExpWdir//'End.dat') 
      do i = 1,xbc_len
         write(4,*) currenttime,xbc(i),hbc(i),hbc_ana(i),
     .            Gbc(i),Gbc_ana(i),ubc(i),ubc_ana(i)
      end do
      close(4)
                 
      
      !Convergence Norms
      normh = sqrt(sum((hbc_ana(:) - hbc(:)) **2)) / 
     . sqrt(sum(hbc_ana(:)**2))       
      normG = sqrt(sum((Gbc_ana(:) - Gbc(:)) **2)) / 
     . sqrt(sum(Gbc_ana(:)**2)) 
      normu = sqrt(sum((ubc_ana(:) - ubc(:)) **2)) / 
     . sqrt(sum(ubc_ana(:)**2)) 
     
      !write out norms
      write(NormFileI,*) dx,normh,normG,normu
      
      open(2, file = IndivExpWdir//'Param.dat') 
      !write out parameters      write(5,*) 'Experiment - Forced Solution, Gaussian Bump'
      write(2,*) 'xstart :',xstart
      write(2,*) 'xend :',xbc(xbc_len - n_GhstCells)
      write(2,*) 'xbc_len :',xbc_len
      write(2,*) 'dx = (x_end - x_start) / (x_len -1)  :' , dx
      write(2,*) 'n_GhstCells :',n_GhstCells
      write(2,*) 'tstart :', tstart
      write(2,*) 'tend :',tend 
      write(2,*) 'actual end time :', currenttime
      write(2,*) 'dt/dx :' , dt/dx
      write(2,*) 'dt = dx*(dt/dx)  :' , dt
      write(2,*) 'gravity :' , ga
      write(2,*) 'theta :' , theta
      write(2,*) 'a0 :' , a0
      write(2,*) 'a1 :' , a1
      write(2,*) 'a2 :' , a2
      write(2,*) 'a3 :' , a3
      write(2,*) 'a4 :' , a4
      write(2,*) 'a5 :' , a5
      write(2,*) 'a6 :' , a6
      write(2,*) 'a7 :' , a7
      close(2)
 
      end     
      



c  ********************************************************************************
c  Part 1 : Initial Conditions
c  Functions that generate Initial Conditions
c
c
c
c
c ********************************************************************************


c==============================
c Function to generate x values at all points: All points in cell 
c for reconstruction + ghost cells
c==============================
      subroutine Generatexbc(x_start,dx,x_len,n_GhstCells,xbc)
          
      implicit none
      DOUBLE PRECISION x_start,dx
      integer x_len,n_GhstCells
      DOUBLE PRECISION xbc(x_len + 2*n_GhstCells)
      
      integer i,xbc_len
      
      xbc_len = x_len + 2*n_GhstCells
      
      !Left boundary
      do i=1,x_len + 2*n_GhstCells  
         xbc(i)  = x_start + (i -1 - n_GhstCells )*dx
      end do  
                    
      end


c ====================
c Function to generate analytic h,u,G at all locations x
c ====================      
      subroutine Initloc(x,x_len,t,a0,a1,a2,a3,a4,h,G,u)
      implicit none
      
      INTEGER x_len
      DOUBLE PRECISION t,a0,a1,a2,a3,a4
      DOUBLE PRECISION h(x_len),G(x_len),u(x_len),x(x_len)
                   
      DOUBLE PRECISION PHI,EXPPHI1,dhdx,dudx,d2udx2
      
      integer i
           
      do i=1,x_len 
         PHI  = x(i) - a2*t
         EXPPHI1 = dexp(-PHI**2 / (2*a3))
         h(i) = EXPPHI1*a1 + a0
         u(i) = EXPPHI1*a4
         dhdx = -EXPPHI1*PHI*a1/a3
         dudx = -EXPPHI1*PHI*a4/a3
         d2udx2 = -a4*(-PHI**2 + a3)*dexp(-PHI**2/(2*a3))/a3**2
                  
         G(i) = h(i)*u(i) -
     .    (d2udx2*h(i)**3 + 3*dhdx*dudx*h(i)**2) / 3d0
      end do 
      
            
      end


c  ********************************************************************************
c  Part 2 : Finite Difference Solver
c  Functions thet solve finite difference matrix Au = G
c  for u. A is tridiagonal, so we solve using Thomas algorithm
c
c
c
c ********************************************************************************

c =====================
c Function that builds the three diagonals of the matrix produced by the finite difference approximation to elliptic equation
c It then uses a matrix solve on the diagonals, and the value of G to solve for u (for gSGN equations)
c ==================      
      subroutine GetufromhG(xbc_len,hbc,Gbc,ubc,dx,n_GhstCells)
      
      implicit none
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION dx
                   
      DOUBLE PRECISION subdiag1(xbc_len),
     . diag(xbc_len),
     . supdiag1(xbc_len),
     . RHS(xbc_len)
          
      DOUBLE PRECISION ht1,ht2,dhc
     
      integer i
                  
      !calculate diagonals in interior
      !set RHS B        
      do i=n_GhstCells+1,xbc_len - n_GhstCells 
      
         call RecondqpmLimGrad(hbc(i-2),hbc(i-1),hbc(i),hbc(i+1),
     . hbc(i+2),dx,dhc)
              
         ht1 = (hbc(i)**3/(dx*dx)) / 3d0
         ht2 = hbc(i)**2/(2.d0*dx)*(dhc)
         
         
         subdiag1(i)  = -ht1 + ht2
         diag(i) = hbc(i) + (2d0*ht1)
         supdiag1(i)  = -ht1 - ht2 
                  
         RHS(i) = Gbc(i)
      end do 
      
      !boundary conditions
      !first and last n_GhstCells  x n_GhstCells in tridiagmatrix
      !Should be identity
      do i = 1, n_GhstCells
         ! left
         subdiag1(i) = 0d0
         diag(i) = 1d0
         supdiag1(i) = 0d0
         RHS(i) = ubc(i)
         
         !right
         subdiag1(xbc_len - n_GhstCells + i) = 0d0
         diag(xbc_len - n_GhstCells + i) = 1d0
         supdiag1(xbc_len - n_GhstCells + i)= 0d0

         RHS(xbc_len - n_GhstCells + i) = 
     .      ubc(xbc_len - n_GhstCells + i)
      end do

      call ThomasTriSolve(xbc_len,subdiag1,diag,supdiag1,RHS,ubc)
      
      end


c ================
c ThomasTriSolve solves a tridiagonal matrix equation Ax = b, where A is tridiagonal
c =================   
      subroutine ThomasTriSolve(n,a,b,c,d,x)
      implicit none
      
      integer n
      DOUBLE PRECISION a(n), c(n), b(n), d(n),x(n)
      
      !  --- Local variables ---
      integer i
      DOUBLE PRECISION q
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
      
      end


c  ********************************************************************************
c  Part 3 : Evolution
c  Functions that evolve h^n, G^n  -> h^{n+1}, G^{n+1}
c  Using u solve, Kurganov method and Rk2 step
c
c
c ********************************************************************************


c ====
c Evolve Step Wrap, is the wrapper function that takes h^n, G^n and produces
c our second order approximation to h^n+1 and G^n+1, using RK time stepping
c ====
      subroutine EvolveStepWrap(xbc_len,hbc,Gbc,ubc,ga,theta,
     . n_GhstCells,dx,dt,
     . xbc,t,a0,a1,a2,a3,a4,a5,a6,a7)
     
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len),
     . xbc(xbc_len)
      DOUBLE PRECISION ga,theta,dx,dt,t,
     . a0,a1,a2,a3,a4,a5,a6,a7
     
      !local variables
      DOUBLE PRECISION hpbc(xbc_len),Gpbc(xbc_len),
     . hppbc(xbc_len),Gppbc(xbc_len)
     
      integer i
         
      !Get ubc, from h and G from FD method (ubc must be initialised with BC properly)
      call GetufromhG(xbc_len,hbc,Gbc,ubc,dx,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+1}, G^{n+1})
      call SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,theta,n_GhstCells,
     . dt,dx,xbc,t,a0,a1,a2,a3,a4,a5,a6,a7,
     . hpbc,Gpbc)
     
      !update BC's, since EulerStep function, assumes its constant
      !BC's
      call Initloc(xbc(1:n_GhstCells),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,hpbc(1:n_GhstCells),
     .   Gpbc(1:n_GhstCells),ubc(1:n_GhstCells))
     
      call Initloc(xbc(xbc_len - n_GhstCells + 1:xbc_len),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,
     .   hpbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   Gpbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   ubc(xbc_len - n_GhstCells + 1:xbc_len))
     
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
      call GetufromhG(xbc_len,hpbc,Gpbc,ubc,dx,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+2}, G^{n+2})     
      call SingleEulerStep(xbc_len,hpbc,Gpbc,ubc,ga,theta,n_GhstCells,
     . dt,dx,xbc,t + dt,a0,a1,a2,a3,a4,a5,a6,a7,
     . hppbc,Gppbc)
      
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      do i= 1,xbc_len
         hbc(i) = ( hbc(i ) + hppbc(i))/2d0
         Gbc(i) = ( Gbc(i ) + Gppbc(i))/2d0
      end do

      !BC's
      call Initloc(xbc(1:n_GhstCells),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,hbc(1:n_GhstCells),
     .   Gbc(1:n_GhstCells),ubc(1:n_GhstCells))
      
      call Initloc(xbc(xbc_len - n_GhstCells + 1:xbc_len),n_GhstCells,
     .   t + dt,a0,a1,a2,a3,a4,
     .   hbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   Gbc(xbc_len - n_GhstCells + 1:xbc_len),
     .   ubc(xbc_len - n_GhstCells + 1:xbc_len))
          
      
      end 
    

c ====
c SingleEulerStep
c produces new cell averages of h,G using forward Euler step,
c with flux approximated using Kurganov's method
c ====
      subroutine SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,
     . theta,n_GhstCells,dt,dx,
     . xbc,t,a0,a1,a2,a3,a4,a5,a6,a7,
     . hpbc,Gpbc)
     
     
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION dt,dx,ga,theta,
     . t,a0,a1,a2,a3,a4,a5,a6,a7,tau
      DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len),
     . hpbc(xbc_len),Gpbc(xbc_len), xbc(xbc_len)
     
     
      DOUBLE PRECISION cdhi,cdGi,cdui, fih,fiG,foh,foG
      
      integer i,ileft,iright
      
      !we assume boundary conditions are constant - will update outside this step
      ileft = 1
      iright = n_GhstCells
      do i = ileft,iright
         hpbc(i) = hbc(i)
         Gpbc(i) = Gbc(i)
         
         hpbc(xbc_len - n_GhstCells + i) = 
     .      hbc(xbc_len - n_GhstCells + i)
         Gpbc(xbc_len - n_GhstCells + i) = 
     .      Gbc(xbc_len - n_GhstCells + i)
      end do
      
      !Now we update interior, first calculating flux across left interior boundary 
      !then loop over cells to get flux in/out and thus updating cell average values
      
      !Do left boundary
      i = n_GhstCells
      
      !initial gradient of h,G across cell i
      call ReconLinLimGrad(hbc(i-1),hbc(i),hbc(i+1),theta,cdhi)
      call ReconLinLimGrad(Gbc(i-1),Gbc(i),Gbc(i+1),theta,cdGi)
      call ReconLinLimGrad(ubc(i-1),ubc(i),ubc(i+1),theta,cdui)
      
      !calculates foh,foG which is flux across x_{i +1/2}
      ! it also updates cdhi,cdGi to be gradient across cell i + 1
      
      tau = a6*(xbc(i) + 0.5*dx - a5*t) + a7
      
      call Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,tau,theta,dx,i,
     . foh,foG,cdhi,cdGi,cdui)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      !loop over interior cells (do not update ghost cells)
      ileft = n_GhstCells + 1
      iright = xbc_len - n_GhstCells
      do i = ileft, iright
      
         !calculates foh,foG which is flux across x_{i +1/2}
         ! it also updates cdhi,cdGi to be gradient across cell i + 1
         tau = a6*(xbc(i) + 0.5*dx - a5*t) + a7
         
         call Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,tau,theta,dx,i,
     . foh,foG,cdhi,cdGi,cdui)
     
         call Forcedht(xbc(i),t,a1,a2,a3,dhdt)
         
         call ForcedGt(xbc(i),t,a0,a1,a2,a3,a4,dGdt)
         
         call Forcedfluxhx(xbc(i),t,a0,a1,a2,a3,a4,dfluxhdx) 
         call ForcedfluxGx(xbc(i),t,ga,a0,a1,a2,a3,a4,a5,
     . a6,a7,dfluxGdx)      
         
         hpbc(i) = hbc(i) + dt*dhdt 
     .    + dt*( dfluxhdx  - (foh - fih)/dx ) 
     
         Gpbc(i) =  Gbc(i)+ dt*dGdt 
     .    + dt*( dfluxGdx  - (foG - fiG)/dx ) 


         !flux out becomes flux in on next cell
         fih = foh
         fiG = foG
         
      end do 
           
      end
      
c subroutine that given arrays, and i calculates flux across x_{i+1/2} for gSGN equations  
c note that - cdhi,cdGi is inout, on the way in cdhi,cdGi are our approximation to the gradient
c of h and G across cell i, on the way out cdhi,cdGi are our approximation to the gradient of h and G
c across cell i +1 (which will be cell i as the method moves to the next cell)    
      subroutine Fluxxiph(xbc_len,hbc,Gbc,ubc,ga,tau,theta,dx,i,
     . foh,foG,cdhi,cdGi,cdui)
     
       integer i,xbc_len
       DOUBLE PRECISION hbc(xbc_len),Gbc(xbc_len),ubc(xbc_len)
       DOUBLE PRECISION ga,tau,theta,dx,
     . cdhi,cdGi,cdui,foh,foG
     
       DOUBLE PRECISION cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l,
     . dhir,ddhir,dhip1l,ddhip1l
                
c     centered approximation
      !use slopes from last calculation         
      hir = hbc(i) + cdhi/2d0
      Gir = Gbc(i) + cdGi/2d0
      uir = ubc(i) + cdui/2d0
      
      call RecondqmLimGrad(ubc(i-2),ubc(i-1),ubc(i),ubc(i+1),
     .  dx,duir)
      call RecondqmLimGrad(hbc(i-2),hbc(i-1),hbc(i),hbc(i+1),
     .  dx,dhir)
      call ReconddqmLimGrad(hbc(i-3),hbc(i-2),hbc(i-1),hbc(i),hbc(i+1),
     . hbc(i+2),dx,ddhir)     

      !reconstruct values on right side of edge x_{i+1/2}
      call ReconLinLimGrad(hbc(i),hbc(i+1),hbc(i+2),theta,cdhip1)
      call ReconLinLimGrad(Gbc(i),Gbc(i+1),Gbc(i+2),theta,cdGip1)
      call ReconLinLimGrad(ubc(i),ubc(i+1),ubc(i+2),theta,cduip1)
      
      hip1l = hbc(i + 1) - cdhip1/2d0
      Gip1l = Gbc(i + 1) - cdGip1/2d0
      uip1l = ubc(i + 1) - cduip1/2d0
      
      call RecondqpLimGrad(ubc(i),ubc(i+1),ubc(i+2),ubc(i+3),
     .  dx,duip1l)
      call RecondqpLimGrad(hbc(i),hbc(i+1),hbc(i+2),hbc(i+3),
     .  dx,dhip1l)
      call ReconddqpLimGrad(hbc(i-1),hbc(i),hbc(i+1),hbc(i+2),
     .  hbc(i+3),hbc(i+4),dx,ddhip1l)
           
      sl  = min(0d0, uir - dsqrt(ga*hir) ,
     . uip1l - dsqrt(ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(ga*hir) ,
     . uip1l + dsqrt(ga*hip1l)  )
      
      !left and right flux
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 
     .      - (2d0/3d0)*hir**3*duir**2
     .      - (tau/2d0)*(2*hir*ddhir - (dhir**2))
     
           
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - (2d0/3d0)*hip1l**3*duip1l**2
     .      - (tau/2d0)*(2*hip1l*ddhip1l - (dhip1l**2))
     
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if 
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
      
      !return gradient of h,G across cell i + 1 (for next iteration)
      cdGi = cdGip1
      cdhi = cdhip1   
      cdui = cduip1 
      end

c ====
c ReconLinLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconLinLimGrad(qim1,qi,qip1,theta,cdq)
      DOUBLE PRECISION qim1,qi,qip1,theta,cdq,
     . fdq,mdq,bdq
      
c      fdq = qip1 - qi
c      mdq = 0.5*(qip1 - qim1)
c      bdq = qi - qim1    
c      call minmod(theta*fdq,mdq,theta*bdq,cdq)

      cdq = 0.5*(qip1 - qim1)      
      end

c ====
c RecondqmLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine RecondqmLimGrad(qim2,qim1,qi,qip1,dx,cdq)
      DOUBLE PRECISION qim2,qim1,qi,qip1,dx,cdq
      DOUBLE PRECISION mdq,bdq
      
c      bdq = (2*qi - 3*qim1 + qim2)  /dx
c      mdq = (qip1 - qi) /dx
c      call minmod(bdq,mdq,bdq,cdq)
      cdq = (qip1 - qi) /dx  
      end

c ====
c RecondqpLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====

      subroutine RecondqpLimGrad(qi,qip1,qip2,qip3,dx,cdq)
      DOUBLE PRECISION qi,qip1,qip2,qip3,dx,cdq
      DOUBLE PRECISION fdq,mdq
      
c      mdq = (qip1 - qi) /dx
c      fdq = (-2*qip1 + 3*qip2 - qip3) /dx      
c      call minmod(fdq,mdq,fdq,cdq)
      cdq = (qip1 - qi) /dx     
      end
           
c ====
c RecondqmLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine RecondqpmLimGrad(qim2,qim1,qi,qip1,
     . qip2,dx,cdq)
      DOUBLE PRECISION qim2,qim1,qi,qip1,qip2,dx,cdq
      DOUBLE PRECISION mdq,bdq,fdq
      
c      bdq = (qim2 - 4*qim1 + 3*qi)  /(2d0*dx)
c      mdq = (qip1 - qim1) /(2d0*dx)
c      fdq = (-3*qi + 4*qip1 - qip2)  /(2d0*dx)            
c      call minmod(bdq,mdq,fdq,cdq)

      cdq = (qip1 - qim1) /(2d0*dx)      
      end

c ====
c ReconddqmLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====

      subroutine ReconddqmLimGrad(qim3,qim2,qim1,qi,qip1,
     . qip2,dx,cddq)
     
      DOUBLE PRECISION qim3,qim2,qim1,qi,qip1,qip2,dx,cddq
      DOUBLE PRECISION mddq,bddq,bpddq
      
c      bddq = (5*qi - 13*qim1 + 11*qim2 - 3*qim3)/(2*dx**2)
c      bpddq = (3*qip1 - 7*qi + 5*qim1 - qim2)/(2*dx**2)  
c      mddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)       
c      call minmod(bpddq,bpddq,mddq,cddq)
      cddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)        
      end

      
c ====
c ReconddqpLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconddqpLimGrad(qim1,qi,qip1,qip2,qip3,
     . qip4,dx,cddq)
      DOUBLE PRECISION qim1,qi,qip1,qip2,qip3,qip4,dx,cddq
      DOUBLE PRECISION mddq,fddq,fpddq
      
c      mddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)
c      fddq  = (5*qip1 - 13*qip2 + 11*qip3 - 3*qip4) /(2*dx**2)
c      fpddq = (3*qi - 7*qip1 + 5*qip2 - qip3) /(2*dx**2)
c      call minmod(fpddq,fpddq,mddq,cddq)


      cddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)      
      end

c ==
c Analytic epxressions for LHS of Forced Solutions
c
c ===
      subroutine Forcedht(x,t,a1,a2,a3,dhdt) 
      DOUBLE PRECISION, intent(out) :: dhdt
      DOUBLE PRECISION, intent(in) :: x,t,a1,a2,a3
      DOUBLE PRECISION :: EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      dhdt = EXPPHI1*PHI*a1*a2/a3
      
      end subroutine Forcedht
      
      subroutine ForcedGt(x,t,a0,a1,a2,a3,a4,dGdt)
      DOUBLE PRECISION, intent(out) :: dGdt 
      DOUBLE PRECISION, intent(in) :: x,t,a0,a1,a2,a3,a4
      DOUBLE PRECISION :: EXPPHI1,PHI,hi,ui
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1

      dGdt = -2*EXPPHI1**2*PHI**3*a1**2*a2*hi*ui/a3**3 
     . + EXPPHI1*PHI*a1*a2*ui/a3 
     . - 3*PHI**3*a1*a2*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**3 
     . - PHI**3*a2*a4*hi**3*dexp(-PHI**2/(2*a3))/(3*a3**3) 
     . + 3*PHI*a1*a2*hi**2*ui*dexp(-PHI**2/(2*a3))/a3**2 
     . + PHI*a2*hi*ui/a3 + PHI*a2*a4*hi**3*dexp(-PHI**2/(2*a3))/a3**2
     
      end subroutine ForcedGt
      
      subroutine Forcedfluxhx(x,t,a0,a1,a2,a3,a4,dfluxhdx) 
      DOUBLE PRECISION, intent(out) :: dfluxhdx
      DOUBLE PRECISION, intent(in) :: x,t,a0,a1,a2,a3,a4
      DOUBLE PRECISION :: EXPPHI1,PHI,hi,ui,dhdx,dudx
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      
      dfluxhdx = ui*dhdx + hi*dudx
      
      end subroutine Forcedfluxhx
      
      subroutine ForcedfluxGx(x,t,ga,a0,a1,a2,a3,a4,a5,
     . a6,a7,dfluxGdx)
      DOUBLE PRECISION, intent(out) :: dfluxGdx
      DOUBLE PRECISION, intent(in) :: x,t,ga,a0,a1,a2,a3,a4,a5,
     . a6,a7
      DOUBLE PRECISION :: EXPPHI1,PHI,hi,ui
      DOUBLE PRECISION :: dhdx,d2hdx2,d3hdx3,dudx,d2udx2,d3udx3,
     . tau,dtaux
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      hi  = a0 + a1*EXPPHI1
      ui  = a4*EXPPHI1
      
      dhdx = -EXPPHI1*PHI*a1/a3
      d2hdx2 = -a1*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3hdx3 = PHI*a1*(-PHI**2 + 3*a3)*EXPPHI1/a3**3
      dudx = -EXPPHI1*PHI*a4/a3
      d2udx2 = -a4*(-PHI**2 + a3)*EXPPHI1/a3**2
      d3udx3 = PHI*a4*(-PHI**2 + 3*a3)*EXPPHI1/a3**3
      
      tau = a6*(x - a5*t) + a7
      dtaux = a6
     
      dfluxGdx = -(4d0/3d0)*d2udx2*dudx*hi**3 
     . - d3hdx3*hi*tau
     . - 2.0*dhdx*dudx**2*hi**2 
     . + dhdx*ga*hi 
     . + dudx*(-d2udx2*(hi**3)/3d0 - dhdx*dudx*hi**2 + hi*ui) 
     . + ui*(-d2hdx2*dudx*hi**2 - 2*d2udx2*dhdx*hi**2 
     . - d3udx3*(hi**3)/3d0 
     . - 2*dhdx**2*dudx*hi + dhdx*ui + dudx*hi) 
     . + (-d2hdx2*hi + (dhdx**2)/2d0)*dtaux
     
      end subroutine ForcedfluxGx



c ===
c minmod is used for limiting and returns smallest size argument if all arguments have same sign,
c otherwise it returns 0.
c ===      
      subroutine minmod(a,b,c,d)
      implicit none
      DOUBLE PRECISION a,b,c, d
      
      if ((a .gt. 0d0) .and. (b .gt. 0d0) .and. (c .gt. 0d0)) then
         d = min(a,b,c)
      else if ((a .lt. 0d0) .and. (b .lt. 0d0) .and. (c .lt. 0d0)) then
         d = max(a,b,c)
      else
         d = 0d0
      end if
      
      end 


c ====================================
c Miscellaneous functions that can be used
c
c ====================================

      subroutine LenTrim(String,n,resultval)
      CHARACTER(len=n) String

      integer n, resultval

      resultval = 1

      do while ((String(resultval:resultval) .NE. ' ' ) .AND.
     . (resultval .LE. n )  )
         resultval = resultval + 1
      end do

      resultval = resultval - 1
      end       


