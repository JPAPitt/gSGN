c ====================================================================================
c Module of fortran subroutines use FDVM2 to solve generalised Serre - Green -Naghdi equations (gSGN)
c it takes initial conditions and returns solutions Q(x,t) where t > tend (Q = [h,G])
c the subroutine also calculates total conserved quantities initially and in the final solution
c ====================================================================================

c=====================================
c Numerical solver program, evolving system through time and producing Q(x,t) where t > tend
c The solver uses arrays xbc,hbc_init,Gbc_init and ubc_init
c to produce: 
c initial energy - Energs_init
c solution for h at final time hbc_fin
c solution for G at final time Gbc_fin
c solution for u at final time ubc_fin
c currenttime -time of produced solution currenttime = tstart + dt*number of timesteps performed
c Notes:
c 1. ubc_init only needs to have u defined at ghost cells
c function will update interior using up to date values of h and G.
c
c 2. Solver assumes beta's are constant and ghost cell values are constant- constant dirichlet boundary conditions
c=====================================
      subroutine NumericalSolveTSPrint(tstart,tend,
     . ga,dx,dt,n_GhstCells,xbc_len,
     . xbc,hbc_init,ubc_init,
     . phbc_init,pubc_init,
     . a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7,
     . currenttime,hbc_fin,ubc_fin,
     . Energ_Init, Energ_Fin,
     . tlist,tlist_len,ExpWdir,expwdirlen)
     
     
      implicit none
      
      integer n_GhstCells,xbc_len,expwdirlen,tlist_len
      CHARACTER(len=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,beta1,beta2,
     . theta,dx,dt,currenttime,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7
      DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
     . ubc_init(xbc_len),hbc_fin(xbc_len),
     . ubc_fin(xbc_len),phbc_init(xbc_len),
     . pubc_init(xbc_len),phbc_fin(xbc_len),
     . pubc_fin(xbc_len)
     
      DOUBLE PRECISION tlist(tlist_len)
      
      DOUBLE PRECISION Energ_Init(3), Energ_Fin(3)
      
      integer i,ileft,iright,filecount
      CHARACTER(len=2) strct
      
      
      
      !initial time
      currenttime  = tstart
      filecount = 1
      
      !loop over and set hbc_fin,Gbc_fin to initial conditions
      ileft = 1
      iright = xbc_len
      do i = ileft,iright
         hbc_fin(i) = hbc_init(i) 
         ubc_fin(i) = ubc_init(i)
         phbc_fin(i) = phbc_init(i) 
         pubc_fin(i) = pubc_init(i)
      end do
      

      !calculate initial Energies     
      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,ga,beta1,beta2,
     . n_GhstCells,dx,Energ_Init)
     
      print *, 1


      !evolve the system through time
      do while (currenttime  .LT. tend )   
      
         if ((currenttime + dt .GE. tlist(filecount)) 
     .      .OR. (filecount .EQ. 1 ))  then
              
            write (strct,'(I2)') filecount
            open(8, file = ExpWdir// strct //'.dat') 
            do i = 1,xbc_len
               write(8,*) currenttime,xbc(i),hbc_fin(i),
     .            ubc_fin(i)
            end do
            close(8)
            
            filecount = filecount + 1
         end if

         call EvolveWrap(hbc_fin,ubc_fin,phbc_fin,pubc_fin,ga,
     .    beta1,beta2,dx,dt,xbc_len,n_GhstCells,
     .    xbc,currenttime,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7) 
     
         currenttime  = currenttime  + dt
         print *, 'Current Time : ', currenttime 
      end do
             
      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,ga,beta1,beta2,
     . n_GhstCells,dx,Energ_Fin)
                 
      end



c  ********************************************************************************
c  Part 2 : Finite Difference Solver
c  Functions thet solve finite difference matrix Au = G
c  for u. A is tridiagonal, so we solve using Thomas algorithm
c
c
c
c ********************************************************************************

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

    
      subroutine EvolveWrap(hbc,ubc,phbc,pubc,g,beta1,beta2,dx,dt,
     . xbc_len,n_GhstCells,
     . xbc,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)   
      integer n_GhstCells,xbc_len
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),phbc(xbc_len),
     . pubc(xbc_len),nhbc(xbc_len),nubc(xbc_len), xbc(xbc_len)
      DOUBLE PRECISION g,dx,dt,beta1,beta2,t,a0,a1,a2,a3,a4,a5
      
      call Evolveh(hbc,ubc,phbc,g,dx,dt,xbc_len,n_GhstCells,nhbc,
     . xbc,t,a0,a1,a2,a3,a4,a5 )
      call Evolveu(hbc,ubc,pubc,beta1,beta2,g,dx,dt,xbc_len,
     . n_GhstCells,nubc,xbc,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
     
     
      !Return updated h,u,ph,pu
      do i = 1, xbc_len
         phbc(i) = hbc(i)
         pubc(i) = ubc(i)
         hbc(i) = nhbc(i)
         ubc(i) = nubc(i)
      end do
      
      end    


      subroutine Evolveh(hbc,ubc,phbc,g,dx,dt,xbc_len,n_GhstCells,nhbc,
     . xbc,t,a0,a1,a2,a3,a4,a5)
      
      integer n_GhstCells,xbc_len,i,j
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),phbc(xbc_len),
     . nhbc(xbc_len), xbc(xbc_len)
      
      DOUBLE PRECISION g,dx,dt,chx,cux,t,a0,a1,a2,a3,a4,a5,
     .  h,u,dhdt,dhdx,dudx
      
      !interior
      do i = n_GhstCells + 1, xbc_len - n_GhstCells
        call ForcedhTerms(xbc(i),t,a0,a1,a2,a3,a4,a5,h,u,dhdt,dhdx,dudx) 
        chx = (hbc(i+1) - hbc(i-1)) / (2*dx)
        cux = (ubc(i+1) - ubc(i-1))/ (2*dx)
        nhbc(i) = phbc(i)  - 2*dt*(hbc(i)*cux + ubc(i)*chx)
     .   + 2*dt*(dhdt + h*dudx + u*dhdx)
      
         
      end do
      
      ! boundaries
      do i = 1,n_GhstCells
        call ForcedhTerms(xbc(i),t+dt,a0,a1,a2,a3,a4,a5,
     .   h,u,dhdt,dhdx,dudx) 
        nhbc(i) = h
        j = xbc_len - n_GhstCells + i
        call ForcedhTerms(xbc(j),t+dt,a0,a1,a2,a3,a4,a5,
     .   h,u,dhdt,dhdx,dudx) 
        nhbc(j) = h
      end do
      
      
      end
      
      subroutine Evolveu(hbc,ubc,pubc,beta1,beta2,g,dx,dt,xbc_len,
     . n_GhstCells,nubc,
     . xbc,t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7)
      integer n_GhstCells,xbc_len,i,j
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),pubc(xbc_len),
     . nubc(xbc_len), subdiag1(xbc_len), diag(xbc_len),
     . supdiag1(xbc_len), RHS(xbc_len),
     . xbc(xbc_len)
      DOUBLE PRECISION beta1,beta2, unterms,
     . chx,cux,chxx,cuxx,chxxx,cuxxx,pux,puxx,t,a0,a1,a2,a3,a4,a5,
     . b1a6,b1a7,b2a6,b2a7,futerm
     
      !print *, t,a0,a1,a2,a3,a4,a5,b1a6,b1a7,b2a6,b2a7
      beta1 = b1a7
      beta2 = b2a7
      do i = n_GhstCells + 1, xbc_len - n_GhstCells


         chx = 0.5*(hbc(i+1) - hbc(i-1))/dx
         cux = 0.5*(ubc(i+1) - ubc(i-1))/dx
         pux = 0.5*(pubc(i+1) - pubc(i-1))/dx
         chxx = (hbc(i+1) - 2*hbc(i) + hbc(i-1))/ (dx**2)
         cuxx = (ubc(i+1) - 2*ubc(i) + ubc(i-1))/ (dx**2)
         puxx = (pubc(i+1) - 2*pubc(i) + pubc(i-1))/ (dx**2)
         chxxx = (hbc(i+2) - 2*hbc(i+1)  + 2*hbc(i-1) - hbc(i-2))/
     .     (2*dx**3)
         cuxxx = (ubc(i+2) - 2*ubc(i+1)  + 2*ubc(i-1) - ubc(i-2))/
     .     (2*dx**3)
     
         unterms = 3*beta1*chx*cux**2*hbc(i)/2 
     .    - 3*beta1*chx*cuxx*hbc(i)*ubc(i)/2 
     .    + beta1*cux*cuxx*hbc(i)**2/2 
     .    - beta1*cuxxx*hbc(i)**2*ubc(i)/2 
     .    - beta2*chx**3*g/2 
     .    - 2*beta2*chx*chxx*g*hbc(i) 
     .    - beta2*chxxx*g*hbc(i)**2/2
     .    + chx*g 
     .    + cux*ubc(i)
         
         
         unm1terms = 3*beta1*chx*hbc(i)*pux/(4*dt) +
     .     beta1*hbc(i)**2*puxx/(4*dt) - pubc(i)/(2*dt)
         
         
         subdiag1(i) = beta1*hbc(i)*(3*chx*dx - 2*hbc(i))/(8*dt*dx**2)
         diag(i) = (beta1*hbc(i)**2 + dx**2)/(2*dt*dx**2)
         supdiag1(i) = -beta1*hbc(i)*(3*chx*dx + 2*hbc(i))/(8*dt*dx**2)
     
         
         call ForceduTerm(xbc(i),t,g,a0,a1,a2,a3,a4,b1a7,b2a7,futerm) 

         RHS(i) = futerm -unm1terms -unterms
         
     
      end do
      
      do i = 1, n_GhstCells
         ! left
         subdiag1(i) = 0d0
         diag(i) = 1d0
         supdiag1(i) = 0d0
         
         call ForcedhTerms(xbc(i),t+dt,a0,a1,a2,a3,a4,a5,
     .   h,u,dhdt,dhdx,dudx) 
     
         RHS(i) = u

         
         !right
         subdiag1(xbc_len - n_GhstCells + i) = 0d0
         diag(xbc_len - n_GhstCells + i) = 1d0
         supdiag1(xbc_len - n_GhstCells + i)= 0d0

         call ForcedhTerms(xbc(xbc_len - n_GhstCells + i),
     .    t+dt,a0,a1,a2,a3,a4,a5,h,u,dhdt,dhdx,dudx) 
         RHS(xbc_len - n_GhstCells + i) = u

      end do
      
      call ThomasTriSolve(xbc_len,subdiag1,diag,supdiag1,RHS,nubc)
      

      end



      subroutine ForcedhTerms(x,t,a0,a1,a2,a3,a4,a5,h,u,dhdt,dhdx,dudx) 
      DOUBLE PRECISION dhdt,dhdx,dudx,h,u
      DOUBLE PRECISION x,t,a0,a1,a2,a3,a4,a5
      DOUBLE PRECISION EXPPHI1,PHI
      
      PHI  = x - a2*t
      EXPPHI1 = dexp(-PHI**2 / (2*a3))
      h  = a0 + a1*EXPPHI1
      u  = a4*EXPPHI1
      dhdt = EXPPHI1*PHI*a1*a2/a3
      dhdx = -EXPPHI1*PHI*a1/a3
      dudx = -EXPPHI1*PHI*a4/a3
      
      end subroutine ForcedhTerms

      subroutine ForceduTerm(x,t,ga,a0,a1,a2,a3,a4,b1a7,b2a7,futerm) 
      DOUBLE PRECISION x,t,ga,a0,a1,a2,a3,a4,b1a7,b2a7
      DOUBLE PRECISION futerm,EXPPHI,PPHI
      
      PPHI  = x - a2*t
      EXPPHI = dexp(-PPHI**2 / (2*a3))
      
      futerm = EXPPHI**3*a1**3*b2a7*ga*(-2*a2*t + 2*x)**3/(16*a3**3) 
     . - 3*EXPPHI**3*a1*a4**2*b1a7*(EXPPHI*a1 + a0)*(-PPHI**2/a3 + 1)
     . *(-2*a2*t + 2*x)/(4*a3**2) 
     . - 3*EXPPHI**3*a1*a4**2*b1a7*(EXPPHI*a1 + a0)
     . *(-2*a2*t + 2*x)**3/(16*a3**3) 
     . - EXPPHI**2*PPHI*a4**2*b1a7*(EXPPHI*a1 + a0)**2
     . *(-PPHI**2/a3 + 3)/(2*a3**2) 
     . + EXPPHI**2*a1**2*b2a7*ga*(EXPPHI*a1 + a0)
     . *(PPHI**2/a3 - 1)*(-2*a2*t + 2*x)/a3**2 
     . + 3*EXPPHI**2*a1*a2*a4*b1a7*(EXPPHI*a1 + a0)
     . *(-PPHI**2/a3 + 1)*(-2*a2*t + 2*x)/(4*a3**2) 
     . - EXPPHI**2*a4**2*(-2*a2*t + 2*x)/(2*a3) 
     . + EXPPHI**2*a4**2*b1a7*(EXPPHI*a1 + a0)**2
     . *(-PPHI**2/a3 + 1)*(-2*a2*t + 2*x)/(4*a3**2) 
     . + EXPPHI*PPHI*a1*b2a7*ga*(EXPPHI*a1 + a0)**2
     . *(PPHI**2/a3 - 3)/(2*a3**2) 
     . + EXPPHI*PPHI*a2*a4/a3 
     . + EXPPHI*PPHI*a2*a4*b1a7*(EXPPHI*a1 + a0)**2
     . *(-PPHI**2/a3 + 3)/(2*a3**2) 
     . - EXPPHI*a1*ga*(-2*a2*t + 2*x)/(2*a3)
      
      end subroutine ForceduTerm
      
 

   

c  ********************************************************************************
c  Part 4 : Analysis Functions
c  Functions that do analyses on solutions
c  1. For Conservation -h, uh, G and Energy  (was at least second order accurate (FD solve, holding it back to second order for G))
c  
c
c
c ********************************************************************************
      
c =====
c Function to generate quartic coefficients using
c f(x_{j-2}),f(x_{j-1}),f(x_{j}),f(x_{j+1}),f(x_{j+2})
c =====      
      subroutine QuarticInterp(qjm2,qjm1,qj,qjp1,qjp2,dx,QuartCoeff)
     
      DOUBLE PRECISION qjm2,qjm1,qj,qjp1,qjp2,dx
      DOUBLE PRECISION QuartCoeff(5)
      
      QuartCoeff(1) = (qjp2 - 4*qjp1 + 6*qj - 4*qjm1 + qjm2) /
     . (24* (dx**4))
      QuartCoeff(2) = (qjp2 - 2*qjp1 + 2*qjm1 - qjm2) /
     . (12* (dx**3))  
      QuartCoeff(3) = (-qjp2 + 16*qjp1 - 30*qj + 16*qjm1 - qjm2) /
     . (24* (dx**2))  
      QuartCoeff(4) = (-qjp2 + 8*qjp1 - 8*qjm1 + qjm2) /
     . (12*dx)
      QuartCoeff(5) = qj
      end


c ========================
c Functions to evaluate Quartic at x, with quartic centered around x_j
c xmxj = x - x_j 
c ====================      
      subroutine QuarticCoeffEvalxj(QuarticCoeff,xmxj,qatxj)
      
      DOUBLE PRECISION QuarticCoeff(5)
      DOUBLE PRECISION xmxj,qatxj
      
      qatxj = QuarticCoeff(1)*xmxj**4 + QuarticCoeff(2)*xmxj**3 +
     . QuarticCoeff(3)*xmxj**2 + QuarticCoeff(4)*xmxj +
     . QuarticCoeff(5)
     
      end
      
      subroutine QuarticCoeffEvalGradxj(QuarticCoeff,xmxj,dqatxj)
      
      DOUBLE PRECISION QuarticCoeff(5)
      DOUBLE PRECISION xmxj,dqatxj
      
      dqatxj = 4*QuarticCoeff(1)*xmxj**3 + 3*QuarticCoeff(2)*xmxj**2 +
     . 2*QuarticCoeff(3)*xmxj + QuarticCoeff(4)
     
      end

c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      subroutine AllEnergiesIntegralCell(xbc_len,h,u,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      integer j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION CellEnergies(3)
      
      integer i
      DOUBLE PRECISION fGPe(3),sGPe(3),tGPe(3)
      
      DOUBLE PRECISION GPmxj,hGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5)
      
      call QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      call QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      fGPe(1) = hGP
      fGPe(2) = hGP*uGP
      fGPe(3) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !second gauss point
      GPmxj = 0.0 
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      sGPe(1) = hGP
      sGPe(2) = hGP*uGP
      sGPe(3) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !third gauss point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      tGPe(1) = hGP
      tGPe(2) = hGP*uGP
      tGPe(3) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !weight the values at gauss points to get approximate integral over cell
      do i = 1,3
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .+ (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      end do
      
      end
      
      !Function to sum all energies
      subroutine TotalEnergy(xbc_len,hbc,ubc,ga,beta1,beta2,
     . n_GhstCells,dx,TotEnergVals)
      
      integer xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION TotEnergVals(3)
      
      DOUBLE PRECISION CellEnergVals(3)
      integer i,j
            
      !running totals for energy values, start at 0
      do i = 1,3
         TotEnergVals(i) = 0.0
      end do
       
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      do j= n_GhstCells + 1, xbc_len - n_GhstCells
         call  AllEnergiesIntegralCell(xbc_len,hbc,ubc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,3
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end

      


