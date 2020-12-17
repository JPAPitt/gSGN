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
     . ga,beta1,beta2,tolcalc,dx,dt,n_GhstCells,xbc_len,cubbc_len,
     . xbc,hcubbc_init,Gcubbc_init,ucubbc_init,
     . hcubbc_fin,Gcubbc_fin,ucubbc_fin,
     . currenttime,Energ_Init, Energ_Fin,
     . tlist,tlist_len,ExpWdir,expwdirlen)
     
     
      implicit none
      
      integer n_GhstCells,xbc_len,cubbc_len,expwdirlen,tlist_len
      CHARACTER(len=expwdirlen) ExpWdir
      DOUBLE PRECISION tstart,tend,ga,beta1,beta2,
     . tolcalc,dx,dt,currenttime
      DOUBLE PRECISION xbc(xbc_len),hcubbc_init(cubbc_len),
     . Gcubbc_init(cubbc_len),ucubbc_init(cubbc_len),
     . hcubbc_fin(cubbc_len),Gcubbc_fin(cubbc_len),
     . ucubbc_fin(cubbc_len), habc(xbc_len),Gabc(xbc_len)

     
      DOUBLE PRECISION tlist(tlist_len)
      
      DOUBLE PRECISION Energ_Init(4), Energ_Fin(4)
      
      integer i,filecount
      CHARACTER(len=2) strct
      
      !initial time
      currenttime  = tstart
      filecount = 1

      !loop over and set hbc_fin,Gbc_fin to initial conditions
      do i = 1,cubbc_len
         hcubbc_fin(i) = hcubbc_init(i) 
         Gcubbc_fin(i) = Gcubbc_init(i) 
         ucubbc_fin(i) = ucubbc_init(i) 
      end do
      
      !calculate initial Energies     
c      call FEMforU(hcubbc_fin,Gcubbc_fin,ucubbc_fin,beta1,dx,xbc_len,
c     .   cubbc_len,3*xbc_len + 1,n_GhstCells)

      call TotalEnergy(xbc_len,cubbc_len,hcubbc_fin,ucubbc_fin,
     . Gcubbc_fin,ga, beta1,beta2,n_GhstCells,dx,Energ_Init)
     
      ! calculate initial cell averages
      call CubicToCellAverage(cubbc_len,xbc_len,hcubbc_init,
     .   dx,habc)
      call CubicToCellAverage(cubbc_len,xbc_len,Gcubbc_init,
     .   dx,Gabc)
     
      !evolve the system through time
      do while (currenttime  .LT. tend )   
      
         if ((currenttime + dt .GE. tlist(filecount)) 
     .      .OR. (filecount .EQ. 1 ))  then
              
            write (strct,'(I2)') filecount
            open(21, file = ExpWdir// strct //'.dat') 
            do i = 1,xbc_len
               write(21,*) currenttime,xbc(i),hcubbc_fin(i),
     .            Gcubbc_fin(i),ucubbc_fin(i)
            end do
            close(21)
            
            filecount = filecount + 1
         end if   
         
         call EvolveStepWrap(xbc_len,cubbc_len,habc,Gabc,
     . ga,beta1,beta2,tolcalc,n_GhstCells,dx,dt,
     . hcubbc_fin,Gcubbc_fin,ucubbc_fin)

     
         currenttime  = currenttime  + dt
         print *, 'Current Time : ', currenttime 
         
      end do  
      
      !calculate cubics
      call CellAverageToCubic(cubbc_len,xbc_len,
     . n_GhstCells,habc,dx,hcubbc_fin,tolcalc)
 
      call CellAverageToCubic(cubbc_len,xbc_len,
     . n_GhstCells,Gabc,dx,Gcubbc_fin,tolcalc)
      
      !calculate end energies  
      call FEMforU(hcubbc_fin,Gcubbc_fin,ucubbc_fin,beta1,dx,xbc_len,
     .   cubbc_len,3*xbc_len + 1,n_GhstCells)
      
      call TotalEnergy(xbc_len,cubbc_len,hcubbc_fin,ucubbc_fin,
     . Gcubbc_fin,ga, beta1,beta2,n_GhstCells,dx,Energ_Fin)
     
      end


c  ********************************************************************************
c  Reconstruction
c  Functions that move between cubic and cell average representation
c
c
c ********************************************************************************

c Cubic to Cell average

      subroutine CubicToCellAverage(cubic_len,cellavg_len,cubic_q,
     .   dx,cellavg_q)
      implicit none
      integer cubic_len,cellavg_len
      DOUBLE PRECISION dx
      DOUBLE PRECISION cubic_q(cubic_len),cellavg_q(cellavg_len)
      
      integer i
      DOUBLE PRECISION  qjmh, qjph,qjms, qjps,pb,pd
      
      do i=1,cellavg_len
           qjmh = cubic_q(4*i  - 3)
           qjms = cubic_q(4*i  - 2)
           qjps = cubic_q(4*i  - 1)
           qjph = cubic_q(4*i)
           
           pb = (9*qjmh - 9*qjms - 9*qjps  + 9*qjph )/ (4*dx**2)
           pd = (-qjmh + 9*qjms + 9*qjps  - qjph )/ 16.0
           
           cellavg_q(i) = 1.0/dx*(2*pb/3.0*(dx/2)**3 +  2*pd*(dx/2) )
      end do
     
      

      end


c Helper Functions - Polynomial Coeffs

      subroutine minmod(a,b,c,d)
      implicit none
      DOUBLE PRECISION a,b,c,d,omin,omax
      omin = min(a,b,c)
      omax = max(a,b,c)
      
      if (omin .gt. 0d0) then
         d = omin
      else if (omax .lt. 0d0) then
         d = omax
      else
         d = 0d0
      end if
c      print *, 'MM I',a,b,c,d
c      print *, 'MM O',e,omin,omax
      end 

      subroutine P1jtojp1AB(qaj,qajp1,dx,A,B)
      implicit none
      DOUBLE PRECISION qaj,qajp1,dx,A,B
      A= (qajp1 - qaj)/ dx
      B = qaj
      end

      subroutine P1jm1tojAB(qajm1,qaj,dx,A,B)
      implicit none
      DOUBLE PRECISION qaj,qajm1,dx,A,B
      A = (qaj - qajm1)/ dx
      B = qaj
      end
      
      subroutine P2jm1tojp1ABC(qajm1,qaj,qajp1,dx,A,B,C)
      implicit none
      DOUBLE PRECISION qaj,qajp1,qajm1,dx,A,B,C
      A = (qajp1 - 2*qaj + qajm1)/ (2*dx**2)
      B = (qajp1 - qajm1)/ (2*dx)
      C =  13*qaj/12 - qajp1/24 - qajm1/24
      end

      
      subroutine P2jtojp2ABC(qaj,qajp1,qajp2,dx,A,B,C)
      implicit none
      DOUBLE PRECISION qaj,qajp1,qajp2,dx,A,B,C      
      A = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
      B = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
      C =  -qajp2/24 + qajp1/12 + 23*qaj/24
      end


      subroutine P2jm2tojABC(qajm2,qajm1,qaj,dx,A,B,C)
      implicit none
      DOUBLE PRECISION qaj,qajm1,qajm2,dx,A,B,C
      A = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
      B = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
      C =  -qajm2/24 + qajm1/12 + 23*qaj/24
      end
      
      subroutine P4jm2tojp2ABCDE(qajm2,qajm1,qaj,qajp1,qajp2,dx,
     . A,B,C,D,E)
      implicit none
      DOUBLE PRECISION qaj,qajm1,qajm2,qajp1,qajp2,dx,A,B,C,D,E
      A = -qajp1/(6*dx**4) + qajp2/(24*dx**4) - qajm1/(6*dx**4) 
     . + qajm2/(24*dx**4) + qaj/(4*dx**4)
      B = -qajp1/(6*dx**3) + qajp2/(12*dx**3) + qajm1/(6*dx**3) 
     . - qajm2/(12*dx**3)
      C = 3*qajp1/(4*dx**2) - qajp2/(16*dx**2) + 3*qajm1/(4*dx**2) 
     . - qajm2/(16*dx**2) - 11*qaj/(8*dx**2)
      D = 17*qajp1/(24*dx) - 5*qajp2/(48*dx) - 17*qajm1/(24*dx) 
     . + 5*qajm2/(48*dx)
      E = -29*qajp1/480 + 3*qajp2/640 - 29*qajm1/480 + 3*qajm2/640 
     . + 1067*qaj/960
      end

      subroutine P3jtojp3ABCD(qaj,qajp1,qajp2,qajp3,dx,A,B,C,D)
      implicit none
      DOUBLE PRECISION qaj,qajp1,qajp2,qajp3,dx,A,B,C,D
      A = (3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
      B = (-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
      C = (69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
      D = 5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
      end

      subroutine P3jm3tojABCD(qajm3,qajm2,qajm1,qaj,dx,A,B,C,D)
      implicit none
      DOUBLE PRECISION qajm3,qajm2,qajm1,qaj,dx,A,B,C,D
      A = (-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
      B = (-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
      C = (-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
      D = 5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
      end
      
      subroutine P6jm3tojp3ABCDEFG(qajm3,qajm2,qajm1,qaj,
     . qajp1,qajp2,qajp3,dx,A,B,C,D,E,F,G)
      implicit none
      DOUBLE PRECISION qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,
     . dx,A,B,C,D,E,F,G
      A = qajp1/(48*dx**6) - qajp2/(120*dx**6) + qajp3/(720*dx**6) 
     . + qajm1/(48*dx**6) - qajm2/(120*dx**6) + qajm3/(720*dx**6) 
     . - qaj/(36*dx**6)
     
      B = qajp1/(48*dx**5) - qajp2/(60*dx**5) + qajp3/(240*dx**5) 
     . - qajm1/(48*dx**5) + qajm2/(60*dx**5) - qajm3/(240*dx**5)
     
      C = -19*qajp1/(64*dx**4) + 3*qajp2/(32*dx**4) 
     . - 5*qajp3/(576*dx**4) - 19*qajm1/(64*dx**4) 
     . + 3*qajm2/(32*dx**4) - 5*qajm3/(576*dx**4) 
     . + 61*qaj/(144*dx**4)
     
      D = -83*qajp1/(288*dx**3) + 13*qajp2/(72*dx**3) 
     . - 7*qajp3/(288*dx**3) + 83*qajm1/(288*dx**3) 
     . - 13*qajm2/(72*dx**3) + 7*qajm3/(288*dx**3)
     
      E = 229*qajp1/(256*dx**2) - 77*qajp2/(640*dx**2) 
     . + 37*qajp3/(3840*dx**2) + 229*qajm1/(256*dx**2) 
     . - 77*qajm2/(640*dx**2) + 37*qajm3/(3840*dx**2) 
     . - 301*qaj/(192*dx**2)
     
      F = 1891*qajp1/(2304*dx) - 559*qajp2/(2880*dx) 
     . + 259*qajp3/(11520*dx) - 1891*qajm1/(2304*dx) 
     . + 559*qajm2/(2880*dx) - 259*qajm3/(11520*dx)
      
      G = -7621*qajp1/107520 + 159*qajp2/17920 
     . - 5*qajp3/7168 - 7621*qajm1/107520 + 159*qajm2/17920 
     . - 5*qajm3/7168 + 30251*qaj/26880
      end

c Cell Average To Cubic
    
     
      subroutine CellAverageToCubic(cubicbc_len,cellavgbc_len,
     . nGhstCells,cellavgbc_q,dx,cubicbc_q,calctol)
      implicit none
      integer cubicbc_len,cellavgbc_len,nGhstCells
      DOUBLE PRECISION dx,calctol
      DOUBLE PRECISION cubicbc_q(cubicbc_len),cellavgbc_q(cellavgbc_len)
      
      integer i
      DOUBLE PRECISION  p3jm3toja3,p3jm3toja2,p3jm3toja1,p3jm3toja0,
     . p3jtojp3a3,p3jtojp3a2,p3jtojp3a1,p3jtojp3a0, 
     . p6jm3tojp3a6,p6jm3tojp3a5,p6jm3tojp3a4,p6jm3tojp3a3,
     . p6jm3tojp3a2,p6jm3tojp3a1,p6jm3tojp3a0,
     . p2jm2toja2,p2jm2toja1,p2jm2toja0,
     . p2jtojp2a2,p2jtojp2a1,p2jtojp2a0,
     . p4jm2tojp2a4,p4jm2tojp2a3,p4jm2tojp2a2,
     . p4jm2tojp2a1,p4jm2tojp2a0,
     . p1jm1toja1,p1jm1toja0,
     . p2jm1tojp1a2,p2jm1tojp1a1,p2jm1tojp1a0,
     . p1jtojp1a1,p1jtojp1a0,
     . p1Wa1,p1Wa0,
     . P1a11w1,P2a22w2,P3a33w3,P3a32w3Corr,P2a21w2Corr,P3a31w3Corr,
     . w1,w2,w3,
     . R3a,R3b,R3c,R3d 
      
      do i=nGhstCells + 1,cellavgbc_len - nGhstCells 
         
         
         call P3jtojp3ABCD(cellavgbc_q(i),cellavgbc_q(i+1),
     .    cellavgbc_q(i+2),cellavgbc_q(i+3),dx,p3jtojp3a3,
     .    p3jtojp3a2,p3jtojp3a1,p3jtojp3a0)  
         
         call P6jm3tojp3ABCDEFG(cellavgbc_q(i-3),cellavgbc_q(i-2),
     .   cellavgbc_q(i-1),cellavgbc_q(i),cellavgbc_q(i+1),
     .   cellavgbc_q(i+2),cellavgbc_q(i+3),dx,
     .   p6jm3tojp3a6,p6jm3tojp3a5,p6jm3tojp3a4,
     .   p6jm3tojp3a3,p6jm3tojp3a2,p6jm3tojp3a1,p6jm3tojp3a0)        

         call P3jm3tojABCD(cellavgbc_q(i-3),cellavgbc_q(i-2),
     .    cellavgbc_q(i-1),cellavgbc_q(i),dx,p3jm3toja3,
     .    p3jm3toja2,p3jm3toja1,p3jm3toja0)  
     
         call P2jtojp2ABC(cellavgbc_q(i),cellavgbc_q(i+1),
     .    cellavgbc_q(i+2),dx,p2jtojp2a2,p2jtojp2a1,p2jtojp2a0)
     
         call P4jm2tojp2ABCDE(cellavgbc_q(i-2),cellavgbc_q(i-1),
     .   cellavgbc_q(i),cellavgbc_q(i+1),cellavgbc_q(i+2),dx,
     .   p4jm2tojp2a4,p4jm2tojp2a3,p4jm2tojp2a2,p4jm2tojp2a1,
     .   p4jm2tojp2a0)  
     
         call P2jm2tojABC(cellavgbc_q(i-2),cellavgbc_q(i-1),
     .    cellavgbc_q(i),dx,p2jm2toja2,p2jm2toja1,p2jm2toja0)   
     
         call P1jtojp1AB(cellavgbc_q(i),cellavgbc_q(i+1),dx,
     .    p1jtojp1a1,p1jtojp1a0)  
         call P1jm1tojAB(cellavgbc_q(i-1),cellavgbc_q(i),dx,
     .    p1jm1toja1,p1jm1toja0)  
         call P2jm1tojp1ABC(cellavgbc_q(i-1),cellavgbc_q(i),
     .    cellavgbc_q(i+1),dx,p2jm1tojp1a2,p2jm1tojp1a1,p2jm1tojp1a0) 

     
         call minmod(p1jtojp1a1,p2jm1tojp1a1,p1jtojp1a1,P1a11w1)
         call minmod(p2jtojp2a2,p4jm2tojp2a2,p2jm2toja2,P2a22w2)
         call minmod(p3jtojp3a3,p6jm3tojp3a3,p3jm3toja3,P3a33w3)

         if (dabs(p6jm3tojp3a3) > 10.0**(-10)) then
             w3 = P3a33w3  / p6jm3tojp3a3
         else
            w3 = 0
         end if

         P3a32w3Corr = w3*(p6jm3tojp3a2 - P2a22w2) + P2a22w2

         if (dabs(p4jm2tojp2a2) > 10.0**(-10)) then
            w2 = P2a22w2  / p4jm2tojp2a2
         else
            w2 = 0
         end if         
         
         P2a21w2Corr = w2*(p4jm2tojp2a1 - P1a11w1) + P1a11w1
         
         P3a31w3Corr = w3*(p6jm3tojp3a1 - P2a21w2Corr) + P2a21w2Corr
 
c         R3a = p3jm1tojp2a3
c         R3b = p3jm1tojp2a2
c         R3c = p3jm1tojp2a1
c         R3d = p3jm1tojp2a0
       
         R3a = P3a33w3
         R3b = P3a32w3Corr
         R3c = P3a31w3Corr
         R3d = cellavgbc_q(i)  - R3b/3*(dx/2)**2 
         
         cubicbc_q(4*i -3) = R3a*(-1.0*dx/2)**3 + R3b*(-1.0*dx/2)**2 +
     .   R3c*(-1.0*dx/2) + R3d
         cubicbc_q(4*i -2) = R3a*(-1.0*dx/6)**3 + R3b*(-1.0*dx/6)**2 
     .   + R3c*(-1.0*dx/6) + R3d
         cubicbc_q(4*i -1 ) = R3a*(dx/6)**3 + R3b*(dx/6)**2 
     .   + R3c*(dx/6) + R3d
         cubicbc_q(4*i) = R3a*(dx/2)**3 + R3b*(dx/2)**2 
     .   + R3c*(dx/2) + R3d
        
     

      end do
     
      

      end
     


c  ********************************************************************************
c  Part 2 : Finite Element Method
c ********************************************************************************

      subroutine FEMforU(hcub_bc,Gcub_bc,uNcub_bc,beta1,dx,xbc_len,
     .   cubbc_len,ncubbc_len,nGhstCell)
      
      !KL = 3
      !KU = 3
      INTEGER cubbc_len,ncubbc_len,nGhstCell,xbc_len
      DOUBLE PRECISION hcub_bc(cubbc_len),uNcub_bc(cubbc_len),
     . Gcub_bc(cubbc_len)
      DOUBLE PRECISION Gelem(4),uhelem(4,4), h3uxelem(4,4)   
      
      DOUBLE PRECISION A(2*3 + 3 + 1,ncubbc_len) , B(ncubbc_len)
      INTEGER IPIV(ncubbc_len),INFO,i,j
      
      if (dabs(beta1) .LE. 10.0**(-10)) then
          do i = 1, cubbc_len
            uNcub_bc(i) = Gcub_bc(i) / hcub_bc(i) 
          end do
      else
      
      ! Zero out
      do i = 1, ncubbc_len
       B(i) = 0.0
       A(4,i) =0.0
       A(5,i) =0.0
       A(6,i) =0.0
       A(7,i) =0.0
       A(8,i) =0.0
       A(9,i) =0.0
       A(10,i) =0.0
      end do
     
      
      !Interior
      do i = nGhstCell + 1,xbc_len - nGhstCell
      
        call  FEMelem(4,Gelem,uhelem,h3uxelem,beta1,dx,
     . hcub_bc(4*i -3),hcub_bc(4*i -2),hcub_bc(4*i -1),hcub_bc(4*i),
     . Gcub_bc(4*i -3),Gcub_bc(4*i -2),Gcub_bc(4*i -1),Gcub_bc(4*i))
              
         B(3*i - 2) = B(3*i - 2) +  Gelem(1)
         B(3*i - 1) = B(3*i - 1) +  Gelem(2)
         B(3*i ) = B(3*i ) +  Gelem(3)
         B(3*i + 1) = B(3*i + 1 ) +  Gelem(4)
 
         !Sup Diag + 2
         A(4,3*i +1) =  A(4,3*i +1) + uhelem(1,4) + h3uxelem(1,4)
         
         !Sup Diag + 2
         A(5,3*i) = A(5,3*i ) +  uhelem(1,3) + h3uxelem(1,3)
         A(5,3*i+1) = A(5,3*i+1) +  uhelem(2,4) + h3uxelem(2,4)

         !Sup Diag + 1
         A(6,3*i-1) = A(6,3*i - 1) +  uhelem(1,2) + h3uxelem(1,2)
         A(6,3*i) = A(6,3*i ) +  uhelem(2,3) + h3uxelem(2,3)
         A(6,3*i+1) = A(6,3*i+1) + uhelem(3,4) + h3uxelem(3,4)
         
         !Diag
         A(7,3*i - 2) = A(7,3*i - 2) +  uhelem(1,1) + h3uxelem(1,1)
         A(7,3*i - 1) = A(7,3*i - 1) +  uhelem(2,2)+ h3uxelem(2,2)
         A(7,3*i) = A(7,3*i) + uhelem(3,3) + h3uxelem(3,3)
         A(7,3*i +1 ) =  A(7,3*i +1) + uhelem(4,4) + h3uxelem(4,4)
  
         !Sub diag
         A(8,3*i - 2) = A(8,3*i - 2) +  uhelem(2,1) + h3uxelem(2,1)
         A(8,3*i -1) = A(8,3*i-1) + uhelem(3,2) + h3uxelem(3,2)
         A(8,3*i ) =  A(8,3*i) + uhelem(4,3) + h3uxelem(4,3)       

         !Sub diag
         A(9,3*i -2) = A(9,3*i -2) + uhelem(3,1) + h3uxelem(3,1)
         A(9,3*i  - 1) =  A(9,3*i - 1 ) + uhelem(4,2) + h3uxelem(4,2)    
         
         !Sub diag
         A(10,3*i -2 ) =  A(10,3*i -2) + uhelem(4,1) + h3uxelem(4,1) 
         
      end do
      

      !Edges
      do i = 1 , nGhstCell
         B(3*i - 2) = uNcub_bc(4*i - 3)
         B(3*i - 1) = uNcub_bc(4*i - 2)
         B(3*i) = uNcub_bc(4*i -1 )
         B(3*i+1) = uNcub_bc(4*i )
 
 
         A(4,3*(i+1) + 1 ) = 0.0

         A(5,3*(i+1) ) = 0.0
         
         A(6,3*(i+1) -1 ) = 0.0
                  
                  
         A(7,3*i - 2) = 1.0
         A(7,3*i - 1) = 1.0
         A(7,3*i ) = 1.0
         A(7,3*i + 1 ) = 1.0

         
         j = (xbc_len + 1 - i)
         B(3*j - 2) = uNcub_bc(4*j - 3)
         B(3*j - 1) = uNcub_bc(4*j -2)
         B(3*j ) = uNcub_bc(4*j-1)
         B(3*j + 1 ) = uNcub_bc(4*j)
         

         A(7,3*j -2 ) = 1.0         
         A(7,3*j -1 ) = 1.0
         A(7,3*j) = 1.0
         A(7,3*j + 1 ) = 1.0

         A(8,3*(j-1) ) = 0.0    

         A(9,3*(j-1) - 1) = 0.0

         A(10,3*(j-1) - 2) = 0.0
         
         
      end do
      

 
      call dgbsv(ncubbc_len,3,3,1,A,10,IPIV, B,ncubbc_len,INFO)
     
      do i = 1, xbc_len
       uNcub_bc(4*i-3) = B(3*i -2) 
       uNcub_bc(4*i-2) = B(3*i -1) 
       uNcub_bc(4*i-1) = B(3*i) 
       uNcub_bc(4*i) = B(3*i + 1) 
      end do
      
      end if

      
      end


      subroutine FEMelem(elem_len,Gelem,uhelem,h3uxelem,beta1,dx,
     . hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph)
      
      INTEGER elem_len,i
      
      DOUBLE PRECISION dx,beta1,hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph
      
      DOUBLE PRECISION Gelem(elem_len),
     .   uhelem(elem_len,elem_len), h3uxelem(elem_len,elem_len)
          
          
      Gelem(1) = dx/2*(16.0/105*Gjmh+ 33.0/280*Gjms 
     .   - 3.0/70*Gjps + 19.0/840*Gjph )
      Gelem(2) = dx/2*(33.0/280*Gjmh + 27.0/35*Gjms
     . - 27.0/280*Gjps - 3.0/70*Gjph) 
      Gelem(3) = dx/2*(-3.0/70*Gjmh- 27.0/280*Gjms
     . + 27.0/35*Gjps + 33.0/280*Gjph)
      Gelem(4) = dx/2*(19.0/840*Gjmh- 3.0/70*Gjms
     . + 33.0/280*Gjps + 16.0/105*Gjph)
      
      uhelem(1,1) = dx/2*(17*hjmh/160 + 9*hjms/140 
     . + hjph/168 - 27*hjps/1120)
      uhelem(1,2) = dx/2*(9*hjmh/140 + 27*hjms/280 
     . + 3*hjph/560 - 27*hjps/560)
      uhelem(1,3) = dx/2*(-27*hjmh/1120 - 27*hjms/560 
     . + 3*hjph/560 + 27*hjps/1120)
      uhelem(1,4) = dx/2*(hjmh/168 + 3*hjms/560 
     . + hjph/168 + 3*hjps/560)
      uhelem(2,1) = dx/2*(9*hjmh/140 + 27*hjms/280 
     . + 3*hjph/560 - 27*hjps/560)
      uhelem(2,2) = dx/2*(27*hjmh/280 + 729*hjms/1120 + 27*hjph/1120)
      uhelem(2,3) = dx/2*(-27*hjmh/560 - 27*hjph/560)
      uhelem(2,4) = dx/2*(3*hjmh/560 + 27*hjms/1120 
     . - 27*hjph/1120 - 27*hjps/560)
      uhelem(3,1) = dx/2*(-27*hjmh/1120 - 27*hjms/560 
     . + 3*hjph/560 + 27*hjps/1120)
      uhelem(3,2) = dx/2*(-27*hjmh/560 - 27*hjph/560)
      uhelem(3,3) = dx/2*(27*hjmh/1120 + 27*hjph/280 
     . + 729*hjps/1120)
      uhelem(3,4) = dx/2*(3*hjmh/560 - 27*hjms/560 
     . + 9*hjph/140 + 27*hjps/280)
      uhelem(4,1) = dx/2*(hjmh/168 + 3*hjms/560 
     . + hjph/168 + 3*hjps/560)
      uhelem(4,2) = dx/2*(3*hjmh/560 + 27*hjms/1120 
     . - 27*hjph/1120 - 27*hjps/560)
      uhelem(4,3) = dx/2*(3*hjmh/560 - 27*hjms/560 
     . + 9*hjph/140 + 27*hjps/280)
      uhelem(4,4) = dx/2*(hjmh/168 - 27*hjms/1120 
     . + 17*hjph/160 + 9*hjps/140)


      h3uxelem(1,1) = (beta1/2.0)*2/dx*(5377*hjmh**3/8960 
     . + 490239*hjmh**2*hjms/640640 + 45223*hjmh**2*hjph/640640 
     . - 7695*hjmh**2*hjps/23296 + 19035*hjmh*hjms**2/23296 
     . + 2601*hjmh*hjms*hjph/20020 - 100521*hjmh*hjms*hjps/160160 
     . + 1163*hjmh*hjph**2/183040 - 8109*hjmh*hjph*hjps/160160 
     . + 14499*hjmh*hjps**2/116480 + 235467*hjms**3/366080 
     . + 297837*hjms**2*hjph/2562560 - 728271*hjms**2*hjps/1281280 
     . + 9963*hjms*hjph**2/2562560 - 1377*hjms*hjph*hjps/16016 
     . + 124659*hjms*hjps**2/640640 + 953*hjph**3/73216 
     . + 8397*hjph**2*hjps/1281280 + 13527*hjph*hjps**2/640640 
     . + 729*hjps**3/1281280)
     
      h3uxelem(1,2) = (beta1/2.0)*2/dx*(-1164861*hjmh**3/1281280 
     . - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 
     . + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 
     . - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 
     . - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 
     . - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 
     . - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 
     . + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 
     . - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 
     . - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 
     . - 150903*hjps**3/1281280)
     
     
      h3uxelem(1,3) = (beta1/2.0)*2/dx*(503469*hjmh**3/1281280 
     . + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 
     . - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 
     . + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 
     . + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 
     . + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 
     . + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 
     . - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 
     . + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 
     . + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 
     . + 6561*hjps**3/116480)
     
     
      h3uxelem(1,4) = (beta1/2.0)*2/dx*(-107519*hjmh**3/1281280 
     . - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 
     . + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 
     . - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 
     . - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 
     . - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 
     . - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 
     . + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 
     . - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 
     . - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 
     . + 78003*hjps**3/1281280)
     
      h3uxelem(2,1) = (beta1/2.0)*2/dx*(-1164861*hjmh**3/1281280 
     . - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 
     . + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 
     . - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 
     . - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 
     . - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 
     . - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 
     . + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 
     . - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 
     . - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 
     . - 150903*hjps**3/1281280)
     
      h3uxelem(2,2) = (beta1/2.0)*2/dx*(323217*hjmh**3/232960 
     . + 1858221*hjmh**2*hjms/1281280 + 13689*hjmh**2*hjph/98560 
     . - 1630773*hjmh**2*hjps/2562560 + 19683*hjmh*hjms**2/18304 
     . + 729*hjmh*hjms*hjph/3640 - 6561*hjmh*hjms*hjps/6160 
     . + 18873*hjmh*hjph**2/640640 - 729*hjmh*hjph*hjps/16016 
     . + 400221*hjmh*hjps**2/2562560 + 85293*hjms**3/98560 
     . + 111537*hjms**2*hjph/1281280 - 19683*hjms**2*hjps/640640 
     . - 21141*hjms*hjph**2/256256 - 6561*hjms*hjph*hjps/20020 
     . + 846369*hjms*hjps**2/1281280 + 50409*hjph**3/183040 
     . + 140697*hjph**2*hjps/640640 + 19683*hjph*hjps**2/183040 
     . + 478953*hjps**3/512512)
     
      h3uxelem(2,3) = (beta1/2.0)*2/dx*(-111429*hjmh**3/183040 
     . - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 
     . + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 
     . - 51759*hjmh*hjms*hjph/1281280 
     . + 465831*hjmh*hjms*hjps/1281280 - 27135*hjmh*hjph**2/512512 
     . - 51759*hjmh*hjph*hjps/1281280 - 6561*hjmh*hjps**2/1281280 
     . - 662661*hjms**3/1281280 - 6561*hjms**2*hjph/1281280 
     . - 1318761*hjms**2*hjps/2562560 + 292329*hjms*hjph**2/1281280 
     . + 465831*hjms*hjph*hjps/1281280 - 1318761*hjms*hjps**2/2562560 
     . - 111429*hjph**3/183040 - 1324593*hjph**2*hjps/2562560 
     . - 98415*hjph*hjps**2/512512 - 662661*hjps**3/1281280)
     
      h3uxelem(2,4) = (beta1/2.0)*2/dx*(47763*hjmh**3/366080 
     . + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 
     . - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 
     . - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 
     . + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 
     . + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 
     . + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 
     . - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 
     . + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 
     . + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 
     . - 59049*hjps**3/197120)
     
      h3uxelem(3,1) = (beta1/2.0)*2/dx*(503469*hjmh**3/1281280 
     . + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 
     . - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 
     . + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 
     . + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 
     . + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 
     . + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 
     . - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 
     . + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 
     . + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 
     . + 6561*hjps**3/116480)
     
      h3uxelem(3,2) = (beta1/2.0)*2/dx*(-111429*hjmh**3/183040 
     . - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 
     . + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 
     . - 51759*hjmh*hjms*hjph/1281280 + 465831*hjmh*hjms*hjps/1281280 
     . - 27135*hjmh*hjph**2/512512 - 51759*hjmh*hjph*hjps/1281280 
     . - 6561*hjmh*hjps**2/1281280 - 662661*hjms**3/1281280 
     . - 6561*hjms**2*hjph/1281280 - 1318761*hjms**2*hjps/2562560 
     . + 292329*hjms*hjph**2/1281280 + 465831*hjms*hjph*hjps/1281280 
     . - 1318761*hjms*hjps**2/2562560 - 111429*hjph**3/183040 
     . - 1324593*hjph**2*hjps/2562560 - 98415*hjph*hjps**2/512512 
     . - 662661*hjps**3/1281280)
     
      h3uxelem(3,3) = (beta1/2.0)*2/dx*(50409*hjmh**3/183040 
     . + 140697*hjmh**2*hjms/640640 + 18873*hjmh**2*hjph/640640 
     . - 21141*hjmh**2*hjps/256256 + 19683*hjmh*hjms**2/183040 
     . - 729*hjmh*hjms*hjph/16016 - 6561*hjmh*hjms*hjps/20020 
     . + 13689*hjmh*hjph**2/98560 + 729*hjmh*hjph*hjps/3640 
     . + 111537*hjmh*hjps**2/1281280 + 478953*hjms**3/512512 
     . + 400221*hjms**2*hjph/2562560 + 846369*hjms**2*hjps/1281280 
     . - 1630773*hjms*hjph**2/2562560 - 6561*hjms*hjph*hjps/6160 
     . - 19683*hjms*hjps**2/640640 + 323217*hjph**3/232960 
     . + 1858221*hjph**2*hjps/1281280 + 19683*hjph*hjps**2/18304 
     . + 85293*hjps**3/98560)
     
      h3uxelem(3,4) = (beta1/2.0)*2/dx*(-6939*hjmh**3/116480 
     . - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 
     . + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 
     . + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 
     . - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 
     . - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 
     . - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 
     . + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 
     . + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280
     . - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 
     . - 518319*hjps**3/1281280)
     
      h3uxelem(4,1) = (beta1/2.0)*2/dx*(-107519*hjmh**3/1281280 
     . - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 
     . + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 
     . - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 
     . - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 
     . - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 
     . - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 
     . + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 
     . - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 
     . - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 
     . + 78003*hjps**3/1281280)
     
      h3uxelem(4,2) = (beta1/2.0)*2/dx*(47763*hjmh**3/366080 
     . + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 
     . - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 
     . - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 
     . + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 
     . + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 
     . + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 
     . - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 
     . + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 
     . + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 
     . - 59049*hjps**3/197120)
     
      h3uxelem(4,3) = (beta1/2.0)*2/dx*(-6939*hjmh**3/116480 
     . - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 
     . + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 
     . + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 
     . - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 
     . - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 
     . - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 
     . + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 
     . + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280 
     . - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 
     . - 518319*hjps**3/1281280)
     
     
      h3uxelem(4,4) = (beta1/2.0)*2/dx*(953*hjmh**3/73216 +
     . 8397*hjmh**2*hjms/1281280 + 1163*hjmh**2*hjph/183040 
     . + 9963*hjmh**2*hjps/2562560 + 13527*hjmh*hjms**2/640640 
     . - 8109*hjmh*hjms*hjph/160160 - 1377*hjmh*hjms*hjps/16016 
     . + 45223*hjmh*hjph**2/640640 + 2601*hjmh*hjph*hjps/20020 
     . + 297837*hjmh*hjps**2/2562560 + 729*hjms**3/1281280 
     . + 14499*hjms**2*hjph/116480 + 124659*hjms**2*hjps/640640 
     . - 7695*hjms*hjph**2/23296 - 100521*hjms*hjph*hjps/160160 
     . - 728271*hjms*hjps**2/1281280 + 5377*hjph**3/8960 
     . + 490239*hjph**2*hjps/640640 + 19035*hjph*hjps**2/23296 
     . + 235467*hjps**3/366080)   

 
     
      
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
      subroutine EvolveStepWrap(xbc_len,cubbc_len,habc,Gabc,
     . ga,beta1,beta2,calctol,nGhstCells,dx,dt, hcub,Gcub,ucub)
          
      integer nGhstCells,xbc_len,cubbc_len
      DOUBLE PRECISION habc(xbc_len),Gabc(xbc_len),uabc(xbc_len)
      DOUBLE PRECISION ga,beta1,beta2,theta,dx,dt
     
      !local variables
      DOUBLE PRECISION hapbc(xbc_len),Gapbc(xbc_len),
     . happbc(xbc_len),Gappbc(xbc_len)

      DOUBLE PRECISION hcub(cubbc_len),Gcub(cubbc_len),
     . ucub(cubbc_len)

     
      integer i,j
     
      !Reconstruct h,G   
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,habc,dx,hcub,calctol)
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,Gabc,dx,Gcub,calctol)
     
      !Calculate u
      call FEMforU(hcub,Gcub,ucub,beta1,dx,xbc_len,
     .   cubbc_len, 3*xbc_len + 1,nGhstCells)
     
      !Euler Step
      call SingleEulerStep(xbc_len,cubbc_len,hcub,Gcub,ucub,ga,
     . beta1,beta2,nGhstCells,dt,dx,habc,Gabc,hapbc,Gapbc)
     
      ! BC
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . hcub(1:4*nGhstCells),
     .   dx,hapbc(1:nGhstCells))
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . Gcub(1:4*nGhstCells),
     .   dx,Gapbc(1:nGhstCells)) 
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . hcub(cubbc_len - 4*nGhstCells + 1 : cubbc_len),
     .   dx,hapbc(xbc_len - nGhstCells + 1:xbc_len))
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . Gcub(cubbc_len - 4*nGhstCells + 1 : cubbc_len),
     .   dx,Gapbc(xbc_len - nGhstCells + 1:xbc_len))   
     
     
c ====== second time step 

      !Reconstruct h,G   
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,hapbc,dx,hcub,calctol)
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,Gapbc,dx,Gcub,calctol)   


      !Calculate u
      call FEMforU(hcub,Gcub,ucub,beta1,dx,xbc_len,
     .   cubbc_len, 3*xbc_len + 1,nGhstCells)
     
      !Euler Step
      call SingleEulerStep(xbc_len,cubbc_len,hcub,Gcub,ucub,ga,
     . beta1,beta2,nGhstCells,dt,dx,hapbc,Gapbc,happbc,Gappbc)
                  
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
      do i= 1,xbc_len
         habc(i) =  ( habc(i ) + happbc(i))/2d0
         Gabc(i) = ( Gabc(i ) + Gappbc(i))/2d0
      end do
      
      ! BC
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . hcub(1:4*nGhstCells),
     .   dx,habc(1:nGhstCells))
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . Gcub(1:4*nGhstCells),
     .   dx,Gabc(1:nGhstCells)) 
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . hcub(cubbc_len - 4*nGhstCells + 1 : cubbc_len),
     .   dx,habc(xbc_len - nGhstCells + 1:xbc_len))
      call CubicToCellAverage(4*nGhstCells,nGhstCells,
     . Gcub(cubbc_len - 4*nGhstCells + 1 : cubbc_len),
     .   dx,Gabc(xbc_len - nGhstCells + 1:xbc_len))   
          
      
      end 
    

c ====
c SingleEulerStep
c produces new cell averages of h,G using forward Euler step,
c with flux approximated using Kurganov's method
c ====
      subroutine SingleEulerStep(xbc_len,cubbc_len,hcub,Gcub,ucub,ga,
     . beta1,beta2,nGhstCells,dt,dx,habc,Gabc,hapbc,Gapbc)
     
     
      integer nGhstCells,xbc_len,cubbc_len
      DOUBLE PRECISION dt,dx,ga,beta1,beta2
      DOUBLE PRECISION hcub(cubbc_len),Gcub(cubbc_len),ucub(cubbc_len)
     . ,hapbc(xbc_len),Gapbc(xbc_len),habc(xbc_len),Gabc(xbc_len)
     
     
      DOUBLE PRECISION cdhi,cdGi, fih,fiG,foh,foG
      
      integer i,ileft,iright

      !Now we update interior, first calculating flux across left interior boundary 
      !then loop over cells to get flux in/out and thus updating cell average values
      
      !Do left boundary
      i = nGhstCells
            
      !calculates foh,foG which is flux across x_{i +1/2}
      ! it also updates cdhi,cdGi to be gradient across cell i + 1
      call Fluxxiph(hcub,Gcub,ucub,cubbc_len,i,ga,beta1,beta2,
     .   dx,foh,foG)
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      do i = nGhstCells + 1, xbc_len - nGhstCells
            
         !calculates foh,foG which is flux across x_{i +1/2}
         ! it also updates cdhi,cdGi to be gradient across cell i + 1
         call Fluxxiph(hcub,Gcub,ucub,cubbc_len,i,ga,beta1,beta2,
     .   dx,foh,foG)
             
         hapbc(i) = habc(i) - dt*(foh - fih)/dx 
     
         Gapbc(i) = Gabc(i) - dt*(foG - fiG)/dx 

         !flux out becomes flux in on next cell
         fih = foh
         fiG = foG
         
      end do 
           
      end
      
c subroutine that cubics, calculates flux across 
      subroutine Fluxxiph(hcub,Gcub,ucub,cubbc_len,i,ga,beta1,beta2,
     .   dx,foh,foG)
     
       INTEGER i,cubbc_len
       DOUBLE PRECISION hcub(cubbc_len),Gcub(cubbc_len),ucub(cubbc_len),
     . ga,beta1,beta2,dx,foh,foG

     
       DOUBLE PRECISION cdGip1,felG,felh,ferG,ferh,sr,sl,isrmsl,
     . hir,Gir,uir,duir,hip1l,Gip1l,uip1l,duip1l,
     . hjmh,hjms,hjps,hjph,
     . hjp1mh,hjp1ms,hjp1ps,hjp1ph,
     . dhir,ddhir,dhip1l,ddhip1l, alpha,
     . pua,pub,puc,pud,pha,phb,phc,phd,
     . phTa,phTb,phTc,phTd,phTe,
     . hT,hp1T,
     . pup1a,pup1b,pup1c,pup1d,php1a,php1b,php1c,php1d
     . phTp1a,phTp1b,phTp1c,phTp1d,phTp1e              

      uir = ucub(4*i)
      hir = hcub(4*i)
      Gir = Gcub(4*i)
      
      pua = (-9*ucub(4*i-3) + 27*ucub(4*i-2) - 27*ucub(4*i-1) 
     . + 9*ucub(4*i) )/ (2*dx**3)
      pub = (9*ucub(4*i-3) - 9*ucub(4*i-2) - 9*ucub(4*i-1) 
     . + 9*ucub(4*i)  )/ (4*dx**2)
      puc = (ucub(4*i-3) - 27*ucub(4*i-2) + 27*ucub(4*i -1) 
     . - ucub(4*i) )/ (8*dx)
      pud = (-ucub(4*i-3) + 9*ucub(4*i-2) + 9*ucub(4*i-1) 
     . - ucub(4*i))/ 16.0

      pha = (-9*hcub(4*i-3) + 27*hcub(4*i-2) - 27*hcub(4*i-1) 
     . + 9*hcub(4*i) )/ (2*dx**3)
      phb = (9*hcub(4*i-3) - 9*hcub(4*i-2) - 9*hcub(4*i-1) 
     . + 9*hcub(4*i)  )/ (4*dx**2)
      phc = (hcub(4*i-3) - 27*hcub(4*i-2) + 27*hcub(4*i -1) 
     . - hcub(4*i) )/ (8*dx)
      phd = (-hcub(4*i-3) + 9*hcub(4*i-2) + 9*hcub(4*i-1) 
     . - hcub(4*i))/ 16.0

c      hjph = 0.5*(hcub(4*i) +hcub(4*i+1)) 
c      pha = (-9*hcub(4*i-3) + 27*hcub(4*i-2) - 27*hcub(4*i-1) 
c     . + 9*hjph )/ (2*dx**3)
c      phb = (9*hcub(4*i-3) - 9*hcub(4*i-2) - 9*hcub(4*i-1) 
c     . + 9*hjph  )/ (4*dx**2)
c      phc = (hcub(4*i-3) - 27*hcub(4*i-2) + 27*hcub(4*i -1) 
c     . - hjph )/ (8*dx)
c      phd = (-hcub(4*i-3) + 9*hcub(4*i-2) + 9*hcub(4*i-1) 
c     . - hjph)/ 16.0

          
      duir = 3*pua*(dx/2)**2 + 2*pub*(dx/2) + puc
      dhir = 3*pha*(dx/2)**2 + 2*phb*(dx/2) + phc
      ddhir = 6*pha*(dx/2) + 2*phb
      
      hjmh = 0.5*(hcub(4*(i-1)) +hcub(4*(i-1)+1)) 
      hjms = hcub(4*i-2)
      hjps = hcub(4*i-1)
      hjph = 0.5*(hcub(4*i) +hcub(4*i+1)) 
      
      phTa = 1701*hjmh**2/(4*dx**6) - 5103*hjmh*hjms/(2*dx**6) 
     . - 1701*hjmh*hjph/(2*dx**6) + 5103*hjmh*hjps/(2*dx**6) 
     . + 15309*hjms**2/(4*dx**6) + 5103*hjms*hjph/(2*dx**6) 
     . - 15309*hjms*hjps/(2*dx**6) + 1701*hjph**2/(4*dx**6) 
     . - 5103*hjph*hjps/(2*dx**6) + 15309*hjps**2/(4*dx**6)
     
      phTb = -567*hjmh**2/(2*dx**5) + 1134*hjmh*hjms/dx**5 
     . - 567*hjmh*hjps/dx**5 - 1701*hjms**2/(2*dx**5) 
     . + 567*hjms*hjph/dx**5 + 567*hjph**2/(2*dx**5) 
     . - 1134*hjph*hjps/dx**5 + 1701*hjps**2/(2*dx**5)
    
      phTc = 243*hjmh**2/(8*dx**4) + 891*hjmh*hjms/(4*dx**4) 
     . + 405*hjmh*hjph/(4*dx**4) - 1539*hjmh*hjps/(4*dx**4) 
     . - 6237*hjms**2/(8*dx**4) - 1539*hjms*hjph/(4*dx**4) 
     . + 6885*hjms*hjps/(4*dx**4) + 243*hjph**2/(8*dx**4) 
     . + 891*hjph*hjps/(4*dx**4) - 6237*hjps**2/(8*dx**4)
     
      phTd = 45*hjmh**2/(8*dx**3) - 207*hjmh*hjms/(2*dx**3) 
     . + 153*hjmh*hjps/(4*dx**3) + 1215*hjms**2/(8*dx**3) 
     . - 153*hjms*hjph/(4*dx**3) - 45*hjph**2/(8*dx**3) 
     . + 207*hjph*hjps/(2*dx**3) - 1215*hjps**2/(8*dx**3)
     
      phTe = -35*hjmh**2/(64*dx**2) + 153*hjmh*hjms/(32*dx**2) 
     . - 37*hjmh*hjph/(32*dx**2) + 207*hjmh*hjps/(32*dx**2) 
     . + 405*hjms**2/(64*dx**2) + 207*hjms*hjph/(32*dx**2) 
     . - 1053*hjms*hjps/(32*dx**2) - 35*hjph**2/(64*dx**2) 
     . + 153*hjph*hjps/(32*dx**2) + 405*hjps**2/(64*dx**2)
     
      hT = phTa*(dx/2)**4 + phTb*(dx/2)**3 + phTc*(dx/2)**2 
     . + phTd*(dx/2) + phTe
   
      uip1l = ucub(4*i + 1)
      hip1l = hcub(4*i + 1)
      Gip1l = Gcub(4*i + 1)
      
      pup1a = (-9*ucub(4*(i+1)-3) + 27*ucub(4*(i+1)-2) 
     . - 27*ucub(4*(i+1)-1)  + 9*ucub(4*(i+1)) )/ (2*dx**3)
      pup1b = (9*ucub(4*(i+1)-3) - 9*ucub(4*(i+1)-2) 
     . - 9*ucub(4*(i+1)-1) + 9*ucub(4*(i+1))  )/ (4*dx**2)
      pup1c = (ucub(4*(i+1)-3) - 27*ucub(4*(i+1)-2) 
     . + 27*ucub(4*(i+1) -1)  - ucub(4*(i+1)) )/ (8*dx)
      pup1d = (-ucub(4*(i+1)-3) + 9*ucub(4*(i+1)-2) 
     . + 9*ucub(4*(i+1)-1)  - ucub(4*(i+1)))/ 16.0

      php1a = (-9*hcub(4*(i+1)-3) + 27*hcub(4*(i+1)-2) 
     . - 27*hcub(4*(i+1)-1)  + 9*hcub(4*(i+1)) )/ (2*dx**3)
      php1b = (9*hcub(4*(i+1)-3) - 9*hcub(4*(i+1)-2) 
     . - 9*hcub(4*(i+1)-1)  + 9*hcub(4*(i+1))  )/ (4*dx**2)
      php1c = (hcub(4*(i+1)-3) - 27*hcub(4*(i+1)-2) 
     . + 27*hcub(4*(i+1) -1)  - hcub(4*(i+1)) )/ (8*dx)
      php1d = (-hcub(4*(i+1)-3) + 9*hcub(4*(i+1)-2) 
     . + 9*hcub(4*(i+1)-1)  - hcub(4*(i+1)))/ 16.0


c      hjp1mh = 0.5*(hcub(4*i) +hcub(4*i+1)) 
c      php1a = (-9*hjp1mh + 27*hcub(4*(i+1)-2) 
c     . - 27*hcub(4*(i+1)-1)  + 9*hcub(4*(i+1)) )/ (2*dx**3)
c      php1b = (9*hjp1mh - 9*hcub(4*(i+1)-2) 
c     . - 9*hcub(4*(i+1)-1)  + 9*hcub(4*(i+1))  )/ (4*dx**2)
c      php1c = (hjp1mh - 27*hcub(4*(i+1)-2) 
c     . + 27*hcub(4*(i+1) -1)  - hcub(4*(i+1)) )/ (8*dx)
c     php1d = (-hjp1mh + 9*hcub(4*(i+1)-2) 
c    . + 9*hcub(4*(i+1)-1)  - hcub(4*(i+1)))/ 16.0
     
      duip1l = 3*pup1a*(-dx/2)**2 + 2*pup1b*(-dx/2) + pup1c
      dhip1l = 3*php1a*(-dx/2)**2 + 2*php1b*(-dx/2) + php1c
      ddhip1l = 6*php1a*(-dx/2) + 2*php1b


      hjp1mh = 0.5*(hcub(4*i) +hcub(4*i+1))
      hjp1ms = hcub(4*(i+1)-2)
      hjp1ps = hcub(4*(i+1)-1)
      hjp1ph = 0.5*(hcub(4*(i+1)) +hcub(4*(i+1)+1))
 
      
    
      phTc = 243*hjmh**2/(8*dx**4) + 891*hjmh*hjms/(4*dx**4) 
     . + 405*hjmh*hjph/(4*dx**4) - 1539*hjmh*hjps/(4*dx**4) 
     . - 6237*hjms**2/(8*dx**4) - 1539*hjms*hjph/(4*dx**4) 
     . + 6885*hjms*hjps/(4*dx**4) + 243*hjph**2/(8*dx**4) 
     . + 891*hjph*hjps/(4*dx**4) - 6237*hjps**2/(8*dx**4)
          
      
      php1Ta = 1701*hjp1mh**2/(4*dx**6) - 5103*hjp1mh*hjp1ms/(2*dx**6) 
     . - 1701*hjp1mh*hjp1ph/(2*dx**6) + 5103*hjp1mh*hjp1ps/(2*dx**6) 
     . + 15309*hjp1ms**2/(4*dx**6) + 5103*hjp1ms*hjp1ph/(2*dx**6) 
     . - 15309*hjp1ms*hjp1ps/(2*dx**6) + 1701*hjp1ph**2/(4*dx**6) 
     . - 5103*hjp1ph*hjp1ps/(2*dx**6) + 15309*hjp1ps**2/(4*dx**6)
    
      php1Tb = -567*hjp1mh**2/(2*dx**5) + 1134*hjp1mh*hjp1ms/dx**5 
     . - 567*hjp1mh*hjp1ps/dx**5 - 1701*hjp1ms**2/(2*dx**5) 
     . + 567*hjp1ms*hjp1ph/dx**5 + 567*hjp1ph**2/(2*dx**5) 
     . - 1134*hjp1ph*hjp1ps/dx**5 + 1701*hjp1ps**2/(2*dx**5)
     
      php1Tc = 243*hjp1mh**2/(8*dx**4) + 891*hjp1mh*hjp1ms/(4*dx**4) 
     . + 405*hjp1mh*hjp1ph/(4*dx**4) - 1539*hjp1mh*hjp1ps/(4*dx**4) 
     . - 6237*hjp1ms**2/(8*dx**4) - 1539*hjp1ms*hjp1ph/(4*dx**4) 
     . + 6885*hjp1ms*hjp1ps/(4*dx**4) + 243*hjp1ph**2/(8*dx**4) 
     . + 891*hjp1ph*hjp1ps/(4*dx**4) - 6237*hjp1ps**2/(8*dx**4)
     
      php1Td = 45*hjp1mh**2/(8*dx**3) - 207*hjp1mh*hjp1ms/(2*dx**3) 
     . + 153*hjp1mh*hjp1ps/(4*dx**3) + 1215*hjp1ms**2/(8*dx**3) 
     . - 153*hjp1ms*hjp1ph/(4*dx**3) - 45*hjp1ph**2/(8*dx**3) 
     . + 207*hjp1ph*hjp1ps/(2*dx**3) - 1215*hjp1ps**2/(8*dx**3)
     
      php1Te = -35*hjp1mh**2/(64*dx**2) + 153*hjp1mh*hjp1ms/(32*dx**2) 
     . - 37*hjp1mh*hjp1ph/(32*dx**2) + 207*hjp1mh*hjp1ps/(32*dx**2) 
     . + 405*hjp1ms**2/(64*dx**2) + 207*hjp1ms*hjp1ph/(32*dx**2) 
     . - 1053*hjp1ms*hjp1ps/(32*dx**2) - 35*hjp1ph**2/(64*dx**2) 
     . + 153*hjp1ph*hjp1ps/(32*dx**2) + 405*hjp1ps**2/(64*dx**2)
     
      hp1T = php1Ta*(-dx/2)**4 + php1Tb*(-dx/2)**3 + php1Tc*(-dx/2)**2 
     . + php1Td*(-dx/2) + php1Te
      
c      duip1l = 0.5*(duip1l + duir)
c      hp1T = 0.5*(hp1T + hT)
c      dhip1l = 0.5*(dhip1l + dhir)
c      ddhip1l = 0.5*(ddhip1l + ddhir)

c      duir = duip1l
c      hT = hp1T
c      dhir = dhip1l
c      ddhir = ddhip1l     
      
      if  (dabs(beta1) < 10d0**(-10))  then
         alpha = 1
      else
         alpha = max(1d0,beta2 / beta1)
      end if
      
      sl  = min(0d0, uir - dsqrt(alpha*ga*hir) ,
     . uip1l - dsqrt(alpha*ga*hip1l)  )
      sr  = max(0d0, uir + dsqrt(alpha*ga*hir) ,
     . uip1l + dsqrt(alpha*ga*hip1l)  )
      
      !left and right flux
      felh = uir*hir
      felG = uir*Gir + ga*(hir**2)/2d0 
     .      - beta1*hir**3*duir**2
     .      - beta2/4d0*ga*(hir**2)*hT
c     .      - beta2/4d0*ga*(hir**2)*(2*hir*ddhir + (dhir**2))
     
           
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - beta1*hip1l**3*duip1l**2
     .      - beta2/4d0*ga*(hip1l**2)*hp1T
c     .      - beta2/4d0*ga*(hip1l**2)*(2*hip1l*ddhip1l + (dhip1l**2))
 
c      print *, '====',i,'===='
c      print *,  felG,ferG
c      print *, beta1*hir**3*duir**2,beta1*hip1l**3*duip1l**2
c      print *,beta2/4d0*ga*(hir**2)*(2*hir*ddhir + (dhir**2)),
c     .  beta2/4d0*ga*(hip1l**2)*(2*hip1l*ddhip1l + (dhip1l**2))
c      print *, '====',i,'===='
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if 
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))

      end



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
      subroutine CubicInterp(qjmh,qjms,qjps,qjph,dx,CubicCoeff)
     
      DOUBLE PRECISION qjmh,qjms,qjps,qjph,dx
      DOUBLE PRECISION CubicCoeff(4)
      
      CubicCoeff(1) = (-9*qjmh + 27*qjms - 27*qjps  + 9*qjph )/ 
     . (2*dx**3)
      CubicCoeff(2) = (9*qjmh - 9*qjms - 9*qjps  + 9*qjph )/ (4*dx**2)
      CubicCoeff(3) = (qjmh - 27*qjms + 27*qjps  - qjph )/ (8*dx)
      CubicCoeff(4) = (-qjmh + 9*qjms + 9*qjps  - qjph )/ 16.0

      end


c ========================
c Functions to evaluate Quartic at x, with quartic centered around x_j
c xmxj = x - x_j 
c ====================      
      subroutine CubicCoeffEvalxj(CubicCoeff,xmxj,qatxj)
      
      DOUBLE PRECISION CubicCoeff(4)
      DOUBLE PRECISION xmxj,qatxj
      
      qatxj = CubicCoeff(1)*xmxj**3 + CubicCoeff(2)*xmxj**2 +
     . CubicCoeff(3)*xmxj + CubicCoeff(4)
     
      end
      
      subroutine CubicCoeffEvalGradxj(CubicCoeff,xmxj,dqatxj)
      
      DOUBLE PRECISION CubicCoeff(4)
      DOUBLE PRECISION xmxj,dqatxj
      
      dqatxj = 3*CubicCoeff(1)*xmxj**2 + 2*CubicCoeff(2)*xmxj +
     . CubicCoeff(3) 
     
      end

c =====
c Functions to get integrals over cell
c ====
      ! Energy function for cell
      subroutine AllEnergiesIntegralCell(xbc_len,cubbc_len,h,u,G,ga,
     . beta1,beta2,j,dx,CellEnergies)
      
      integer j,xbc_len,cubbc_len
      DOUBLE PRECISION h(cubbc_len),u(cubbc_len),G(cubbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION CellEnergies(4)
      
      integer i
      DOUBLE PRECISION fGPe(4),sGPe(4),tGPe(4)
      
      DOUBLE PRECISION GPmxj,hGP,GGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(4), uCoeff(4), GCoeff(4)
      
      call CubicInterp(h(4*j-3),h(4*j-2),h(4*j-1),h(4*j),dx,hCoeff)
      call CubicInterp(u(4*j-3),u(4*j-2),u(4*j-1),u(4*j),dx,uCoeff)
      call CubicInterp(G(4*j-3),G(4*j-2),G(4*j-1),G(4*j),dx,GCoeff)

      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2

      call CubicCoeffEvalxj(hCoeff,GPmxj,hGP)
      call CubicCoeffEvalxj(uCoeff,GPmxj,uGP)
      call CubicCoeffEvalxj(GCoeff,GPmxj,GGP)
      
      call CubicCoeffEvalGradxj(hCoeff,xmxj,hxGP)
      call CubicCoeffEvalGradxj(uCoeff,xmxj,uxGP)
      
      
      fGPe(1) = hGP
      fGPe(2) = GGP
      fGPe(3) = hGP*uGP
      fGPe(4) = (hGP*uGP**2 + 0.5*beta1*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !second gauss point
      GPmxj = 0.0 
      call CubicCoeffEvalxj(hCoeff,GPmxj,hGP)
      call CubicCoeffEvalxj(uCoeff,GPmxj,uGP)
      call CubicCoeffEvalxj(GCoeff,GPmxj,GGP)
      
      call CubicCoeffEvalGradxj(hCoeff,xmxj,hxGP)
      call CubicCoeffEvalGradxj(uCoeff,xmxj,uxGP)
      
      sGPe(1) = hGP
      sGPe(2) = GGP
      sGPe(3) = hGP*uGP
      sGPe(4) = (hGP*uGP**2 + 0.5*beta1*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !third gauss point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      call CubicCoeffEvalxj(hCoeff,GPmxj,hGP)
      call CubicCoeffEvalxj(uCoeff,GPmxj,uGP)
      call CubicCoeffEvalxj(GCoeff,GPmxj,GGP)
      
      call CubicCoeffEvalGradxj(hCoeff,xmxj,hxGP)
      call CubicCoeffEvalGradxj(uCoeff,xmxj,uxGP)
      
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + 0.5*beta1*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !weight the values at gauss points to get approximate integral over cell
      do i = 1,4
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .+ (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      end do
      
      end
      
      !Function to sum all energies
      subroutine TotalEnergy(xbc_len,cubbc_len,hcubbc,ucubbc,Gcubbc,ga,
     . beta1,beta2,n_GhstCells,dx,TotEnergVals)
      
      integer xbc_len,n_GhstCells,cubbc_len
      DOUBLE PRECISION hcubbc(cubbc_len),ucubbc(cubbc_len),
     . Gcubbc(cubbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION TotEnergVals(4)
      
      DOUBLE PRECISION CellEnergVals(4)
      integer i,j

      !running totals for energy values, start at 0
      do i = 1,4
         TotEnergVals(i) = 0.0
      end do
      
      !just loop over interior of hbc,Gbc, ubc which have interior values + ghost cell values
      do j= n_GhstCells + 1, xbc_len - n_GhstCells
      
         call  AllEnergiesIntegralCell(xbc_len,cubbc_len,hcubbc,
     .    ucubbc,Gcubbc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,4
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end

      


