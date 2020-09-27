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
c      subroutine NumericalSolveTSPrint(tstart,tend,
c     . ga,beta1,beta2,theta,dx,dt,n_GhstCells,xbc_len,
c     . xbc,hbc_init,Gbc_init,ubc_init,
c     . currenttime,hbc_fin,Gbc_fin,ubc_fin,
c     . Energ_Init, Energ_Fin,
c     . tlist,tlist_len,ExpWdir,expwdirlen)
c     
c     
c      implicit none
c      
c      integer n_GhstCells,xbc_len,expwdirlen,tlist_len
cc      CHARACTER(len=expwdirlen) ExpWdir
c     DOUBLE PRECISION tstart,tend,ga,beta1,beta2,
c    . theta,dx,dt,currenttime
c     DOUBLE PRECISION xbc(xbc_len),hbc_init(xbc_len),
c     . Gbc_init(xbc_len),
c     . ubc_init(xbc_len),hbc_fin(xbc_len),
c     . Gbc_fin(xbc_len), ubc_fin(xbc_len)
c     
c      DOUBLE PRECISION tlist(tlist_len)
c      
c      DOUBLE PRECISION Energ_Init(4), Energ_Fin(4)
c      
c      integer i,ileft,iright,filecount
c      CHARACTER(len=2) strct
c      
c      
c      
c      !initial time
c      currenttime  = tstart
c      filecount = 1
c     
c     !loop over and set hbc_fin,Gbc_fin to initial conditions
c     ileft = 1
c     iright = xbc_len
cc      do i = ileft,iright
c         hbc_fin(i) = hbc_init(i) 
c         Gbc_fin(i) = Gbc_init(i) 
c         ubc_fin(i) = ubc_init(i)
c      end do
c      
c      
c      !calculate initial Energies
c      call GetufromhG(xbc_len,hbc_fin,Gbc_fin,ubc_fin,
c     . beta1,dx,n_GhstCells)
c     
c      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,Gbc_fin,ga,beta1,beta2,
c     . n_GhstCells,dx,Energ_Init)
c
c      !evolve the system through time
c      do while (currenttime  .LT. tend )   
c      
c         if ((currenttime + dt .GE. tlist(filecount)) 
c     .      .OR. (filecount .EQ. 1 ))  then
c              
c            write (strct,'(I2)') filecount
c            open(8, file = ExpWdir// strct //'.dat') 
c            do i = 1,xbc_len
c               write(8,*) currenttime,xbc(i),hbc_fin(i),
c     .            Gbc_fin(i),ubc_fin(i)
c            end do
c            close(8)
c            
c            filecount = filecount + 1
c         end if

         
c         call EvolveStepWrap(xbc_len,hbc_fin,Gbc_fin,ubc_fin,ga,
c     .   beta1,beta2,theta,n_GhstCells,dx,dt)
     
c         currenttime  = currenttime  + dt
c         print *, 'Current Time : ', currenttime 
c      end do
       
      !calculate end energies  
c      call GetufromhG(xbc_len,hbc_fin,Gbc_fin,ubc_fin,
c     . beta1,dx,n_GhstCells)
      
c      call TotalEnergy(xbc_len,hbc_fin,ubc_fin,Gbc_fin,ga,beta1,beta2,
c     . n_GhstCells,dx,Energ_Fin)
c                 
c      end


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


c Cell Average To Cubic

      subroutine CellAverageToCubic(cubicbc_len,cellavgbc_len,
     . nGhstCells,cellavgbc_q,dx,cubicbc_q,calctol)
      implicit none
      integer cubicbc_len,cellavgbc_len,nGhstCells
      DOUBLE PRECISION dx,calctol
      DOUBLE PRECISION cubicbc_q(cubicbc_len),cellavgbc_q(cellavgbc_len)
      
      integer i
      DOUBLE PRECISION  pjm3toja,pjm3tojb,pjm3tojc,pjm3tojd,
     . pjm2tojp1a,pjm2tojp1b,pjm2tojp1c,pjm2tojp1d,
     . pjm1tojp2a,pjm1tojp2b,pjm1tojp2c,pjm1tojp2d,
     .  pjtojp3a,pjtojp3b,pjtojp3c,pjtojp3d
       
      DOUBLE PRECISION Bjm3toj,Bjm2tojp1,Bjm1tojp2,Bjtojp3,
     . iw1,iw2,iw3,iw4,w1,w2,w3,w4,qa,qb,qc,qd
      
      
      do i=nGhstCells + 1,cellavgbc_len - nGhstCells 
           
         pjm3toja  =(cellavgbc_q(i) - 3*cellavgbc_q(i-1) 
     .  + 3*cellavgbc_q(i-2) - cellavgbc_q(i-3))/(6*dx**3)    
         pjm3tojb  =(2*cellavgbc_q(i) - 5*cellavgbc_q(i-1)
     . + 4*cellavgbc_q(i-2) - cellavgbc_q(i-3))/(2*dx**2)
         pjm3tojc  =(43*cellavgbc_q(i) - 69*cellavgbc_q(i-1) 
     .   + 33*cellavgbc_q(i-2) - 7*cellavgbc_q(i-3))/(24*dx)
         pjm3tojd  = 11*cellavgbc_q(i)/12 + 5*cellavgbc_q(i-1)/24
     .    - cellavgbc_q(i-2)/6 + cellavgbc_q(i-3)/24
     
         !interpolating j-2 to j+1    
         pjm2tojp1a  = (-3*cellavgbc_q(i) + cellavgbc_q(i+1) 
     .   + 3*cellavgbc_q(i-1) - cellavgbc_q(i-2))/(6*dx**3)
         pjm2tojp1b  =(-2*cellavgbc_q(i) + cellavgbc_q(i+1) 
     .   + cellavgbc_q(i-1))/(2*dx**2)
         pjm2tojp1c  = (15*cellavgbc_q(i) + 7*cellavgbc_q(i+1) 
     .   - 27*cellavgbc_q(i-1) + 5*cellavgbc_q(i-2))/(24*dx)
         pjm2tojp1d  =13*cellavgbc_q(i)/12 - cellavgbc_q(i+1)/24 
     .   - cellavgbc_q(i-1)/24
     
         !interpolating j-1 to j+2  
         pjm1tojp2a  = (3*cellavgbc_q(i) - 3*cellavgbc_q(i+1) 
     .   + cellavgbc_q(i+2) - cellavgbc_q(i-1))/(6*dx**3)
         pjm1tojp2b  = (-2*cellavgbc_q(i) + cellavgbc_q(i+1) 
     .   + cellavgbc_q(i-1))/(2*dx**2)
         pjm1tojp2c  = (-15*cellavgbc_q(i) + 27*cellavgbc_q(i+1) 
     .   - 5*cellavgbc_q(i+2) - 7*cellavgbc_q(i-1))/(24*dx)
         pjm1tojp2d  = 13*cellavgbc_q(i)/12 
     .   - cellavgbc_q(i+1)/24 - cellavgbc_q(i-1)/24
     
     
         !interpolating j to j+3
         pjtojp3a  = (-cellavgbc_q(i) + 3*cellavgbc_q(i+1) 
     .   - 3*cellavgbc_q(i+2) + cellavgbc_q(i+3))/(6*dx**3)
         pjtojp3b  = (2*cellavgbc_q(i) - 5*cellavgbc_q(i+1) 
     .   + 4*cellavgbc_q(i+2) - cellavgbc_q(i+3))/(2*dx**2)
         pjtojp3c  = (-43*cellavgbc_q(i) + 69*cellavgbc_q(i+1) 
     .   - 33*cellavgbc_q(i+2) + 7*cellavgbc_q(i+3))/(24*dx)
         pjtojp3d  = 11*cellavgbc_q(i)/12 + 5*cellavgbc_q(i+1)/24 
     .   - cellavgbc_q(i+2)/6 + cellavgbc_q(i+3)/24 
     
         !Smoothness indicators
         Bjm3toj = 2107*cellavgbc_q(i)**2/240 
     .    - 1567*cellavgbc_q(i)*cellavgbc_q(i-1)/40 
     .    + 3521*cellavgbc_q(i)*cellavgbc_q(i-2)/120 
     .    - 309*cellavgbc_q(i)*cellavgbc_q(i-3)/40 
     .    + 11003*cellavgbc_q(i-1)**2/240
     .    - 8623*cellavgbc_q(i-1)*cellavgbc_q(i-2)/120 
     .    + 2321*cellavgbc_q(i-1)*cellavgbc_q(i-3)/120 
     .    + 7043*cellavgbc_q(i-2)**2/240 
     .    - 647*cellavgbc_q(i-2)*cellavgbc_q(i-3)/40 
     .    + 547*cellavgbc_q(i-3)**2/240
         
         Bjm2tojp1 = 3443*cellavgbc_q(i)**2/240 
     .    - 1261*cellavgbc_q(i)*cellavgbc_q(i+1)/120 
     .    - 2983*cellavgbc_q(i)*cellavgbc_q(i-1)/120 
     .    + 267*cellavgbc_q(i)*cellavgbc_q(i-2)/40 
     .    + 547*cellavgbc_q(i+1)**2/240 
     .    + 961*cellavgbc_q(i+1)*cellavgbc_q(i-1)/120 
     .    - 247*cellavgbc_q(i+1)*cellavgbc_q(i-2)/120 
     .    + 2843*cellavgbc_q(i-1)**2/240 
     .    - 821*cellavgbc_q(i-1)*cellavgbc_q(i-2)/120 
     .    + 89*cellavgbc_q(i-2)**2/80
         
         Bjm1tojp2 = 3443*cellavgbc_q(i)**2/240 
     .    - 2983*cellavgbc_q(i)*cellavgbc_q(i+1)/120 
     .    + 267*cellavgbc_q(i)*cellavgbc_q(i+2)/40 
     .    - 1261*cellavgbc_q(i)*cellavgbc_q(i-1)/120
     .    + 2843*cellavgbc_q(i+1)**2/240 
     .    - 821*cellavgbc_q(i+1)*cellavgbc_q(i+2)/120 
     .    + 961*cellavgbc_q(i+1)*cellavgbc_q(i-1)/120 
     .    + 89*cellavgbc_q(i+2)**2/80 
     .    - 247*cellavgbc_q(i+2)*cellavgbc_q(i-1)/120 
     .    + 547*cellavgbc_q(i-1)**2/240

         Bjtojp3 =2107*cellavgbc_q(i)**2/240 
     .    - 1567*cellavgbc_q(i)*cellavgbc_q(i+1)/40 
     .    + 3521*cellavgbc_q(i)*cellavgbc_q(i+2)/120 
     .    - 309*cellavgbc_q(i)*cellavgbc_q(i+3)/40 
     .    + 11003*cellavgbc_q(i+1)**2/240 
     .    - 8623*cellavgbc_q(i+1)*cellavgbc_q(i+2)/120 
     .    + 2321*cellavgbc_q(i+1)*cellavgbc_q(i+3)/120 
     .    + 7043*cellavgbc_q(i+2)**2/240 
     .    - 647*cellavgbc_q(i+2)*cellavgbc_q(i+3)/40 
     .    + 547*cellavgbc_q(i+3)**2/240
     
         iw1 = (1.0/35.0) *( 1.0 / (calctol +Bjm3toj )**2)
         iw2 = (12.0/35.0) * ( 1.0/ (calctol +Bjm2tojp1 )**2)
         iw3 = (18.0/35.0) * (1.0 / (calctol +Bjm1tojp2 )**2)
         iw4 = (4.0/35.0) * (1.0 / (calctol +Bjtojp3 )**2)
         
         w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
         w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
         w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
         w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
   
         
         qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
         qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
         qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
         qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
                  
         cubicbc_q(4*i -3) = qa*(-1.0*dx/2)**3 + qb*(-1.0*dx/2)**2 +
     .   qc*(-1.0*dx/2) + qd
         cubicbc_q(4*i -2) = qa*(-1.0*dx/6)**3 + qb*(-1.0*dx/6)**2 
     .   + qc*(-1.0*dx/6) + qd
         cubicbc_q(4*i -1 ) = qa*(dx/6)**3 + qb*(dx/6)**2 
     .   + qc*(dx/6) + qd
         cubicbc_q(4*i) = qa*(dx/2)**3 + qb*(dx/2)**2 
     .   + qc*(dx/2) + qd
     

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
      
      
      !print *, 1
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
      
      print *, 1
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,habc,dx,hcub,calctol)
     
      print *, 2
      call CellAverageToCubic(cubbc_len,xbc_len,
     . nGhstCells,Gabc,dx,Gcub,calctol)
     
      print *, 3
      !Calculate u
      call FEMforU(hcub,Gcub,ucub,beta1,dx,xbc_len,
     .   cubbc_len, 3*xbc_len + 1,nGhstCells)
     
      print *, 4
      call SingleEulerStep(xbc_len,cubbc_len,hcub,Gcub,ucub,ga,
     . beta1,beta2,nGhstCells,dt,dx,habc,Gabc,hapbc,Gapbc)
      print *, 5
            
      !Update cell averages (first order approximation to h^{n+1}, G^{n+1})
!      call SingleEulerStep(xbc_len,hbc,Gbc,ubc,ga,beta1,beta2
!     . ,theta,n_GhstCells,dt,dx,hpbc,Gpbc)
     
      !Get ubc, from hp and Gp from FD method (upbc must be initialised with BC properly)
!      call GetufromhG(xbc_len,hpbc,Gpbc,ubc,beta1,dx,n_GhstCells)
      
      !Update cell averages (first order approximation to h^{n+2}, G^{n+2})
!      call SingleEulerStep(xbc_len,hpbc,Gpbc,ubc,ga,beta1,beta2
!     . ,theta,n_GhstCells,dt,dx,hppbc,Gppbc)
      
      !use RK timestepping to convert h^n,G^n and first order approximation to h^{n+2}, G^{n+2}
      ! to second order approximation to approximation to h^{n+1}, G^{n+1}
      !since boundary conditions are constant, the average will be the initial value
!      do i= 1,xbc_len
!         hbc(i) = ( hbc(i ) + hppbc(i))/2d0
!         Gbc(i) = ( Gbc(i ) + Gppbc(i))/2d0
!      end do
          
      
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
      
      !print *,i, xbc_len,cubbc_len,ga,beta1,beta2,nGhstCells,dt,dx
            
      !calculates foh,foG which is flux across x_{i +1/2}
      ! it also updates cdhi,cdGi to be gradient across cell i + 1
      call Fluxxiph(hcub,Gcub,ucub,cubbc_len,i,ga,beta1,beta2,
     .   dx,foh,foG)
     
     
      !print *, i
     
      !flux out becomes flux in on next cell
      fih = foh
      fiG = foG
      do i = nGhstCells + 1, xbc_len - nGhstCells
      
         print *, i
      
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
     . dhir,ddhir,dhip1l,ddhip1l, alpha,
     . pua,pub,puc,pud,pha,phb,phc,phd,
     . pup1a,pup1b,pup1c,pup1d,php1a,php1b,php1c,php1d
     
      print *, cubbc_len,i,ga,beta1,beta2,dx,foh,foG
                

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
      phd = (-ucub(4*i-3) + 9*hcub(4*i-2) + 9*hcub(4*i-1) 
     . - hcub(4*i))/ 16.0
          
      duir = 3*pua*(dx/2)**2 + 2*pub*(dx/2) + puc
      dhir = 3*pha*(dx/2)**2 + 2*phb*(dx/2) + phc
      ddhir = 6*pha*(dx/2) + 2*phb
      
      print *, 'edge mid'
      
      uip1l = ucub(4*(i+1))
      hip1l = hcub(4*(i+1))
      Gip1l = Gcub(4*(i+1))
      
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
      php1d = (-ucub(4*(i+1)-3) + 9*hcub(4*(i+1)-2) 
     . + 9*hcub(4*(i+1)-1)  - hcub(4*(i+1)))/ 16.0
     
      duip1l = 3*pup1a*(-dx/2)**2 + 2*pup1b*(-dx/2) + pup1c
      dhip1l = 3*php1a*(-dx/2)**2 + 2*php1b*(-dx/2) + php1c
      ddhip1l = 6*php1a*(-dx/2) + 2*php1b
      
      print *, 'edge end'
      
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
     .      - beta2/4d0*ga*(hir**2)*(2*hir*ddhir + (dhir**2))
     
           
      ferh = uip1l*hip1l
      ferG = uip1l*Gip1l + ga*(hip1l**2)/2d0 
     .      - beta1*hip1l**3*duip1l**2
     .      - beta2/4d0*ga*(hip1l**2)*(2*hip1l*ddhip1l + (dhip1l**2))
     
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if 
      
      !calculate flux from cell i to cell i + 1 (Kurganov)
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))

      end

c ====
c ReconLinLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconLinLimGrad(qim1,qi,qip1,theta,cdq)
      DOUBLE PRECISION qim1,qi,qip1,theta,cdq,
     . fdq,mdq,bdq
      
      fdq = qip1 - qi
      mdq = 0.5*(qip1 - qim1)
      bdq = qi - qim1
      
      call minmod(theta*fdq,mdq,theta*bdq,cdq)
      
      end

c ====
c RecondqmLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine RecondqmLimGrad(qim2,qim1,qi,qip1,dx,cdq)
      DOUBLE PRECISION qim2,qim1,qi,qip1,dx,cdq
      DOUBLE PRECISION mdq,bdq
      
c      bdq = (2*qi - 3*qim1 + qim2)  /dx
      cdq = (qip1 - qi) /dx
      
c      call minmod(bdq,mdq,bdq,cdq)
      
      end

c ====
c RecondqpLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====

      subroutine RecondqpLimGrad(qi,qip1,qip2,qip3,dx,cdq)
      DOUBLE PRECISION qi,qip1,qip2,qip3,dx,cdq
      DOUBLE PRECISION fdq,mdq
      
      cdq = (qip1 - qi) /dx
c      fdq = (-2*qip1 + 3*qip2 - qip3) /dx
c      
c      call minmod(fdq,mdq,fdq,cdq)
      
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
      cdq = (qip1 - qim1) /(2d0*dx)
c      fdq = (-3*qi + 4*qip1 - qip2)  /(2d0*dx)
c      
c      call minmod(bdq,mdq,fdq,cdq)
      
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
     
      cddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)   
c      
c      call minmod(bpddq,bpddq,mddq,cddq)
      
      end

      
c ====
c ReconddqpLimGrad - produces the gradient across cell i using neighbouring values
c and a limiting function, with limiting parameter theta
c ====
      subroutine ReconddqpLimGrad(qim1,qi,qip1,qip2,qip3,
     . qip4,dx,cddq)
      DOUBLE PRECISION qim1,qi,qip1,qip2,qip3,qip4,dx,cddq
      DOUBLE PRECISION mddq,fddq,fpddq
      
      cddq = (qip2  - qip1 - qi + qim1 )/(2*dx**2)

c      fddq  = (5*qip1 - 13*qip2 + 11*qip3 - 3*qip4) /(2*dx**2)
c      fpddq = (3*qi - 7*qip1 + 5*qip2 - qip3) /(2*dx**2)
c      
c      call minmod(fpddq,fpddq,mddq,cddq)
      
      end
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
      subroutine AllEnergiesIntegralCell(xbc_len,h,u,G,ga,beta1,beta2,j
     . ,dx,CellEnergies)
      
      integer j,xbc_len
      DOUBLE PRECISION h(xbc_len),u(xbc_len),G(xbc_len)
      DOUBLE PRECISION dx,ga,beta1,beta2
      DOUBLE PRECISION CellEnergies(4)
      
      integer i
      DOUBLE PRECISION fGPe(4),sGPe(4),tGPe(4)
      
      DOUBLE PRECISION GPmxj,hGP,GGP,uGP,uxGP,hxGP
      DOUBLE PRECISION hCoeff(5), uCoeff(5), GCoeff(5)
      
      call QuarticInterp(h(j-2),h(j-1),h(j),h(j+1),h(j+2),dx,hCoeff)
      call QuarticInterp(u(j-2),u(j-1),u(j),u(j+1),u(j+2),dx,uCoeff)
      call QuarticInterp(G(j-2),G(j-1),G(j),G(j+1),G(j+2),dx,GCoeff)
      
      !first gauss point
      GPmxj = -dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      fGPe(1) = hGP
      fGPe(2) = GGP
      fGPe(3) = hGP*uGP
      fGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0

      !second gauss point
      GPmxj = 0.0 
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      sGPe(1) = hGP
      sGPe(2) = GGP
      sGPe(3) = hGP*uGP
      sGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !third gauss point
      GPmxj = dx*DSQRT(3.0d0/5.0d0)/2
      call QuarticCoeffEvalxj(hCoeff,GPmxj,hGP)
      call QuarticCoeffEvalGradxj(hCoeff,GPmxj,hxGP)
      call QuarticCoeffEvalxj(GCoeff,GPmxj,GGP)
      call QuarticCoeffEvalxj(uCoeff,GPmxj,uGP)
      call QuarticCoeffEvalGradxj(uCoeff,GPmxj,uxGP)
      
      tGPe(1) = hGP
      tGPe(2) = GGP
      tGPe(3) = hGP*uGP
      tGPe(4) = (hGP*uGP**2 + (1d0/3d0 + beta1/2d0)*(uxGP**2)*(hGP**3)
     . + ga*hGP**2*(1d0 + beta2/2d0*hxGP**2 )   )/2d0
      
      !weight the values at gauss points to get approximate integral over cell
      do i = 1,4
         CellEnergies(i) = (dx /2d0)*( (5.0/9.0)*fgpe(i) 
     .+ (8.0/9.0)*sgpe(i) + (5.0/9.0)*tgpe(i))
      end do
      
      end
      
      !Function to sum all energies
      subroutine TotalEnergy(xbc_len,hbc,ubc,Gbc,ga,beta1,beta2,
     . n_GhstCells,dx,TotEnergVals)
      
      integer xbc_len,n_GhstCells
      DOUBLE PRECISION hbc(xbc_len),ubc(xbc_len),Gbc(xbc_len)
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
         call  AllEnergiesIntegralCell(xbc_len,hbc,ubc,Gbc,ga,
     .      beta1,beta2,j,dx,CellEnergVals)
     
         !add cell energy value to running total
         do i = 1,4
            TotEnergVals(i) = TotEnergVals(i) + CellEnergVals(i)
         end do
      end do
      
      
      end

      


