      
c     Sine recont
      subroutine Simpleq(a0,a1,a2,x_bc,qa_bc,qcubic_bc,
     . dx,cubic_len,xbc_len)
      
      integer xbc_len,cubic_len
      DOUBLE PRECISION dx,a0,a1,a2
      DOUBLE PRECISION x_bc(xbc_len),qa_bc(xbc_len),
     .   qcubic_bc(cubic_len)
      
      integer i
      
      DOUBLE PRECISION xjmh, xjms, xjps,xjph
      
      
      do i=1,xbc_len
           xjmh = x_bc(i) - 0.5*dx
           xjms = x_bc(i) - dx/6.0
           xjps = x_bc(i) + dx/6.0
           xjph = x_bc(i) + 0.5*dx
           
           qcubic_bc(4*i - 3) = a0 + a1*dsin(a2*xjmh)
           qcubic_bc(4*i - 2) = a0 + a1*dsin(a2*xjms)
           qcubic_bc(4*i - 1) = a0 + a1*dsin(a2*xjps)
           qcubic_bc(4*i ) = a0 + a1*dsin(a2*xjph)

           qa_bc(i) = 1.0/dx*( (a0*xjph - a1/a2*(dcos(a2*xjph))) 
     .            - (a0*xjmh - a1/a2*(dcos(a2*xjmh))) )
           
      end do
      
      
      end
     
      subroutine SimpleqhuG(a0,a1,a2,a3,a4,beta1,x_bc,
     . hcub_bc,ucub_bc,Gcub_bc,
     . dx,cubbc_len,xbc_len)
      
      integer xbc_len,cubbc_len
      DOUBLE PRECISION dx,a0,a1,a2,a3,a4,beta1
      DOUBLE PRECISION x_bc(xbc_len),
     .   hcub_bc(cubbc_len), ucub_bc(cubbc_len), Gcub_bc(cubbc_len)
      
      integer i
      
      DOUBLE PRECISION xjmh, xjms, xjps,xjph,
     . hxmh,hxms,hxps,hxph,uxmh,uxms,uxps,uxph,
     . uxxmh,uxxms,uxxps,uxxph
      
      do i=1,xbc_len
           xjmh = x_bc(i) - 0.5*dx
           xjms = x_bc(i) - dx/6.0
           xjps = x_bc(i) + dx/6.0
           xjph = x_bc(i) + 0.5*dx
           
           hcub_bc(4*i - 3) = a0 + a1*dsin(a2*xjmh)
           hcub_bc(4*i - 2) = a0 + a1*dsin(a2*xjms)
           hcub_bc(4*i - 1) = a0 + a1*dsin(a2*xjps)
           hcub_bc(4*i ) = a0 + a1*dsin(a2*xjph)
           
           hxmh = a1*a2*dcos(a2*xjmh)
           hxms = a1*a2*dcos(a2*xjms)
           hxps = a1*a2*dcos(a2*xjps)
           hxph = a1*a2*dcos(a2*xjph)
           
           ucub_bc(4*i - 3) = a3*dcos(a4*xjmh)
           ucub_bc(4*i - 2) = a3*dcos(a4*xjms)
           ucub_bc(4*i - 1) = a3*dcos(a4*xjps)
           ucub_bc(4*i ) = a3*dcos(a4*xjph)
           
           uxmh = -a3*a4*dsin(a4*xjmh)
           uxms = -a3*a4*dsin(a4*xjms)
           uxps = -a3*a4*dsin(a4*xjps)
           uxph = -a3*a4*dsin(a4*xjph)
           
           uxxmh = -a3*a4*a4*dcos(a4*xjmh)
           uxxms = -a3*a4*a4*dcos(a4*xjms)
           uxxps = -a3*a4*a4*dcos(a4*xjps)
           uxxph = -a3*a4*a4*dcos(a4*xjph)
                      
           Gcub_bc(4*i - 3) = ucub_bc(4*i - 3)*hcub_bc(4*i - 3) - 
     .      beta1/2.0*(3*(hcub_bc(4*i - 3)**2)*hxmh*uxmh +
     .      (hcub_bc(4*i - 3)**3)*uxxmh)
     
           Gcub_bc(4*i - 2) = ucub_bc(4*i - 2)*hcub_bc(4*i - 2) - 
     .      beta1/2.0*(3*(hcub_bc(4*i - 2)**2)*hxms*uxms +
     .      (hcub_bc(4*i - 2)**3)*uxxms)

           Gcub_bc(4*i - 1) = ucub_bc(4*i - 1)*hcub_bc(4*i - 1) - 
     .      beta1/2.0*(3*(hcub_bc(4*i - 1)**2)*hxps*uxps +
     .      (hcub_bc(4*i - 1)**3)*uxxps)

           Gcub_bc(4*i) = ucub_bc(4*i)*hcub_bc(4*i) - 
     .      beta1/2.0*(3*(hcub_bc(4*i)**2)*hxph*uxph +
     .      (hcub_bc(4*i)**3)*uxxph)
           
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



c    
      subroutine L2Norm(xbc_len,qbc_a,qbc_f,L2n) 
      
      integer xbc_len
      DOUBLE PRECISION qbc_a(xbc_len),qbc_f(xbc_len)
      integer i
      
      DOUBLE PRECISION L2num,L2den,L2n
      
      L2num = 0.0
      L2den = 0.0
      
      do i = 1, xbc_len
         L2num = L2num + (qbc_a(i) - qbc_f(i))**2
         L2den = L2den+  (qbc_a(i))**2
      end do  
      
      if (L2den .LT. 10.0**(-12)) then
         L2n = dsqrt(L2num)
      
      else 
         L2n = dsqrt(L2num /L2den)
      end if 
      
      end    
      
      
      subroutine Experiment(a0,a1,a2,a3,a4,beta1,xstart,dx,xbc_len,
     .   nGhstCells,cubbc_len,L2eG,L2eh,L2eu,wdir,effeclenwdir)
      
         Integer xbc_len,cubbc_len,nGhstCells,effeclenwdir,i
         DOUBLE PRECISION a0,a1,a2,a3,a4,beta1,xstart,dx,L2eG,L2eh
     .   ,L2eu
         
         CHARACTER(len=effeclenwdir) wdir
         
         DOUBLE PRECISION x_bc(xbc_len),ha_bc(xbc_len),Ga_bc(xbc_len),
     .   ucub_bc(cubbc_len),Gcub_bc(cubbc_len),hcub_bc(cubbc_len),
     .   hcubN_bc(cubbc_len),ucubN_bc(cubbc_len),GcubN_bc(cubbc_len)
         
         call Generatexbc(xstart,dx,xbc_len,nGhstCells,x_bc)
         
     
         call SimpleqhuG(a0,a1,a2,a3,a4,beta1,x_bc,
     . hcub_bc,ucub_bc,Gcub_bc,dx,cubbc_len,xbc_len)
     
         do i = 1,cubbc_len
           hcubN_bc(i)= hcub_bc(i)
           ucubN_bc(i)= ucub_bc(i)
           GcubN_bc(i)= Gcub_bc(i)
         end do
         
         open(1, file = wdir//'xhuGInit.dat') 
         do i = 1,xbc_len
            write(1,*) x_bc(i) -0.5*dx, hcub_bc(4*i -3), ucub_bc(4*i -3)
     .        , Gcub_bc(4*i -3)
            write(1,*) x_bc(i) - dx/6, hcub_bc(4*i -2), ucub_bc(4*i -2)
     .       , Gcub_bc(4*i -2)
            write(1,*) x_bc(i) + dx/6, hcub_bc(4*i -1), ucub_bc(4*i -1)
     .       , Gcub_bc(4*i -1)
            write(1,*) x_bc(i) + 0.5*dx, hcub_bc(4*i), ucub_bc(4*i)
     .       , Gcub_bc(4*i )
         end do
         
         close(1)

         call CubicToCellAverage(cubbc_len,xbc_len,hcub_bc,dx,ha_bc)     
         call CubicToCellAverage(cubbc_len,xbc_len,Gcub_bc,dx,Ga_bc)
     
         open(1, file = wdir//'xhaGa.dat') 
         do i = 1,xbc_len
            write(1,*) x_bc(i), ha_bc(i), Ga_bc(i)
         end do
         
         close(1)
         print *, xbc_len,cubbc_len,beta1,beta2,
     .      nGhstCells,dx,dt
         call EvolveStepWrap(xbc_len,cubbc_len,ha_bc,Ga_bc,10.0,
     . beta1,beta2,10.0**(-12),nGhstCells,dx,dx, 
     . hcubN_bc,GcubN_bc,ucubN_bc)
     
c        print *, 'Recon',cubic_len,xbc_len, nGhstCells,dx
c        call CellAverageToCubic(cubbc_len,xbc_len,
c     . nGhstCells,ha_bc,dx,hcubN_bc,10.0**(-12))

c        call CellAverageToCubic(cubbc_len,xbc_len,
c     . nGhstCells,Ga_bc,dx,GcubN_bc,10.0**(-12))     
        ! test reconstruction
     
c        print *, 'Norm 1'
        call L2Norm(cubbc_len,Gcub_bc,GcubN_bc,L2eG) 
        
c        print *, 'Norm 2'
        call L2Norm(cubbc_len,hcub_bc,hcubN_bc,L2eh) 
        
        
c        call FEM(hcub_bc,Gcub_bc,ucubN_bc,beta1,dx,xbc_len,
c     .   cubbc_len, 3*xbc_len + 1,nGhstCells)
     
     
     
        call L2Norm(cubbc_len,ucub_bc,ucubN_bc,L2eu) 
     
         open(1, file = wdir//'xhuGNum.dat') 
         do i = 1,xbc_len
            write(1,*) x_bc(i) -0.5*dx, hcubN_bc(4*i -3),
     .        ucubN_bc(4*i -3), Gcubn_bc(4*i -3)
            write(1,*) x_bc(i) - dx/6, hcubN_bc(4*i -2),
     .      ucubN_bc(4*i -2) , GcubN_bc(4*i -2)
            write(1,*) x_bc(i) + dx/6, hcubN_bc(4*i -1),
     .       ucubN_bc(4*i -1), GcubN_bc(4*i -1)
            write(1,*) x_bc(i) + 0.5*dx, hcubN_bc(4*i),
     .       ucubN_bc(4*i), GcubN_bc(4*i )
         end do   
         
         close(1)  
        
      
      end
      
      
      
      program main
         
      implicit none
   
      Integer wdirlen
      PARAMETER(wdirlen= 300)
      
      INTEGER effeclenwdir,n_GhstCells,lowestresx,x_len,xbc_len,
     . cubic_len, expi
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) strdiri
      
      DOUBLE PRECISION a0,a1,a2,a3,a4,beta1,xstart,xend,dx,
     . L2eG,L2eh,L2eu
      
      wdir = "../Results/Validation/Recon/Test1/"
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      CALL SYSTEM('rm -rf '//wdir(1:effeclenwdir))
      CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir)) 
      
      a0 = 1.0
      a1 = 0.5
      a2 = 0.1
      a3 = 2.0
      a4 = 0.1
      !beta1 = 0.0
      beta1 = 2.0/3.0
      
      xstart = -50d0
      xend = 50d0
      
      n_GhstCells = 3
      lowestresx = 3
      
      open(2, file = wdir(1:effeclenwdir)//'L2err.dat') 
      
      do expi = 0,0
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells
         cubic_len = 4*xbc_len 
         
         dx = (xend - xstart) / (x_len -1)
         !dx = (xend - xstart) / (x_len)
         
         print *,'++++++++++++ Experiment : ',expi ,' || ', '# Cells :',
     .    x_len , '++++++++++++'
         
         call Experiment(a0,a1,a2,a3,a4,beta1,xstart,dx,xbc_len,
     .   n_GhstCells,cubic_len,L2eG,L2eh,L2eu,
     . wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)
     
     
         write(2,*) dx,L2eG,L2eh,L2eu

      end do
      
      close(2)
     
      end
 
