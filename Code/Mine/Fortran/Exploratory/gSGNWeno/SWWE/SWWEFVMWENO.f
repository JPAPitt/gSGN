      subroutine Reconjmh(qajm2,qajm1,qaj,qajp1,qajp2,calctol,dx,qjmh)
      implicit none
      DOUBLE PRECISION qajm2,qajm1,qaj,qajp1,qajp2,calctol,
     . pjm2toja,pjm2tojb,pjm2tojc,pjm1tojp1a,pjm1tojp1b,pjm1tojp1c,
     . pjtojp2a,pjtojp2b,pjtojp2c,Bjm2toj,Bjm1tojp1,Bjtojp2,
     . djm2toj,djm1tojp1,djtojp2,
     . iwjm2toj,iwjm1tojp1,iwjtojp2,wjm2toj,wjm1tojp1,wjtojp2,
     . qa,qb,qc,sumw,dx,qjmh
      
      djm2toj = 3.0/10.0
      djm1tojp1 = 3.0/5.0
      djtojp2 = 1.0/10.0
      
      pjm2toja  =(qaj - 2*qajm1 + qajm2)/(2*dx**2)
      pjm2tojb  =(3*qaj - 4*qajm1 + qajm2)/(2*dx)
      pjm2tojc  =23*qaj/24 + qajm1/12 - qajm2/24
       
      pjm1tojp1a  = (-2*qaj + qajm1 + qajp1)/(2*dx**2)
      pjm1tojp1b  = (-qajm1 + qajp1)/(2*dx)
      pjm1tojp1c  = 13*qaj/12 - qajm1/24 - qajp1/24
      
      pjtojp2a  =(qaj - 2*qajp1 + qajp2)/(2*dx**2)
      pjtojp2b  =(-3*qaj + 4*qajp1 - qajp2)/(2*dx)
      pjtojp2c  =23*qaj/24 + qajp1/12 - qajp2/24      


      Bjm2toj = 10*qaj**2/3 - 31*qaj*qajm1/3 + 11*qaj*qajm2/3
     . + 25*qajm1**2/3 - 19*qajm1*qajm2/3 + 4*qajm2**2/3
     
      Bjm1tojp1 =qaj**2/3 - qaj*qajm1/3 - qaj*qajp1/3 + qajm1**2/3 
     . - qajm1*qajp1/3 + qajp1**2/3 + (-2*qaj + qajm1 + qajp1)**2
     
      Bjtojp2 = 10*qaj**2/3 - 31*qaj*qajp1/3 + 11*qaj*qajp2/3 
     . + 25*qajp1**2/3 - 19*qajp1*qajp2/3 + 4*qajp2**2/3
      
      iwjm2toj = (djm2toj / (calctol +Bjm2toj )**2)
      iwjm1tojp1 = (djm1tojp1 / (calctol +Bjm1tojp1 )**2)
      iwjtojp2 = (djtojp2 / (calctol +Bjtojp2 )**2)
      
      sumw = iwjm2toj  + iwjm1tojp1 + iwjtojp2
      wjm2toj = iwjm2toj  / sumw
      wjm1tojp1 = iwjm1tojp1 / sumw
      wjtojp2 = iwjtojp2  / sumw
      
      qa = wjm2toj*pjm2toja + wjm1tojp1*pjm1tojp1a + wjtojp2 *pjtojp2a
      qb = wjm2toj*pjm2tojb + wjm1tojp1*pjm1tojp1b + wjtojp2 *pjtojp2b
      qc = wjm2toj*pjm2tojc + wjm1tojp1*pjm1tojp1c + wjtojp2 *pjtojp2c
      
      qjmh = qa*(-dx/2)**2 + qb*(-dx/2) + qc
      
      end
 
 
       subroutine Reconjph(qajm2,qajm1,qaj,qajp1,qajp2,calctol,dx,qjph)
      implicit none
      DOUBLE PRECISION qajm2,qajm1,qaj,qajp1,qajp2,calctol,
     . pjm2toja,pjm2tojb,pjm2tojc,pjm1tojp1a,pjm1tojp1b,pjm1tojp1c,
     . pjtojp2a,pjtojp2b,pjtojp2c,Bjm2toj,Bjm1tojp1,Bjtojp2,
     . djm2toj,djm1tojp1,djtojp2,
     . iwjm2toj,iwjm1tojp1,iwjtojp2,wjm2toj,wjm1tojp1,wjtojp2,
     . qa,qb,qc,sumw,dx,qjph
      
      djm2toj = 1.0/10.0
      djm1tojp1 = 3.0/5.0
      djtojp2 = 3.0/10.0
      
      pjm2toja  =(qaj - 2*qajm1 + qajm2)/(2*dx**2)
      pjm2tojb  =(3*qaj - 4*qajm1 + qajm2)/(2*dx)
      pjm2tojc  =23*qaj/24 + qajm1/12 - qajm2/24
       
      pjm1tojp1a  = (-2*qaj + qajm1 + qajp1)/(2*dx**2)
      pjm1tojp1b  = (-qajm1 + qajp1)/(2*dx)
      pjm1tojp1c  = 13*qaj/12 - qajm1/24 - qajp1/24
      
      pjtojp2a  =(qaj - 2*qajp1 + qajp2)/(2*dx**2)
      pjtojp2b  =(-3*qaj + 4*qajp1 - qajp2)/(2*dx)
      pjtojp2c  =23*qaj/24 + qajp1/12 - qajp2/24      


      Bjm2toj = 10*qaj**2/3 - 31*qaj*qajm1/3 + 11*qaj*qajm2/3
     . + 25*qajm1**2/3 - 19*qajm1*qajm2/3 + 4*qajm2**2/3
     
      Bjm1tojp1 =qaj**2/3 - qaj*qajm1/3 - qaj*qajp1/3 + qajm1**2/3 
     . - qajm1*qajp1/3 + qajp1**2/3 + (-2*qaj + qajm1 + qajp1)**2
     
      Bjtojp2 = 10*qaj**2/3 - 31*qaj*qajp1/3 + 11*qaj*qajp2/3 
     . + 25*qajp1**2/3 - 19*qajp1*qajp2/3 + 4*qajp2**2/3
      
      iwjm2toj = (djm2toj / (calctol +Bjm2toj )**2)
      iwjm1tojp1 = (djm1tojp1 / (calctol +Bjm1tojp1 )**2)
      iwjtojp2 = (djtojp2 / (calctol +Bjtojp2 )**2)
      
      sumw = iwjm2toj  + iwjm1tojp1 + iwjtojp2
      wjm2toj = iwjm2toj  / sumw
      wjm1tojp1 = iwjm1tojp1 / sumw
      wjtojp2 = iwjtojp2  / sumw
      
      qa = wjm2toj*pjm2toja + wjm1tojp1*pjm1tojp1a + wjtojp2 *pjtojp2a
      qb = wjm2toj*pjm2tojb + wjm1tojp1*pjm1tojp1b + wjtojp2 *pjtojp2b
      qc = wjm2toj*pjm2tojc + wjm1tojp1*pjm1tojp1c + wjtojp2 *pjtojp2c
      
      qjph = qa*(dx/2)**2 + qb*(dx/2) + qc
      
      end

      
      subroutine FVMSWWE(hA,GA,hAp,GAp,xbc_len,n_GhstCells,grav,
     . dx,dt,calctol)
      implicit none
      INTEGER xbc_len,n_GhstCells,j
      DOUBLE PRECISION grav,dx,dt,calctol
      DOUBLE PRECISION hA(xbc_len),GA(xbc_len),hAp(xbc_len),
     . GAp(xbc_len)
     
      DOUBLE PRECISION hir,hip1l,Gir,Gip1l,uir,uip1l,
     . sl,sr,felh,ferh,felG,ferG,isrmsl,foh,foG,fih,fiG
     
      !boundaries
      do j = 1, xbc_len
         hAp(j) = hA(j)
         GAp(j) = GA(j)
      end do
      
      !interior
     
      j = n_GhstCells - 1
      call Reconjph(hA(j-2),hA(j-1),hA(j),hA(j+1),hA(j+2),
     . calctol,dx,hir)
      call Reconjmh(hA(j-1),hA(j),hA(j+1),hA(j+2),hA(j+3),
     . calctol,dx,hip1l)

      call Reconjph(GA(j-2),GA(j-1),GA(j),GA(j+1),GA(j+2),
     . calctol,dx,Gir)
      call Reconjmh(GA(j-1),GA(j),GA(j+1),GA(j+2),GA(j+3),
     . calctol,dx,Gip1l)
      
      uir = Gir/hir 
      uip1l = Gip1l/hip1l
      
      sl = min(0.0,uir - dsqrt(grav*hir), uip1l  - dsqrt(grav*hip1l))
      sr = max(0.0,uir + dsqrt(grav*hir), uip1l + dsqrt(grav*hip1l))
      
      felh = Gir
      ferh = Gip1l

      felG = uir*Gir + grav*hir*hir/2.0
      ferG = uip1l*Gip1l + grav*hip1l*hip1l/2.0
      
      if (sr == sl) then
         isrmsl = 0.0
      else
         isrmsl = 1.0 / (sr - sl)
      end if      
      
      foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
      foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir)) 
      
      fih = foh
      fiG = foG
      
      do j = n_GhstCells, xbc_len- n_GhstCells
      
         call Reconjph(hA(j-2),hA(j-1),hA(j),hA(j+1),hA(j+2),
     . calctol,dx,hir)
         call Reconjmh(hA(j-1),hA(j),hA(j+1),hA(j+2),hA(j+3),
     . calctol,dx,hip1l)

         call Reconjph(GA(j-2),GA(j-1),GA(j),GA(j+1),GA(j+2),
     . calctol,dx,Gir)
         call Reconjmh(GA(j-1),GA(j),GA(j+1),GA(j+2),GA(j+3),
     . calctol,dx,Gip1l)
         
         uir = Gir/hir 
         uip1l = Gip1l/hip1l
         
         sl = min(0.0,uir - dsqrt(grav*hir), uip1l  - dsqrt(grav*hip1l))
         sr = max(0.0,uir + dsqrt(grav*hir), uip1l + dsqrt(grav*hip1l))
         
         felh = Gir
         ferh = Gip1l

         felG = uir*Gir + grav*hir*hir/2.0
         ferG = uip1l*Gip1l + grav*hip1l*hip1l/2.0
         
         if (sr == sl) then
            isrmsl = 0.0
         else
            isrmsl = 1.0 / (sr - sl)
         end if      
         
         foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
         foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
         
         
         hAp(j) = hA(j) - dt*(foh - fih)/dx
         GAp(j) = GA(j) - dt*(foG - fiG)/dx 
         
         fih = foh
         fiG = foG
      
      end do
            
      end 
      
      
      subroutine RKStep(hA,GA,xbc_len,n_GhstCells,grav,dx,dt,calctol)
      implicit none
      INTEGER xbc_len,n_GhstCells,j
      DOUBLE PRECISION grav,dx,dt,calctol
      DOUBLE PRECISION hA(xbc_len),GA(xbc_len),hAp(xbc_len),
     . GAp(xbc_len),hApp(xbc_len),GApp(xbc_len)
     
      call FVMSWWE(hA,GA,hAp,GAp,xbc_len,n_GhstCells,grav,
     . dx,dt,calctol)

      call FVMSWWE(hAp,GAp,hApp,GApp,xbc_len,n_GhstCells,grav,
     . dx,dt,calctol)
     
      do j = 1,xbc_len  
         hA(j) = 0.5*(hA(j) + hApp(j))
         GA(j) = 0.5*(GA(j) + GApp(j))
      end do
     
     
      end 
      
      subroutine SWWESolver(hA,GA,xbc_len,n_GhstCells,grav,dx,dt,
     . calctol,startt,endt,currentt)
      INTEGER xbc_len,n_GhstCells
      DOUBLE PRECISION grav,dx,dt,calctol,startt,endt,currentt
      DOUBLE PRECISION hA(xbc_len),GA(xbc_len)
      
      currentt = startt
      do while (currentt  .LT. endt ) 
      
         call RKStep(hA,GA,xbc_len,n_GhstCells,grav,dx,dt,calctol)
         currentt  = currentt  + dt
         print *, 'Current Time : ', currentt 
      
      end do
      
      
      end
