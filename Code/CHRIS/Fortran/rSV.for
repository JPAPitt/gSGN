      program rSV

c============================================================
c This program solves the Serre equations written in conservative 
c form using second-order TVD Runga-Kutta scheme. 
c
c Horizontal bed.
c============================================================

      implicit none

      integer md, nxmax, nxd, mn, i
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)
      
      real*8 h(nxd), Gg(nxd), z2(nxd), x(nxd), u(nxd)
      integer nx, nxx, version, isteps

      real*8 epsilon

      common /ivariables/ nx, nxx, version, isteps, epsilon

      real*8 toll,  theta,     g,     dx,      dt,    hobs,  
     .        uobs,    u0,    tf,     tc,      hl,      hr,
     .          h0,    x0,    u1,     h1,      u2,      h2,
     .          x1,    x2,    x3,    tau,  linf_h, linf1_h,
     .       l11_h,  l2_h, l21_h, linf_u, linf1_u,    l1_u,
     .       l11_u,  l2_u, l21_u,    cfl,      Gz,     Gzz,
     .        Gzzz,  l1_h


      common /rvariables/ toll, theta,  g,  dx,  dt, tf, tc, hl, hr, 
     .                       h0,   x0, tau,  h1

      call util_open_files

      epsilon=1.0d0
      
c-----------------------------------
c Parameters
c-----------------------------------
      nx = 10001
      tf = 10.d0 ! 48.d0

      nxx = nx + 2*md
      if (nx.gt.nxmax) write(*,*)'nx too large'

      toll = 1.0d-9
      theta = 1.d-1
      g = 1.d0

      hl  = 0.d0
      hr  = 0.d0

      h0 = 1.d0
      u0 = 0.d0
      h1 = 2.d0
      u1 = 0.d0

      x0 = 50.d0
      dx = 100.d0/dble(nx-1)
      cfl = 0.2d0

c-----------------------------------
c Boussinesq terms set to 0 to turn off
c or to 1 to include Boussinesq terms
c-----------------------------------
      Gz   = 1.d0*epsilon
      Gzz  = 1.d0*epsilon
      Gzzz = 1.d0*epsilon
	
      dt  = dx/(dsqrt(g*h1))*cfl
c      dt = 0.01d0
      tau = 0.0d0
      tc  = 0.d0

c-----------------------------------
c Setup Terrain
c-----------------------------------
      do i = 1, nxx
        x(i)  = dble(i-md-1)*dx
c
c     Horizontal bed
c----------------------------------
        z2(i) = 0.d0
      end do

c-----------------------------------
c Initial condition
c-----------------------------------
      do i = 1, nxx
        call initial condition(x(i), x0, h0, h1, h(i), u(i), Gg(i))
      end do

c-----------------------------------
c Boundary condition
c-----------------------------------
      call boundary_apply(h0, h1, h, u, Gg)

c-----------------------------------
c Setup for Time looping
c-----------------------------------
      isteps = 0
      call util_dump_solution(h, u, x)

c      call timer(ibegin)

c-----------------------------------
c Evolve
c-----------------------------------
      call evolver_rk_2(h, u, Gg, z2, Gz, Gzz, Gzzz)

      call util_dump_solution(h, u, x)

c      call timer(ifinish)

c-----------------------------------
c  Calculate the L1, L2 and L_infinity norm
c-----------------------------------
      linf_h  = 0.d0
      linf1_h = 0.d0
      l1_h    = 0.d0
      l11_h   = 0.d0
      l2_h    = 0.d0
      l21_h   = 0.d0
      linf_u  = 0.d0
      linf1_u = 0.d0
      l1_u    = 0.d0
      l11_u   = 0.d0
      l2_u    = 0.d0
      l21_u   = 0.d0

c-----------------------------------
c Analytical solution
c-----------------------------------
      call analytical solution(h0, h1, g, x0, x1, x2, x3, h2, u2, tc)

      write(22,*)tc
      do i = 1, nxx
        if(x(i).le.x1)then
          hobs = h1
          uobs = u1
        else if(x(i).lt.x3)then
          hobs = (((x0 - x(i))/tc + 2.d0*dsqrt(g*h1))/3.d0)**2/g
          uobs = 2.d0*dsqrt(g)*(dsqrt(h1) - dsqrt(hobs))
        else if(x(i).le.x2)then
          hobs = h2
          uobs = u2
        else
          hobs = h0
          uobs = u0
        end if
        
        write(20,*)x(i), hobs, uobs
        
        linf_h  = dmax1(linf_h,dabs(h(i) - hobs))
        linf1_h = dmax1(linf1_h,dabs(hobs))
        l1_h    = l1_h  + dabs(h(i) - hobs)
        l11_h   = l11_h + dabs(hobs)
        l2_h    = l2_h  + (h(i) - hobs)**2
        l21_h   = l21_h + hobs**2
        linf_u  = dmax1(linf_u,dabs(u(i) - uobs))
        linf1_u = dmax1(linf1_u,dabs(uobs))
        l1_u    = l1_u  + dabs(u(i) - uobs)
        l11_u   = l11_u + dabs(uobs)
        l2_u    = l2_u  + (u(i) - uobs)**2
        l21_u   = l21_u + uobs**2
      end do

      write(*,*)dx, dt
      linf1_h = linf_h/linf1_h
      l11_h   = l1_h/l11_h
      l2_h    = dsqrt(l2_h)
      l21_h   = l2_h/dsqrt(l21_h)

      linf1_u = linf_u/linf1_u
      l11_u   = l1_u/l11_u
      l2_u    = dsqrt(l2_u)
      l21_u   = l2_u/dsqrt(l21_u)

      write(*,*)l11_h, l21_h, linf1_h
      write(*,*)l11_u, l21_u, linf1_u
c      write(23,100)dx, l11_h, linf1_h, l11_u, linf1_u
c  100 format(5f15.8)

c-----------------------------------
c Close down
c-----------------------------------
      call util_close_files

      stop
      end

c===========================================================
c Boussinesq Flux
c===========================================================
      subroutine flux_hG(f1, f2,  h,  Gg, uh, em_x, l1, l2, g, toll)

      implicit none

      real*8  f1, f2, l1, l2, em_x, g, toll, uh, u, h, Gg, cvel 

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)     

      if(h.le.toll)then
        u  = 0.d0
        h  = 0.d0
        Gg = 0.d0
      else
        u = uh/h   
      endif
      
      cvel = dsqrt(g*h)
      l1   = u - cvel
      l2   = u + cvel
      em_x = dabs(u) + cvel

c-----------------------------------
c Flux
c-----------------------------------
      f1 = u*h
      f2 = u*Gg + 0.5d0*g*h**2

      return
      end

c===============================================================
c array limiter
c===============================================================
      subroutine array_limit(u, u_l, u_r)

      implicit none

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)
      integer nx, nxx, version, isteps
      real*8 epsilon
      common /ivariables/ nx, nxx, version, isteps, epsilon

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       tau, h1

      common /rvariables/ toll, theta,   g, dx, dt, tf, tc, hl, hr,
     .                       h0,    x0, tau, h1

      integer i

      real*8 a, b, t, xmin, xmic, u_x, d1, d2, u(nxd), u_l(nxd),
     .       u_r(nxd)

      xmin(a,b) = 0.5d0*(dsign(1.d0,a) + dsign(1.d0,b))*
     .            dmin1(dabs(a),dabs(b))
      xmic(t,a,b) = xmin(t*xmin(a,b), 0.5d0*(a + b) )

c
c Second-order reconstruction with General minimod limiter
c------------------------------------------------------------
      do i = md-1, nx + md+1
        d1     = u(i)   - u(i-1)
        d2     = u(i+1) - u(i)
        u_x    = xmic( theta, d1, d2 )
        u_l(i) = u(i) - 0.5d0*u_x
        u_r(i) = u(i) + 0.5d0*u_x
      end do

      return
      end

c===========================================================
c Output Routine
c===========================================================
      subroutine util_dump_solution(h, u, x)

      implicit none

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)

      real*8 h(nxd), u(nxd), x(nxd)

      integer nx, nxx, version, isteps
      real*8 epsilon
      common /ivariables/ nx, nxx, version, isteps, epsilon
      integer i

      isteps = isteps + 1

c      do i = md + 1,nx + md
      do i = 1,nxx
        write(21,100) x(i), h(i), u(i)*h(i)
      end do
  100 format(4f20.15)

      return
      end

c===============================================================
c The evolver using TVD Second-order Runga-Kutta method
c===============================================================
      subroutine evolver_rk_2(h, u, Gg, z2, Gz, Gzz, Gzzz)

      implicit none

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)

      integer nx, nxx, version, isteps
      real*8 epsilon
      common /ivariables/ nx, nxx, version, isteps, epsilon
      
      real*8 toll, theta, g, dx, dt, tf, tc, hl, hr, uhl, uhr, h0, 
     .         x0,    tau,  em_x, h1
      common /rvariables/ toll, theta,   g, dx, dt, tf, tc, hl, hr, 
     .                       h0,    x0, tau, h1

      real*8     h(nxd),    uh(nxd),     u(nxd),    Gg(nxd), 
     .          h_1(nxd),  uh_1(nxd),  u_1(nxd),  Gg_1(nxd),
     .          h_2(nxd),  uh_2(nxd),  u_2(nxd),  Gg_2(nxd),
     .          f_h(nxd),     v(nxd),       TVD,
     .          f_G(nxd),   h_l(nxd),  h_r(nxd),   uh_l(nxd),
     .          u_l(nxd),   u_r(nxd),  
     .         uh_r(nxd),  Gg_l(nxd), Gg_r(nxd),     z2(nxd),
     .           a1(nxd),    a2(nxd),   a3(nxd),          hh,
     .               hhh,         Gz,       Gzz,        Gzzz,
     .           AA(nxd),   AAA(nxd),    r(nxd)     

      integer i, j, nt, ip1, im1

      nt = 0

      do while (.true.)
      
c        epsilon = dmin1(1.d0/tc**3,1.d0)
        write(*,*) epsilon
        do i = 1, nxx
          if (h(i).lt.toll) then
            write(*,*)'i = ',i,' h(i) = ', h(i), tc
          endif
          f_h(i)  = 0.0d0
          f_G(i)  = 0.0d0
          h_l(i)  = h(i)
          h_r(i)  = h(i)
           uh(i)  = u(i)*h(i)
          uh_l(i) = uh(i)
          uh_r(i) = uh(i)
        end do

        nt = nt + 1

c-----------------------------------
c First step
c-----------------------------------    
      
c
c     Calculate the conserved quantity Gg
c-----------------------------------
          do i = md, nx+md 
            ip1 = i + 1
            im1 = i - 1
            AA(i)  = Gzz*h(i)*h(i)*(u(ip1) - u(im1))/2.d0/dx*
     .            (h(ip1) - h(im1))/2.d0/dx
            AAA(i) = Gz/3.d0*h(i)**3*(u(ip1) - 2.d0*u(i) + u(im1))/dx/dx
            Gg(i) = uh(i) - AAA(i) - AA(i)
            if(h(i).le.toll)then
              u(i)  = 0.d0
              h(i)  = 0.d0
              Gg(i) = 0.d0
            end if
          end do

c          
c     Boundary condition
c-----------------------------------
          call boundary_apply(h0, h1, h, u, Gg)

c
c     Calculate the momentum
c-----------------------------------
          do i = 1, nxx
            uh(i) = h(i)*u(i)     
          end do

c
c     Limit h, uh and G
c-----------------------------------         
          call array_limit  ( h,  h_l,  h_r)
          call array_limit  (uh, uh_l, uh_r)
          call array_limit  (Gg, Gg_l, Gg_r)

c
c     Set up ghost nodes
c-----------------------------------
          do i = 1, md - 1
            h_l(i)   = h_l(md)
            h_r(i)   = h_r(md)
            uh_l(i)  = uh_l(md)
            uh_r(i)  = uh_r(md)
            Gg_l(i)  = Gg_l(md)
            Gg_r(i)  = Gg_r(md)
          end do
          do i = nx + md + 1, nxx
            h_l(i)   = h_l(nx+md)
            h_r(i)   = h_r(nx+md)
            uh_l(i)  = uh_l(nx+md)
            uh_r(i)  = uh_r(nx+md)
            Gg_l(i)  = Gg_l(nx+md)
            Gg_r(i)  = Gg_r(nx+md)
          end do   
  
c          
c     Calculated the Flux
c-----------------------------------          
          call numerical_flux_central(uh_l, uh_r,  h_l,  h_r, Gg_l, 
     .                                 Gg_r,  f_h,  f_G, em_x, Gzzz, h)

c
c     Update the conserved quantities
c-----------------------------------
          do i = md, nx + md
             h_1(i)  =  h(i)  - dt/dx*(f_h(i) - f_h(i-1))
             Gg_1(i) =  Gg(i) - dt/dx*(f_G(i) - f_G(i-1))

c
c     Bed Slope
c-----------------------------------
             Gg_1(i) = Gg_1(i) - dt/dx*g*(z2(i) - z2(i-1))*(h(i))

c
c     Friction
c-----------------------------------
             if(h_1(i).gt.toll)then
               Gg_1(i) = Gg_1(i) - dt*tau*uh(i)
             else
               h_1(i)  = 0.d0
               Gg_1(i) = 0.d0
             end if

             if (h_1(i).lt.-toll) then
                write(*,*)' h_1(',i,') = ', h_1(i)
             endif
             if (Gg_1(i).lt.-toll) then
c                write(*,*)' Gg_1(',i,') = ', Gg_1(i)
c                write(*,*) tc, uh(i), aa(i), aaa(i)
             endif
          end do         
c         
c     Solve for u_1 by solving a Laplace equation which
c     results in  tridiagonal system of equations.
c
c     Set up coefficient matrix
c-------------------------------------------------
          do i = md + 1, nx + md - 1
            ip1 = i + 1
            im1 = i - 1
            j = i - md + 1
            hh  = Gz*h_1(i)**3/(3.d0*dx*dx)   
            hhh = Gzz*h_1(i)**2/(4.d0*dx*dx)*(h_1(ip1) - h_1(im1))
            a1(j)  =    -hh + hhh
            a2(j)  = h_1(i) + 2.d0*hh
            a3(j)  =    -hh - hhh
             r(j)  = Gg_1(i)
          end do
          
c
c           Upstream
c-------------------------------------------------
          a2(1) = 1.d0
          a3(1) = 0.d0
           r(1) = Gg_1(md)*0.d0

c
c           Downstream
c-------------------------------------------------
          a1(nx) = 0.d0
          a2(nx) = 1.0d0
           r(nx) = Gg_1(nx+md-1)*0.d0

c
c     Solve tridiagonal system for u_1
c-----------------------------------
          call tridag(a1, a2, a3, r, v, nx)  
     
          do j = 1, nx
            i = j + md - 1
            u_1(i) = v(j)
          end do
          do i = 1, md - 1
            u_1(i) = 0.d0
          end do
          do i = nx + md + 1, nxx
            u_1(i) = 0.d0
          end do
          
c-----------------------------------
c Second step
c-----------------------------------
c
c     Boundary condition
c-----------------------------------
          call boundary_apply(h0, h1, h_1, u_1, Gg_1)

c
c     Calculate the momentum
c-----------------------------------
          do i = 1, nxx
            uh_1(i) = h_1(i)*u_1(i)     
          end do

c
c     Set up ghost nodes
c-----------------------------------
          do i = 1, md - 1
            h_l(i)  = h_l(md)
            h_r(i)  = h_r(md)
            uh_l(i) = uh_l(md)
            uh_r(i) = uh_r(md)
            Gg_l(i) = Gg_l(md)
            Gg_r(i) = Gg_r(md)
          end do
          do i = nx + md + 1, nxx
            h_l(i)  = h_l(nx+md)
            h_r(i)  = h_r(nx+md)
            uh_l(i) = uh_l(nx+md)
            uh_r(i) = uh_r(nx+md)
            Gg_l(i) = Gg_l(nx+md)
            Gg_r(i) = Gg_r(nx+md)
          end do  
   
c
c     Limit h_1, u_1 and Gg_1
c-----------------------------------     
          call array_limit  ( h_1,   h_l,  h_r)
          call array_limit  (uh_1,  uh_l, uh_r)
          call array_limit  (Gg_1,  Gg_l, Gg_r)

c
c     Calculated the Flux
c-----------------------------------     
          call numerical_flux_central(uh_l, uh_r, h_l,  h_r, Gg_l, 
     .                                 Gg_r,  f_h, f_G, em_x, Gzzz, h)

c
c     Second step
c-----------------------------------
          do i = md, md + nx
             h_2(i)  =  h_1(i)  - dt/dx*(f_h(i) - f_h(i-1))
             Gg_2(i) =  Gg_1(i) - dt/dx*(f_G(i) - f_G(i-1))

c
c     Bed Slope
c-----------------------------------
             Gg_2(i) = Gg_2(i) - dt/dx*g*(z2(i) - z2(i-1))*(h_1(i))

c
c     Friction
c-----------------------------------
             if(h_2(i).gt.toll)then
               Gg_2(i) = Gg_2(i) - dt*tau*Gg_2(i)
             else
               h_2(i)  = 0.d0
               Gg_2(i) = 0.d0
             end if
             
             if (h_2(i).lt.-toll) then
c                write(*,*)' h_2(',i,') = ',h_2(i)
             endif
             if (Gg_2(i).lt.-toll) then
c                write(*,*)' Gg_2(',i,') = ',Gg_2(i)
             endif
                          
          end do

c-------------------------------------------------     
c Third step
c-------------------------------------------------     
          do i = md, nx + md

c    
c     Second-order Runga-Kutta
c-------------------------------------------------     
             h(i)  = (h_2(i)  +  h(i))/2.d0
             Gg(i) = (Gg_2(i) + Gg(i))/2.d0
          end do
          
c-------------------------------------------------          
c     Solve for u by solving a Laplace equation which
c     results in  tridiagonal system of equations.
c
c     Set up coefficient matrix
c-------------------------------------------------
          do i = md + 1, nx + md - 1
            ip1 = i + 1
            im1 = i - 1
            j = i - md + 1
            hh  = Gz*h(i)**3/(3.d0*dx*dx)   
            hhh = Gzz*h(i)**2/(4.d0*dx*dx)*(h(ip1) - h(im1))
            a1(j)  =  -hh + hhh
            a2(j)  = h(i) + 2.d0*hh
            a3(j)  =  -hh - hhh
             r(j)  = Gg(i)
          end do

c
c     Upstream
c-------------------------------------------------
          a2(1) = 1.d0
          a3(1) = 0.d0
           r(1) = Gg(md)*0.d0

c
c     Downstream
c-------------------------------------------------
          a1(nx) = 0.d0
          a2(nx) = 1.d0
           r(nx) = Gg(nx+md-1)*0.d0

c
c     Solve tridiagonal system for u_star
c-----------------------------------
          call tridag(a1, a2, a3, r, v, nx)

          do j = 1, nx
            i = j + md - 1
            u(i) = v(j)
          end do
          do i = 1, md - 1
            u(i) = 0.d0
          end do
          do i = nx + md + 1, nxx
            u(i) = 0.d0
          end do
          
          do i = 1, nxx
            uh(i) = h(i)*u(i)            
          end do
        
          tc = tc + dt
          write(*,*)tc
          isteps = isteps + 1

c
c     Boundary condition
c-----------------------------------
          call boundary_apply(h0, h1, h, u, Gg)
          
c
c     Calculate the momentum
c-----------------------------------
          do i = 1, nxx
            uh(i) = h(i)*u(i)
          end do

          if(tc.ge.tf)go to 1001

c-----------------------------------
c End timestep
c-----------------------------------
      end do

 1001 write(*,*)isteps

c-----------------------------------
c     Calculate TVD new
c-----------------------------------
      TVD = 0.d0
      do i = md, nx + md - 1
        TVD = TVD + abs(h(i+1) - h(i))
      end do
      write(25,*)dx,TVD
      
      return
      end

c===========================================================
c Numerical Central Flux
c===========================================================
      subroutine numerical_flux_central(uh_l, uh_r, h_l,  h_r, Gg_l,
     .                                   Gg_r,  f_h, f_G, em_x, Gzzz, h)

      implicit none

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)
      integer nx, nxx, version, isteps
      real*8 epsilon
      common /ivariables/ nx, nxx, version, isteps, epsilon
      
      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       tau, h1
      common /rvariables/ toll,  theta, g, dx, dt,
     .                    tf, tc, hl, hr, 
     .                    h0, x0, tau, h1

      real*8 f_h(nxd),   f_G(nxd),   h_l(nxd),      em_x,
     .        h_r(nxd),  uh_l(nxd), uh_r(nxd),  Gg_l(nxd), Gg_r(nxd),
     .       fh_l(nxd),  fh_r(nxd), fG_l(nxd),
     .       fG_r(nxd),  l1_l(nxd), l2_l(nxd),  l1_r(nxd), l2_r(nxd),
     .           a_max,      a_min,        aa,        axa,      Gzzz

      integer i, ip1, im1, ip2

      real *8   dh1,    dh2,    dh3, Bond, h(nxd)
            

      Bond = epsilon*1.d0/3.d0*g*h0**2
      
      
      em_x = 1.0d-15

      do i = 1, nxx - 1
        fh_l(i)  = 0.d0
        fh_r(i)  = 0.d0
        fG_l(i)  = 0.d0
        fG_r(i)  = 0.d0
        
c
c     Calculate flux
c-----------------------------------
         call flux_hG(fh_l(i), fG_l(i),  h_r(i),    Gg_r(i), uh_r(i),
     .                    em_x, l1_l(i), l2_l(i),          g, toll)
         call flux_hG(fh_r(i), fG_r(i), h_l(i+1), Gg_l(i+1), uh_l(i+1),
     .                    em_x, l1_r(i), l2_r(i),          g, toll)
      end do

      do i = 2, nxx - 2
        ip1 = i + 1
        im1 = i - 1
        ip2 = i + 2    
    
c
c     Additional Boussinesq terms
c-----------------------------------
         fG_l(i) = fG_l(i) - Gzzz*2.d0/3.d0*h_l(i)**3
     .            *((uh_l(i)/h_l(i) - uh_l(im1)/h_l(im1))/dx)**2
         fG_r(i) = fG_r(i) - Gzzz*2.d0/3.d0*h_r(i)**3
     .            *((uh_r(ip1)/h_r(ip1) - uh_r(i)/h_r(i))/dx)**2
    
         
c
c     Calculate max wave speed
c-----------------------------------
         a_max = dmax1(l2_l(i), l2_r(i), 0.d0)
         a_min = dmin1(l1_l(i), l1_r(i), 0.d0)

c
c     Compute the numerical flux
c-----------------------------------
         aa  = a_max - a_min
         axa = a_max*a_min
         if ( aa. gt. dx*dt/1000.0d0) then
            f_h(i)  = (a_max*fh_l(i) - a_min*fh_r(i) +
     .                  axa*(h_l(i+1) - h_r(i)))/aa
            f_G(i) = (a_max*fG_l(i)  - a_min*fG_r(i) +
     .                 axa*(Gg_l(i+1) - Gg_r(i)))/aa
            dh1 = (h(i+1)-h(i-1))/dx/2.d0
            dh2 = (h(i+1)-2.d0*h(i)+h(i-1))/dx/dx
            f_G(i) = f_G(i) - Bond*(h(i)**3*dh2+0.5d0*(dh1*h(i))**2)
         else
            f_h(i) = 0.d0
            f_G(i) = 0.d0
         endif

      end do

      return
      end

c===========================================================
c Open Files
c===========================================================
      subroutine util_open_files

      open(20,file='Serre_rk_2_40_Surface_Tension.obs')
      open(21,file='Serre_rk_2_40_Surface_Tension.r')
      open(22,file='Serre_rk_2_40_Surface_Tension.t')
      open(24,file='Serre_rk_2_40_Surface_Tension.out')
      open(23,file='Serre_rk_2_40_Surface_Tension.e',access='append')
      open(25,file='Serre_rk_2_40_Surface_Tension.tvd',access='append')

      return
      end

c===========================================================
c Close Files
c===========================================================
      subroutine util_close_files

      close(20)
      close(21)
      close(22)
      close(23)
      close(24)

      return
      end

c===========================================================
c Initial condition
c===========================================================
      subroutine initial condition(x, x0, h0, h1, water_depth, u, Gg)

      implicit none

      real*8 x, x0, h0, h1, water_depth, u, Gg
      real*8 hr, hl, delta

      water_depth = h0
      u           = 0.d0
      Gg          = 0.d0

      
      hr = 0.0
      hl = h1-h0
      delta = 0.5
      
!      if (x.lt.x0)then
        water_depth = h0+hl+0.5D0*(hr-hl)*(1.D0+dtanh(delta*(x-x0)))
!      end if


      return
      end

c===========================================================
c Boundary Conditions
c===========================================================
      subroutine boundary_apply(h0, h1, h, u, Gg)

      implicit none

      integer md, nxmax, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)
      real*8 h(nxd), u(nxd), Gg(nxd), h0, h1
      integer nx, nxx, version, isteps
      real*8 epsilon
      common /ivariables/ nx, nxx, version, isteps, epsilon
      
      integer i

      do i = 1, md

c
c Left Boundary
c-----------------------------------
        h(i) = h1
        u(i) = 0.d0
       Gg(i) = 0.d0
      end do

      do i = nx + md, nx + 2*md

c
c Right Boundary
c-----------------------------------
        h(i) = h0
        u(i) = 0.d0
       Gg(i) = 0.d0
      end do

      return
      end

c===========================================================
c Analytical solution
c===========================================================
      subroutine analytical solution(h0, h1, g, x0, x1, x2, x3,
     .                               h2, u2, tc)

      implicit none

      integer i
      real*8 h0, h1, g, x0, x1, x2, x3, h2, u2, zmin, zmax,
     .        z, c0, c1, c2, tc, func

      c0 = dsqrt(g*h0)
      c1 = dsqrt(g*h1)

c----------------------------------------------
c     Derive the analytical solution using the bisection method
c----------------------------------------------
      zmin = -100.d0
      zmax =  101.d0
      do i = 1,100
        z    = (zmin + zmax)/2.d0
        u2   = z - c0*c0/4.d0/z*(1.d0 + dsqrt(1.d0 + 8.d0*z*z/c0/c0))
        c2   = c0*dsqrt(0.5d0*(dsqrt(1.d0 + 8.d0*z*z/c0/c0) - 1.d0))
        func = 2.d0*c1/c0 - u2/c0 - 2.d0*c2/c0
        if(func.gt.0.d0)then
          zmin = z
        else
          zmax = z
        end if
      end do

      if(dabs(z).gt.99.d0)stop 'no convergence'

      h2 = h0/(1.d0 - u2/z)
      x3 = (u2 - c2)*tc + x0
      x2 = z*tc + x0
      x1 = -c1*tc + x0

      return
      end

c===========================================================
c Tridiagonal solver
c===========================================================
      subroutine tridag(a, b, c, r, u, n)

c
c Tridiagonal matrix solver from Numerical Recipes
c     Coefficient Matrix
c          b diagonal elements
c           a below diagonal elements
c           c above diagonal elements
c     r vector of right hand side
c     u vector of unknowns (solution)
c---------------------------------------------------------------------------
      implicit none

      integer j, n, nxmax, md, nxd, mn
      parameter (md = 3, nxmax = 150000, nxd = nxmax + 2*md, mn = 2)

      real*8 bet, gam(nxd), a(nxd), b(nxd), c(nxd), r(nxd), u(nxd)
      
      if(b(1).eq.0.d0) pause 'tridag: rewrite equations'
      bet  = b(1)
      u(1) = r(1)/bet
      do j = 2,n
        gam(j) = c(j-1)/bet
        bet    = b(j) - a(j)*gam(j)
        if(bet.eq.0.d0) pause 'tridag failed'
        u(j) = (r(j) - a(j)*u(j-1))/bet
      end do

      do j = n - 1, 1, -1
        u(j) = u(j) - gam(j+1)*u(j+1)
      end do

      return     
      end
