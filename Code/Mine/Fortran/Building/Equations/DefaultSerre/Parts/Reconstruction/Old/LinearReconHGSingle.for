c=====================================
c Program that given array of h calculates h* the reconstructed array
c=====================================
      program LinearRecon
      
      implicit none

      integer x_len, xR_len , n_PtsCell, n_GhstCells,nGP, i
      
      !xR_len = n_PtsCell*(x_len + 2*n_GhstCells)
      !nGP = n_PtsCell * n_GhstCells 
      parameter(x_len = 1000, n_PtsCell = 3, n_GhstCells = 1,
     .      xR_len= 3006, nGP=3)
                 
      real*8 ha(x_len),hc(x_len),xc(x_len)
      real*8 xReconGhost(xR_len),hR(xR_len),hRec(xR_len)
      real*8 hbeg(nGP),hend(nGP)


      !Initial data
      real*8 theta, x_start, x_end, dx,
     .   PI,Amp,Freq
      
      theta = 1.2
      x_start = 0
      x_end = 100
      dx = (x_end - x_start) / (x_len) 
      Amp = 1.5
      PI = 2*DASIN(1.0_8)
      Freq = 2*PI/ 10.0
      
      print *, "Linear Reconstruction"
      print *, "Left X Value" , x_start
      print *, "Right X Value" , x_end
      print *, "Delta x" , dx
      print *, "Length of x", x_len
      print *, "Length of xR", xR_len
      print *, "Theta", theta
      print *, "Sine Amp", Amp
      print *, "Sine Freq", Freq
      
      !Get locations of reconstruction plus Ghost Cells
      xReconGhost = GenerateAllLoc(x_start,dx,x_len,xR_len,n_PtsCell,
     . n_GhstCells)
     
      xc = xReconGhost(n_PtsCell*n_GhstCells + 2 :
     . xR_len-n_PtsCell*n_GhstCells : 3)
            
      !define ha - cell averages
      ha(:)  = (-Amp/Freq*(DCOS(Freq* (xc(:) + 0.5*dx)) 
     .           - DCOS(Freq* (xc(:) - 0.5*dx))))/dx

      ! define hc - cell node values
      hc(:)  = Amp*DSIN(Freq*xc(:) )     
      
      ! define h - analytic value of h at reconstructed points
      hR(:)  = Amp*DSIN(Freq*xReconGhost(:) ) 
      
      !Write Out input data
      open(1, file = 'InitialNodal.dat', status = 'unknown')  
      do i=1,x_len  
         write(1,*) xc(i),ha(i), hc(i)
      end do  
      close(1)  
      
      !Do Recontruction
      !Reconstruct each cell to giver x^+_j-1/2, x_j, x^-_j+1/2
      
      !get hbeg and hend
      hbeg(:) = hR(1:nGP)
      hend(:) = hR(xR_len- (nGP-1):xR_len)
      
      !Reconstruct h
      hRec =PerformReconstruction(ha,xc,dx,theta,hbeg,hend,
     . n_PtsCell, n_GhstCells)
      

      !Write Out At Reconstructed Points
      open(1, file = 'InitialRecon.dat', status = 'unknown')  
      do i=1,xR_len
         write(1,*) xReconGhost(i),hR(i), hRec(i)
      end do  
      close(1)   
      
      !Write out all, and error
         
      contains



c==============================
c Function to generate x values at all points: All points in cell 
c for reconstruction + ghost cells
c==============================
      function GenerateAllLoc(x_start,dx,x_len,xR_len,n_PtsCell,
     . n_GhstCells) result (xR)
      
      implicit none
      real*8,intent(in) :: x_start,dx
      integer,intent(in) :: x_len,xR_len,n_PtsCell,n_GhstCells
      real*8 xR(xR_len)
      
      integer i
      
      !Left boundary
      do i=1,n_GhstCells 
         xR(n_PtsCell*(i) -2)  =  x_start
     .   - (n_GhstCells- i +1)*dx
     
         xR(n_PtsCell*(i) -1)  =  x_start
     .   - (n_GhstCells- i +1)*dx + 0.5*dx
     
         xR(n_PtsCell*(i))  = x_start 
     .   - (n_GhstCells- i + 1)*dx + dx
     
      end do 

      !Interior
      do i=1,x_len  
         xR(n_PtsCell*(i +n_GhstCells )- 2)  = x_start
     .      + (i-1)*dx
         xR(n_PtsCell*(i +n_GhstCells )- 1)  = x_start
     .      + (i-1)*dx + 0.5*dx
         xR(n_PtsCell*(i +n_GhstCells ))  = x_start
     .      + (i)*dx
      end do  
      
      !Right boundary
      do i=1,n_GhstCells 
         xR(n_PtsCell*(x_len +n_GhstCells + i) -2)  = x_start
     .      + (x_len + i -1)*dx
         xR(n_PtsCell*(x_len +n_GhstCells + i) -1)  = x_start
     .      + (x_len + i -1)*dx + 0.5*dx
         xR(n_PtsCell*(x_len +n_GhstCells + i))  = x_start
     .      + (x_len + i)*dx
      end do 
         
      end function

c==============================
c Function to do Reconstruction - get h at x^+_j-1/2, x_j, x^-_j+1/2
c==============================
      function PerformReconstruction(h,x,dx,theta,hbeg,hend,
     . n_PtsCell, n_GhstCells) result (hR)
     
      implicit none
      
      integer, intent(in) :: n_PtsCell, n_GhstCells
      
      real*8,intent(in) :: theta,dx
      real*8, dimension(:), intent(in)   :: h,x
      real*8,intent(in) :: hbeg(n_PtsCell*n_GhstCells),
     . hend(n_PtsCell*n_GhstCells)
      
      real*8, dimension( n_PtsCell*(SIZE(x)+ 2*n_GhstCells) ) :: hR
      
      integer i,x_len
      real*8 fdh,mdh,bdh,cdh
      x_len = SIZE(x)
      
      !Left boundary
      hR(1: n_PtsCell*n_GhstCells) = hbeg(:)
      
      !Left most cell - needs hbeg to calculate derivative
      !cell 1 by i count, left cell point is n_PtsCell*n_GhstCells
      i = 1
      fdh = h(2) - h(1)
      mdh = 0.5*(h(2) - hbeg(2))
      bdh = h(1) - hbeg(2)
      
      cdh = minmod(theta*fdh,mdh,theta*bdh)
      
      hR(n_PtsCell*(i+n_GhstCells)-2) = h(1) - 0.5*cdh
      hR(n_PtsCell*(i+n_GhstCells)-1) = h(1)
      hR(n_PtsCell*(i+n_GhstCells)) = h(1) + 0.5*cdh
      
      !Interior
      do i=2,x_len-1
         fdh = h(i+1) - h(i)
         mdh = 0.5*(h(i+1) - h(i-1))
         bdh = h(i) - h(i-1)
         
         cdh = minmod(theta*fdh,mdh,theta*bdh)
         
         hR(n_PtsCell*(i+n_GhstCells)-2) = h(i) - 0.5*cdh
         hR(n_PtsCell*(i+n_GhstCells)-1) = h(i)
         hR(n_PtsCell*(i+n_GhstCells)) = h(i) + 0.5*cdh
      end do
      
      !Right most cell - needs hend to calculate derivative
      !cell x_len by i count, left cell point is n_PtsCell*n_GhstCells
      i=x_len
      fdh = hend(2) - h(x_len)
      mdh = 0.5*(hend(2) - h(x_len-1))
      bdh = h(x_len) - h(x_len-1)
      
      cdh = minmod(theta*fdh,mdh,theta*bdh)
      
      hR(n_PtsCell*(i+n_GhstCells)-2) = h(i) - 0.5*cdh
      hR(n_PtsCell*(i+n_GhstCells)-1) = h(i)
      hR(n_PtsCell*(i+n_GhstCells)) = h(i) + 0.5*cdh
         
      !Right Boundary
      hR(xR_len- (n_PtsCell*n_GhstCells-1):xR_len) = hend(:)
      
      end function
      
      function minmod(a,b,c) result (d)
      implicit none
      real*8, intent(in) :: a,b,c
      real*8 d
      
      if ((a > 0) .and. (b>0) .and. (c>0)) then
         d = min(a,b,c)
      else if ((a < 0) .and. (b<0) .and. (c<0)) then
         d = max(a,b,c)
      else
         d = 0.0
      end if
      
      end function
   
      
      end program LinearRecon
      
