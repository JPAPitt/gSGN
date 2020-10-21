      subroutine ReconTestExp(a0,a1,a2,xstart,xbc_len,n_GhstCells,dx,
     . tolcalc,NormFile,ExpWdir,ExpWdir_len)
      
     
      implicit none
      integer xbc_len,n_GhstCells,i,ExpWdir_len,NormFile
      CHARACTER(len=ExpWdir_len) ExpWdir
      
      DOUBLE PRECISION xbc(xbc_len),hA(xbc_len),
     . hjphNp(xbc_len),hjphA(xbc_len),hjphNm(xbc_len),
     . a0,a1,a2,xstart,tolcalc,hjph,L2Nm,L2Np,dx
     
     
      !generate cell nodes
      call Generatexbc(xstart,dx,xbc_len,n_GhstCells,xbc)
      
      !Get values
      call ReconTests(xbc,xbc_len,dx,hA,hjphA,a0,a1,a2) 
      
      !boundary
      do i = 1,xbc_len
         hjphNp(i) = hjphA(i)
         hjphNm(i) = hjphA(i)
      end do
      
      !interior 
      do i = n_GhstCells,xbc_len-n_GhstCells
         call Reconjph(hA(i-2),hA(i-1),hA(i),hA(i+1),hA(i+2),
     .    tolcalc,dx,hjph)
         hjphNm(i) = hjph
         call Reconjmh(hA(i-1),hA(i),hA(i+1),hA(i+2),hA(i+3),
     .    tolcalc,dx,hjph)
         hjphNp(i) = hjph
      end do
      
      
      call L2Norm(xbc_len,hjphA,hjphNm,L2Nm) 
      call L2Norm(xbc_len,hjphA,hjphNp,L2Np) 

      write(NormFile,*)  dx, L2Nm, L2Np
      
      open(11, file = ExpWdir//'Recons.dat') 
      do i = 1,xbc_len
         write(11,*) xbc(i), hjphNm(i), hjphNp(i), hjphA(i)
      end do
      close(11)

      
      end
      
      
      
      program main
         
      implicit none
   
      Integer wdirlen,NormFile
      PARAMETER(wdirlen= 300,NormFile=1)
      
      INTEGER effeclenwdir,n_GhstCells,lowestresx,x_len,xbc_len,
     . expi
      
      CHARACTER(len =wdirlen) wdir
      CHARACTER(len =2) strdiri
      
      DOUBLE PRECISION a0,a1,a2,xstart,xend,dx,
     . tolcalc
     
      
      wdir = "../Results/Validation/Recon/SINE/"
      call LenTrim(wdir,wdirlen,effeclenwdir)
      
      CALL SYSTEM('rm -rf '//wdir(1:effeclenwdir))
      CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir)) 
      
      open(NormFile, file = wdir(1:effeclenwdir)//'Norms.dat') 
      
      tolcalc = 10d0**(-12)
      
      a0 = 1.0
      a1 = 0.5
      a2 = 0.1
      
      xstart = -100.0
      xend = 100.0
      
      n_GhstCells = 3
      lowestresx = 100
      
      
      do expi = 1,10
      
         write (strdiri,'(I2.2)') expi
         
         CALL SYSTEM('mkdir -p '// wdir(1:effeclenwdir) //strdiri)
         
         x_len = lowestresx *2**(expi)
         xbc_len = x_len + 2 *n_GhstCells

         dx = (xend - xstart) / (x_len -1)

         print *, 'Experiment', expi
         call ReconTestExp(a0,a1,a2,xstart,xbc_len,n_GhstCells,dx,
     . tolcalc,NormFile,
     . wdir(1:effeclenwdir)//strdiri//'/',effeclenwdir+3)
         


      end do

      close(NormFile)

      
     
      end
 
