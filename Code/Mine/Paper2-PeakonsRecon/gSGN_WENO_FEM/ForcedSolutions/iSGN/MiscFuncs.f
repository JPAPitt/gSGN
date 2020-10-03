c ====================================
c Miscellaneous functions that can be used
c
c ====================================

c get length of string without trailing whitespace
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
      
      
c fill array with equallyspaced points

      subroutine EqualSpaced(startv,endv,num,stepsize, res)
      DOUBLE PRECISION startv,endv
      INTEGER num
      DOUBLE PRECISION res(num)
      
      DOUBLE PRECISION stepsize
      INTEGER i
            
      do i = 1,num
         res(i) = startv + (i-1)*stepsize
      end do
      
   
   
      end
      
c calculate L2 norm
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
   
   
