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
   
   
