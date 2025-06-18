      function amag3(a)
c
c computes the magnitude of a vector using a dot product
c
      implicit none
      real*8 amag3,a(3)

      amag3 = dsqrt(a(1)**2 +a(2)**2 + a(3)**2)
      return
      end
