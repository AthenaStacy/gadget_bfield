      real*8 function chebyshev(n,x)
      integer n
      real*8  x

      if (n .eq. 0) then
        chebyshev = 1
        return
      else if (n .eq. 1) then
        chebyshev = x
        return
      else
        chebyshev = 2d0 * x * chebyshev(n-1,x) - chebyshev(n-2,x)
        return
      endif

      end
