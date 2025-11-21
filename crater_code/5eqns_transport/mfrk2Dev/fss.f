c
c ---------------------------------------------------
c
      function fss(s)
      implicit double precision (a-h,o-z)
      common /fsscoef/ a00,a10,b00,b10
c   
      x = a00+a10*s
      y = b00+b10*s
      fss = fdisc_local(x,y)
      return
c
      end
