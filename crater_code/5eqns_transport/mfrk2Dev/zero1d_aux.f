c
c -----------------------------------------------------------------
c
      double precision function zero1d_aux(x0,f,df,tol,iflag,
     &                              q0,ubc,H0,gs0,phibc,geos0)

      implicit double precision (a-h,o-z)
      parameter(itmax=1200,c0=0.5d0,alf0=1.0d-10)
      external f,df
c
c     # root-finding routine
c     # hybrid safeguard Newton and secant method
c
      iflag = 1
c
      alf  = 1.0d0
      iter = 0
      EPS  = epsilon(1.0d0)
c
      x1  = x0
      x2  = x1
      f1  = f(x1,q0,ubc,H0,gs0,phibc,geos0)
      df1 = df(x1,q0,ubc,H0,gs0,phibc,geos0)
c
 19   continue
c     # stopping criterion
      if ((dabs(f1)  .lt. tol) .or.
     &    (dabs(df1) .lt. EPS)) then     
          zero1d_aux = x1
          return
      endif
      if (iter .gt. itmax) go to 89
c
 29   continue
c     # safeguard Newton's iteration
      xn = x1-alf*f1/df1
c
c      write(66,*) 'x1=',x1,xn,f1,df1
c
      if (xn .lt. EPS) then
c         # safeguard step size reduction
c
          alf = c0*alf
c
          if (alf .lt. alf0) then
              write(66,*) 'error in zero1d_aux: alf=',alf
c
              iflag = -1
              return
          endif
          go to 29
      else
          fn = f(xn,q0,ubc,H0,gs0,phibc,geos0)
          if (dabs(fn) .lt. tol .or.
     &        dabs((x1-xn)/x1) .lt. tol .or.
     &        dabs(f1-fn) .lt. tol) then
              zero1d_aux = xn
              return
          endif
c
          if (f1*fn .le. 0.d0) then
              x2 = xn
              go to 99
          else
              x1 = xn
              x2 = xn
              f1  = fn
              df1 = df(x2,q0,ubc,H0,gs0,phibc,geos0)
c
              alf  = 1.d0
              iter = iter+1
              go to 19
          endif
      endif
c
 89   continue
c      write(66,*) 'zero1d_aux'
c      write(66,*) 'iter=',iter,', alf=',alf
c      write(66,*) 'x0=',x0,', x1=',x1,', x2=',x2
c      write(66,*) 'f1=',f1,', df1=',df1
c
      iflag = -2
c      stop
      return
c
 99   continue
      zero1d_aux = zeroin_local(x1,x2,f,tol,q0,ubc,H0,gs0,phibc,geos0)
      return
c
      end
