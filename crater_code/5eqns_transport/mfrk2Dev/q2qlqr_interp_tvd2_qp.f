c
c--------------------------------------------------------------------
c

      subroutine q2qlqr_interp(ixy,maxm,meqn,mbc,maux,mx,
     &           dt,dx,q1d,aux1d,ql,qr,auxl,auxr,s,wave,qp,dq)
!                uu,hh,evl,evr)

      use amr_module, only : mthlim, mwaves

      implicit double precision (a-h,o-z)
      dimension q1d(meqn,1-mbc:maxm+mbc)
      dimension aux1d(maux,1-mbc:maxm+mbc)
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
      dimension qp(meqn,1-mbc:maxm+mbc)
      dimension dq(meqn,1-mbc:maxm+mbc)
      dimension uu(meqn,2,1-mbc:maxm+mbc)

      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension s(1-mbc:maxm+mbc,mwaves)

      ! didnt change the orderings below since 
      ! dont appear to be used
      dimension hh(1-mbc:maxm+mbc,-2:2)
      dimension evl(1-mbc:maxm+mbc,meqn,mwaves)
      dimension evr(1-mbc:maxm+mbc,meqn,mwaves)
c
c     # Second order TVD reconstruction 
c     # based on limited cell-average   
c
c     # initialization: piecewise constant
      do i=2-mbc,mx+mbc
      do m=1,meqn
         qr(m,i-1) = q1d(m,i-1)
         ql(m,i)   = q1d(m,i)
      enddo
      enddo
c
      do i=2-mbc,mx+mbc
      do ma=1,maux
         auxr(ma,i-1) = aux1d(ma,i-1)
         auxl(ma,i)   = aux1d(ma,i)
      enddo
      enddo
c
c     # change conservative to primitive 
      call c2prm_1d(maxm,meqn,mbc,maux,mx,q1d,qp,aux1d)

c     # compute limited slope
c     do i=1-mbc,mx+mbc
      do i=2-mbc,mx+mbc-1
      do m=1,meqn
c        if (i .eq. 1-mbc) then
c            # one-sided difference
c
c            dq(m,i) = qp(m,i+1)-qp(m,i)
c        else if (i .eq. mx+mbc) then
c            # one-sided difference
c
c            dq(m,i) = qp(m,i)-qp(m,i-1)
c            else
c            # central difference

             dqL = qp(m,i)-qp(m,i-1)
             dqR = qp(m,i+1)-qp(m,i)
c
             if (dabs(dqL) .lt. 1.0d-16) then
                 dq(m,i) = 0.0d0
             else
c                dq(m,i) = philim(dqL,dqR,mthlim(1))*dqL
                 dq(m,i) = max(0.d0,min(1.d0,dqR/dqL))*dqL
             endif
c        endif
      enddo
      enddo
c
c     # piecewise linear reconstruction
c     do i=1-mbc,mx+mbc
      do i=2-mbc,mx+mbc-1
      do m=1,meqn
         uu(m,1,i) = qp(m,i)-0.5d0*dq(m,i)
         uu(m,2,i) = qp(m,i)+0.5d0*dq(m,i)
      enddo    
      enddo    
c
c     # change primitive to conservative
      call prm2c_edge(maxm,meqn,mbc,maux,mx,
     &      ql,qr,uu,q1d,aux1d)
      return
c
      end
