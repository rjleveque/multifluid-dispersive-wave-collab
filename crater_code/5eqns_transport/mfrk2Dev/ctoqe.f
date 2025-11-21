c
c ----------------------------------------------------------
c
      subroutine c2prm_1d(maxm,meqn,mbc,maux,mx,q,qp,aux)

      implicit double precision (a-h,o-z)
      !dimension q(meqn,1-mbc:mx+mbc)
      dimension q(meqn,1-mbc:maxm+mbc)
      !dimension aux(maux,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:maxm+mbc)
      !dimension qp(meqn,1-mbc:mx+mbc)
      dimension qp(meqn,1-mbc:maxm+mbc)
c
c     # 5-equation transport model
c     # change conservative to primitive
c     # 2-phase case
c
      do i=1-mbc,mx+mbc
         qp(1,i) = q(1,i)
         qp(2,i) = q(2,i)
         qp(3,i) = q(3,i)/(q(1,i)+q(2,i))
         qp(4,i) = q(4,i)/(q(1,i)+q(2,i))
         qp(5,i) = q(5,i)-0.5d0*(q(3,i)**2+q(4,i)**2)/
     &             (q(1,i)+q(2,i))
         qp(6,i) = q(6,i)
         enddo
      return
c
      end
c
c --------------------------------------------------------------
c
      subroutine prm2c_edge(maxm,meqn,mbc,maux,mx,
     &                      ql,qr,uu,q,aux)

      implicit double precision (a-h,o-z)
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      dimension q(meqn,1-mbc:maxm+mbc)
      dimension aux(maux,1-mbc:maxm+mbc)
      dimension uu(meqn,2,1-mbc:maxm+mbc)
c
c     # 5-equation transport model
c     # change primitive variables to conservative variables
c     # 2-phase case
c
      do i=0,mx+1
c        # left cell-edge 
c
         ql(1,i) = uu(1,1,i)
         ql(2,i) = uu(2,1,i)
         ql(3,i) = (ql(1,i)+ql(2,i))*uu(3,1,i)
         ql(4,i) = (ql(1,i)+ql(2,i))*uu(4,1,i)
         ql(5,i) = uu(5,1,i)+
     &             0.5d0*(ql(3,i)**2+ql(4,i)**2)/
     &             (ql(1,i)+ql(2,i))
         ql(6,i) = uu(6,1,i)
c
c        # right cell-edge
         qr(1,i) = uu(1,2,i)
         qr(2,i) = uu(2,2,i)
         qr(3,i) = (qr(1,i)+qr(2,i))*uu(3,2,i)
         qr(4,i) = (qr(1,i)+qr(2,i))*uu(4,2,i)
         qr(5,i) = uu(5,2,i)+
     &             0.5d0*(qr(3,i)**2+qr(4,i)**2)/
     &             (qr(1,i)+qr(2,i))
         qr(6,i) = uu(6,2,i)
         enddo
c
c     # positivity-preserving check
      do i=0,mx+1      
         if ((ql(1,i) .lt. 0.d0) .or. 
     &       (ql(2,i) .lt. 0.d0) .or.
     &       (qr(1,i) .lt. 0.d0) .or.
     &       (qr(2,i) .lt. 0.d0)) then
             write(*,*) 'i=',i
             write(*,*) 'left-state'
             write(*,*) (qr(m,i-1),m=1,meqn)
             write(*,*) 'right-state'
             write(*,*) (ql(m,i),m=1,meqn)
             write(*,*) 'orig-state'
             write(*,*) (q(m,i),m=1,meqn)
             stop
         endif
      enddo
c
      return
      end
