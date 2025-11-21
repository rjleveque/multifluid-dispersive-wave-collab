c
c ---------------------------------------------------------------
c
      subroutine rp_solver(ql_state,qr_state,aux1,aux2,
     &           delta,s_local,wave_local,meqn)

      implicit double precision (a-h,o-z)
      dimension ql_state(20),qr_state(20)
      dimension aux1(30),aux2(30)
      dimension delta(20)
      dimension s_local(20),wave_local(20,20)
      dimension flux1(20),flux2(20)
      dimension qhll(20),work(20)
      dimension qml(20),qmr(20)
c
c     # 5-equation transport model
c     # HLLC approximate Riemann solver
c     # Davis' wave speed estimate
c
      call ctomfd(ql_state,rhoa1,rhob1,u1,v1,rhoh1,p1,
     &     grue1,pref1,zfa1,zfb1,c21)
      call ctomfd(qr_state,rhoa2,rhob2,u2,v2,rhoh2,p2,
     &     grue2,pref2,zfa2,zfb2,c22)
c
      if (c21  .lt. 0.d0 .or. c22  .lt. 0.d0) then
          write(6,*) 'error in hllc_solver:'
          write(6,*) 'c21=',c21,', c22=',c22
          write(6,*) 'rhoa1=',rhoa1,', rhob1=',rhob1
          write(6,*) 'rhoa2=',rhoa2,', rhob2=',rhob2
          write(6,*) 'p1=',p1,', p2=',p2
          write(6,*) 'zfa1=',zfa1,', zfa2=',zfa2
          write(6,*) 'zfb1=',zfb1,', zfb2=',zfb2
          write(6,*) 'ql_state=',(ql_state(m),m=1,meqn)
          write(6,*) 'qr_state=',(qr_state(m),m=1,meqn)
          stop
      endif
c
      s_local(1) = dmin1(u1-dsqrt(c21),u2-dsqrt(c22))
      s_local(3) = dmax1(u1+dsqrt(c21),u2+dsqrt(c22))
c
      flux1(1) = ql_state(1)*u1
      flux1(2) = ql_state(2)*u1
      flux1(3) = ql_state(3)*u1+p1
      flux1(4) = ql_state(4)*u1
      flux1(5) = (ql_state(5)+p1)*u1
c
      flux2(1) = qr_state(1)*u2
      flux2(2) = qr_state(2)*u2
      flux2(3) = qr_state(3)*u2+p2
      flux2(4) = qr_state(4)*u2
      flux2(5) = (qr_state(5)+p2)*u2
c
c     # compute qhll state     
      do m=1,3
         qhll(m) = ((s_local(3)*qr_state(m)-
     &      s_local(1)*ql_state(m))-
     &      (flux2(m)-flux1(m)))/(s_local(3)-s_local(1))
         enddo
c            
c     # set mid-state velocity
      s_local(2) = qhll(3)/(qhll(1)+qhll(2))
c
c     # set middle state (left)
      do m=1,5
         work(m) = s_local(1)*ql_state(m)-flux1(m)
         enddo
      qml(1) = work(1)/(s_local(1)-s_local(2))
      qml(2) = work(2)/(s_local(1)-s_local(2))
      qml(3) = (qml(1)+qml(2))*s_local(2)
      qml(4) = (qml(1)+qml(2))*v1
      pm1    = s_local(2)*(work(1)+work(2))-work(3)
      qml(5) = (work(5)+s_local(2)*pm1)/(s_local(1)-s_local(2))
      qml(6) = ql_state(6)
c
c     # set middle state (right)
      do m=1,5
         work(m) = s_local(3)*qr_state(m)-flux2(m)
         enddo
      qmr(1) = work(1)/(s_local(3)-s_local(2))
      qmr(2) = work(2)/(s_local(3)-s_local(2))
      qmr(3) = (qmr(1)+qmr(2))*s_local(2)
      qmr(4) = (qmr(1)+qmr(2))*v2
      pm2    = s_local(2)*(work(1)+work(2))-work(3)
      qmr(5) = (work(5)+s_local(2)*pm2)/(s_local(3)-s_local(2))
      qmr(6) = qr_state(6)
c
      do m=1,meqn
         wave_local(m,1) = qml(m)-ql_state(m)
         wave_local(m,2) = qmr(m)-qml(m)
         wave_local(m,3) = qr_state(m)-qmr(m)
         enddo
      return
c
      end
