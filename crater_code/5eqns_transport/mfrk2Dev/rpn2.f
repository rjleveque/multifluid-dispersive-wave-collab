c
c ---------------------------------------------------------------
c
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &           auxl,auxr,wave,s,amdq,apdq)
      implicit double precision (a-h,o-z)
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      dimension s(mwaves,1-mbc:maxm+mbc)
      dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
      dimension ql_state(20),qr_state(20)
      dimension aux1(30),aux2(30)
      dimension delta(20)
      dimension wave_local(20,20)
      dimension s_local(20)
      dimension speeds(2,20)
c
c     # Riemann solver in the normal direction for the Euler
c     # equations on a cartesian grid
c
      if (ixy .eq. 1) then
          mu = 3
          mv = 4
      else
          mu = 4
          mv = 3
      endif
c
      do i=2-mbc,mx+mbc-1
         do m=1,meqn
            ql_state(m) = qr(m,i-1)
            qr_state(m) = ql(m,i)
            enddo
c
         ql_state(3) = qr(mu,i-1)
         ql_state(4) = qr(mv,i-1)
         qr_state(3) = ql(mu,i)
         qr_state(4) = ql(mv,i)
c
         do ma=1,maux
            aux1(ma) = auxr(ma,i-1)
            aux2(ma) = auxl(ma,i)
            enddo
c
         do m=1,meqn
            delta(m) = qr_state(m)-ql_state(m)
            enddo
c
         call rp_solver(ql_state,qr_state,aux1,aux2,
     &        delta,s_local,wave_local,meqn)
c
         area = 1.0d0               
         do mw=1,mwaves
            wave_3loc = wave_local(3,mw)
            wave_4loc = wave_local(4,mw)
c
            wave_local(mu,mw) = wave_3loc
            wave_local(mv,mw) = wave_4loc
c
            speeds(1,mw) = area*dmin1(s_local(mw),0.d0)
            speeds(2,mw) = area*dmax1(s_local(mw),0.d0)
            s(mw,i) = speeds(1,mw)+speeds(2,mw)
            enddo

         do m=1,meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
               wave(m,mw,i) = wave_local(m,mw)
               amdq(m,i) = amdq(m,i)+speeds(1,mw)*wave(m,mw,i)
               apdq(m,i) = apdq(m,i)+speeds(2,mw)*wave(m,mw,i)
               enddo
            enddo
         enddo
      return
c
      end
