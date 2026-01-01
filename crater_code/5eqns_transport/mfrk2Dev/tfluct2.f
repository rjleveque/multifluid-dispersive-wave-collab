c
c -------------------------------------------------------------
c
      subroutine tfluct2(ixy,maxm,meqn,mwaves,mbc,maux,mx,
     &           ql,qr,auxl,auxr,s,adq)

      implicit double precision (a-h,o-z)
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)

      dimension s(mwaves,1-mbc:maxm+mbc)
      dimension adq(meqn,1-mbc:maxm+mbc)
      dimension ql_state(20),qr_state(20)
      dimension aux1(30),aux2(30)
      dimension delta(20)
      dimension wave_local(20,20),s_local(20),speeds(20,2)
c
c     # compute fluctuation within cell
c
      if (ixy .eq. 1) then
          mu = 3
          mv = 4
      else
          mu = 4
          mv = 3
      endif
c
      do i=2-mbc,mx+mbc
         do ma=1,maux
            aux1(ma) = auxl(ma,i)
            aux2(ma) = auxr(ma,i)
            enddo
c
         do m=1,meqn
            ql_state(m) = ql(m,i)
            qr_state(m) = qr(m,i)
            enddo
c
         ql_state(3) = ql(mu,i)
         ql_state(4) = ql(mv,i)
         qr_state(3) = qr(mu,i)
         qr_state(4) = qr(mv,i)
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
            speeds(mw,1) = area*dmin1(s_local(mw),0.d0)
            speeds(mw,2) = area*dmax1(s_local(mw),0.d0)
            s(mw,i) = speeds(mw,1)+speeds(mw,2)
            enddo
c
         do m=1,meqn
            amdq = 0.d0
            apdq = 0.d0
            do mw=1,mwaves
               amdq = amdq+speeds(mw,1)*wave_local(m,mw)
               apdq = apdq+speeds(mw,2)*wave_local(m,mw)
               enddo
c
            adq(i,m) = amdq+apdq
            enddo
         enddo
      return
c
      end
