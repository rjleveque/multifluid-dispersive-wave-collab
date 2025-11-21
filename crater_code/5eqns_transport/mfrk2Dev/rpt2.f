c
c ---------------------------------------------------------------
c           
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
     &           aux1,aux2,aux3,asdq,bmasdq,bpasdq)

      implicit double precision (a-h,o-z)
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      !dimension s(1-mbc:maxm+mbc,mwaves)
      !dimension wave(1-mbc:maxm+mbc,meqn,mwaves)
      dimension asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc)
      dimension aux2(maux,1-mbc:maxm+mbc)
      dimension aux3(maux,1-mbc:maxm+mbc)
      dimension speeds(2,20),s_local(20),wave_local(20,20)
      dimension ql_state(20),qr_state(20)
      dimension delta(20)
      dimension rot(4)
c
c     # Riemann solver in the transverse direction for the Euler
c     # equations on a curvilinear grid
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Use the same idea as in rpn2 but now rotate into the direction
c     # normal to the cell edge above or below this cell.
c
      if (ixy .eq. 1) then
          rot(1) = 0.0d0
          rot(2) = 1.0d0
          rot(3) = -1.0d0
          rot(4) = 0.0d0
      else
          rot(1) = 1.0d0
          rot(2) = 0.0d0
          rot(3) = 0.0d0
          rot(4) = 1.0d0
      endif
c
      do i=2-mbc,mx+mbc
         do m=1,meqn
            ql_state(m) = qr(m,i-1)   
            qr_state(m) = ql(m,i)    
c
            delta(m) = asdq(m,i)
            enddo
c
         call state2_rotate(rot,ql_state)
         call state2_rotate(rot,qr_state)
c
         call state2_rotate(rot,delta)
c
c        # solve for minus side
         call rpt_solver(ql_state,qr_state,delta,
     &        s_local,wave_local,meqn)
c
         area = 1.0d0
         do mw=1,mwaves
            call state2_rotate_tr(rot,wave_local(1,mw))
            speeds(1,mw) = area*dmin1(s_local(mw),0.d0)
            enddo
c
         do m=1,meqn
            bmasdq(m,i) = 0.d0
            do mw=1,mwaves
               bmasdq(m,i) = bmasdq(m,i)+
     &                       speeds(1,mw)*wave_local(m,mw)
               enddo
            enddo
c
c        # solve for plus side
         call rpt_solver(ql_state,qr_state,delta,
     &        s_local,wave_local,meqn)
c
         area = 1.0d0              
         do mw=1,mwaves
            call state2_rotate_tr(rot,wave_local(1,mw))
            speeds(2,mw) = area*dmax1(s_local(mw),0.d0)
            enddo
c
         do m=1,meqn
            bpasdq(m,i) = 0.d0
            do mw=1,mwaves
               bpasdq(m,i) = bpasdq(m,i)+ 
     &                       speeds(2,mw)*wave_local(m,mw)
               enddo
            enddo
         enddo                         
      return
c      
      end
