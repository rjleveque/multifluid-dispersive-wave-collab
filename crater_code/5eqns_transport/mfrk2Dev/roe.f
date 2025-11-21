c
c-----------------------------------------------------------------
c
      subroutine rpt_solver(ql_state,qr_state,
     &           delta,s_local,wave_local,meqn)
      implicit double precision (a-h,o-z)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),em(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      dimension ql_state(20),qr_state(20)
      dimension s_local(20),wave_local(20,20)
      dimension delta(20)
c
c     # Roe approximate Riemann solver for 5-equation model
c     # Ref: G. Allaire, S. Clerc, & S. Kokh,
c     #      A five-equation model for the simulation of
c     #      interfaces between compressible fluids,
c     #      JCP 181, 577-616 (2002)    
c
c     # Compute Roe-averaged quantities:
c
      call ctomfd(ql_state,rhoa1,rhob1,u1,v1,rhoh1,p1,      
     &     grue1,pref1,zfa1,zfb1,c21)
      call ctomfd(qr_state,rhoa2,rhob2,u2,v2,rhoh2,p2,      
     &     grue2,pref2,zfa2,zfb2,c22)
c
c     # compute Roe averaged values
      rhol   = rhoa1+rhob1
      rhor   = rhoa2+rhob2
      rsqtl  = dsqrt(rhol)
      rsqtr  = dsqrt(rhor)
      rsqtlr = rsqtl+rsqtr
c
      rho = (rhor*rsqtl+rhol*rsqtr)/rsqtlr
      u   = (u1*rsqtl+u2*rsqtr)/rsqtlr
      v   = (v1*rsqtl+v2*rsqtr)/rsqtlr
      zy1 = (rhoa1/rsqtl+rhoa2/rsqtr)/rsqtlr
      zy2 = (rhob1/rsqtl+rhob2/rsqtr)/rsqtlr
c
      enth = (rhoh1/rsqtl+rhoh2/rsqtr)/rsqtlr
      engk = 0.5d0*(u**2+v**2)
      grue = rsqtlr/(rsqtl/grue1+rsqtr/grue2)
      pr   = grue*(p1/grue1*rsqtl+
     &             p2/grue2*rsqtr)/rsqtlr
c
      rhoe1 = (pr+geos(1)*bn(1))/(geos(1)-1.0d0)
      rhoe2 = (pr+geos(2)*bn(2))/(geos(2)-1.0d0)
      drhoe = rhoe1-rhoe2
c
      c2 = grue*(enth-engk)
c
      if (c2 .le. 0.d0) then
          write(6,*) 'error in c2=',c2
c
          write(6,*) 'c21=',c21,', c22=',c22
          write(6,*) 'rhoa1=',rhoa1,', rhob1=',rhob1
          write(6,*) 'rhoa2=',rhoa2,', rhob2=',rhob2
          write(6,*) 'p1=',p1,', p2=',p2
          write(6,*) 'u1=',u1,', u2=',u2
          write(6,*) 'v1=',v1,', v2=',v2
          write(6,*) 'zfa1=',zfa1,', zfa2=',zfa2
          flush(6)
          stop
      else
          c = dsqrt(c2)
      endif
c
c     # wave strength
      a0 = -grue/c2*(-engk*(delta(1)+delta(2))+
     &      u*delta(3)+v*delta(4)-delta(5)+drhoe*delta(6))
      a2 = delta(1)-zy1*a0
      a3 = delta(2)-zy2*a0
      a4 = -v*(delta(1)+delta(2))+delta(4)
      a6 = ((c-u)*(delta(1)+delta(2))+delta(3)-
     &      c*(a2+a3))/(2.d0*c)
      a1 = a0-a6
      a5 = delta(6)
c
c     # waves
      s_local(1) = u-c
      wave_local(1,1) = a1*zy1
      wave_local(2,1) = a1*zy2
      wave_local(3,1) = a1*(u-c)
      wave_local(4,1) = a1*v
      wave_local(5,1) = a1*(enth-u*c)
      wave_local(6,1) = 0.d0             
c
      s_local(2) = u
      wave_local(1,2) = a2
      wave_local(2,2) = a3   
      wave_local(3,2) = (a2+a3)*u
      wave_local(4,2) = (a2+a3)*v+a4
      wave_local(5,2) = (a2+a3)*engk+a4*v+a5*drhoe
      wave_local(6,2) = a5     
c
      s_local(3) = u+c
      wave_local(1,3) = a6*zy1
      wave_local(2,3) = a6*zy2
      wave_local(3,3) = a6*(u+c)
      wave_local(4,3) = a6*v
      wave_local(5,3) = a6*(enth+u*c)
      wave_local(6,3) = 0.d0
      return    
c
      end
