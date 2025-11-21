c
c ------------------------------------------------------------
c
      subroutine ctomfd(q,rhoa0,rhob0,vx0,vy0,rhoh0,p0,
     &           grue0,pref0,zfa0,zfb0,c20)
      implicit double precision (a-h,o-z)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        acav(10),rhosat(10),psat(10),acsat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      dimension q(20)           
c
c     # 5-equation transport model
c     # change variables from conservative to primitive
c     # 2-phase case
c
      rhoa0 = q(1)
      rhob0 = q(2)
      vx0   = q(3)/(q(1)+q(2))
      vy0   = q(4)/(q(1)+q(2))
      rhoe0 = q(5)-0.5d0*(q(3)**2+q(4)**2)/(q(1)+q(2))
      zfa0  = q(6)
      zfb0  = 1.0d0-q(6)     
c
c     # phase1: ideal gas EOS
      ph1   = 0.d0
      eh1   = 0.d0              
      grue1 = geos(1)-1.d0
c
c     # phase 2: stiffened gas EOS
      ph2   = -geos(2)*bn(2)
      eh2   = 0.0d0  
      grue2 = geos(2)-1.d0
c
      grue0 = 1.d0/(q(6)/grue1+(1.0d0-q(6))/grue2)
      eref0 = q(1)*eh1+q(2)*eh2
      pref0 = q(6)*ph1/grue1+(1.0d0-q(6))*ph2/grue2
c
      p0 = grue0*(rhoe0-eref0+pref0)
c
c     # Saurel etal: frozen speed of sound
      zc1 = geos(1)*p0
      zc2 = geos(2)*(p0+bn(2))
      c20 = (zfa0*zc1+zfb0*zc2)/(rhoa0+rhob0)
c
      rhoh0 = q(5)+p0
      return
c
      end
