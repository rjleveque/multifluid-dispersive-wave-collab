c
c --------------------------------------------------------------
c
      subroutine prmtoc_local(q,rho1,rho2,vx,vy,p1,p2,zfa0,zfb0)
      implicit double precision (a-h,o-z)

      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      dimension q(20)
c
c     # 5-equation transport model 
c     # change primitive variables to conservative variables
c     # 2-phase case
c
c     # phase 1: ideal gas EOS
      rhoe1 = rhoe_prho_ideal(p1,geos(1))
      rhok1 = 0.5d0*rho1*(vx**2+vy**2)
c
c     # phase 2: stiffened gas EOS
      rhoe2 = rhoe_prho_sg(p2,geos(2),bn(2))
      rhok2 = 0.5d0*rho2*(vx**2+vy**2)
c
c     # partial density
      q(1) = zfa0*rho1
c
c     # partial density
      q(2) = zfb0*rho2
c
c     # mixture momentum in x-dir
      q(3) = (zfa0*rho1+zfb0*rho2)*vx
c
c     # mixture momentum in y-dir
      q(4) = (zfa0*rho1+zfb0*rho2)*vy
c
c     # total energy
      q(5) = zfa0*(rhoe1+rhok1)+zfb0*(rhoe2+rhok2)
c
      q(6) = zfa0
      return
c
      end
