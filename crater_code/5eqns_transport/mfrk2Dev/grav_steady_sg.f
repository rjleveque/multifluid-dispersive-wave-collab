c
c ------------------------------------------------------------
c
      double precision function grav_steady_sg(rho,
     &                          q0,u0,H0,gs0,phi,geos0)

      implicit double precision (a-h,o-z)
      !common /grav_steady_info/ q0,u0,H0,gS0,phi,geos0
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      rhoK = 0.5d0*(q0**2/rho+rho*u0**2)
      grav_steady_sg = rho*(H0-phi)-
     &             (gS0*rho**geos0+rhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', grav_steady_sg=',grav_steady_sg
      return
c
      end
