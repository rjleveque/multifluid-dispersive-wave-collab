c
c ------------------------------------------------------------
c
      double precision function dgrav_steady_sg(rho,
     &                       q0,u0,h0,gs0,phi,geos0)

      implicit double precision (a-h,o-z)
      !common /grav_steady_info/ q0,u0,H0,gS0,phi,geos0
c
c     # gravitational source terms
c     # steady equilibrium condition
c     # stiffened gas EOS case
c
      drhoK = 0.5d0*(-q0**2/rho**2+u0**2)
      dgrav_steady_sg = (H0-phi)-
     &       (geos0*gS0*rho**(geos0-1.0d0)+drhoK)
c
c      write(66,*) 'rho=',rho,
c     &      ', dgrav_steady_sg=',dgrav_steady_sg
      return
      end
