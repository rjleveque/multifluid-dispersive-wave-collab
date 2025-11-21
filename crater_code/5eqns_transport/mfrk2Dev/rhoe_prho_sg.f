c
c ------------------------------------------------------------
c
      double precision function rhoe_prho_sg(p,geos0,bn0)
      implicit double precision (a-h,o-z)
c
c     # compute specific internal energy as function of 
c     # density and pressure     
c     # stiffened gas EOS
c
      rhoe_prho_sg = (p+geos0*bn0)/(geos0-1.d0)
      return
c
      end
