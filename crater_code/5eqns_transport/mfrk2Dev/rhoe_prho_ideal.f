c
c ------------------------------------------------------------
c
      double precision function rhoe_prho_ideal(p,geos0)  
      implicit double precision (a-h,o-z)
c
c     # compute phasic internal energy
c     # ideal gas EOS
c
      rhoe_prho_ideal = p/(geos0-1.d0)
      return
c
      end
