c
c ---------------------------------------------------------------------
c
      subroutine state2_rotate(rot,velcomps)
      implicit none
      double precision rot(4),velcomps(20),v1,v2
c
      v1 = velcomps(3)
      v2 = velcomps(4)
c
      velcomps(3) = rot(1)*v1+rot(2)*v2
      velcomps(4) = rot(3)*v1+rot(4)*v2
      return
c
      end
c
c ----------------------------------------------------------
c
      subroutine state2_rotate_tr(rot,velcomps)
      implicit none
      double precision velcomps(20),rot(4), v1, v2
c
      v1 = velcomps(3)
      v2 = velcomps(4)
c
      velcomps(3) = rot(1)*v1+rot(3)*v2
      velcomps(4) = rot(2)*v1+rot(4)*v2
      return
c
      end
