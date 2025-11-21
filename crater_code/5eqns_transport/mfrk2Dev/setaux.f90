!
! -----------------------------------------------------------
!
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
!
!     # radial distance
      do j=1-mbc,my+mbc
         do i=1-mbc,mx+mbc
            xc = xlower+dble(i-1)*dx+0.5d0*dx
            aux(1,i,j) = xc
         enddo
      enddo

!     # set auxiliary state
!     # gravitational potential
!     # psi = g y
!
      grav = 9.81d0
      do j=1-mbc,my+mbc
         yc = ylower+dble(j-1)*dy+0.5d0*dy
         do i=1-mbc,mx+mbc
            aux(2,i,j) = grav*yc
         enddo
      enddo
      return
!
      end
