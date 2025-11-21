subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)


!> Called to update q by solving source term equation 
!! $q_t = \psi(q)$ over time dt starting at time t.
!!
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8) :: qloc(meqn), psi(meqn)
    real(kind=8) :: rhoa0,rhob0,vx0,vy0,rhoh0,p0,grue0
    real(kind=8) :: pref0,zeta0,zfa0,zfb0,c20,gaxis,xc,dt2
    integer :: i,j,m,ndim
    real(kind=8) :: grav
    common /cgrav/ grav

! modified KehMings routine in file tsrc2_gaxis.f for claw5 format
!     # source terms for cylindrical symmetry in 2d Euler equations
!     # about y-axis (so x=radius)
!
!     # aux(1,i,j) stores x coordinate of cell center
!
      ndim = 2
      dt2  = dt/2.0d0

      do j = 1,my
      do i = 1,mx
         xc = xlower+dble(i-1)*dx+0.5d0*dx
         qloc(1:meqn) = q(:,i,j)
!
         call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,         &
                 rhoh0,p0,grue0,pref0,zfa0,zfb0,c20)
!
         gaxis  = -dble(ndim-1)/xc
         psi(1) = gaxis*qloc(1)*vx0
         psi(2) = gaxis*qloc(2)*vx0
         psi(3) = gaxis*qloc(3)*vx0
         psi(4) = gaxis*qloc(4)*vx0 - grav*(qloc(1)+qloc(2))
         psi(5) = gaxis*rhoh0*vx0 - grav*qloc(4)
         psi(6) = 0.d0
!
         do m = 1, meqn
            qloc(m) = qloc(m)+dt2*psi(m)
         end do
!
         call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,rhoh0,p0,   &
                    grue0,pref0,zfa0,zfb0,c20)

         psi(1) = gaxis*qloc(1)*vx0
         psi(2) = gaxis*qloc(2)*vx0
         psi(3) = gaxis*qloc(3)*vx0
         psi(4) = gaxis*qloc(4)*vx0 - grav*(qloc(1)+qloc(2))
         psi(5) = gaxis*rhoh0*vx0 - grav*qloc(4)
         psi(6) = 0.0d0

         do m=1,meqn
!           # midpoint rule
            q(m,i,j) = q(m,i,j)+dt*psi(m)
         end do

      end do
      end do

     return
end subroutine src2
