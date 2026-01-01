!> Called before each call to step2.
!! Use to set time-dependent aux arrays or perform other tasks.
!!
!! This default version does nothing. 
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,mptr,level,istage)

 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux,mptr,istage,level
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer :: i,j,m

    logical  IS_GHOST

    IS_GHOST(i,j) = (i .lt. 1 .or. i .gt. mx .or.   &
                     j .lt. 1 .or. j .gt. my)

    do i = 1-mbc, mx+mbc
    do j = 1-mbc, my+mbc
      if (q(1,i,j) .lt. 0.d0 .or. q(2,i,j) .lt. 0.d0) then
         write(*,*)"Found negative density on grid ",mptr," level ",level,   &
                   " time ",t," istage ",istage
         write(*,*)"Grid of size mx,my = ",mx,my
         write(*,900) i,j,(q(m,i,j),m=1,meqn)
 900     format(2i5,6e25.15)
         if (IS_GHOST(i,j)) then
           write(*,*)"This is a ghost cell"
         endif
      endif
    end do
    end do

end subroutine b4step2
