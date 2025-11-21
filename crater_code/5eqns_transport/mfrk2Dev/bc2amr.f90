! :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
!> \callgraph
!! \callergraph
!!  Take a grid patch with mesh widths **hx**,**hy**, of dimensions **nrow** by
!!  **ncol**,  and set the values of any piece of
!!  of the patch which extends outside the physical domain 
!!  using the boundary conditions. 
!!
!!  ### Standard boundary condition choices for amr2ez in clawpack
!!
!!  At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
!!
!!  mthbc(k) =  
!!  * 0  for user-supplied BC's (must be inserted!)
!!  * 1  for zero-order extrapolation
!!  * 2  for periodic boundary conditions
!!  * 3  for solid walls, assuming this can be implemented
!!                   by reflecting the data about the boundary and then
!!                   negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!!                   component of q.
!!  * 4  for sphere bcs (left half maps to right half of same side, and vice versa), as if domain folded in half
!!
!!  The corners of the grid patch are at 
!!     (xlo_patch,ylo_patch)  --  lower left corner
!!     (xhi_patch,yhi_patch) --  upper right corner
!!
!!  The physical domain itself is a rectangle bounded by
!!     (xlower,ylower)  -- lower left corner
!!     (xupper,yupper)  -- upper right corner
!!  
!   This figure below does not work with doxygen
!   the picture is the following: 
!  ____________________________________________________
! 
!                _____________________ (xupper,yupper)
!               |                     |  
!           ____|____ (xhi_patch,yhi_patch)   
!           |   |    |                |
!           |   |    |                |
!           |   |    |                |
!           |___|____|                |
!  (xlo_patch,ylo_patch) |            |
!               |                     |
!               |_____________________|   
!    (xlower,ylower)
!  ____________________________________________________
!!
!!
!>  Any cells that lie outside the physical domain are ghost cells whose
!!  values should be set in this routine.  This is tested for by comparing
!!  xlo_patch with xlower to see if values need to be set at the left
!   as in the figure above, 
!
!>  and similarly at the other boundaries.
!!  Patches are guaranteed to have at least 1 row of cells filled
!!  with interior values so it is possible to extrapolate. 
!!  Fix [trimbd()](@ref trimbd) if you want more than 1 row pre-set.
!!
!!  Make sure the order the boundaries are specified is correct
!!  so that diagonal corner cells are also properly taken care of.
!!
!!  Periodic boundaries are set before calling this routine, so if the
!!  domain is periodic in one direction only you
!!  can safely extrapolate in the other direction. 
!!
!!  Don't overwrite ghost cells in periodic directions!
!!
!! \param val data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables 
!! \param nrow number of cells in *i* direction on this grid
!! \param ncol number of cells in *j* direction on this grid
!! \param meqn number of equations for the system
!! \param naux number of auxiliary variables
!! \param hx spacing (mesh size) in *i* direction
!! \param hy spacing (mesh size) in *j* direction
!! \param level AMR level of this grid
!! \param time setting ghost cell values at time **time**
!! \param xlo_patch left bound of the input grid
!! \param xhi_patch right bound of the input grid 
!! \param ylo_patch lower bound of the input grid 
!! \param yhi_patch upper bound of the input grid 
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

subroutine bc2amr(val,aux,nrow,ncol,meqn,naux, hx, hy, level, time,   &
                  xlo_patch, xhi_patch, ylo_patch, yhi_patch) 

    use amr_module, only: mthbc, xlower, ylower, xupper, yupper
    use amr_module, only: xperdom,yperdom,spheredom

    implicit none

    ! Input/Output
    integer, intent(in) :: nrow, ncol, meqn, naux, level
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in) :: ylo_patch, yhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow, ncol)
    real(kind=8), intent(in out) :: aux(naux, nrow, ncol)
    
    ! Local storage
    integer :: i, j, ibeg, jbeg, nxl, nxr, nyb, nyt
    real(kind=8) :: hxmarg, hymarg
    real(kind=8) :: qbc(20),zfa0,zfb0,Y20,Y10,w0,vbc,q0,v0
    real(kind=8) :: u0,ubc,s0,rhobc,rhob0,rhoa0,rho0,rho2_bc,rho1_bc
    real(kind=8) :: rho00,phibc,phi0,pbc,p0,H0,gS0,geos0,bn0
    real(kind=8) :: zero1d_aux, grue0,rhoh0,pref0,c20
    integer :: jbc,m,iflag
    real(kind=8) :: rhoeps
    parameter(rhoeps=1.0d-8)

    !common /grav_steady_info/ q0,ubc,H0,gS0,phibc,geos0

    external grav_steady_sg, dgrav_steady_sg

    hxmarg = hx * .01d0
    hymarg = hy * .01d0

    ! Use periodic boundary condition specialized code only, if only one 
    ! boundary is periodic we still proceed below
    if (xperdom .and. (yperdom .or. spheredom)) then
        return
    end if

    ! Each check has an initial check to ensure that the boundary is a real
    ! boundary condition and otherwise skips the code.  Otherwise 
    !-------------------------------------------------------
    ! Left boundary:
    !-------------------------------------------------------
    if (xlo_patch < xlower-hxmarg) then
        ! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch) / hx)

        select case(mthbc(1))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do j = 1, ncol
                    do i=1, nxl
                        val(:, i, j) = val(:, nxl + 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = 1, ncol
                    do i=1, nxl
                        val(:, i, j) = val(:, 2 * nxl + 1 - i, j)
                    end do
                end do
                ! negate the normal velocity: 
                ! THIS IS component 3 in mfluid
                do j = 1, ncol
                    do i=1, nxl
                        val(3, i, j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop
        end select
    end if

    !-------------------------------------------------------
    ! Right boundary:
    !-------------------------------------------------------
    if (xhi_patch > xupper+hxmarg) then

        ! number of grid cells lying outside physical domain:
        nxr = int((xhi_patch - xupper + hxmarg) / hx)
        ibeg = max(nrow - nxr + 1, 1)

        select case(mthbc(2))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do i = ibeg, nrow
                    do j = 1, ncol
                        val(:, i, j) = val(:, ibeg - 1, j)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do i=ibeg, nrow
                    do j = 1, ncol
                        val(:, i, j) = val(:, 2 * ibeg - 1 - i, j)
                    end do
                end do
                ! negate the normal velocity:
                ! THIS IS Component 3 in Mfluid
                do i = ibeg, nrow
                    do j = 1, ncol
                        val(3, i, j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Bottom boundary:
    !-------------------------------------------------------
    if (ylo_patch < ylower - hymarg) then

        ! number of grid cells lying outside physical domain:
        nyb = int((ylower + hymarg - ylo_patch) / hy)

        select case(mthbc(3))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                write(*,*)" Copy bottom bc from MFLUID after all"
                stop "A user defined boundary condition was not provided."
            
            case(1) ! Zero-order extrapolation
                do j = 1, nyb
                    do i = 1, nrow
                        val(:,i,j) = val(:, i, nyb + 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = 1, nyb
                    do i = 1, nrow
                        val(:,i,j) = val(:, i, 2 * nyb + 1 - j)
                    end do
                end do
                ! negate the normal velocity:
                ! THIS IS Component 4 for MFLUID
                do j = 1, nyb
                    do i = 1, nrow
                        val(4,i,j) = -val(4, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Top boundary:
    !-------------------------------------------------------
    if (yhi_patch > yupper + hymarg) then

        ! number of grid cells lying outside physical domain:
        nyt = int((yhi_patch - yupper + hymarg) / hy)
        jbeg = max(ncol - nyt + 1, 1)

        select case(mthbc(4))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                !stop "A user defined boundary condition was not provided."
                ! from Keh-Ming code, 
                ! # gravitational source terms
                ! # steady-state equilibrium conditions
            ! loop to  put 3 last j ghost cells at top  into qbc for that i
            do j = jbeg, ncol
            do i = 1, nrow
              qbc(1:meqn) = val(:,i,2*jbeg-1-j)

              call ctomfd(qbc,rhoa0,rhob0,u0,v0,rhoh0,p0,       &
                           grue0,pref0,zfa0,zfb0,c20)

              ! EOS constants
              geos0 = grue0 + 1.0d0
              bn0 = -pref0*grue0/(grue0+1.0d0)

!             # total density
              rho0 = rhoa0+rhob0
!
!             # mass fraction
              Y10  = rhoa0/rho0
              Y20  = rhob0/rho0
!
!             # constant momentum
              q0 = rho0*v0
              w0 = rho0*u0*v0
!
!             # constant entropy
              S0 = (p0+bn0)/rho0**geos0
!
!             # constant modified enthalpy
              phi0 = aux(2,i,2*jbeg-1-j)  
              !H0   = (qbc(5)+p0)/rho0+phi0
              H0   = rhoh0/rho0 + phi0
!
!             # transverse velocity
              ubc = u0
!
              !gS0   = geos0/(geos0-1.0d0)*S0
              gS0   = geos0/grue0 * S0
              phibc = aux(2,i,2*jbeg-1-j+1)
!
!             # solve for density
              rho00 = rho0
              rhobc = zero1d_aux(rho00,grav_steady_sg,     &
                          dgrav_steady_sg,rhoeps,iflag,    &
                          q0,ubc,H0,gs0,phibc,geos0)
              if (iflag .ne. 1) then
!                  write(66,*) 'error in top bc: time=',t,
!     &                         ', iflag=',iflag
!                  write(66,*) 'qbc=',(qbc(m),m=1,meqn)
!                  write(66,*) 'rho0=',rho0,u0,v0,p0
!
!                  stop
              else
                 vbc = q0/rhobc
                 pbc = S0*rhobc**geos0-bn0
!
                 rho1_bc = Y10*rhobc/zfa0
                 rho2_bc = Y20*rhobc/zfb0
!
                 call prmtoc_local(qbc,rho1_bc,rho2_bc,ubc,vbc,    &
                      pbc,pbc,zfa0,zfb0)
              endif
!
              val(:,i,j) = qbc(1:meqn)
            end do
            end do

            case(1) ! Zero-order extrapolation
                do j = jbeg, ncol
                    do i = 1, nrow
                        val(:, i, j) = val(:, i, jbeg - 1)
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do j = jbeg, ncol 
                    do i = 1, nrow
                        val(:, i, j) = val(:, i, 2 * jbeg - 1 - j)
                    end do
                end do
                ! negate the normal velocity:
                do j = jbeg, ncol
                    do i = 1, nrow
                        val(3, i, j) = -val(3, i, j)
                    end do
                end do

            case(4) ! Spherical domain
                continue

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

end subroutine bc2amr

