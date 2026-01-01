!> Compute all fluxes at cell edges 
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge

subroutine step2(maxm,meqn,maux,mbc,mx,my,q,aux,dx,dy,dt,    &
                 cflgrid,dq,rpn2,rpt2,mptr)

!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!        q has initial data for this step
!        and is unchanged in this version.
!    
    
    use amr_module

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my,mptr
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: dq(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    
    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) ::  dq1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux1(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux3(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)
    
    !real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    !real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    !real(kind=8) ::  amdq(meqn,1-mbc:maxm + mbc)
    !real(kind=8) ::  apdq(meqn,1-mbc:maxm + mbc)
    !real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
    !real(kind=8) :: bmadq(meqn,1-mbc:maxm + mbc)
    !real(kind=8) :: bpadq(meqn,1-mbc:maxm + mbc)
    
    ! Looping scalar storage
    integer :: i,j,m
    real(kind=8) :: dtdx,dtdy,cfl1d
    
    ! Common block storage
    integer :: icom,jcom
    real(kind=8) :: dtcom,dxcom,dycom,tcom
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
    
    ! Store mesh parameters in common block
    dxcom = dx
    dycom = dy
    dtcom = dt
    
    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

    ! dq will accumulate src term plus stage updates for q
    ! call to src2 moved to stepgrid where it usually is
    !dq = 0  
    !if (method(5) .eq. 1) then
    !   call src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,       &
    !        maux,aux,dt,dq)
    !endif
    
    
    ! ============================================================================
    ! Perform X-Sweeps
    do j = 0,my+1
        ! Copy old q into 1d slice
        q1d(:,1-mbc:mx+mbc) = q(:,1-mbc:mx+mbc,j)
        
        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif
        
        ! Copy aux array into slices
        if (maux > 0) then
            aux1(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j-1)
            aux2(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j  )
            aux3(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j+1)
        endif
        
        ! Store value of j along the slice into common block
        jcom = j

        ! Compute modifications fadd and gadd to fluxes along this slice:
        call flux2(1,maxm,meqn,maux,mbc,mx,q1d,dtdx1d,aux1,aux2,aux3,     &
                   dq1d,cfl1d,rpn2,rpt2)       
                   
        cflgrid = max(cflgrid,cfl1d)

        ! put update into dq as in KM
        do i = 1,mx
        do m = 1, meqn
           dq(m,i,j) = dq(m,i,j) + dq1d(m,i)
        end do
        end do

        
    enddo

    ! ============================================================================
    !  y-sweeps    
    !
    do i = 0,mx+1
        
        ! Copy data along a slice into 1d arrays:
        q1d(:,1-mbc:my+mbc) = q(:,i,1-mbc:my+mbc)

        ! Set dt/dy ratio in slice
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! Copy aux slices
        if (maux .gt. 0)  then
            aux1(:,1-mbc:my+mbc) = aux(:,i-1,1-mbc:my+mbc)
            aux2(:,1-mbc:my+mbc) = aux(:,i,1-mbc:my+mbc)
            aux3(:,1-mbc:my+mbc) = aux(:,i+1,1-mbc:my+mbc)
        endif
        
        ! Store the value of i along this slice in the common block
        icom = i
        
        ! Compute modifications fadd and gadd to fluxes along this slice
        call flux2(2,maxm,meqn,maux,mbc,my,q1d,dtdy1d,aux1,aux2,aux3, &
                   dq1d,cfl1d,rpn2,rpt2)

        cflgrid = max(cflgrid,cfl1d)

        ! put update into dq as in KM
        do j = 1,my
        do m = 1, meqn
           dq(m,i,j) = dq(m,i,j) + dq1d(m,j)
        end do
        end do

    enddo
    return
end subroutine step2
