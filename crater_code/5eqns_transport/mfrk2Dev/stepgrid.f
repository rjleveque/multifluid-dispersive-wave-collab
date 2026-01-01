c
c -------------------------------------------------------------
c
!> Take a time step on a single grid **mptr** and overwrite solution array **q**. 
!! A modified version of the clawpack routine step2 is used.
!! 
!! Return new solution **q** as well as fluxes in fm,fp and gm,gp.
!! Patch has room for ghost cells (mbc of them) around the grid.
!! Everything is the enlarged size (**mitot** by **mjtot**).
!! 
!! \param[in] mbc number of ghost cells  (= lwidth)
!! \param[in] mptr grid number (for debugging)
!! \param[in] xlow left edge of enlarged grid (including ghost cells).
!! \param[in] ylow lower edge of enlarged grid (including ghost cells).
!! \param[in] dt incoming time step
!! \param[in] dx mesh size in x direction for this grid
!! \param[in] dx mesh size in y direction for this grid
!! \param[in,out] q solution array
!! \param[out] dtnew  return suggested new time step for this grid's soln.
!! \param[out] fm fluxes on the left side of each vertical edge
!! \param[out] fp fluxes on the right side of each vertical edge
!! \param[out] gm fluxes on the lower side of each horizontal edge
!! \param[out] gp fluxes on the upper side of each horizontal edge
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,
     &                  dx,dy,nvar,xlow,ylow,time,mptr,maux,aux,istage)
c
c          
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)
      external rpn2,rpt2

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
      common /RKmethod/ coeff(5),mstage


      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)

      dimension dq(nvar,mitot,mjtot)

      logical    debug,  dump
      data       debug/.false./
      data     dump/.false./
      !data    dump/.true./

c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      tcom = time

      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i-mbc,j-mbc,(q(ivar,i,j),ivar=1,nvar) 
            !xi = xlow+i*dx
            !yj = ylow+j*dy
            !write(outunit,545) xi,yj,i,j,(q(ivar,i,j),ivar=1,nvar) 
c    .                  ,(aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,6e15.7)
! 545        format(2e14.6,2i4,6e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c      This call has been moved out to advanc
C       call b4step2(mbc,mx,my,nvar,q,
C      &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
c
      dq = 0  
      if (method(5) .eq. 1) then
         call src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q, 
     &        maux,aux,dt,dq)
      endif
c
c     # take one step on the conservation law:
c
c     CHECK ARGS TO STEP2
      call step2(mbig,nvar,maux,mbc,mx,my,
     &           q, aux,dx,dy,dt,cflgrid,
     &           dq,rpn2,rpt2,mptr)
c
c
        !write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
        !write(*,1001) mptr, node(nestlevel,mptr),cflgrid
 1001   format(' Courant # of grid', i4,
     &        ' on level', i3, ' is  ', e10.3)
c

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy

c     # following Keh-Ming storing
c     # (dt)*L(q^n) = (dt)*d/dt(q^n(t^n) in qrk
c     # KM calls it dq in step2, qrk in calling prog
c     # NB that KM also updates ghost cells, stepgrid 
c     # usually doesnt
        do  j=1,mjtot
        do  i=1,mitot
        do  m=1,nvar
           q(m,i,j) = q(m,i,j) + dq(m,i,j)
        end do 
        end do 
        end do 
c
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," after stepgrid"
         do i = mbc+1, mitot-mbc
         do j = mbc+1, mjtot-mbc
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
            xi = xlow+i*dx
            yj = ylow+j*dy
            !write(outunit,545) xi,yj,i,j,(q(ivar,i,j),ivar=1,nvar) 
         end do
         end do
      endif

c     test for negative densities after each stage
      if (.true.) then
         do i = mbc+1, mitot-mbc
         do j = mbc+1, mjtot-mbc
            if (q(1,i,j) .lt. 0.d0 .or. q(2,i,j) .lt. 0.d0) then 
               write(outunit,*) "Grid ",mptr," level ",level,
     &                " after stepgrid istage ",istage
               write(*,*) "Grid ",mptr," level ",level,
     &                " after stepgrid istage ",istage
               write(outunit,*)"Found negative density on grid ",mptr,
     &                " time ",time," END of stage ",istage
               write(*,*)"Found negative density on grid ",mptr,
     &                " time ",time," END of stage ",istage
               write(outunit,900) i,j,(q(m,i,j),m=1,meqn)
               write(*,900) i,j,(q(m,i,j),m=1,meqn)
 900          format(2i5,6e25.15)
            endif
         end do
         end do
      endif

      return
      end


