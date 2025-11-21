      subroutine setprob

      use amr_module
      implicit real*8 (a-h,o-z)
      character(len=25) fname
      logical nolimiter,quad
c
      common /RKmethod/ coeff(5), mstage
      common  /order2/ ssw,quad,nolimiter
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) mstage
      write(*,*) "will use RK scheme with ", mstage," stages"

      coeff = 0.d0 
      if (mstage .eq. 1) then
         coeff(1) = 1.d0
      else if (mstage .eq. 2) then
         coeff(1) = 0.5d0
         coeff(2) = 1.0d0
      endif

      gamma = 1.4d0
      gamma1 = gamma - 1.d0
      xprob = xupper
      yprob = yupper
      cflcart = cfl  ! have to go through and only use one

      call setprob_local()


      return
      end
