c
c--------------------------------------------------------
c
      subroutine cellave(xlow,ylow,dx,dy,wl)
      implicit double precision (a-h,o-z)
      parameter(flagno=1234567890.0)
      common /fsscoef/ a00,a10,b00,b10
      common /intar/ area,xc,yc
      logical fl(5),alll,allr
      dimension x(10),y(10),xx(5),yy(5)
      dimension po(20,2)
      external fss
c   
c     # compute wl, fraction of cell that lies in left state.
c     # For initial data with two states ql and qr separated by a 
c     # discontinuity. The curve along which the discontinuity lies is
c     # specified by the function fdisc, which should return a value that
c     # is negative on the side where ql lies and positive on the qr side.
c
c     # xlow,ylow is the coordinate of the lower left corner of the cell.
c     # dx, dy are grid spacing in x and y.
c
      call mapc2p(xlow,ylow,xx(1),yy(1))
      call mapc2p(xlow,ylow+dy,xx(2),yy(2))
      call mapc2p(xlow+dx,ylow+dy,xx(3),yy(3))
      call mapc2p(xlow+dx,ylow,xx(4),yy(4))
      xx(5) = xx(1)
      yy(5) = yy(1)
c
      area0 = 0.d0
      do i=1,4
         area0 = area0+
     &           0.5d0*(yy(i)+yy(i+1))*(xx(i+1)-xx(i))
         enddo
c
      alll = .true.
      allr = .true.
      do i=1,4
         fl(i) = fdisc_local(xx(i),yy(i)) .lt. 0.d0
         alll = alll .and. fl(i)
         allr = allr .and. (.not. fl(i))
         enddo    
      fl(5) = fl(1)
c
      if (alll) then
          wl = 1.d0
          return
      endif
c
      if (allr) then
          wl = 0.d0
          return
      endif
c
      iv = 0
      do i=1,4
         if (fl(i)) then
             iv = iv+1
             x(iv) = xx(i)
             y(iv) = yy(i)
         endif
c
         if (fl(i) .neqv. fl(i+1)) then
             iv = iv+1
             a00 = xx(i)
             a10 = xx(i+1)-xx(i)
             b00 = yy(i)
             b10 = yy(i+1)-yy(i)
             z   = zeroin(0.d0,1.d0,fss,1.d-10)
             x(iv) = a00+a10*z
             y(iv) = b00+b10*z
             endif
         enddo           
c
c     # compute area:
c
      if (iv .eq. 0) then
          wl = 0.d0
          return
      endif
c
c      x(iv+1) = x(1)
c      y(iv+1) = y(1)
c      area = 0.d0
c      do i=1,iv
c         area = area+.5d0*(y(i)+y(i+1))*(x(i+1)-x(i))
c         enddo
c
      do ii=1,iv
         po(ii,1) = x(ii)
         po(ii,2) = y(ii)
         enddo
c
      po(iv+1,1) = x(1)
      po(iv+1,2) = y(1)
      po(iv+2,1) = -flagno
      po(iv+2,2) = 0.d0
c     
c     # compute area & cell center
      call cent2d(po,xc,yc,area)
c      
      wl = area/area0   
      return
c
      end
c
c --------------------------------------------------------------
c
      subroutine cent2d(p,xc,yc,area)
      implicit double precision (a-h,o-z)
      parameter(flagno=1234567890.0)
      dimension p(20,2)
      dimension subxc(20),subyc(20),subar(20)
c
c     # Compute the centroid (xc,yc) of the polygon p.
c     # Also returns the area.
c
      do i=3,20
         if (p(i+1,1) .eq. -flagno) then
             if (i .eq. 3) then
                 write(6,*) 'degenerate polygon in cent2d'
                 go to 109
             endif
             go to 200
         endif
c
c        # Compute centroid of triangle with vertices 1, i-1, and i
         subxc(i) = (p(i-1,1)+p(i,1)+p(1,1))/3.d0
         subyc(i) = (p(i-1,2)+p(i,2)+p(1,2))/3.d0
c
c        # compute area
         subar(i) = 0.5d0*((p(1,2)+p(i-1,2))*(p(i-1,1)-p(1,1))+
     &                     (p(i-1,2)+p(i,2))*(p(i,1)-p(i-1,1))+
     &                     (p(i,2)+p(1,2))*(p(1,1)-p(i,1)))
         enddo
c
      write(6,*) 'error in cent2d: too many vertices'
 109  continue
      do ii=1,10
         write(6,*) p(ii,1),p(ii,2)
         enddo
      stop
c
 200  continue
c     # take weighted average of triangle centroids to find 
c     # centroid of polygon:
      area = 0.d0
      xc   = 0.d0
      yc   = 0.d0
      do j=3,i-1
         area = area+subar(j)
         xc   = xc+subar(j)*subxc(j)
         yc   = yc+subar(j)*subyc(j)
         enddo
c
      if (area .eq. 0.d0) then
          xc = p(1,1)
          yc = p(1,2)
      else
          xc   = xc/area
          yc   = yc/area
      endif
      return
c
      end
