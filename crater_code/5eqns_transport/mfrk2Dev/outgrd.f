c
c ------------------------------------------------------------
c
      subroutine outgrd(mbc,mx,my,xlower,ylower,
     &           dx,dy,iframe,aux,maux)

      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      character*10 fname
c
c     # output physical grid system
c
c      fname = 'fort.pxxxx'
c      nstp = iframe
c      do ipos=10,7,-1
c         idigit = mod(nstp,10)
c         fname(ipos:ipos) = char(ichar('0') + idigit)
c         nstp = nstp/10
c         enddo
c
c      open(unit=50,file=fname,status='unknown',form='formatted')
c
c     # cell vertices
c
c      write(50,*) mx+1,my+1
c
c      do j=1,my+1
c         do i=1,mx+1
c            write(50,1001) aux(i,j,10),aux(i,j,11)
c            enddo
c         enddo
c
c      close(unit=50)
c
c     # cell center
c
      fname = 'fort.gxxxx'
      nstp = iframe
      do ipos=10,7,-1
         idigit = mod(nstp,10)
         fname(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp/10
         enddo
c
      open(unit=50,file=fname,status='unknown',form='formatted')
c
      write(50,*) mx,my
c
      do j=1,my
         yc = ylower+dble(j-1)*dy+0.5d0*dy
         do i=1,mx
            xc = xlower+dble(i-1)*dx+0.5d0*dx
c
            write(50,1001) xc,yc
            enddo
         enddo
c
      close(unit=50)
      return
c
 1001 format(2f26.12)
      end
