c
c ----------------------------------------------------------
c
      subroutine bc2_lr(ixy,maxm,meqn,mbc,mx,xlower,dx,
     &           ql,qr,maux,auxl,auxr,t,dt,mthbc)
      implicit double precision (a-h,o-z)
      dimension ql(1-mbc:maxm+mbc,meqn)
      dimension qr(1-mbc:maxm+mbc,meqn)
      dimension auxl(1-mbc:maxm+mbc,maux)
      dimension auxr(1-mbc:maxm+mbc,maux)
      dimension mthbc(4)
c
c     # set cell-edge boundary condition
c
      go to (1,2) ixy
c
 1    continue
c     # x-sweep
c
c     # left boundary:
      go to (100,110,120,130) mthbc(1)+1
      goto 199
c
 100  continue
      go to 199
c
 110  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            ql(1-ibc,m) = ql(1,m)
            qr(1-ibc,m) = qr(1,m)
            enddo
         enddo
      go to 199

 120  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            ql(1-ibc,m) = ql(mx+1-ibc,m)
            qr(1-ibc,m) = qr(mx+1-ibc,m)
            enddo
         enddo
      go to 199
c
 130  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            ql(1-ibc,m) = ql(ibc,m)
            qr(1-ibc,m) = qr(ibc,m)
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         ql(1-ibc,3) = -ql(1-ibc,3)
         qr(1-ibc,3) = -qr(1-ibc,3)
         enddo
      go to 199
c
 199  continue
c     # right boundary:
      go to (200,210,220,230) mthbc(2)+1
      goto 299
c
 200  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

 210  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            ql(mx+ibc,m) = ql(mx,m)
            qr(mx+ibc,m) = qr(mx,m)
            enddo
         enddo
      go to 299

 220  continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            ql(mx+ibc,m) = ql(ibc,m)
            qr(mx+ibc,m) = qr(ibc,m)
            enddo
         enddo
      go to 299

 230  continue
c     # solid wall
      do m=1,meqn
         do ibc=1,mbc
            ql(mx+ibc,m) = ql(mx+1-ibc,m)
            qr(mx+ibc,m) = qr(mx+1-ibc,m)
            enddo
         enddo
c
c     # negate the normal velocity:
      do ibc=1,mbc
         ql(mx+ibc,3) = -ql(mx+ibc,3)
         qr(mx+ibc,3) = -qr(mx+ibc,3)
         enddo
      go to 299

 299  continue
      return
c
 2    continue
c     # y-sweep
c     # bottom boundary:
      go to (300,310,320,330) mthbc(3)+1
      goto 399
c
 300  continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
 310  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            ql(1-jbc,m) = ql(1,m)
            qr(1-jbc,m) = qr(1,m)
            enddo
         enddo
      go to 399
c
 320  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            ql(1-jbc,m) = ql(mx+1-jbc,m)
            qr(1-jbc,m) = qr(mx+1-jbc,m)
            enddo
         enddo
      go to 399

 330  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            ql(1-jbc,m) = ql(jbc,m)
            qr(1-jbc,m) = qr(jbc,m)
            enddo
         enddo
c
c     # negate the normal velocity:
      do jbc=1,mbc
         ql(1-jbc,4) = -ql(1-jbc,4)
         qr(1-jbc,4) = -qr(1-jbc,4)
         enddo    
      go to 399

 399  continue
c     # top boundary:
      go to (400,410,420,430) mthbc(4)+1
      goto 499
c
 400  continue
      go to 499

 410  continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            ql(mx+jbc,m) = ql(mx,m)
            qr(mx+jbc,m) = qr(mx,m)
            enddo
         enddo
      go to 499
c
 420  continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            ql(mx+jbc,m) = ql(jbc,m)
            qr(mx+jbc,m) = qr(jbc,m)
            enddo
         enddo
      go to 499
c
 430  continue
c     # solid wall
      do m=1,meqn
         do jbc=1,mbc
            ql(mx+jbc,m) = ql(mx+1-jbc,m)
            qr(mx+jbc,m) = qr(mx+1-jbc,m)
            enddo
         enddo
c
c     # negate the normal velocity:
      do jbc=1,mbc
         ql(mx+jbc,4) = -ql(mx+jbc,4)
         qr(mx+jbc,4) = -qr(mx+jbc,4)
         enddo
      go to 499
c
 499  continue
      return
c
      end
