c
c called tsrc2_gaxis in KM code
c
c----------------------------------------------------------------
c
      subroutine src2(meqn,mbc,mx,my,
     &           xlower,ylower,dx,dy,q,maux,aux,
     &           dt,dq)
      implicit double precision(a-h,o-z)
      common /cgrav/ grav
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension dq(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension qloc(20)
      dimension psi(20)
c
c     # source terms for cylindrical symmetry in 2d Euler equations
c     # about y-axis (so x=radius)
c
c     # aux(1,i,j) stores x coordinate of cell center
c
      ndim = 2

      do j=1,my
      do i=1,mx
         do m=1,meqn
            qloc(m) = q(m,i,j)
         enddo
c
         call ctomfd(qloc,rhoa0,rhob0,vx0,vy0,
     &          rhoh0,p0,grue0,pref0,zfa0,zfb0,c20)
c
         gaxis  = -dble(ndim-1)/aux(1,i,j)
         psi(1) = dt*gaxis*qloc(1)*vx0
         psi(2) = dt*gaxis*qloc(2)*vx0
         psi(3) = dt*gaxis*(qloc(1)+qloc(2))*vx0**2
         psi(4) = dt*gaxis*(qloc(1)+qloc(2))*vx0*vy0
         psi(5) = dt*gaxis*rhoh0*vx0
         psi(6) = 0.d0
c
c        # vertical momentum
         psi(4) = psi(4) - dt*grav*(q(1,i,j)+q(2,i,j))
c        # total energy
         psi(5) = psi(5) - dt*grav*q(4,i,j)

         do m=1,meqn
            dq(m,i,j) = dq(m,i,j)+psi(m)
         enddo   
      enddo   
      enddo   


c
      go to 99
c     # source terms for gravitational force
      do j=1,my
      do i=1,mx
c        # vertical momentum
         psi(4) = -dt*grav*(q(1,i,j)+q(2,i,j))
c
c        # total energy
         psi(5) = -dt*grav*q(4,i,j)
c     
         do m=4,5
            dq(m,i,j) = dq(m,i,j)+psi(m)
         enddo
      enddo
      enddo
c
 99   return
      end
