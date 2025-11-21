c
c-----------------------------------------------------------------------
c
      subroutine estdt(val,mitot,mjtot,nvar,dx,dy,dt,nghost,aux,
     &                 naux,cfl)
c
c :::::::::::::::::::::::: ESTDT :::::::::::::::::::::::::::;
c  estimate the initial time step for the given values
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

       implicit double precision (a-h, o-z)
       dimension val(nvar,mitot,mjtot)
       dimension aux(naux,mitot,mjtot)
       dimension qloc(nvar)
c
c     # set initial variable time step
c
      cmax = 0.d0
c
      !write(234,*)"from estdt"
      do j=1, mjtot
      do i=1, mitot
            qloc(:) = val(:,i,j)
c            
            call ctomfd(qloc,rhoa0,rhob0,u0,v0,
     &           rhoh0,p0,grue0,zeta0,zfa0,zfb0,c20)
c
            cmax  = dmax1(cmax,dabs(u0)+dsqrt(c20),
     &                         dabs(v0)+dsqrt(c20))
      !write(47,234) i-nghost,j-nghost,c20,u0,v0,cmax
 234  format(2i5,4e15.7)
      enddo
      enddo
      !stop
      go to 999  ! this is probably leftover from looping over grids
c      
 999  continue
      if (cmax .ne. 0.d0) dt = cfl*dmin1(dx,dy)/cmax
      return
c
      end
