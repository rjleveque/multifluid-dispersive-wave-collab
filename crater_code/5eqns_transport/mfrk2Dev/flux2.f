c
c ---------------------------------------------------------------
c
      subroutine flux2(ixy,maxm,meqn,maux,mbc,mx,
     &           q1d,dtdx1d,aux1,aux2,aux3,
     &           dq1d,cfl1d,rpn2,rpt2)
     
      use amr_module, only: mthlim, method,mwaves,mthbc

      implicit double precision (a-h,o-z)

      dimension q1d(meqn,1-mbc:maxm+mbc)
      dimension dq1d(meqn,1-mbc:maxm+mbc)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc)
      dimension aux2(maux,1-mbc:maxm+mbc)
      dimension aux3(maux,1-mbc:maxm+mbc)

      ! below moved to local storage, removed from parameters
      dimension ql(meqn,1-mbc:maxm+mbc)
      dimension qr(meqn,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension amdq2(meqn,1-mbc:maxm+mbc)
      dimension apdq2(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
      ! not changing order of dimensions below for
      ! now since only talks to KM routines
      !dimension uu(1-mbc:maxm+mbc,meqn,2)
      !dimension hh(1-mbc:maxm+mbc,-2:2)
      !dimension evl(1-mbc:maxm+mbc,meqn,mwaves)
      !dimension evr(1-mbc:maxm+mbc,meqn,mwaves)
      dimension s(1-mbc:maxm+mbc,mwaves)
      dimension wave(1-mbc:maxm+mbc,meqn,mwaves)

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c     # Evaluate (delta t) * dq/dt along a 1D slice
c     #    ixy = 1  if it is a slice in x
c     #          2  if it is a slice in y
c     # This value is passed into the Riemann solvers. The scaled fluctuations
c     # go into the array dq1d.  The notation is written assuming
c     # we are solving along a 1D slice in the x-direction.
c
c     # dq1d(i,.) modifies Q of cell i
c
c     # No transverse waves or genuinely multi-d reconstruction
c
c     Note that if method(6)=1 then the capa array comes into the second 
c     order correction terms, and is already included in dtdx1d:
c     If ixy = 1 then
c        dtdx1d(i) = dt/dx                 if method(6) = 0
c                  = dt/(dx*capa(i,jcom))  if method(6) = 1
c     If ixy = 2 then
c        dtdx1d(j) = dt/dy                 if method(6) = 0
c                  = dt/(dy*capa(icom,j))  if method(6) = 1
c
c     Notation:
c        The jump in q (q1d(i,:)-q1d(i-1,:))  is split by rpn2 into
c            amdq =  the left-going flux difference  A^- Delta q  
c            apdq = the right-going flux difference  A^+ Delta q  
c
c     # initialize dq,apdq/amdq
      do i=1-mbc,mx+mbc
      do m=1,meqn
         dq1d(m,i) = 0.d0
      enddo    
      enddo    
c
      go to (9,10,20) method(9)+1
c
 9    continue
!c     # piecewise constant reconstruction
!c
!      call q2qlqr(ixy,maxm,meqn,mbc,maux,mx,dt,dx,
!     &     q1d,aux1,aux2,aux3,ql,qr,auxl,auxr,amdq,1)
      go to 90
c
 10   continue
c     # higher order reconstruction
c
      ! amdq, apdq called qp,dq, in interp routine
      call q2qlqr_interp(ixy,maxm,meqn,mbc,maux,mx,
     &     dt,dx,q1d,aux2,ql,qr,auxl,auxr,s,wave,
     &     amdq,apdq)
!    !&     uu,hh,evl,evr)
c
c     # set boundary condition 

!!!!! dont call bc routine since should have enough ghost cells 
!!!! already set
!     call bc2_lr(ixy,maxm,meqn,mbc,mx,xlower,dx,ql,qr,
!    &     maux,auxl,auxr,t,dt)
c
!! routine below just returned, did nothing, didnt copy
!!      call auxbc2_lr(ixy,maxm,mbc,mx,xlower,dx,
!!     &     maux,auxl,auxr,t,dt,mthbc)
      go to 90
c
 20   continue
!c     # hybrid THINC & other reconstruction
!c
!c     # reconstruction scheme (basis)
!      call q2qlqr_interp(ixy,maxm,meqn,mwaves,mbc,maux,mx,
!     &     dt,dx,q1d,aux2,ql,qr,auxl,auxr,s,wave,amdq,apdq,
!     &     uu,hh,evl,evr,mthlim)
!c
!c     # THINC reconstruction (interface-sharpening)
!      call q2qlqr(ixy,maxm,meqn,mbc,maux,mx,dt,dx,
!     &     q1d,aux1,aux2,aux3,ql,qr,auxl,auxr,amdq,0)
!c
!c     # set boundary condition 
!      call bc2_lr(ixy,maxm,meqn,mbc,mx,xlower,dx,ql,qr,
!     &     maux,auxl,auxr,t,dt,mthbc)
!c
!      call auxbc2_lr(ixy,maxm,mbc,mx,xlower,dx,
!     &     maux,auxl,auxr,t,dt,mthbc)
!      go to 90 
!c
 90   continue
c     # solve Riemann problem at each interface and compute updates
      call rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &     auxl,auxr,wave,s,amdq,apdq)
c
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.0d0
      do mw=1,mwaves
         do i=1,mx+1
c           # if s>0 use dtdx1d(i) to compute CFL,
c           # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(i,mw), 
     &                          -dtdx1d(i-1)*s(i,mw))
            enddo
         enddo 
c
c     # find total fluctuation within each cell
      if (method(8) .eq. 0) then
c         # ignore fluctutation within cells
          do m=1,meqn
             do i=1,mx
                dq1d(m,i) = dq1d(m,i)-dtdx1d(i)*(apdq(m,i)+amdq(m,i+1))
                enddo
             enddo
      else if (method(8) .eq. 1) then
c         # tfluct2 should be a special solver that uses the parameters aux(i)
c         # to solve a Riemann problem with left state ql(i)
c         # and right state qr(i), and returns a total fluctuation in amdq2
c         # note that amdq2 is a total fluctuation and should be called
c         # adq; we do it this way to save on storage
c
          call tfluct2(ixy,maxm,meqn,mwaves,mbc,maux,mx,
     &         ql,qr,aux2,aux2,wave,s,amdq2)
c
          do m=1,meqn
             do i=1,mx
                dq1d(m,i) = dq1d(m,i)-dtdx1d(i)*(apdq(m,i)+amdq2(m,i)+
     &                      amdq(m,i+1))
                enddo
             enddo
      elseif (method(8) .eq. 2) then
c         # swap things around and use the usual Riemann solver
c         # This may be more convenient, but is less efficient.
c
          do i=2-mbc,mx+mbc
             do m=1,meqn
                qr(m,i-1) = ql(m,i)
                ql(m,i)   = qr(m,i)
                enddo
c
             do ma=1,maux
                auxr(ma,i-1) = aux2(ma,i)
                auxl(ma,i)   = aux2(ma,i)
                enddo
             enddo
c
          call rpn2(ixy,maxm,meqn,mwaves,mbc,maux,mx,ql,qr,
     &         auxl,auxr,wave,s,amdq2,apdq2)
c
          do m=1,meqn
             do i=1,mx+1
                dq1d(m,i) = dq1d(m,i)-
     &                dtdx1d(i)*(apdq(m,i)+apdq2(m,i)+
     &                           amdq(m,i+1)+amdq2(m,i))
                enddo
             enddo
      endif
c
      return
      end
