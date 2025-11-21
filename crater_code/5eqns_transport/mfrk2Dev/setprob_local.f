c
c -----------------------------------------------------------------
c
      subroutine setprob_local

      use amr_module, only : method

      implicit double precision (a-h,o-z)
      common /loops/ xloop(10),yloop(10),zloop(10),rloop(10),
     &        axloop(10),bxloop(10),inbloop(2,2,10)
      common /mfluid/ geos(10),an(10),bn(10),cn(10),
     &        dn(10),en(10),rref(10),c2ref(10),tref(10),
     &        cv(10),emu(10),rhocav(10),pcav(10),
     &        ccav(10),rhosat(10),psat(10),csat(10),
     &        rhovm1(10),rhovm2(10),psi(10),zf(10),zv(10),
     &        mwoods,mphase
      common /phypar/ den0,velx0,vely0,velz0,pr0
      common /compi/ pi
      common /cgrav/ grav
      common /crater_info/ RC,RD,DC,r1,rdecay,linfactor
      common /crater_reg/ smoothed
      common /thinc0/ thinc_beta,thinc_eta,gp1,gp2,vofmin,
     &        intf_test1,intf_test2
      logical smoothed
c
      pi = 4.d0*datan(1.d0)
c
      open(unit=7,file='setprob_local.txt',status='old',
     &     form='formatted')
c
      read(7,*) nlayer        
      read(7,*) (geos(i),i=1,nlayer)
      read(7,*) (an(i),i=1,nlayer)
      read(7,*) (bn(i),i=1,nlayer)
      read(7,*) (cn(i),i=1,nlayer)
      read(7,*) (dn(i),i=1,nlayer)
      read(7,*) (en(i),i=1,nlayer)
      read(7,*) (cv(i),i=1,nlayer)
      read(7,*) (rref(i),i=1,nlayer)
      read(7,*) (tref(i),i=1,nlayer)
      read(7,*) (pcav(i),i=1,nlayer)
      read(7,*) (zf(i),i=1,nlayer)
      read(7,*) (zv(i),i=1,nlayer)       
      read(7,*) den0,velx0,vely0,pr0  
      read(7,*) xloop(1),yloop(1),zloop(1),rloop(1)
      read(7,*) xloop(2),yloop(2),zloop(2),rloop(2)
      read(7,*) grav
c
c     # crater info
c     # RC = inner radius
      read(7,*) RC
c     
c     # RD = outer radius
      read(7,*) RD    
c    
c     # DC used in Figure 5, from Table 1 not 2*RC/3
      read(7,*) DC    
c    
c     # smoothed -- use Gaussian smoothing?
      read(7,*) smoothed
c    
c     # r1 used for Gaussian smoothing of lip
      read(7,*) r1     
c     
c     # rdecay used for Gaussian smoothing of lip
      read(7,*) rdecay      
c     
c     # linfactor to scale DC for linearization
      read(7,*) linfactor
c
      read(7,*) thinc_beta
      read(7,*) vofmin  
      close(7)

c     from KM code, since clawpack has dimension 7 for method
      method(8) =  2
      method(9) =  1
c
      return
      end
