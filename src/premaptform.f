c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,true1,true2,iup)
c
c     This routine calculates constants for routine maptform and puts
c     them in common block mptf.
c
c      real*8 true1r,true2r,xlatcr,cotrue1r
      real*4 true1r,true2r,xlatcr,cotrue1r
c
      include "commptf"
c
      pi_mptf=4.*atan(1.)     ! 3.1415...
      rpd_mptf=pi_mptf/180.    ! radians per degree
      rearth1_mptf=6370.949  ! true mean radius of planet, in km
      rearth2_mptf=6366.2  ! radius of planet as 20,000 km / pi
      dskmc_mptf=dskmc
      xloncr_mptf=xlonc*rpd_mptf  ! longitude parallel to y-axis in rad
      nproj_mptf=nproj
      ciy_mptf=.5*(1.+miycors)
      cjx_mptf=.5*(1.+mjxcors)
c
      true1r=true1*rpd_mptf
      true2r=true2*rpd_mptf
      xlatcr=xlatc*rpd_mptf
c
      if (nproj_mptf.eq.1) then   ! Lambert Comformal
c
c   Make sure xlatcr, true1r, and true2r are all in same hemisphere,
c      and calculate ihm_mptf.
c
      if (xlatcr.gt.0..and.true1r.gt.0..and.true2r.gt.0.) then
         ihm_mptf=1
      elseif (xlatcr.lt.0..and.true1r.lt.0..and.true2r.lt.0.) then
         ihm_mptf=-1
      else
         write(iup,*)'Invalid latitude parameters for LC map.'
         stop
      endif
c
c   Calculate cone factor
c
      if (true1r.ne.true2r) then
         cone_mptf=log10(cos(true1r)/cos(true2r))/
     &        log10(tan(.25*pi_mptf-ihm_mptf*.5*true1r)/
     &              tan(.25*pi_mptf-ihm_mptf*.5*true2r))
      else
         cone_mptf=cos(0.5*pi_mptf-ihm_mptf*true1r)
      endif
c
c   Calculate other constants
c
      conei_mptf=1./cone_mptf
      cotrue1r=ihm_mptf*0.5*pi_mptf-true1r
      c1_mptf=rearth1_mptf*sin(cotrue1r)/
     &   (cone_mptf*(ihm_mptf*tan(.5*cotrue1r))**cone_mptf)
      c2_mptf=tan(.5*cotrue1r)*(cone_mptf/
     &   (ihm_mptf*rearth1_mptf*sin(cotrue1r)))**conei_mptf
      yc_mptf=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &   2.*xlatcr)))**cone_mptf
c
      elseif (nproj_mptf.eq.2) then   ! Polar Stereographic
c
c   Make sure xlatcr and true1r are both in same hemisphere,
c      and calculate ihm_mptf.
c
      if (xlatcr.gt.0..and.true1r.gt.0.) then
         ihm_mptf=1
      elseif (xlatcr.lt.0..and.true1r.lt.0.) then
         ihm_mptf=-1
      else
         write(iup,*)'Invalid latitude parameters for PS map.'
         stop
      endif
c
c   Calculate cone factor
c
      cone_mptf=1.
c
c   Calculate other constants
c
      conei_mptf=1./cone_mptf
      cotrue1r=ihm_mptf*0.5*pi_mptf-true1r
      c1_mptf=1.+cos(cotrue1r)
      c2_mptf=1. !not used
      yc_mptf=-rearth1_mptf*sin(.5*ihm_mptf*pi_mptf-xlatcr)*
     &   c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-xlatcr))
c
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then   ! Mercator
c
c       (Note: nproj=0 means "idealized" (i.e. no map), but we'll treat
c        it as Mercator so nothing in RIP goes haywire.)
c
      ihm_mptf=1
      cone_mptf=1.
      conei_mptf=1.
      c1_mptf=cos(true1r)
      c2_mptf=1.  !not used
      yc_mptf=rearth1_mptf*c1_mptf*log((1.+sin(xlatcr))/
     &   cos(xlatcr))
c
      elseif (nproj_mptf.eq.4) then
c
c     stretched rotated cylindrical equidistant (SRCE), i.e., the map
c     projection on which an Eta or NMM grid has straight
c     perpendicular grid lines and square grid boxes
c
      ihm_mptf=1   ! not used
      cone_mptf=true1  ! this is neither a true latitude nor a cone factor;
                       ! it is actually the "stretch factor" (dxdeg/dydeg)
                       ! that was set up for the NMM grid (it is
                       ! normally slightly > 1.0)
      conei_mptf=1.! not used
      c1_mptf=cos(xlatcr)
      c2_mptf=sin(xlatcr)
      yc_mptf=0.
c
      endif
c
      return
      end
