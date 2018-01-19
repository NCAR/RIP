c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine maptform(riy,rjx,rlat,rlon,idir)
c
c     This routine converts a coarse domain dot grid point, <riy,rjx>,
c     into a lat/lon point <rlat,rlon> if idir=1, or vice versa if
c     idir=-1. It works for Lambert Conformal (LC,1), Polar
c     Stereographic (ST,2), Mercator (ME,3), or stretched rotated
c     cylindrical equidistant (SRCE,4) projections.  It is assumed that
c     premaptform has been called prior to this so that the proper
c     constants have been placed in the common block called mptf.
c
c      real*8 ypoint,xpoint,rlatr,rlonr,rlatrrot,rlonrrot,dlon,dlonr
      real*4 ypoint,xpoint,rlatr,rlonr,rlatrrot,rlonrrot,dlon,dlonr
c
      include "commptf"
c
      if (idir.eq.1) then   ! First, deal with idir=1
c
      ypoint=(riy-ciy_mptf)*dskmc_mptf+yc_mptf
      xpoint=(rjx-cjx_mptf)*dskmc_mptf
c
      if (nproj_mptf.eq.1) then ! Lambert Conformal
         rlatr=.5*ihm_mptf*pi_mptf-2.*atan(c2_mptf*
     &      (sqrt(xpoint**2+ypoint**2))**conei_mptf)
         rlonr=xloncr_mptf+(conei_mptf*atan2(xpoint,
     &      -ihm_mptf*ypoint))
      elseif (nproj_mptf.eq.2) then ! Polar Stereographic
         rlatr=(.5*ihm_mptf*pi_mptf-ihm_mptf*2.*atan(sqrt(xpoint**2+
     &      ypoint**2)/(rearth1_mptf*c1_mptf)))
         if(xpoint.eq.0..and.ypoint.eq.0.) then
            rlonr=xloncr_mptf
         else
            rlonr=xloncr_mptf+atan2(xpoint,-ihm_mptf*ypoint)
         endif
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then ! Mercator
c
c       (Note: nproj=0 means "idealized" (i.e. no map), but we'll treat
c        it as Mercator so nothing in RIP goes haywire.)
c
         rlatr=2.*atan(exp(ypoint/(rearth1_mptf*c1_mptf)))-.5*pi_mptf
         rlonr=xloncr_mptf+xpoint/(rearth1_mptf*c1_mptf)
      elseif (nproj_mptf.eq.4) then ! stretch-rot-cyl-equid (SRCE)
         rlatrrot = ypoint          /rearth2_mptf
         rlonrrot = xpoint*cone_mptf/rearth2_mptf
c
c     Note, if the equator is rotated, then the above lat/lon values are
c     in the rotated lat/lon system.  To get the natural lat/lon values,
c     need to do the following transformation:
c
         if (c1_mptf.ne.1.0) then
            rlatr=asin ( c2_mptf*cos(rlatrrot)*cos(rlonrrot) +
     &         c1_mptf*sin(rlatrrot) )
            rlonr=xloncr_mptf+asin(sin(rlonrrot)*cos(rlatrrot)/
     &         cos(rlatr))
         else
            rlatr=rlatrrot
            rlonr=rlonrrot
         endif
      endif
      rlat=rlatr/rpd_mptf
      rlon=rlonr/rpd_mptf
      rlon=(mod(rlon+900.,360.)-180.)
c
      else   ! Otherwise, deal with idir=-1
c
      dlon=rlon-xloncr_mptf/rpd_mptf
      dlon=(mod(dlon+900.,360.)-180.)
      dlonr=dlon*rpd_mptf
      rlatr=rlat*rpd_mptf
      if (nproj_mptf.eq.1) then ! Lambert Conformal
         ypoint=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rlatr)))**cone_mptf*cos(cone_mptf*dlonr)
         xpoint=ihm_mptf*c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rlatr)))**cone_mptf*sin(cone_mptf*dlonr)
      elseif (nproj_mptf.eq.2) then ! Polar Stereographic
         ypoint=-rearth1_mptf*sin(.5*ihm_mptf*pi_mptf-rlatr)*
     &      c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rlatr))*
     &      cos(dlonr)
         xpoint=ihm_mptf*rearth1_mptf*sin(.5*ihm_mptf*pi_mptf-
     &      rlatr)*c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rlatr))*
     &      sin(dlonr)
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then ! Mercator
         ypoint=rearth1_mptf*c1_mptf*log((1.+sin(rlatr))/
     &      cos(rlatr))
         xpoint=dlonr*rearth1_mptf*c1_mptf
      elseif (nproj_mptf.eq.4) then ! stretch-rot-cyl-equid (SRCE)
c
c     Note, the input lat/lon values are in the natural (unrotated)
c     lat/lon system.  If the equator is rotated, these must be
c     converted to lat/lon values in the rotated system with the
c     following transformation:
c
         if (c1_mptf.ne.1.0) then
            rlatrrot=asin(c1_mptf*sin(rlatr)-
     &         c2_mptf*cos(rlatr)*cos(dlonr))
            rlonrrot=atan2(cos(rlatr)*sin(dlonr),
     &         c1_mptf*cos(rlatr)*cos(dlonr)+c2_mptf*sin(rlatr) )
         else
            rlatrrot=rlatr
            rlonrrot=rlonr
         endif
         ypoint=rlatrrot*rearth2_mptf
         xpoint=rlonrrot*rearth2_mptf/cone_mptf
      endif
      riy=(ypoint-yc_mptf)/dskmc_mptf+ciy_mptf
      rjx=xpoint/dskmc_mptf+cjx_mptf
c
      endif
c
      return
      end
