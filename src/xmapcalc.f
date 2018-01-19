c                                                                     c
c*********************************************************************c
c                                                                     c
      function xmapcalc(riy,rjx)
c
c     This function calculates the map scale factor at a coarse domain
c     dot grid point, <riy,rjx>.  It works for Lambert Conformal (LC,1), Polar
c     Stereographic (ST,2), Mercator (ME,3), or stretched rotated
c     cylindrical equidistant (SRCE,4) projections.  It is assumed that
c     premaptform has been called prior to this so that the proper
c     constants have been placed in the common block called mptf.
c
      include "commptf"
c
      if (nproj_mptf.ne.4) then
         call maptform(riy,rjx,rlat,rlon,1)
         rlatr=rpd_mptf*rlat
      else
         ypoint=(riy-ciy_mptf)*dskmc_mptf+yc_mptf
         rlatrrot = ypoint/rearth2_mptf
      endif
c
      if (nproj_mptf.eq.1) then   ! Lambert Conformal
         psi=ihm_mptf*pi_mptf*0.5-rlatr
         xmapcalc=cone_mptf*c1_mptf/rearth1_mptf*
     &      (ihm_mptf*tan(0.5*psi))**cone_mptf/sin(psi)
      elseif (nproj_mptf.eq.2) then   ! Polar Stereographic
         psi=ihm_mptf*pi_mptf*0.5-rlatr
         xmapcalc=c1_mptf/(1.+cos(psi))
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then ! Mercator
c
c       (Note: nproj=0 means "idealized" (i.e. no map), but we'll treat
c        it as Mercator so nothing in RIP goes haywire.)
c
         xmapcalc=c1_mptf/cos(rlatr)
      elseif (nproj_mptf.eq.4) then ! stretch-rot-cyl-equid (SRCE)
         xmapcalc=1./(cone_mptf*cos(rlatrrot))
      endif
c
      return
      end
