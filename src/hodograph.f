C_____________________________________________________________________
C
      SUBROUTINE HODOGRAPH(UW,VW,PRES,NLEVELS,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,icomax,barbfac,iwkidcgm)
C
c  subroutine to plot hodgraph in rip. Pres is in hPa, U and V come
c  in as whatever was requested in svecdraw (so it's the same units as
c  the wind barbs).
c
      PARAMETER(NANGLES=360)
      DIMENSION PRES(*)
      REAL UW(*),VW(*),XCIRC(NANGLES),YCIRC(NANGLES)
      real uu(nlevels), vv(nlevels)
      INTEGER IRING(4)
      CHARACTER CTIT*2, units*3
C
C  Convert from dir,speed to u,v and find max speed for nice speed rings
C
      PI=3.14159
      nl = 0
      XMAXSPD  = 0.
      k500 = 0
      DO I=NLEVELS,1,-1
            if(pres(i).lt.225.0) goto 25
            nl = nl + 1
c    uw(i) = 1.94 * uw(i)
c    vw(i) = 1.94 * vw(i)
	    spd = sqrt(uw(i)*uw(i) + vw(i)*vw(i))
c    write(6,*) 'i = ',i,' pres = ',pres(i),' spd = ',spd
            if(pres(i).gt.500.0) k500=nl
            XMAXSPD = MAX(XMAXSPD,SPD)
      ENDDO
 25   continue
c     write(6,*) 'max spd = ',XMAXSPD,' k500 = ',k500
c     write(6,*) 'nl = ',nl
      IF(XMAXSPD.GT.50.0) THEN
         IRING(1)=25
         IRING(2)=50
         IRING(3)=75
         IRING(4)=5
      ELSEIF(XMAXSPD.LE.50.0 .AND. XMAXSPD.GT.25.0) THEN
         IRING(1)=10
         IRING(2)=25
         IRING(3)=50
         IRING(4)=5
      ELSEIF(XMAXSPD.LE.25.0) THEN
         IRING(1)=5
         IRING(2)=10
         IRING(3)=15
         IRING(4)=1
      ELSE
         STOP 'IRING'
      ENDIF
C
C  Set up window for hodograph in top left corner of skewt
C
c     CALL SET(.05,.95,.05,.95,-19.0,27.1,-.9346217,44.061,1)
c     CALL SET(.05,.95,.05,.95,-19.0,27.1,-.9346217,44.061,1)
c     CALL SET(.0,1.,.0,.9,-19.0,27.1,-.9346217,44.061,1)
      space=3.
c     CALL SET(.05,.95,.05,.95,-19.0-space,27.1+space,
      CALL SET(flminsou,frmaxsou,fbminsou,ftmaxsou,
     &   -19.0-space,27.1+space,-.9346217-.2*space,44.061+space,1)
      CALL GETSET(SL,SR,SB,ST,UL,UR,UB,UT,LL)
c     write(6,*) 'sl = ',sl,' sr = ',sr,' sb = ',sb,' st = ',st
c     write(6,*) 'ul = ',ul,' ur = ',ur,' ub = ',ub,' ut = ',ut
c     UTN=UT
c     ULN=UL
      UTB=30.796  ! y coord of the 200 hPa level; 44.061 = 100 hPa
c     URN=-4.5
c     SLN=SL
c     STN=ST
c     SBN = ST - ( (ST-SB)*(UTN-UTB) )/(UT-UB)
c     SRN = SL + ( (SR-SL)*(URN-ULN) )/(UR-UL)
      sln = cufx(-19.)
      srn = cufx(-4.5)
      srn = .31
      sbn = cufy(utb)
      stn = cufy(44.061)
c     write(6,*) 'sln = ',sln,' srn = ',srn,' sbn = ',sbn,
c    & ' stn = ',stn
      CALL SET(SLN,SRN,SBN,STN
     +        ,-XMAXSPD,XMAXSPD,-XMAXSPD,XMAXSPD,0)
c
c blank out the hodo...
      call ngdots(0.,0.,1,1000.,0)
C
C  Draw 18 and 28 m/s anulus around 9 km AGL wind
C
c  find some important wind levels for plotting.
c  change end point from 9km to 10km
      z_0 = std_atmos(pres(nlevels)*100.)
      k1k = -99
      k3k = -99
      DO K=NLEVELS,1,-1
         z_k = std_atmos(pres(k)*100.)
	 if (z_k .gt. z_0+1000. .and. k1k .lt. 0) k1k = k
	 if (z_k .gt. z_0+3000. .and. k3k .lt. 0) k3k = k
         if (z_k .gt. z_0+10000.) goto 35
      ENDDO
 35   continue
      k_9km = k+1
      k_9km = min(k_9km,nlevels)
      go to 987
      call storm_motion_rs(uw,vw,pres,nl,u_storm,v_storm)
      xoff=uw(k_9km) 
      yoff=vw(k_9km)
c
c   Define special colors for hodograph
c
      icolgray=icomax+1
      icolgreen=icomax+2
      icolblue=icomax+3
      call gscr (iwkidcgm,icolgray, 0.75,0.75,0.75)   ! light gray
      call gscr (iwkidcgm,icolgreen,0.50,1.00,0.50)   ! light green
      call gscr (iwkidcgm,icolblue, 0.50,0.50,1.00)   ! light blue
      icomax=icomax+3
      
      x1 = 28.0*1.94
      y1 = 0.0
      call ngdots(xoff,yoff,1,x1*2.0,icolgray)
c                  x,   y, num, size, icolor

      x1 = 18.0*1.94
      y1 = 0.0
      call ngdots(xoff,yoff,1,x1*2.0,icolgreen)

      x1 = 12.0*1.94
      y1 = 0.0
      call ngdots(xoff,yoff,1,x1*2.0,0)

      xoff=u_storm
      yoff=v_storm
      x1 = 4.0*1.94
      y1 = 0.0
      call ngdots(xoff,yoff,1,x1*2.0,icolblue)
c     force = sqrt(u_storm*u_storm + v_storm*v_storm)
c     force = min(force,99.0)
c     direc = atan2(-u_storm,-v_storm)*180./3.14159
c     if (direc.lt.0.0) direc = direc+360.0
c     m_storm = nint(direc)*100 + nint(force)
  987 continue
C
C  Draw hodograph rings
C
      CALL GSPLCI(1)
      xoff=0.
      yoff=0.
c     CALL PLCHMQ(XOFF,YOFF,'+',20.,0.,0.)
      call gsplci(3)
      call line (-XMAXSPD,0.,XMAXSPD,0.)
      call line (0.,-XMAXSPD,0.,XMAXSPD)
      call gsplci(1)
      DO N=1,3
         X1 = FLOAT(IRING(N))
         Y1 = 0.
         if (n.eq.1) then
            n_beg = 15
         elseif (n.eq.2) then
            n_beg = 8
         elseif (n.eq.3) then
            n_beg = 6
         endif
         DO IANGLE=1,NANGLES
            XANGLE=(135.+FLOAT(IANGLE))*PI/180.
            XCIRC(IANGLE)=(COS(XANGLE)*X1 - SIN(XANGLE)*Y1)+ XOFF
            YCIRC(IANGLE)=(SIN(XANGLE)*X1 - COS(XANGLE)*Y1)+ YOFF
         enddo
         CALL CURVE(XCIRC(n_beg),YCIRC(n_beg),NANGLES-n_beg*2+1)
         WRITE(CTIT,'(I2)') IRING(N)
         CALL PLCHMQ(XCIRC(NANGLES),YCIRC(NANGLES),CTIT(1:2),9.,-45.,0.)
      enddo
C
C  Draw the hodograph itself
C
      do nn = nlevels,1,-1     ! invert wind arrays, find level above 9km
        n = nlevels-nn + 1
        uu(nn) = uw(n)
        vv(nn) = vw(n)
        if (nn .eq. k_9km) nlevelscurve = min(n+1,nlevels)
        if (nn .eq. k1k) n1k = min(n,nlevels)
        if (nn .eq. k3k) n3k = min(n,nlevels)
      enddo
      CALL GSLWSC(2.)
c     CALL CURVE(UW,VW,NLevels)
      call frstpt (uu(1),vv(1))
      xlw = 1.5
      CALL GSLWSC(xlw)
      call gsplci(13)  !red
      do nn = 2, nlevelscurve
c       CALL GSLWSC(xlw)
        if ( nn .eq. n1k ) call gsplci(38)
        if ( nn .eq. n3k ) call gsplci(1)
	call vector(uu(nn),vv(nn))
	CALL PLOTIF (0.,0.,2)
c       xlw = xlw+.1
      enddo
c     CALL CURVE(Uu,Vv,nlevelscurve)
      CALL GSLWSC(1.)
c     angl=atan2( (vw(nlevels-1)-vw(nlevels)),(uw(nlevels-1)-
c    &   uw(nlevels)) )*180./pi
c     CALL PLCHMQ(UW(nlevels),VW(nlevels),'V',8.,angl+90.,0.)
c     if (k500.gt.0 .and. k500.le.NL) then
c        angl=atan2( (vw(k500)-vw(k500-1)),(uw(k500)-uw(k500-1)) )
c    +         *180./pi
c        CALL PLCHMQ(UW(k500),VW(k500),'V',8.,angl+90.,0.)
c     endif
c restore u, v to m/s
c     do n = nlevels,1,-1
cif ( pres(n) .lt. 225. ) return
cuw(n) = uw(n) / 1.94
cvw(n) = vw(n) / 1.94
c     enddo
C
      call line (-XMAXSPD,-XMAXSPD,-XMAXSPD,XMAXSPD)
      call line (-XMAXSPD,XMAXSPD,XMAXSPD,XMAXSPD)
      call line (XMAXSPD,XMAXSPD,XMAXSPD,-XMAXSPD)
      call line (XMAXSPD,-XMAXSPD,-XMAXSPD,-XMAXSPD)
      if (barbfac .eq. 1.) then
	units = 'm/s'
      else
	units = 'kts'
      endif
      call gsplci(1)        ! set color to def.foreground
c plot units in lower-left corner
      call plchmq(-.94*XMAXSPD,-.92*XMAXSPD,units,11.,0.,-1.)  
      RETURN
      END
C_____________________________________________________________________
C
      subroutine storm_motion_rs(u,v,p,nlevs,u_storm,v_storm)
      dimension u(*),v(*),p(*)
      parameter (pi=3.14159265, dtr=pi/180., dpr=180./pi)
c
      u_storm = u(nlevs)
      v_storm = v(nlevs)
      u_pbl = u(nlevs)
      v_pbl = v(nlevs)
      k_pbl = 1
c     write(6,*) 'u sfc = ',u_pbl,' v sfc = ',v_pbl
      z_1k = std_atmos(p(nlevs)*100.)
c     write(6,*) 'z_1k = ',z_1k
      do k = nlevs-1,1,-1
         if (p(k).ge.p(nlevs)-100.0) then
            u_pbl = u_pbl + u(k)
            v_pbl = v_pbl + v(k)
            k_pbl = k_pbl + 1
         endif
         z_k = std_atmos(p(k)*100.)
         if (z_k.gt.z_1k+4000.) goto 45
      enddo
 45   continue
c     write(6,*) 'k_pbl = ',k_pbl
      u_pbl = u_pbl/k_pbl
      v_pbl = v_pbl/k_pbl
      k_4km = min(k,nlevs)
c     write(6,*) 'u_pbl = ',u_pbl,' v_pbl = ',v_pbl
c     write(6,*) 'k_4km = ',k_4km
      if (k_4km .lt. nlevs) then
         u_4km = u(k_4km)
         v_4km = v(k_4km)
         ub = .6 * (u_4km - u_pbl)
         vb = .6 * (v_4km - v_pbl)
	 if (ub .eq. 0. .and. vb .eq. 0.) then
	   adr = 0.
	 else
	   adr = dpr * (pi + atan2(ub,vb))
	 endif
c        write(6,*) 'adr = ',adr
c Now find the point that's 8 m/s orthogonal to the BL-4km shear
          bsp = 8.*1.94
          bdr = adr + 90.
          if (bdr .gt. 360.) bdr = bdr-360.
c  write(6,*) 'bdr = ',bdr
          cu = -bsp * sin(bdr*dtr)
          cv = -bsp * cos(bdr*dtr)
c  write(6,*) 'cu = ',cu,' cv = ',cv
C Here's our motion vector
          cu = u_pbl + ub + cu
          cv = v_pbl + vb + cv
c  write(6,*) 'motion u = cu',cu,' cv = ',cv
	  u_storm = cu
	  v_storm = cv
c        shr0_4u = u_4km - u_pbl
c        shr0_4v = v_4km - v_pbl
c        shr0_4m = sqrt( (shr0_4u*0.6)**2 + (shr0_4v*0.6)**2)
c        theta = atan(8.0*1.94/shr0_4m)
c        hypot = sqrt(shr0_4m**2 + (8.0*1.94)**2)
c        u_storm = hypot*cos(-1.*theta) + u_pbl
c        v_storm = hypot*sin(-1.*theta) + v_pbl
      endif
C
      RETURN
      END
c+---+-----------------------------------------------------------------+
c
      function std_atmos(pres_Pa)
c..
c..            pres_Pa = Pressure in Pascals
c..            standard atmos height in meters is returned for given p
c..
      pr = pres_Pa*0.01
      height = 44307.692 * (1.0 - (pr/1013.25)**0.190)
      std_atmos=height
      return
      end
