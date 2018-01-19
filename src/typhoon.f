c this file contains subroutines to support the typhoon track
c plotting.
c
c Developed under sponsorship from Taiwan's CAA, it is designed to track
c typhoons in domains with grid spacings that range from 135 km to 5 km.
c It is possible to tune the algorithm to perform better for a particular grid
c spacing.
c
c sample plspec entry:
c
c  feld=map; ptyp=hb; ouds=solid; oulw=2; cint=10.; mllm=land;>
c    mfco=light.cerulean,peach
c  feld=map; ptyp=hb; ouds=solid; oulw=2; cint=10.; mllm=land
c  feld=track; ptyp=hb; colr=dark.red; tslb=.014; tynt=6.
c  feld=tic; ptyp=hb; axlg=0.; axtg=0.
c
c note that ptimes should be set to a single time (usually max forecast hour)
c
c If a file named 'tcdat' exists in the current directory, the routine will
c read information from that file to get the initial typhoon location instead
c of using the model's zero hour output. See the format statement for formatting 
c information.

c The track algorithm uses multiple criteria to determine typhoon location and
c to eliminate transients and extratropical cyclones. The location is determined
c by a weighted combination of minimum central pressure, minimum surface winds,
c maximum 650 hPa vorticity, and maximum 850 hPa vorticity. Thresholds for vorticity,
c surface temperature, and 700 hPa temperature are also used.

c In general, the performance of the algorithm depends on the model grid spacing
c and strength of the cyclone. For small grid spacings (say less than 8 km), where
c an eyewall can be resolved, the max vorticity and max winds may be found in the
c eyewall away from the typhoon center. For large spacings, an eye cannot be resolved.
c For weak cyclones, the computed track sometimes zig-zags as convection wraps around
c the cyclone center. Unsteady, or noisy, tracks are often due to the model forecast
c rather than deficiencies in the algorithm. It is useful to compare model forecast
c fields with the track plot before modifying the algorithm. Note that mean SLPs in
c the West Pacific are lower than in the West Atlantic. Weak storms in the Atlantic
c may require an adjustment of the SLP threshold for accurate track plotting.

c  Developed by Jim Bresch NCAR/MMM. 
c  Track plotting and labelling by Michael Duda NCAR/MMM
c  Additional contributions from Wei Wang.
c
      subroutine typhoon (uuu,vvv,tmk,qvp,www,prs,ght,prs_tsf,ter,
     &     ixwin,iywin,irota,dmap,xmap,cor,u10,v10,k700,
     &     miy,mjx,mkzh,maxptimes, itime, tloc)

c routine to compute typhoon locations

      parameter (MAX_STORMS = 5)      
      parameter (MAX_PARAMETERS = 7)  

      dimension uuu(miy,mjx,mkzh),vvv(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),www(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),ght(miy,mjx,mkzh),prs_tsf(miy,mjx),
     &   dmap(miy,mjx),xmap(miy,mjx),cor(miy,mjx),
     &   u10(miy,mjx),v10(miy,mjx),k700(miy,mjx),
     &   ter(miy,mjx), tloc(maxptimes,max_parameters,max_storms)
      character tnames(max_storms)*20

      parameter (sealpmin = 1004.)
      real vor(miy,mjx), slp(miy,mjx)
      integer ltyph(max_storms)
      logical debug, ltcdat
      character line*80, tname*8, strng(2)*20
c
      include 'comconst'

      debug = .false.
c     debug = .true.
c     if ( itime .eq. 13 ) debug = .true.
c
      if (debug) write(6,*) 'Begin typhoon, time = ',itime
      expon=rgas*ussalr/grav
      exponi=1./expon

c typhoons must be separated by 850 km or 8 grid spaces
      sepmin = max(8.*dskm,  850.)    
c     sepmin = max(8.*dskm, 450.)     ! wei uses 450
c set the vorticity threshold based on the grid spacing
      vorlimit = max(min(1000.*(1./dskm),120.),10.)
      vorlimit = 10.
c     if (debug) write(6,*) 'vorlimit = ',vorlimit

      if (itime .eq. 1) then
c initialize the array at the first time
	do i = 1, maxptimes
	do j = 1, max_parameters
	do k = 1, max_storms
	  tloc(i,j,k) = -1.
	enddo
	enddo
	enddo
        inquire (FILE='tcdat', EXIST=ltcdat)
        if (ltcdat) then
c   read the tcbogus file to find the location of the typhoon at initial time
        open(47,file='tcdat',form='formatted',status='unknown')
        read(47,'(i8)',err=601,end=602) mdate
        read(47,'(i1)',err=601,end=602) nlines
        do id = 1, nlines
          read(47,'(a)') line
          read (line,'(a8,1x,f4.1,3x,f5.1,3x,f4.1,3x,f5.1,3x,i4,4x,i3,
     &       4x,i2)') tname, xlatc, ylonc, xlat6, ylon6, ip, irad, isp
          if (line(15:15) .eq. 'S' ) xlatc = -1. * xlatc
          if (line(23:23) .eq. 'W' ) ylonc = -1. * ylonc
          write(6,*) 'name = ',tname
          write(6,*) 'xlat6 = ',xlat6,' ylon6 = ',ylon6
          write(6,*) 'xlatc = ',xlatc,' ylonc = ',ylonc
          write(6,*) 'central pressure = ',ip,' radius = ',irad
	  call maptform (riy,rjx,xlatc,xlonc,-1)
	  write(strng(1),'(f8.3,a3)') xlatc,'lat'
	  write(strng(2),'(f8.3,a3)') ylonc,'lon'
	  write(6,*) 'strng(1) = ',strng(1)
	  write(6,*) 'strng(2) = ',strng(2)
          call locinterp(strng,gridx,gridy,
     &              rlat,rlon,iwmo,icaoid,stelv,tnames(id),' ')

	  rjx = 1. + (gridx - xjcorn) * refrat
	  riy = 1. + (gridy - yicorn) * refrat
          write(6,*) 'i = ',riy,' j = ',rjx
	  if (riy .lt. miy .and. riy .ge. 1 .and. 
     &        rjx .lt. mjx .and. rjx .ge. 1) then
	    tnames(id) = tname
	    tloc(1,1,id) = riy
	    tloc(1,2,id) = rjx
	    tloc(1,3,id) = 0
	    tloc(1,4,id) = float(ip)
	    if (debug) then
	      write(6,*) 'y = ',tloc(1,1,id),' x = ',tloc(1,2,id) 
	    endif
	  endif
        enddo
  601  continue
  602  continue
      endif
      endif
c if there's ever been a typhoon reported, track it
      do l = 1, max_storms
	ltyph(l) = 0
	do i = 1, itime
	  if (tloc(i,1,l) .ne. -1. ) ltyph(l)=1
	enddo
      enddo
c compute the slp array
      slpmean = 0.
      do j=2,mjx-1
      do i=2,miy-1
	tvlhsl=virtual(tmk(i,j,mkzh),qvp(i,j,mkzh))
	prslhsl=prs(i,j,mkzh)
	psurf = prs_tsf(i,j)
	ghtlhsl=ter(i,j)+tvlhsl/ussalr*
     &     ((psurf/prslhsl)**expon - 1.)
	slp(i,j)=prslhsl*(1.+ussalr/tvlhsl*ghtlhsl)**exponi
	slpmean = slpmean + slp(i,j)
      enddo
      enddo
      slpmean = slpmean / ((mjx-3) * (miy-3))
      if (debug) write(6,*) 'slpmean = ',slpmean
      do nty = 1, max_storms
        if (debug) write(6,*) '        Begin loop, nty = ',nty
	iy = -1
	jx = -1
	slpmin = sealpmin
	vorm = 0.
	vorn = 0.
	vors = 0.
	k850 = 0
	do j = 1, mjx-1
	do i = 1, miy-1
	k7 = 0
	do k = mkzh, 1, -1
c k700 actually contains the 650 mb level
	  if (prs(i,j,k) .lt. 650. .and. k7 .eq. 0) then
	    k7 = k
	  endif
	enddo
	k700(i,j) = k7
	enddo
	enddo
	do k = mkzh, 1, -1
	  if (prs(miy/2,mjx/2,k) .lt. 850. .and. k850 .eq. 0) then
	    k850 = k
	  endif
	enddo
        if (debug) then
	do k = mkzh, 1, -1
	  write(6,'(i4,f20.8)') k,prs(55,73,k)
	enddo
        write(6,*) 'k700 = ',k700(miy/2,mjx/2)
        write(6,*) 'ltyph(nty) = ltyph(',nty,') = ',ltyph(nty)
	endif
	if ( ltyph(nty) .eq. 0 ) then
	  do 190 j=2,mjx-1
	  do 190 i=2,miy-1
c allow a tc to form only over water, except at initial time
            if ( itime .eq. 1 .or. ter(i,j) .le. 0. ) then
c    &          nint(xlus(i,j)) .eq. iwater) ) then   xlus n/a in rip4
              vor7 = cvor (uuu, vvv, i, j, k700(i,j), ixwin, iywin,
     &                irota, miy,mjx,mkzh)
c check to see if it's a closed low
              ifail = 0
              if ( slp(i,j) .gt. 1004.) then
                irad = max (1, nint (120./dskm))
                imin = max(2,i-irad)
                imax = min(miy-1,i+irad)
                jmin = max(2,j-irad)
                jmax = min(mjx-1,j+irad)
c               if (debug) then
c               write(6,*) 'i = ',i,' j = ',j
c               write(6,*) 'slp(i,j) = ',slp(i,j),
c    &   '        slp(imax,j) = ',slp(imax,j)
c               endif
                if (slp(i,j) - slp(imax,j) .gt. -2.) ifail = 1
c               if (debug) then
c               write(6,*) 'slp(i,j) = ',slp(i,j),' slp(imin,j) = ',
c    &           slp(imin,j)
c               endif
                if (slp(i,j) - slp(i,jmax) .gt. -2.) ifail = 1
c               if (debug) then
c               write(6,*) 'slp(i,j) = ',slp(i,j),' slp(i,jmax) = ',
c    &           slp(i,jmax)
c               write(6,*) 'slp(i,j) = ',slp(i,j),' slp(i,jmin) = ',
c    &           slp(i,jmin)
c               endif
                if (slp(i,j) - slp(imin,j) .gt. -2.) ifail = 1
                if (slp(i,j) - slp(i,jmin) .gt. -2.) ifail = 1
c               if ( debug .and. ifail .eq. 1) then
c                 write(6,*) 'closed low failure'
c               endif
              endif
              if (slp(i,j) .lt. slpmin .and. tmk(i,j,k700(i,j)) .gt.270.
     &          .and. vor7 .gt. vorlimit .and. ifail .eq. 0) then
            if (debug) then
            write(6,*) 'slpmin = ',slpmin,' slp = ',slp(i,j),' t700 = ',
     &        tmk(i,j,k700(i,j)), ' x = ',j,' y = ',i
            endif
		lcheck = 0
		do lch = 1, max_storms
                  if (ltyph(lch) .ne. 0. ) lcheck = 1
		enddo
                if (lcheck .ne. 0 ) then
c               if (ltyph(1) .ne. 0. .or. ltyph(2) .ne. 0 .or. 
c    &              ltyph(3) .ne. 0 ) then
                  if ( nty .eq. 1) then
                    l = 2
                  else
                    l = 1
                  endif
                  do 187 l = 1, max_storms
                    if ( l .eq. nty) go to 187
c                if (debug) write(6,*) 'checking for proximity, l = ',l
                    do it = itime, 1, -1
                      if (tloc(it,1,l) .ne. -1) then
                        di = abs(i - tloc(it,1,l)) * dskm
                        dj = abs(j - tloc(it,2,l)) * dskm
                        dis = sqrt (dj**2 + di**2)
c                       if (debug) then
c                       write(6,*) 'inew = ',i,' iold = ',tloc(it,1,l)
c                       write(6,*) 'jnew = ',j,' iold = ',tloc(it,2,l)
c                       write(6,*) 'dis = ',dis
c                       endif
		        if ( dis .le. sepmin ) go to 190
		      endif
                    enddo
  187             continue
c                 if (debug) write(6,*) 'finished proximity loop'
                endif
                slpmin = slp(i,j)
                iy = i
                jx = j
              endif      ! end of limits tests
	    endif             ! end of water test
 190      continue
c       if (debug) then
c        write(6,*) 'min of ',slpmin,' x = ',jx,' y = ',iy,' itime = ',
c    &     itime
c       endif
        else
c         do it = itime-1, 1, -1
	  do it = itime, 1, -1
	    if (tloc(it,1,nty) .ne. -1) then
	      iy = tloc(it,1,nty)
	      jx = tloc(it,2,nty)
	      goto 195
	    endif
	  enddo
          if (debug) write(6,*) 'why am I here?'
	endif
  195   continue
	if ( iy .ne. -1) then
c now search within 200 km of slpmin for a typhoon center
c We assume the typhoon hasn't moved, jumped, or reformed more than
c 200 km away from it's previous position.
          irad = max (3, nint (200./dskm))
          imin = max(2,iy-irad)
          imax = min(miy-1,iy+irad)
          jmin = max(2,jx-irad)
          jmax = min(mjx-1,jx+irad)
	  if (debug) then
          write(6,*) 'checking for jx = ',jx,' iy = ',iy,' irad = ',irad
          write(6,*) 'imin = ',imin,' imax = ',imax
          write(6,*) 'jmin = ',jmin,' jmax = ',jmax
	  endif
          slpmin = sealpmin
	  do 10 j = jmin, jmax
	  do 10 i = imin, imax
c find min slp in the box
            if (slp(i,j) .lt. slpmin) then
              slpmin = slp(i,j)
              iy = i
              jx = j
            endif
   10     continue
          if(debug) then
	    write(6,*) 'slpmin is at x = ',jx,' y = ',iy
	  endif
c extrapolated center location
	  if (itime-2 .gt. 0) then
	    if (tloc(itime-1,1,nty) .gt. 0 .and. 
     &          tloc(itime-2,1,nty) .gt. 0) then
	      iye = nint(tloc(itime-1,1,nty) - tloc(itime-2,1,nty) +
     &                  tloc(itime-1,1,nty))
	      jxe = nint(tloc(itime-1,2,nty) - tloc(itime-2,2,nty) +
     &                  tloc(itime-1,2,nty))
	    endif
	  endif
          if(debug) then
	    write(6,*) 'extrapolated center is at x = ',jxe,' y = ',iye
	  endif
c now check near the new center for the vorticity and wind speed 
          irad = max (3, min(nint (200./dskm),8))
          imin = max(2,iy-irad)
          imax = min(miy-1,iy+irad)
          jmin = max(2,jx-irad)
          jmax = min(mjx-1,jx+irad)
          smin = 99999.
          smax = -9999.
          smin1 = 99999.
          smax1 = -9999.
	  vor7sum = 0.
	  vorssum = 0.
	  icnt = 0 
          do 11 j = jmin, jmax
          do 11 i = imin, imax
c an estimate of the 700 and 850 mb vorticity
            vor7 = cvor (uuu, vvv, i, j, k700(i,j), ixwin, iywin, 
     &                irota, miy,mjx,mkzh)
            vor8 = cvor (uuu, vvv, i, j, k850, ixwin, iywin, 
     &                irota, miy,mjx,mkzh)
            vor0 = cvor (uuu, vvv, i, j, mkzh, ixwin, iywin, 
     &                irota, miy,mjx,mkzh)
c      if (debug) write(6,*)' x = ',j,' y = ',i,' vor7 = ',vor7,
c    &    ' vorm = ',vorm
	    vorssum = vorssum + vor0
	    vor7sum = vor7sum + vor7
	    icnt = icnt + 1
            if (vor0 .gt. vors) then
              vors = vor0
              iy6 = i
              jx6 = j
            endif
            if (vor7 .gt. vorm) then
              vorm = vor7
              iy2 = i
              jx2 = j
            endif
            if (vor8 .gt. vorn) then
              vorn = vor8
              iy5 = i
              jx5 = j
            endif
            wsp1 = sqrt(vvv(i,j,mkzh)**2 + uuu(i,j,mkzh)**2)
c use u10 and v10
            wsp = sqrt(u10(i,j)**2 + v10(i,j)**2)
            if (itime.eq.1) wsp = wsp1
c  minimum wind speed should only be used over water
            if (wsp .lt. smin .and. ter(i,j) .le. 0.) then
              smin = wsp
              iy3 = i
              jx3 = j
c        if (debug) then
c     write(6,*) 'setting smin at i = ',i,' j = ',j,' ter = ',ter(i,j),
c    &  ' wsp = ',wsp
c        endif
            endif
            if (wsp1 .lt. smin1 .and. ter(i,j) .le. 0.) then
              smin1 = wsp1
              iy8 = i
              jx8 = j
c      if (debug) then
c     write(6,*) 'setting smin1 at i = ',i,' j = ',j,' ter = ',ter(i,j),
c    &  ' wsp1 = ',wsp1
c      endif
            endif
c           if (wsp .gt. smax) then
c             smax = wsp
c             iy4 = i
c             jx4 = j
c           endif
   11     continue
c wei: use a different radius to find max sfc wind (~160 km)
          if (dskm.lt.15.) irad = max (3, min(nint (200./dskm),13))
          if (dskm.lt. 5.) irad = max (3, min(nint (200./dskm),40))
c         if(itime.eq.1) print *, 'irad for max wind speed is ',irad
          imin = max(2,iy-irad)
          imax = min(miy-1,iy+irad)
          jmin = max(2,jx-irad)
          jmax = min(mjx-1,jx+irad)
          smax = -9999.
          smax1 = -9999.
      do i = imin,imax
      do j = jmin,jmax
            wsp = sqrt(u10(i,j)**2 + v10(i,j)**2)
            wsp1 = sqrt(vvv(i,j,mkzh)**2 + uuu(i,j,mkzh)**2)
            if (itime.eq.1) wsp=wsp1
            if (wsp .gt. smax) then
              smax = wsp
              iy4 = i
              jx4 = j
            endif
            if (wsp1 .gt. smax1) then
              smax1 = wsp1
              iy7 = i
              jx7 = j
            endif
      end do
      end do
        if (debug) then
        write(6,*) 'min slp of ',slpmin,' x = ',jx,' y = ',iy
        write(6,*) 'max 700 vor of ',vorm,' x = ',jx2,' y = ',iy2
        write(6,*) 'min wsp of ',smin,' x = ',jx3,' y = ',iy3
        write(6,*) 'max 10m wsp of ',smax,' x = ',jx4,' y = ',iy4
        write(6,*) 'max 1st level wsp of ',smax1,' x = ',jx7,' y = ',iy7
        write(6,*) 'max 850 vor of ',vorn,' x = ',jx5,' y = ',iy5
        write(6,*) 'max sfc vor of ',vors,' x = ',jx6,' y = ',iy6
        endif
	  avg7 = 0.
	  avg0 = 0.
	  if ( icnt .gt. 0) then
	    avg7 = vor7sum/float(icnt)
	    avg0 = vorssum/float(icnt)
	    if (debug) then
            write(6,*) 'avg sfc vor is ',avg0
            write(6,*) 'avg 700 vor is ',avg7
	    endif
	  endif
	  if (debug) then
          write(6,*) 'official position',' x = ',jx,' y = ',iy,
     &  ' itime = ',itime,' nty = ',nty
          endif
   12     continue
          a = tmk(iy,jx,k700(iy,jx))-tmk(iy,jmin,k700(iy,jmin))
          b = tmk(iy,jx,k700(iy,jx))-tmk(imax,jx,k700(imax,jx))
          c = tmk(iy,jx,k700(iy,jx))-tmk(iy,jmax,k700(iy,jmax))
          d = tmk(iy,jx,k700(iy,jx))-tmk(imin,jx,k700(imin,jx))
          if (debug) then
          write(6,*) 'tdif at A = ', a 
          write(6,*) 'tdif at B = ', b 
          write(6,*) 'tdif at C = ', c 
          write(6,*) 'tdif at D = ', d 
	  endif
          tmean = (a + b + c + d ) * .25
          t700 = tmk(iy,jx,k700(iy,jx))-273.16
	  if (debug) then
            write(6,*) 'tmean = ',tmean,' 700 mb T = ',t700
          endif
c         if ( tmean .gt. 0. .and. vorm .gt. 10. .and. t700 .gt. 0. ) then
c         if ( tmk(iy,jx,mkzh) .gt. 280. .and. vorm .gt. vorlimit .and.
c    &         t700 .gt. 1. .and. smin .lt. 200. ) then
          if ( tmk(iy,jx,mkzh) .gt. 280. .and. vorm .gt. vorlimit .and.
     &         t700 .gt. 1. .and. smin .lt. 200. .and. 
     &         ( dskm .ge. 90. .or. (avg7 .gt. vorlimit*.33333 .and. 
     &         avg0 .gt. vorlimit*.33333)) ) then
c take a mean of the 4 positions for an "official" position
c if the grid spacing is greater than 90 km, then it can't resolve an
c eye and so don't use min wind speed.
c           iy = nint ((iy + iy2 + iy3) / 3.)
c           jx = nint ((jx + jx2 + jx3) / 3.)
            if ( dskm .le. 90.) then
              x = (jx + jx2-.5 + jx3-.5 + jx5-.5) / 4.
              y = (iy + iy2-.5 + iy3-.5 + iy5-.5) / 4.
c             x = jx
c             y = iy
            else
              x = (jx + jx2-.5 + jx5-.5) / 3.
              y = (iy + iy2-.5 + iy5-.5) / 3.
c             x = jx2-.5
c             y = iy2-.5
            endif
            if (debug) then
            write(6,*) 'x = ',x
            write(6,*) 'y = ',y
            endif
c wei: compute the distance of max sfc wind to storm center
c           disw = sqrt( (float(iy4)-y)**2 + (float(jx4)-x)**2 )*dskm
            if(debug) then
             write(6,*)' itime = ',itime,' tloc y = ',tloc(itime,1,nty)
            endif
	    if (itime .ne. 1 .or. tloc(itime,1,nty) .lt. 0.) then
c  the first time's location may have been read in from an observation file
              tloc(itime,1,nty) = y
              tloc(itime,2,nty) = x
c             tloc(itime,3,nty) = mkzh
              tloc(itime,3,nty) = itime-1
              tloc(itime,4,nty) = slpmin
	    endif
            tloc(itime,5,nty) = smax
            tloc(itime,6,nty) = jx4
            tloc(itime,7,nty) = iy4
            ltyph(nty) = 1
          else
c it's not a typhoon
	    if(debug) then
	      write(6,*) "I don't think it's a typhoon"
	      write(6,*) 'tmk = ',tmk(iy,jx,mkzh)
	      write(6,*) 'vorm = ',vorm,' t700 = ',t700
	      write(6,*) 'smin = ',smin,' avg7 = ',avg7,' avg0 = ',avg0
	    endif
            iy = -1
            jx = -1
          endif
        endif
      enddo   ! end of nty loop
c now clean up any transient depressions
      if ( itime .ge. 3 ) then
	do nty = 1, max_storms
	  if (tloc(itime,1,nty) .eq. -1) then      ! no typhoon at present
	    if (tloc(itime-1,4,nty) .ne. -1 .and. 
     &          tloc(itime-1,4,nty) .gt. 1004.) then ! if a depression existed at the previous time
	      if (tloc(itime-2,1,nty) .eq. -1) then  ! but none existed at the time before that, it's transient
                do j = 1,max_parameters
                  tloc(itime-1,j,nty) = -1          ! remove the transient
		enddo
                if (debug) then
 	          write(6,*) 'removed transient vortex, nty = ',nty,
     &              ' itime = ',itime
                endif
                isthere = 0              ! if the removal means no typhoon at all, reset ltyph
                do i = 1, itime
                  if (tloc(i,1,nty) .ne. -1) isthere = 1
                enddo
                if (isthere .eq. 0) ltyph(nty) = 0
              endif
	    endif
	  endif
	enddo     ! end short nty loop
	  !  check to see if typhoons have merged   9/26/08
	do nty = 1, max_storms
	  do 179 nty2 = nty+1, max_storms
	    if (nty .ne. nty2 .and. tloc(itime,1,nty) .ne. -1 .and. 
     &                              tloc(itime,1,nty2) .ne. -1 ) then
              di = abs(tloc(itime,1,nty) - tloc(itime,1,nty2)) * dskm
              dj = abs(tloc(itime,2,nty) - tloc(itime,2,nty2)) * dskm
              dis = sqrt (dj**2 + di**2)
c             write(6,*) 'nty = ',nty,' nty2 = ',nty2,' dis = ',dis
	      if ( dis .lt. sepmin/3.) then
	        do j = 1, max_parameters
		  tloc(itime,j,nty2) = -1          ! remove the duplicate
		enddo
	        if (tloc(itime-1,4,nty2) .ne. -1 .and. 
     &              tloc(itime-1,4,nty2) .gt. 1004.) then ! if a depression existed at the previous time
	          if (tloc(itime-2,1,nty2) .eq. -1) then  ! but none existed at the time before that, it's transient
                    do j = 1,max_parameters
                      tloc(itime-1,j,nty2) = -1          ! remove the transient
		    enddo
                    if (debug) then
 	              write(6,*) 'removed transient vortex, nty2 = ',
     &                  nty2,' itime = ',itime
                    endif
                    isthere = 0              ! if the removal means no typhoon at all, reset ltyph
                    do i = 1, itime
                      if (tloc(i,1,nty2) .ne. -1) isthere = 1
                    enddo
                    if (isthere .eq. 0) ltyph(nty2) = 0
                  endif
                endif
	      endif
	    endif
  179     continue  ! end nty2 loop
	enddo     ! end nty loop for distance check
      endif
      return
      end
c---------------------------------------------------------------
      function cvor (uuu, vvv, i, j, k, ixwin, iywin, irota,
     &      miy,mjx,mkzh)
      dimension uuu(miy,mjx,mkzh), vvv(miy,mjx,mkzh)
      include 'comconst'
c a crude estimate of the vorticity at a point
      ib = max (i-1,1)
      it = min (i+1,miy)
      jl = max (j-1,1)
      jr = min (j+1,mjx)
      cvor = 1.e5 * (((vvv(i,jr,k)-vvv(i,jl,k))/(2.*ds)) -
     &               ((uuu(it,j,k)-uuu(ib,j,k))/(2.*ds)))
      riy=fy(float(j),float(i))
      if (riy .eq. float(i)) then
        rix=fx(float(j),float(i))
	xpp = xjcorn + (rix+(ixwin-1) -1.)/refrat
	ypp = yicorn + (riy+(iywin-1) -1.)/refrat
	call maptform(ypp,xpp,rlat,rlon,1)
	if (rlat .lt. 0.) cvor = -1. * cvor
      endif
      end
c---------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c The do_track subroutine calculates appropriate locations to label points
c  on typhoon tracks (given in storm_val) and draws the points on the
c  track along with their associated labels.
c If the dimensions of storm_val(:::) are changed, the parameters MAX_STORMS
c  and MAX_PARAMETERS must be set in the following subroutines in this file:
c  do_track(), calc_initial(), iterate(), and uncross().
c Points on a typhoon track are labeled if the forecast hour of the point
c  is a multiple of LABEL_INT, if the point is adjacent to a "gap" in the
c  track, or if the point is the last point on the track. If a value in the
c  storm_val(:::) array is not to be drawn (and hence not labeled) because
c  it contains no valid typhoon position or otherwise, the value of the
c  y and x coordinate should be set to a negative value.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine do_track(storm_val,rcrag,rcrbg,nxtavl,storm,rtslb,
     &   idash, ixwin,iywin,icolr,ipl,irota,ilinw,storm_cutoff1,
     &   storm_cutoff2,init_date,cxtimeavl,maxpl,maxptimes,maxtavl,
     &   rtynt)

      ! Constants
      parameter (MAX_STORMS = 5)      ! The 3rd dimension of storm_val()
      parameter (MAX_PARAMETERS = 7)  ! The 2nd dimension of storm_val()
      parameter (NUM_ITERATIONS = 60) ! Number of iterations of algorithm
      parameter (LABEL_INT = 12)      ! The interval (hours) to label points

      ! Arguments
      integer maxpl,    ! The maximum number of plots
     &  maxptimes,      ! The dimension of ptimes (and first dim of storm_val) 
     &  maxtavl         ! The dimension of cxtimeavl
      real storm_val(maxptimes,MAX_PARAMETERS,MAX_STORMS), 
                                     ! storm_val( ) is the array holding typhoon info:
                                     ! 1st dimension: ptime
                                     ! 2nd dimension: y, x, mkzh, slpmin, smax
                                     ! 3rd dimension: typhoon index
     &  storm_cutoff1,        ! smax threshold for using 2nd symbol
     &  storm_cutoff2         ! smax threshold for using 3rd symbol
      real rcrag(2,maxpl),    ! Used for calling hlinedraw and hbulldraw
     &  rcrbg(2,maxpl),       ! Used for calling hlinedraw and hbulldraw
     &  rtslb(maxpl),         ! Text size for labels
     &  rtynt(maxpl)          ! Typhoon symbol interval, hours
      integer nxtavl,         ! Actual number of values in cxtimeavl
     &  ipl,                  ! Current plot index
     &  init_date             ! Forecast initialization date
      integer ixwin(2,maxpl), ! For plot windowing, the min/max x-coordinate
     &  iywin(2,maxpl),       ! For plot windowing, the min/max y-coordinate
     &  icolr(maxpl),         ! Track color
     &  ilinw(maxpl),         ! Track line width
     &  irota(maxpl),         ! Plot rotation (for POLSTR)
     &  idash(maxpl)          ! dash pattern
      character storm(maxpl)*82      !
      character cxtimeavl(maxtavl)*10 !

      ! Local variables
      real prev_label(maxptimes,2,MAX_STORMS), prev_label2(2), dist3
      real current_fchr, xdif, ydif, theta, label_dist
      real temp, inst, out_interval, old_tslb
      real downset, insetx
      real dist1, dist2, prev_text(maxptimes,2), track_center(2)
      integer old_linw
      integer times_through,icurrent_fchr,
     &   init_hr,new_date,iscp, comp_point
      integer s_width, s_height, old_font, old_color
      integer x_off, y_off, center_count
      integer y, z, iii, i, j, k, kk
      logical cross, track_read, inwin(maxptimes, MAX_STORMS)
      logical do_label(maxptimes, MAX_STORMS)
      logical go_ul

      include 'comconst'

      ! Get the width and height, in gridpoints, of the window
      s_width=(ixwin(2,ipl)-ixwin(1,ipl))
      s_height=(iywin(2,ipl)-iywin(1,ipl)) 

      ! inst is the distance in from the edge that we want labels
      inst = 5.0*rtslb(ipl)*s_height

      ! Save values of line width, font, and color
      old_linw = ilinw(ipl)
      call pcgeti('FN',old_font)      
      call pcgeti('CC',old_color)

      ! first determine if there are any storms to plot, and if not, write
      !   a mesage in the plot saying that there are no storms
      do k=1,MAX_STORMS
      do i=1,maxptimes
        if(storm_val(i,1,k).ge.iywin(1,ipl).and.
     &    storm_val(i,1,k).le.iywin(2,ipl)) then
        if(storm_val(i,2,k).ge.ixwin(1,ipl).and.
     &    storm_val(i,2,k).le.ixwin(2,ipl)) then
            goto 20
        endif
        endif
      enddo
      enddo

      ! if we get to this point, then we have found no storms, so write message
      !   and jump to end of code
      rcrbg(1,ipl)=(iywin(2,ipl)-s_height*4.*
     &  rtslb(ipl))/refrat+yicorn
      rcrbg(2,ipl)=(ixwin(1,ipl)+s_width/2.)
     &  /refrat+xjcorn
      write(storm(ipl),'(a30)') 'No tropical storms forecasted.'
      write(6,*) 'No tropical storms forecasted.'
      call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrbg,
     &  maxpl,ipl,irota)
      goto 900 
   20 continue

      ! if the storm is out of the window, but in the domain, we need to avoid
      !   labeling the track points that are out of the window
      do k=1,MAX_STORMS
      do i=1,maxptimes
        inwin(i,k)=.true.
        if(storm_val(i,1,k).lt.iywin(1,ipl).or.
     &    storm_val(i,1,k).gt.iywin(2,ipl)) inwin(i,k)=.false.
        if(storm_val(i,2,k).lt.ixwin(1,ipl).or.
     &    storm_val(i,2,k).gt.ixwin(2,ipl)) inwin(i,k)=.false.
      enddo
      enddo

      ! Calculate the ideal distance a label should be from the pt. it labels
      label_dist = 0.
      do i=1,nxtavl-1
        do j=1,MAX_STORMS
          if (inwin(i,j) .and. inwin(i+1,j)) then
            if (fdist(storm_val(i,1,j), storm_val(i,2,j), 
     &                storm_val(i+1,1,j), storm_val(i+1,2,j)) .gt. 
     &          label_dist) then
                label_dist = fdist(storm_val(i,1,j), storm_val(i,2,j),
     &                       storm_val(i+1,1,j), storm_val(i+1,2,j))
            end if
          end if
        end do
      end do

      if (label_dist .gt. s_width/9.) then
        label_dist = s_width/9.
      end if
      if (label_dist .lt. 
     &    (rtslb(ipl)*100.*s_width/10.)) then
        label_dist = rtslb(ipl)*100.*s_width/10. 
      end if

      ! set up a value for typhoon symbol interval (if not defined by user, set 
      !   interval to be the time between model output times)
      if(rtynt(ipl).lt..5) then
        if(nxtavl.gt.1) then
          read(cxtimeavl(2),'(f6.1)') current_fchr
          rtynt(ipl) = current_fchr
          read(cxtimeavl(1),'(f6.1)') current_fchr
          rtynt(ipl) = nint(rtynt(ipl)) - nint(current_fchr)
        else
          rtynt(ipl) = 1.
        endif
      endif
      read(cxtimeavl(2),'(f6.1)') current_fchr
      out_interval = current_fchr
      read(cxtimeavl(1),'(f6.1)') current_fchr
      out_interval = nint(out_interval) - nint(current_fchr)
      ! this should not happen...
      if(out_interval.eq.0) then
        write(6,*) 'Warning: Model output is every 0 hours or only one',
     &    ' time exists'
        out_interval = 1.
      endif
      if(rtynt(ipl)/out_interval .ne.
     &   nint(rtynt(ipl)/out_interval)) then 
        write(6,*) 'Warning: rtynt is not a multiple of ',
     &              'model output interval'
      endif

      call pcseti('FN',137)
      call pcseti('CC',icolr(ipl))

      ! we would like to keep the symbols a little smaller than label text
      old_tslb = rtslb(ipl)
      rtslb(ipl) = rtslb(ipl) - .01
      if(rtslb(ipl).lt..015) rtslb(ipl) = .015

      ! This do loop draws all of the lines, open circles, and open storms
      do kk=1, MAX_STORMS 
        times_through = 1

        do iii=1, nxtavl 

          if(storm_val(iii,1,kk).ge.0 .and. 
     &       storm_val(iii,2,kk).ge.0) then 
            rcrbg(1,ipl)=rcrag(1,ipl)
            rcrbg(2,ipl)=rcrag(2,ipl)
            rcrag(1,ipl)=(storm_val(iii,1,kk)-1.)/refrat+yicorn
            rcrag(2,ipl)=(storm_val(iii,2,kk)-1.)/refrat+xjcorn

            read(cxtimeavl(iii),'(f6.1)') current_fchr
            if(nint(current_fchr)/rtynt(ipl).eq.
     &         nint(nint(current_fchr)/rtynt(ipl))) then 

            if(storm_val(iii,5,kk).lt.storm_cutoff1) then 
            ! storm(ipl)(1:1)=char(109) 
              storm(ipl)(1:1)='m'
              call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrag,
     &           maxpl,ipl,irota)
            elseif(storm_val(iii,5,kk).ge.storm_cutoff1 .and. 
     &           storm_val(iii,5,kk).lt.storm_cutoff2) then
            ! storm(ipl)(1:1)=char(112) 
              storm(ipl)(1:1)='p'
              call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrag,
     &          maxpl,ipl,irota)
	    endif 
            endif 

            if(times_through .gt. 1) then 
              call hlinedraw(ilinw,idash,ixwin,iywin,icolr,rcrag,rcrbg,
     &           maxpl,ipl,irota)
            endif 
            times_through = times_through + 1
	  else
              times_through=1
          endif 
                
        enddo 
      enddo 

      ! This do loop draws any closed storm symbols
      ! We do this separately so that we don't have to change fonts more
      !   than twice
      ! This loop is only entered if it is necessary to draw the
      !   closed storm symbols, to save a possible font change
      do kk=1,MAX_STORMS 
        track_read=.false.
        do iii=1, nxtavl 
          if(storm_val(iii,5,kk).ge.storm_cutoff2)
     &       track_read=.true. 
        enddo 
        if(track_read) then  
          call pcseti('FN',37)

          do iii=1, nxtavl 

            read(cxtimeavl(iii),'(f6.1)') current_fchr
            if(nint(current_fchr)/rtynt(ipl).eq.
     &         nint(nint(current_fchr)/rtynt(ipl))) then  

              if(storm_val(iii,1,kk).ge.0 .and. 
     &          storm_val(iii,2,kk).ge.0) then 
                rcrag(1,ipl)=(storm_val(iii,1,kk)-1.)/refrat+yicorn
                rcrag(2,ipl)=(storm_val(iii,2,kk)-1.)/refrat+xjcorn

                if(storm_val(iii,5,kk).ge.storm_cutoff2) then 
                ! storm(ipl)(1:1)=char(112) 
                  storm(ipl)(1:1)='p'
                  call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrag,
     &              maxpl,ipl,irota)
                endif	

              endif 

            endif 

          enddo 
        endif 
      enddo 
      call pcseti('FN',old_font)
      rtslb(ipl) = old_tslb

      ! Get an array that will tell labelling algorithm which labels will be drawn
      do kk = 1, MAX_STORMS 
        do i = 1, nxtavl
          read(cxtimeavl(i),'(f6.1)') current_fchr
          icurrent_fchr=nint(current_fchr)

          ! Only label symbols every LABEL_INT hours or if there is a break in the track
          if(inwin(i,kk) .and. ((mod(icurrent_fchr,LABEL_INT).eq.0) .or. 
     &                          (i.gt.1 .and. .not. inwin(i-1,kk)) .or. 
     &                          (i.lt.nxtavl .and. .not. inwin(i+1,kk)) 
     &      )) then
            do_label(i,kk) = .true.
          else
            do_label(i,kk) = .false.
          end if
        end do
      end do
    
      ! Get an initial position for each label
      call calc_initial(storm_val, prev_label, maxptimes, nxtavl, 
     &                  s_width, s_height, label_dist, ixwin(1,ipl), 
     &                  ixwin(2,ipl), iywin(1,ipl), iywin(2,ipl))

      ! Adjust label positions for a set number of iterations
      do kk=1,NUM_ITERATIONS
        call iterate(storm_val, prev_label, do_label, maxptimes, nxtavl,
     &               s_width, s_height, label_dist, ixwin(1,ipl), 
     &               ixwin(2,ipl), iywin(1,ipl), iywin(2,ipl), inwin)
      enddo

      ! finally, draw the labels and lines
      do kk=1, MAX_STORMS
      do iii=1, nxtavl
      if(do_label(iii,kk)) then

        ydif=prev_label(iii,1,kk)-storm_val(iii,1,kk)
        xdif=prev_label(iii,2,kk)-storm_val(iii,2,kk)
        if(abs(ydif).lt..0001) ydif = sign(.0001,ydif)
        theta=atan(ydif/xdif)

        if(xdif.gt.0) then
         if(ydif.gt.0) then
          prev_text(iii,1) = prev_label(iii,1,kk)+
     &     inst*abs(sin(theta))/2.7
          prev_text(iii,2) = prev_label(iii,2,kk)+
     &     inst*abs(cos(theta))/2.7
         else
          prev_text(iii,1) = prev_label(iii,1,kk)-
     &     inst*abs(sin(theta))/2.7
          prev_text(iii,2) = prev_label(iii,2,kk)+
     &     inst*abs(cos(theta))/2.7
         endif
        else
         if(ydif.gt.0) then
          prev_text(iii,1) = prev_label(iii,1,kk)+
     &     inst*abs(sin(theta))/2.7
          prev_text(iii,2) = prev_label(iii,2,kk)-
     &     inst*abs(cos(theta))/2.7
         else
          prev_text(iii,1) = prev_label(iii,1,kk)-
     &     inst*abs(sin(theta))/2.7
          prev_text(iii,2) = prev_label(iii,2,kk)-
     &     inst*abs(cos(theta))/2.7
         endif
        endif

        read(cxtimeavl(iii),'(f6.1)') current_fchr
        icurrent_fchr=nint(current_fchr)

cwrite(6,*) 'current_fchr = ',current_fchr,' icurrent_fchr = ',
c    & icurrent_fchr
cwrite(6,*) 'calling mconvert, init_hr = ',init_hr
cwrite(6,*) 'calling mconvert, init_date = ',init_date
        call mconvert(init_date,init_hr,1,1970)
        init_hr = init_hr + icurrent_fchr
cwrite(6,*) 'init_hr = ',init_hr
        call mconvert(new_date,init_hr,-1,1970)
        write(storm(ipl)(1:5),'(i2.2,a1,i2.2)')
     &    (new_date-(new_date/10000)*10000)/100,
     &    '/',mod(new_date,100)
cwrite(6,*) 'ipl = ',ipl,' new_date = ',new_date
cwrite(6,*) 'day = ',(new_date-(new_date/10000)*10000)/100
cwrite(6,*) 'init_hr = ',init_hr
cwrite(6,*) 'init_date = ',init_date
        iscp=nint(storm_val(iii,4,kk))

        if(storm_val(iii,1,kk).ge.0 .and.
     &     storm_val(iii,2,kk).ge.0) then
          rcrag(1,ipl)=(storm_val(iii,1,kk)-1.)/refrat+yicorn
          rcrag(2,ipl)=(storm_val(iii,2,kk)-1.)/refrat+xjcorn
          rcrbg(1,ipl)=(prev_label(iii,1,kk)-1.)/refrat+yicorn
          rcrbg(2,ipl)=(prev_label(iii,2,kk)-1.)/refrat+xjcorn
          call hlinedraw(ilinw,idash,ixwin,iywin,icolr,rcrag,rcrbg,
     &       maxpl,ipl,irota)
          rcrbg(1,ipl)=(prev_text(iii,1)-1.)/refrat+yicorn
          rcrbg(2,ipl)=(prev_text(iii,2)-1.)/refrat+xjcorn
          call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrbg,
     &       maxpl,ipl,irota)
          rcrbg(1,ipl) = rcrbg(1,ipl)-iywin(2,ipl)
     &                   *rtslb(ipl)*2.1/refrat
          if(iscp.ge.1000) then
            write(storm(ipl)(1:5),'(i4.4,a1)') iscp,' '
          else
            write(storm(ipl)(1:5),'(a1,i3.3,a1)') 
     &       ' ',iscp,' '
          endif
          call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrbg,
     &       maxpl,ipl,irota)
        endif

      endif
      enddo
      enddo

      ! Also, we should make a legend in one of the left corners of the
      !   plot window describing what the text of each label represents
      old_tslb = rtslb(ipl)
      rtslb(ipl) = .010
      downset = real(s_height)*rtslb(ipl)*3./refrat

      ! First determine if we can safely put the info in the lower left.
      !   If we can't, then the info will be put in the upper left.
      if(s_width.gt.s_height) then
        insetx = real(s_width)/6.5
      else
        insetx = real(s_height)/6.5
      endif

      go_ul = .false.
      do i=1,nxtavl
        do j=1,MAX_STORMS
          if(storm_val(i,1,j).le.(iywin(1,ipl)+insetx) .and.
     &       storm_val(i,1,j).ge.(iywin(1,ipl)-1.)) then
            if(storm_val(i,2,j).le.(ixwin(1,ipl)+insetx) .and.
     &         storm_val(i,2,j).ge.(ixwin(1,ipl)-1.)) then
              go_ul = .true.
            endif
          endif
          if(prev_label(i,1,j).le.(iywin(1,ipl)+insetx) .and.
     &       prev_label(i,1,j).ge.(iywin(1,ipl)-1.) .and.
     &       storm_val(i,1,j).gt.0 .and. storm_val(i,2,j).gt.0) then
            if(prev_label(i,2,j).le.(ixwin(1,ipl)+insetx) .and.
     &         prev_label(i,2,j).ge.(ixwin(1,ipl)-1.)) then
              go_ul = .true.
            endif
          endif
        enddo
      enddo

      if(go_ul) then
        rcrag(1,ipl)=(iywin(2,ipl)-insetx/2.-1.)/refrat+yicorn+
     &                downset
        rcrag(2,ipl)=(ixwin(1,ipl)+insetx/2.-1.)/refrat+xjcorn
        rcrbg(1,ipl)=rcrag(1,ipl)-downset
        rcrbg(2,ipl)=rcrag(2,ipl)
      else
        rcrag(1,ipl)=(iywin(1,ipl)+insetx/2.-1.)/refrat+yicorn
        rcrag(2,ipl)=(ixwin(1,ipl)+insetx/2.-1.)/refrat+xjcorn
        rcrbg(1,ipl)=rcrag(1,ipl)-downset
        rcrbg(2,ipl)=rcrag(2,ipl)
      endif

      write(storm(ipl)(1:8),'(a8)') 'Day/Hour'
      call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrag,
     &         maxpl,ipl,irota)

      write(storm(ipl)(1:8),'(a8)') 'Min SLP '
      call hbulldraw(rtslb,storm,ixwin,iywin,icolr,rcrbg,
     &         maxpl,ipl,irota)
      rtslb(ipl) = old_tslb

      rtslb(ipl) = rtslb(ipl) * 2
      ilinw(ipl) = old_linw
      storm(ipl)(1:13)='               '
      call pcseti('CC',old_color)

  900 continue

      return
      end	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uncross(points, labels, do_label, maxptimes, nxtavl)

      implicit none

      integer MAX_STORMS
      integer MAX_PARAMETERS
      parameter (MAX_STORMS = 5)     ! The 3rd dimension of points()
      parameter (MAX_PARAMETERS = 7) ! The 2nd dimension of points()

      ! Arguments
      integer maxptimes, nxtavl
      real points(maxptimes,MAX_PARAMETERS,MAX_STORMS)
      real labels(maxptimes,2,MAX_STORMS)
      logical do_label(maxptimes,MAX_STORMS)

      ! local variables
      integer i, j, k, m, n
      real m1, m2, b1, b2, x, temp

      do k=1, MAX_STORMS
        do i=1, nxtavl
          do m=k, MAX_STORMS
            do j=1, nxtavl
            if ((i .ne. j) .and. do_label(j,m) .and. do_label(i,k)) then

            ! y=mx+b
            if (points(i,2,k) .ne. labels(i,2,k)) then
              m1 = (points(i,1,k)-labels(i,1,k))/
     &             (points(i,2,k)-labels(i,2,k)) 
            else
              m1 = 10000.
            end if

            if (points(j,2,m) .ne. labels(j,2,m)) then
              m2 = (points(j,1,m)-labels(j,1,m))/
     &             (points(j,2,m)-labels(j,2,m)) 
            else
              m2 = 10000.
            end if

            n = 0

            ! Check for intersection if segments are not parallel
            if (m1 .ne. m2) then
              b1 = points(i,1,k) - (m1 * points(i,2,k))
              b2 = points(j,1,m) - (m2 * points(j,2,m))
              x = (b2-b1)/(m1-m2)

              ! If x-coordinate of intersection is between the x-coordinates of 
              !   endpoints of first segment
              if (points(i,2,k) .gt. labels(i,2,k)) then
                if ((x .ge. labels(i,2,k)) .and. (x .le. points(i,2,k)))
     &            n = n + 1
              else
                if ((x .le. labels(i,2,k)) .and. (x .ge. points(i,2,k)))
     &            n = n + 1
              end if

              ! If x-coordinate of intersection is between the x-coordinates of 
              !   endpoints of second segment
              if (points(j,2,m) .gt. labels(j,2,m)) then
                if ((x .ge. labels(j,2,m)) .and. (x .le. points(j,2,m)))
     &            n = n + 1
              else
                if ((x .le. labels(j,2,m)) .and. (x .ge. points(j,2,m)))
     &            n = n + 1
              end if

              ! Houston, we have an intersection...
              if (n .eq. 2) then
                temp = labels(j,2,m)
                labels(j,2,m) = labels(i,2,k)
                labels(i,2,k) = temp
                temp = labels(j,1,m)
                labels(j,1,m) = labels(i,1,k)
                labels(i,1,k) = temp
              end if
            end if
            end if
            end do
          end do
        end do
      end do

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_initial(points, labels, maxptimes, nxtavl, width, 
     &                        height, label_dist, xmin, xmax, ymin, 
     &                        ymax)

      integer MAX_STORMS
      integer MAX_PARAMETERS
      real MASS
      real DT
      real LABEL_REPULSION
      real POINT_REPULSION
      real EDGE_REPULSION
      real LABEL_ATTRACTION
      real POINT_ATTRACTION
      parameter (MAX_STORMS = 5)         ! The 3rd dimension of points()
      parameter (MAX_PARAMETERS = 7)     ! The 2nd dimension of points()
      parameter (MASS = 1.0)             ! The "mass" of each label
      parameter (DT = 0.12)              ! The squared time step for each iteration
      parameter (LABEL_REPULSION = 4.0)  ! How much label is repelled by other labels
      parameter (POINT_REPULSION = 2.0)  ! How much label is repelled by track points 
      parameter (EDGE_REPULSION = 4.0)   ! How much label is repelled by window edges 
      parameter (LABEL_ATTRACTION = 0.5) ! How much label is attracted to neighbor labels 
      parameter (POINT_ATTRACTION = 0.5) ! How much label is attracted to point it labels

      ! Arguments
      integer maxptimes, nxtavl, width, height
      integer xmin, xmax, ymin, ymax
      real label_dist
      real points(maxptimes,MAX_PARAMETERS,MAX_STORMS)
      real labels(maxptimes,2,MAX_STORMS)

      ! Local variables
      integer i, j, k, m, n
      real theta, slope, fx, fy, ax, ay, ave_x, ave_y, distance

      ! Initial position of labels, before adjustment, should be just beside the point they label
      do j=1,MAX_STORMS 
        do i=1,nxtavl 
          labels(i,2,j) = points(i,2,j)+0.5
          labels(i,1,j) = points(i,1,j)+0.5

          fx = 0.0
          fy = 0.0

          ! Calculate repulsion from edges of plot
          if (labels(i,2,j) .lt. (xmin + label_dist))
     &      fx = fx + EDGE_REPULSION*(xmin + label_dist - labels(i,2,j))
          if (labels(i,2,j) .gt. (xmax - label_dist))
     &      fx = fx - EDGE_REPULSION*(label_dist - 
     &           (xmax - labels(i,2,j)))
          if (labels(i,1,j) .lt. (ymin + label_dist))
     &      fy = fy + EDGE_REPULSION*(ymin + label_dist - labels(i,1,j))
          if (labels(i,1,j) .gt. (ymax - label_dist))
     &      fy = fy - EDGE_REPULSION*(label_dist - 
     &           (ymax - labels(i,1,j)))

          ! F=ma
          ax = fx / MASS
          ay = fy / MASS

          ! Find new position after DT seconds
          labels(i,2,j) = 0.5*ax*DT + labels(i,2,j)
          labels(i,1,j) = 0.5*ay*DT + labels(i,1,j)
        end do 
      end do 

      ! Need initial positions for labels on all typhoons
      do j=1,MAX_STORMS 
        ! All other storms affect the initial placement of current typhoon's labels
        do m=1,MAX_STORMS 
          if(m .ne. j) then
            ! Get average location of points on the other typhoon being considered
            ave_x = 0.0 
            ave_y = 0.0
            n = 0
            do k=1,nxtavl 
              ave_x = ave_x + points(k,2,m)
              ave_y = ave_y + points(k,1,m)
              n = n + 1
            end do 
            if (n .gt. 0) then
              ave_x = ave_x / n
              ave_y = ave_y / n
            end if
  
            ! Now see how a point is affected by the other typhoon
            do i=1,nxtavl 
  
              fx = 0.0
              fy = 0.0
              distance = fdist(points(i,2,j),points(i,1,j),ave_x,ave_y)
              if (distance .lt. 4.0*label_dist) then 
                if (points(i,2,j) .ne. ave_x) then
                  theta = atan((points(i,1,j) - ave_y)/
     &                         (points(i,2,j) - ave_x))
                else
                  theta = atan(1.0/0.001)
                end if
  
                ! Get x and y components of force. We want to nudge a label in a direction away from other
                !  typhoons
                if (points(i,2,j) .ge. ave_x) then
                  fx = fx + POINT_REPULSION/3.0 * 
     &                      (3.0*label_dist - distance) * cos(theta)
                  fy = fy + POINT_REPULSION/3.0 * 
     &                      (3.0*label_dist - distance) * sin(theta)
                else
                  fx = fx + POINT_REPULSION/3.0 * 
     &                      (distance - 3.0*label_dist) * cos(theta)
                  fy = fy + POINT_REPULSION/3.0 * 
     &                      (distance - 3.0*label_dist) * sin(theta)
                end if
              end if 
  
              ! F=ma
              ax = fx / MASS
              ay = fy / MASS
  
              ! Find new position after DT seconds
              labels(i,2,j) = 0.5*ax*DT + labels(i,2,j)
              labels(i,1,j) = 0.5*ay*DT + labels(i,1,j)
            end do 
          end if
        end do 
      end do 

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine iterate(points, labels, do_label, maxptimes, nxtavl, 
     &  width, height, label_dist, xmin, xmax, ymin, ymax, inwin)
 
      integer MAX_STORMS
      integer MAX_PARAMETERS
      real MASS
      real DT
      real LABEL_REPULSION
      real POINT_REPULSION
      real EDGE_REPULSION
      real LABEL_ATTRACTION
      real POINT_ATTRACTION
      parameter (MAX_STORMS = 5)         ! The 3rd dimension of points()
      parameter (MAX_PARAMETERS = 7)     ! The 2nd dimension of points()
      parameter (MASS = 1.0)             ! The "mass" of each label
      parameter (DT = 0.12)              ! The squared time step for each iteration
      parameter (LABEL_REPULSION = 4.0)  ! How much label is repelled by other labels
      parameter (POINT_REPULSION = 2.0)  ! How much label is repelled by track points 
      parameter (EDGE_REPULSION = 4.0)   ! How much label is repelled by window edges 
      parameter (LABEL_ATTRACTION = 0.5) ! How much label is attracted to neighbor labels 
      parameter (POINT_ATTRACTION = 0.5) ! How much label is attracted to point it labels

      ! Arguments
      integer maxptimes, nxtavl, width, height
      integer xmin, xmax, ymin, ymax
      real label_dist
      real points(maxptimes,MAX_PARAMETERS,MAX_STORMS)
      real labels(maxptimes,2,MAX_STORMS)
      logical do_label(maxptimes,MAX_STORMS)
      logical inwin(maxptimes, MAX_STORMS)

      ! Local variables
      integer i, j, k, m
      real distance, point_distance, fx, fy, force, theta, ax, ay
      real tfx, tfy

      do k=1,MAX_STORMS !{
      do i=1,nxtavl !{
        if (inwin(i,k)) then
        fx = 0.0
        fy = 0.0

        ! Calculate repulsion from points on the typhoon and other labels
        do j=1,nxtavl !{
          do m=1,MAX_STORMS !{
            ! Repulsion from points on the typhoon track
              distance = fdist(labels(i,2,k),labels(i,1,k),
     &                        points(j,2,m),points(j,1,m))
              if (distance .lt. label_dist) then !{
                if (labels(i,2,k) .ne. points(j,2,m)) then
                  theta = atan((labels(i,1,k) - points(j,1,m))/
     &                            (labels(i,2,k) - points(j,2,m)))
                else
                  theta = atan(1.0/0.001)
                end if

                ! Get x and y components of force
                if (labels(i,2,k) .ge. points(j,2,m)) then !{
                  tfx = POINT_REPULSION * (label_dist - distance) 
     &                  * cos(theta)
                  tfy = POINT_REPULSION * (label_dist - distance)
     &                 * sin(theta)
                  fx = fx + tfx
                  fy = fy + tfy
c                 write(6,*) i,k,'Force ',tfx,tfy,' from point ',j,m
                !}
                else !{
                  tfx = POINT_REPULSION * (distance - label_dist) 
     &                  * cos(theta)
                  tfy = POINT_REPULSION * (distance - label_dist)
     &                 * sin(theta)
                  fx = fx + tfx
                  fy = fy + tfy
c                 write(6,*) i,k,'Force ',tfx,tfy,' from point ',j,m
                end if !}
              end if !}

            ! Repulsion from other labels
            if ( ((i .ne. j) .or. (m .ne. k)) .and. do_label(i,k) 
     &           .and. do_label(j,m) ) then !{
              distance = fdist(labels(i,2,k),labels(i,1,k),
     &                        labels(j,2,m),labels(j,1,m))
              if (distance .lt. label_dist) then !{
                if (labels(i,2,k) .ne. labels(j,2,m)) then
                  theta = atan((labels(i,1,k) - labels(j,1,m))/
     &                              (labels(i,2,k) - labels(j,2,m)))
                else
                  theta = atan(1.0/0.001)
                end if
  
                ! Get x and y components of force
                if (labels(i,2,k) .ge. labels(j,2,m)) then !{
                  tfx = LABEL_REPULSION * (label_dist - distance) 
     &                  * cos(theta)
                  tfy = LABEL_REPULSION * (label_dist - distance)
     &                 * sin(theta)
                  fx = fx + tfx
                  fy = fy + tfy
c                 write(6,*) i,k,'Force ',tfx,tfy,' from label ',j,m
                !}
                else !{
                  tfx = LABEL_REPULSION * (distance - label_dist) 
     &                  * cos(theta)
                  tfy = LABEL_REPULSION * (distance - label_dist)
     &                 * sin(theta)
                  fx = fx + tfx
                  fy = fy + tfy
c                 write(6,*) i,k,'Force ',tfx,tfy,' from label ',j,m
                end if !}
              end if !}
            end if !}
          end do !}
        end do !}

        ! Calculate repulsion from edges of plot
        if (labels(i,2,k) .lt. (xmin + label_dist))
     &    fx = fx + EDGE_REPULSION*(xmin + label_dist - labels(i,2,k))
        if (labels(i,2,k) .gt. (xmax - label_dist))
     &    fx = fx - EDGE_REPULSION*(label_dist - 
     &         (xmax - labels(i,2,k)))
        if (labels(i,1,k) .lt. (ymin + label_dist))
     &    fy = fy + EDGE_REPULSION*(ymin + label_dist - labels(i,1,k))
        if (labels(i,1,k) .gt. (ymax - label_dist))
     &    fy = fy - EDGE_REPULSION*(label_dist - 
     &         (ymax - labels(i,1,k)))

        ! Calculate attractive force from a label's corresponding point on the typhoon
        distance = fdist(labels(i,2,k),labels(i,1,k),
     &                   points(i,2,k),points(i,1,k))
        if (distance .gt. label_dist) then !{
          if (labels(i,2,k) .ne. points(i,2,k)) then
            theta = atan((labels(i,1,k) - points(i,1,k))/
     &                        (labels(i,2,k) - points(i,2,k)))
          else
            theta = atan(1.0/0.001)
          end if

          ! Get x and y components of force
          if (labels(i,2,k) .ge. points(i,2,k)) then !{
                  tfx = POINT_ATTRACTION * (label_dist - distance) 
     &                  * cos(theta)
                  tfy = POINT_ATTRACTION * (label_dist - distance)
     &                 * sin(theta)
            fx = fx + tfx
            fy = fy + tfy
c           write(6,*) i,k,'Force ',tfx,tfy,' from own label '
          !}
          else !{
                  tfx = POINT_ATTRACTION * (distance - label_dist) 
     &                  * cos(theta)
                  tfy = POINT_ATTRACTION * (distance - label_dist)
     &                 * sin(theta)
            fx = fx + tfx
            fy = fy + tfy
c           write(6,*) i,k,'Force ',tfx,tfy,' from own label '
          end if !}
        end if !}

        ! Calculate force from neighboring labels
        if (i .gt. 1) then !{
          point_distance = fdist(points(i,2,k),points(i,1,k),
     &                          points(i-1,2,k),points(i-1,1,k))
          distance = fdist(labels(i,2,k),labels(i,1,k),
     &                    labels(i-1,2,k),labels(i-1,1,k))
          if (distance .gt. point_distance) then !{
            if (labels(i,2,k) .ne. labels(i-1,2,k)) then
              theta = atan((labels(i,1,k) - labels(i-1,1,k))/
     &                          (labels(i,2,k) - labels(i-1,2,k)))
            else
              theta = atan(1.0/0.001)
            end if

            ! Get x and y components of force
            if (labels(i,2,k) .ge. labels(i-1,2,k)) 
     &      then !{
                  tfx = LABEL_ATTRACTION * (point_distance - distance) 
     &                  * cos(theta)
                  tfy = LABEL_ATTRACTION * (point_distance - distance)
     &                 * sin(theta)
              fx = fx + tfx
              fy = fy + tfy
c         write(6,*) i,k,'Force ',tfx,tfy,' from neighboring label ',
c    &          i-1,k
            !}
            else !{
                  tfx = LABEL_ATTRACTION * (distance - point_distance) 
     &                  * cos(theta)
                  tfy = LABEL_ATTRACTION * (distance - point_distance)
     &                 * sin(theta)
              fx = fx + tfx
              fy = fy + tfy
c         write(6,*) i,k,'Force ',tfx,tfy,' from neighboring label ',
c    &          i-1,k
            end if !}
          end if !}
        end if !}
        if (i .lt. nxtavl) then !{
          point_distance = fdist(points(i,2,k),points(i,1,k),
     &                          points(i+1,2,k),points(i+1,1,k))
          distance = fdist(labels(i,2,k),labels(i,1,k),
     &                    labels(i+1,2,k),labels(i+1,1,k))
          if (distance .gt. point_distance) then !{
            if (labels(i,2,k) .ne. labels(i+1,2,k)) then
              theta = atan((labels(i,1,k) - labels(i+1,1,k))/
     &                          (labels(i,2,k) - labels(i+1,2,k)))
            else
              theta = atan(1.0/0.001)
            end if

            ! Get x and y components of force
            if (labels(i,2,k) .ge. labels(i+1,2,k)) then !{
                  tfx = LABEL_ATTRACTION * (point_distance - distance) 
     &                  * cos(theta)
                  tfy = LABEL_ATTRACTION * (point_distance - distance)
     &                 * sin(theta)
              fx = fx + tfx
              fy = fy + tfy
c         write(6,*) i,k,'Force ',tfx,tfy,' from neighboring label ',
c    &          i+1,k
            !}
            else !{
                  tfx = LABEL_ATTRACTION * (distance - point_distance) 
     &                  * cos(theta)
                  tfy = LABEL_ATTRACTION * (distance - point_distance)
     &                 * sin(theta)
              fx = fx + tfx
              fy = fy + tfy
c         write(6,*) i,k,'Force ',tfx,tfy,' from neighboring label ',
c    &          i+1,k
            end if !}
          end if !}
        end if !}

        ! F=ma
c       write(6,*) 'FX=',fx,' FY=',fy
        if (abs(fx) .gt. 700.) then
          if (fx .lt. 0) then
            fx = -700.
          else
            fx = 700.
          endif
        endif
        if (abs(fy) .gt. 700.) then
          if (fy .lt. 0) then
            fy = -700.
          else
            fy = 700.
          endif
        endif
        ax = fx / MASS
        ay = fy / MASS

        ! Find new position after DT seconds
        labels(i,2,k) = labels(i,2,k) + 0.5*ax*DT
        labels(i,1,k) = labels(i,1,k) + 0.5*ay*DT

        call uncross(points, labels, do_label, maxptimes, nxtavl)
        end if
      end do !}
      end do !}
 
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fdist(x1, y1, x2, y2)

      real x1, y1, x2, y2
      real fdist
    
      fdist = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

      return
      end
