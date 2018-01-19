c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine profvelcalc(comp,fname,mdate,rhour,xtime,
     &   cosa,sina,strmvel,rip_root,unorth,vnorth,ght,vel,
     &   miy,mjx,mkzh)
c
c     This routine reads in profiler data, and analyzes it to a strip of
c     model grid points by means of time-to-space translation of the
c     profiler data.  The strip of model grid points straddle a straight
c     line defined in the usual manner by crsa and crsb (which
c     subsequently determine the values of cosa and sina).  It is assumed that
c     there is some fixed 2D structure passing by the station at a
c     constant speed of strmvel (allowing time-space translation).  The
c     cross section should be chosen to be perpendicular to the 2D
c     feature, passing through the profiler station.
c
      dimension unorth(miy,mjx),vnorth(miy,mjx),ght(miy,mjx,mkzh),
     &   vel(miy,mjx,mkzh)
      character comp*1, fname*256, rip_root*256
      character string(2)*20,icaoid*4,locdesc*44
      dimension distpr(300),nz(300),ht(100,300),
     &   uwly(100,300),usly(100,300)
c
      parameter (rinft=600.,rinfz=500.) ! rad infl for time(s) and ht(m)
c
      include 'comconst'
c
      rinfgp=rinft*strmvel/ds
c
      open(unit=iutraj,file=fname,form='formatted',status='old')
      read(iutraj,'(7x,a4)')string(1)
      string(2)='missing             '
      call locinterp (string,gridx,gridy,rlat,rlon,iwmo,icaoid,
     &   stelev,locdesc,rip_root)
      gridx=1.+(gridx-xjcorn)*refrat
      gridy=1.+(gridy-yicorn)*refrat
c
      read(iutraj,*)
c
      is=0
 100  continue   ! top of time (or distance) loop for reading in data
c
      is=is+1
      read(iutraj,'(21x,3(1x,i2),1x,2(1x,i2))',end=400)
     &   imo,idy,iyr,ihr,imn
      mdatepr=1000000*iyr+10000*imo+
     &   100*idy+ihr
      rhourpr=imn/60.
      call mconvert(mdatepr,mhourpr,1,1940)
      xtimepr=float(mhourpr-mhourb)+rhourpr-rhourb
c
c   Assign a "distance" (in grid intervals) for this time using time-to-space
c   correction.  Note distpr=0. refers to the station location.
c
      distpr(is)=strmvel*(xtime-xtimepr)*3600./ds
      read(iutraj,*)
      read(iutraj,*)
c
      iz=0
 200  continue   ! top of height loop for reading in data
c
      iz=iz+1
 250  read(iutraj,'(3x,f4.0,2f8.0)')htftagl,uwly(iz,is),usly(iz,is)
      if (htftagl.eq.9e9) goto 300   ! end of this time period
      if (htftagl.ge.900.) goto 250  ! skip entry edited for bad quality
      ht(iz,is)=htftagl*30.48 + stelev ! convert from 100s ft AGL to m AMSL
      goto 200
 300  nz(is)=iz-1
c
      goto 100
c
 400  ns=is-1
c
c     OK, we have the profiler data, now we loop through all model grid
c     points to analyze the profiler data to the strip of model grid points.
c
      phi=rpd*(90.-strmdir)
c
      do k=1,mkzh
      do j=1,mjx
      do i=1,miy
c
      bx=float(j)-gridx
      by=float(i)-gridy
c
c   Calculate distance from the station (in gr.points) of the projection
c   of this grid point onto the cross section
c
      distgp=cosa*bx+sina*by
c
c   Calculate the normal distance from this grid point to the cross section
c
      doffline=abs(cosa*by-sina*bx)
c
c   Since we're only interested in a cross sectiion through the
c   profiler site that runs parallel to storm motion, we will only
c   analyze to points that are within 3 grid points of that line.
c
      if (doffline.gt.3) then
         vel(i,j,k)=rmsg
         goto 500
      endif
c
c   We now have height of the grid point (ght) and distance (distgp)
c   of the grid point.  We have the same for
c   all the profiler data points (ht and distpr), so we can do an "analysis".
c
      totwly=0.
      totsly=0.
      totwt=0.
      do is=1,ns
c
c      Calculate hor distance from grid point to profiler data point,
c      normalized by radius of influence.
c
         ddhor=(distpr(is)-distgp)/rinfgp
         if (ddhor.le.1.) then   ! within rad. of inf, so include in anal.
            do iz=1,nz(is)
c
c            Calculate ver distance from grid point to profiler data point,
c            normalized by radius of influence.
c
               ddver=(ht(iz,is)-ght(i,j,k))/rinfz
               if (ddver.le.1.) then
                  dd2d=sqrt(ddhor*ddhor+ddver*ddver)
                  wt=max(0.,(1.-dd2d)/(1.+dd2d))
                  totwly=totwly+wt*uwly(iz,is)
                  totsly=totsly+wt*usly(iz,is)
                  totwt=totwt+wt
               endif
            enddo
         endif
      enddo
c
 33   if (totwt.gt.0.) then
         uwlygp=totwly/totwt
         uslygp=totsly/totwt
         if (comp.eq.'u') then
            vel(i,j,k)=vnorth(i,j)*uwlygp+unorth(i,j)*uslygp
         elseif (comp.eq.'v') then
            vel(i,j,k)=vnorth(i,j)*uslygp-unorth(i,j)*uwlygp
         endif
      else
         vel(i,j,k)=rmsg
      endif
c
 500  continue
c
      enddo
      enddo
      enddo
c
      close (iutraj)
      return
      end
