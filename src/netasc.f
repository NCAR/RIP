      subroutine netasc(casename,iendc,cxtimeavl,xtimeavl,
     &   ncxc,nxtavl,maxtavl,prs1,prs2,eh1,eh2,
     &   ctjfl,rtjst,rtjen,tacch,asc,miy,mjx,mkzh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c   This subroutine reads in a gridswarm trajectory position file and  c
c   calculates the neta ascent for a specified time period at all      c
c   domain grid points that are within the gridswarm grid.             c

c The program assumes:
c
c 1) The entire trajectory file is one 3D gridswarm, i.e. trajectory
c    points lie on a 3D rectangular grid at the initial time.  The routine
c    will work for either forward or backward trajectories.
c 2) The gridswarm points lie on model cross points.  This means the
c    initial trajectory x and y locations should all be an
c    integer + 0.5., and z locations should be model levels (not some
c    other vertical coordinate).
c 3) The grid has a horizontal grid spacing that is some integer
c    multiple of the model grid space.  Interpolation will be used
c    to get net ascent at the model points in between the trajectory
c    grid points.
c 4) The gridswarm vertical levels must be consecutive model levels.
c 5) The trajectories are arranged such that the x-index of the
c    gridswarm varies most rapidly (left to right),
c    then the y (bottom to top), and then the z (top to bottom).
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10, casename*(*)
c
      parameter (maxtraj=7000,maxtrajtime=200)
      dimension prs1(miy,mjx,mkzh),prs2(miy,mjx,mkzh),
     &   eh1(miy,mjx,mkzh),eh2(miy,mjx,mkzh),
     &   asc(miy,mjx,mkzh)
      dimension stortr(maxtrajtime,maxtraj,3),asctr(maxtraj)
      character ctjfl*256
c
      include 'comconst'
c
c   Read in data from trajectory file
c
      open (unit=iutrajin,file=ctjfl,form='unformatted',status='old')
      read (iutrajin)
      read (iutrajin) rtim,ctim,dttraj,ntraj
      ntrajtime=nint(abs(rtim-ctim)/dttraj*3600) + 1
      if (rtim.lt.ctim) then
         itm1=1
         itm2=ntrajtime
         itmi=1
      else
         itm1=ntrajtime
         itm2=1
         itmi=-1
      endif
      if (rtjst.eq.rmsg) then
         starttim=min(rtim,ctim)
      else
         starttim=max(rtjst,min(rtim,ctim))
      endif
      if (rtjen.eq.rmsg) then
         endtim=max(rtim,ctim)
      else
         endtim=min(rtjen,max(rtim,ctim))
      endif
c
c   Check if starttim, endtim, and release time are available
c
      nxt_start=-1
      nxt_end=-1
      nxt_rel=-1
      do i=1,nxtavl
         if (abs(xtimeavl(i)-starttim).le.tacch) then
            nxt_start=i
         endif
         if (abs(xtimeavl(i)-endtim).le.tacch) then
            nxt_end=i
         endif
         if (abs(xtimeavl(i)-rtim).le.tacch) then
            nxt_rtim=i
         endif
      enddo
      if (nxt_start.eq.-1.or.nxt_end.eq.-1) then
         write(iup,*)'In netasc, either start or end time requested',
     &      ' is not available. Stopping.'
         stop
      endif
c
c   Trj. data is read in so that it is always in chronologically forward
c   direction, regarless of whether trajectories were calculated
c   backward or forward.  itm1 will be the index of the release time
c   in the stortr array, and itm2 the index of the complete time.
c
      do itm=itm1,itm2,itmi
         read(iutrajin) (stortr(itm,itr,1),itr=1,ntraj),
     &       (stortr(itm,itr,2),itr=1,ntraj),
     &       (stortr(itm,itr,3),itr=1,ntraj)
      enddo
      close (iutrajin)
      itma1=1+nint((starttim-min(rtim,ctim))/dttraj*3600)
      itma2=itma1+nint((endtim-starttim)/dttraj*3600)
c
c   Read in yicorn, xjcorn, pressure, and height data at earliest time
c   in net ascent window.
c
      call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt_start,maxtavl,ncxc,nproj,
     &   miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn_start,xjcorn_start,
     &   rhourb,dsc,ds,refrat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_start,ncxc,
     &     'prs       ',miy,mjx,mkzh,maxtavl,3,1,prs1,
     &     istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_start,ncxc,
     &     'ght       ',miy,mjx,mkzh,maxtavl,3,1,eh1,
     &     istat)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         eh1(i,j,k)=exp(-eh1(i,j,k)/sclht)
      enddo
      enddo
      enddo
c
c   Read in yicorn, xjcorn, pressure, and height data at latest time
c   in net ascent window.
c
      call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt_end,maxtavl,ncxc,nproj,
     &   miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn_end,xjcorn_end,
     &   rhourb,dsc,ds,refrat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_end,ncxc,
     &     'prs       ',miy,mjx,mkzh,maxtavl,3,1,prs2,
     &     istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_end,ncxc,
     &     'ght       ',miy,mjx,mkzh,maxtavl,3,1,eh2,
     &     istat)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         eh2(i,j,k)=exp(-eh2(i,j,k)/sclht)
      enddo
      enddo
      enddo
c
c   Read in yicorn and xjcorn at release time.
c
      call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt_rtim,maxtavl,ncxc,nproj,
     &   miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn_rtim,xjcorn_rtim,
     &   rhourb,dsc,ds,refrat)
c
c   Get pressure at trajectory locations.
c   Also, convert x and y values (which are in coarse dom. grid points)
c   to grid points in the current domain.
c
      do itr=1,ntraj
         if (stortr(itm1,itr,1).eq.rmsg.or.
     &       stortr(itm1,itr,2).eq.rmsg.or.
     &       stortr(itm1,itr,3).eq.rmsg) then
            write(iup,*)'problem with position of a',
     &         ' trajectory at rtim.'
            write(iup,*)'itr,st1,st2,st3=',itr,stortr(itm1,itr,1),
     &         stortr(itm1,itr,2),stortr(itm1,itr,3)
            stop
         elseif (itm1.ne.itma1.and.itm1.ne.itma2) then
            stortr(itm1,itr,1)=
     &         1.+refrat*(stortr(itm1,itr,1)-yicorn_rtim)
            stortr(itm1,itr,2)=
     &         1.+refrat*(stortr(itm1,itr,2)-xjcorn_rtim)
         endif
         if (stortr(itma1,itr,1).ne.rmsg) then
            pr1=finterp(eh1,prs1,1,miy,mjx,mkzh,
     &         stortr(itma1,itr,1),stortr(itma1,itr,2),
     &         stortr(itma1,itr,3),refrat,yicorn_start,xjcorn_start,
     &         rmsg,iup)
            if (pr1.eq.rmsg) then
               stortr(itma1,itr,1)=rmsg
               stortr(itma1,itr,2)=rmsg
               stortr(itma1,itr,3)=rmsg
            else
               stortr(itma1,itr,1)=
     &            1.+refrat*(stortr(itma1,itr,1)-yicorn_start)
               stortr(itma1,itr,2)=
     &            1.+refrat*(stortr(itma1,itr,2)-xjcorn_start)
            endif
         endif
         if (stortr(itma2,itr,1).ne.rmsg) then
            pr2=finterp(eh2,prs2,1,miy,mjx,mkzh,
     &         stortr(itma2,itr,1),stortr(itma2,itr,2),
     &         stortr(itma2,itr,3),refrat,yicorn_end,xjcorn_end,
     &         rmsg,iup)
            if (pr2.eq.rmsg) then
               stortr(itma2,itr,1)=rmsg
               stortr(itma2,itr,2)=rmsg
               stortr(itma2,itr,3)=rmsg
            else
               stortr(itma2,itr,1)=
     &            1.+refrat*(stortr(itma2,itr,1)-yicorn_end)
               stortr(itma2,itr,2)=
     &            1.+refrat*(stortr(itma2,itr,2)-xjcorn_end)
            endif
         endif
         if (stortr(itma1,itr,1).ne.rmsg.and.
     &       stortr(itma2,itr,1).ne.rmsg) then
            asctr(itr)=pr2-pr1
         else
            asctr(itr)=rmsg
         endif
      enddo
c
c   Figure out horizontal grid space and dimensions of gridswarm.
c
      iytest=abs(nint(100.*(stortr(itm1,1,1)-nint(stortr(itm1,1,1)))))
      jxtest=abs(nint(100.*(stortr(itm1,1,2)-nint(stortr(itm1,1,2)))))
      if (iytest.ne.50.or.jxtest.ne.50) then
         write(iup,*)'In netasc: trj. swarm corner not on cross point.'
         write(iup,*)'x,y vals. = ',stortr(itm1,1,2),stortr(itm1,1,1)
         write(iup,*)'jxtest,iytest=',jxtest,iytest
         stop
      endif
      iorientest=nint(100.*(stortr(itm1,2,1)-stortr(itm1,1,1)))
      if (iorientest.ne.0) then
         write(iup,*)'In netasc: trj. swarm not properly oriented.'
         write(iup,*)'First two y vals. = ',stortr(itm1,1,1),
     &      stortr(itm1,2,1)
         stop
      endif
      ids=nint(stortr(itm1,2,2)-stortr(itm1,1,2))
      iycorn=nint(stortr(itm1,1,1)-.5)
      jxcorn=nint(stortr(itm1,1,2)-.5)
      igotnx=0
      do itr=2,ntraj
         if (igotnx.eq.0.and.
     &       stortr(itm1,itr,2).lt.stortr(itm1,itr-1,2)) then
            igotnx=1
            nx=itr-1
         endif
         if (stortr(itm1,itr,1).lt.stortr(itm1,itr-1,1)-.01) then
            ny=(itr-1)/nx !prob
            goto 33
         endif
      enddo
 33   continue
c
c   Read in height data at release time.
c
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_rtim,ncxc,
     &     'ght       ',miy,mjx,mkzh,maxtavl,3,1,eh2,
     &     istat)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         eh2(i,j,k)=exp(-eh2(i,j,k)/sclht)
      enddo
      enddo
      enddo
c
c   Figure out vertical grid space and dimensions of gridswarm.
c
      ehcornertop=stortr(itm1,1,3)
      ehcornerbot=stortr(itm1,ntraj-nx*ny+1,3)
      itopct=0
      ibotct=0
      do k=1,mkzh
         if (abs(ehcornertop-eh2(iycorn,jxcorn,k)).lt..00005) then ! about 0.5 m
            ktop=k
            itopct=itopct+1
         endif
         if (abs(ehcornerbot-eh2(iycorn,jxcorn,k)).lt..00005) then ! about 0.5 m
            kbot=k
            ibotct=ibotct+1
         endif
      enddo
      if (itopct.ne.1.or.ibotct.ne.1) then
         write(iup,*)'Problem with finding vertical levels'
         write(iup,*)'in routine netasc.'
      endif
c
c   Fill the model array with the values from the gridswarm.
c      
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            asc(i,j,k)=rmsg
         enddo
         enddo
      enddo
c
      itr=0
      ibeg=iycorn
      iend=ibeg+(ny-1)*ids
      jbeg=jxcorn
      jend=jbeg+(nx-1)*ids
c
      do k=ktop,kbot
         do i=ibeg,iend,ids
         do j=jbeg,jend,ids
            itr=itr+1
            asc(i,j,k)=asctr(itr)
         enddo
         enddo
      enddo
c
      return
c
 170  write(iup,*)'In netasc: The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  write(iup,*)'In netasc: Unexpected EOF reached when trying'
      write(iup,*)'to read model data header.  Stopping.'
      stop
c
 190  write(iup,*)'In netasc: Error in reading the model data array.'//
     &   ' Stopping.'
      stop
c
 195  write(iup,*)'In netasc: Unexpected EOF encountered while'
      write(iup,*)'reading the model data array. Stopping.'
      stop
c
      end
