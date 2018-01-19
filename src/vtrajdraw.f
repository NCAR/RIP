c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vtrajdraw(ilinw,vc3d,eh,idash,
     &         cnohl,lnolb,cvcor,icolr,ilcll,ilchl,rtslb,rtshl,
     &         rtjsp,itjns,itjid,itjni,rtjar,cfeld,rtjst,
     &         rtjen,rtjti,rstrm,igdir,set1,set2,xdist,ydist,rcrag,
     &         xseclen,rtim,ctim,dttraj,
     &         ntraj,ntrajtime,xtime,maxpl,miy,mjx,mkzh,ipl)
c
      dimension vc3d(miy,mjx,mkzh),eh(miy,mjx,mkzh),
     &   idash(maxpl),rcrag(2,maxpl),
     &   ilinw(maxpl),icolr(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),rtjsp(3,50,maxpl),
     &   itjns(maxpl),itjid(30,maxpl),itjni(maxpl),
     &   rtjar(2,maxpl),rtjst(maxpl),rtjen(maxpl),
     &   rtjti(maxpl),rstrm(2,maxpl),igdir(maxpl)
      logical lnolb(maxpl)
      character cvcor(maxpl)*1,cfeld(3,maxpl)*10,cnohl(maxpl)*1
c
      parameter (maxtraj=7000,maxtrajtime=2800)
c
      dimension stortr(maxtrajtime,maxtraj,3),stormpx(maxtrajtime),
     &   stormpy(maxtrajtime),tmd(maxtrajtime),
     &   tmv(maxtrajtime),itrchoos(maxtraj),tld(maxtrajtime),
     &   tlv(maxtrajtime),trd(maxtrajtime),trv(maxtrajtime),
     &   ichanged(maxtraj)
c      dimension ftalt(maxtrajtime)
      character trlab*4
c
      include 'comconst'
c
c      xtimenow=rtjen(ipl)
      xtimenow=xtime
c
      if (ntraj.gt.maxtraj) then
         write(iup,*)'In vtrajdraw, ntraj,maxtraj=',ntraj,maxtraj
         write(iup,*)'Increase maxtraj parameter to exceed ntraj.'
         stop
      elseif (ntrajtime.gt.maxtrajtime) then
         write(iup,*)'In vtrajdraw, ntrajtime,maxtrajtime=',
     &      ntrajtime,maxtrajtime
         write(iup,*)'Increase maxtrajtime parameter to exceed',
     &      ' ntrajtime.'
         stop
      endif
c
      fb=fbmin
      ft=ftmax
      fl=flmin
      fr=frmax
c
c   Make set call for conrec.
c
      call set(fl,fr,fb,ft,0.,xseclen,set1,set2,1)
      ftoud=xseclen/(fr-fl)
      ftouv=(set2-set1)/(ft-fb)
      xsecleni=1./xseclen
      if (igdir(ipl).ne.362) then
         tanb=tan(float(igdir(ipl))*rpd)
      else
         tanb=0.
      endif
c
c   Convert plspecs to usable values for line width and dash pattern.
c
      lwidth=ilinw(ipl)*1000
      call getdash(idash(ipl),ndot)
c
c   Determine forward/backward.
c
      if (rtim.lt.ctim) then
         trendtime=ctim
         trbegtime=rtim
         itm1=1
         itm2=ntrajtime
         itmi=1
      else
         trendtime=rtim
         trbegtime=ctim
         itm1=ntrajtime
         itm2=1
         itmi=-1
      endif
c
c   Get the trajectory position info.
c   Note, at this point, stortr(n,m,1) is y-coordinate in coarse domain
c   dot-point grid, stortr(n,m,2) is x-coordinate, and stortr(n,m,3) is
c   exponential height [i.e., exp(-z/H), z in meters].
c      
      do itm=itm1,itm2,itmi
         read(iutrajin) (stortr(itm,itr,1),itr=1,ntraj),
     &       (stortr(itm,itr,2),itr=1,ntraj),
     &       (stortr(itm,itr,3),itr=1,ntraj)
      enddo
c   
c   Make storm position arrays, and change storm position to be
c       relative to the position at the current time.  Also,
c       convert storm position from coarse domain grid values
c       (as specified in tjsp) to current domain grid values.
c   
      if (itjns(ipl).gt.0) then
         do itm=1,ntrajtime
            xtimetraj=trbegtime+(itm-1)*dttraj/3600.
            do j=1,itjns(ipl)-1
               if (xtimetraj.ge.rtjsp(1,j,ipl).and.
     &             xtimetraj.le.rtjsp(1,j+1,ipl)) then
                  fac1=(xtimetraj-rtjsp(1,j,ipl))/
     &                 (rtjsp(1,j+1,ipl) -rtjsp(1,j,ipl))
                  fac2=(rtjsp(1,j+1,ipl)-xtimetraj)/
     &                 (rtjsp(1,j+1,ipl) -rtjsp(1,j,ipl))
                  stormpx(itm)=fac2*rtjsp(3,j,ipl)+
     &                       fac1*rtjsp(3,j+1,ipl)
                  stormpx(itm)=1.+(stormpx(itm)-xjcorn)*refrat
                  stormpy(itm)=fac2*rtjsp(2,j,ipl)+
     &                       fac1*rtjsp(2,j+1,ipl)
                  stormpy(itm)=1.+(stormpy(itm)-yicorn)*refrat
                  goto 117
               endif
            enddo
            write(iup,*)'trajectory time not within storm position'
            write(iup,*)'time range.'
            stop
 117        continue
         enddo
         xposnow=stormpx( nint((xtimenow-trbegtime)/dttraj*3600.)+1 )
         yposnow=stormpy( nint((xtimenow-trbegtime)/dttraj*3600.)+1 )
         do itm=1,ntrajtime
            stormpx(itm)=stormpx(itm)-xposnow
            stormpy(itm)=stormpy(itm)-yposnow
         enddo
      else
         do itm=1,ntrajtime
            xtimetraj=trbegtime+(itm-1)*dttraj/3600.
            stormpx(itm)=rstrm(2,ipl)*(xtimetraj-xtimenow)*3600./ds
            stormpy(itm)=rstrm(1,ipl)*(xtimetraj-xtimenow)*3600./ds
         enddo
      endif
c
c   Figure out important time steps
c
      if (rtjst(ipl).eq.rmsg) then
         plbegtime=trbegtime
      else
         plbegtime=rtjst(ipl)
      endif
      if (rtjen(ipl).eq.rmsg) then
         plendtime=trendtime
      else
         plendtime=rtjen(ipl)
      endif
      iplbt=nint((plbegtime-trbegtime)*3600./dttraj) + 1
      iplet=nint((plendtime-trbegtime)*3600./dttraj) + 1
      if (rtjti(ipl).ge.0.) iait=nint(rtjti(ipl)*3600./dttraj)
      izero=nint((-trbegtime)*3600./dttraj) + 1
      inow=nint((xtimenow-trbegtime)*3600./dttraj) + 1
c
c   Determine which trajectories to use.
c
      if (cfeld(1,ipl).ne.'gridswarm ') then
         if (itjni(ipl).eq.0) then
            do ich=1,ntraj
               itrchoos(ich)=ich
            enddo
            ntrchoos=ntraj
         else
            ii=0
            itrajid=0
  100       ii=ii+1
            if (ii.gt.itjni(ipl)) goto 120
            if (itjid(ii,ipl).ge.0) then
               itrajid=itrajid+1
               itrchoos(itrajid)=itjid(ii,ipl)
            else
               ii=ii+1
               if (itjid(ii,ipl).gt.0) then
                  istart=itjid(ii-2,ipl)
                  iend=-itjid(ii-1,ipl)
                  iinc=itjid(ii,ipl)
                  idist=iend-istart
                  isign=idist/abs(idist)
                  nlseries=abs(idist)/iinc + 1
                  do 110 i=2,nlseries
                     itrajid=itrajid+1
                     itrchoos(itrajid)=itrchoos(itrajid-1)+isign*iinc
  110             continue
               else
                  write(iup,*)'Error in traj. id series.'
                  stop
               endif
            endif
            goto 100
  120       ntrchoos=itrajid
         endif
      else
         if (itjni(ipl).eq.0) then
            write(iup,*)'In vtrajdraw, you specified feld=gridswarm,'
            write(iup,*)'but didn''t give any values for tjid.'
            write(iup,*)'You must supply 3 values in tjid.'
            stop
         endif
         icorner=itjid(1,ipl)
         ngrapid=itjid(2,ipl)
         ngslow=itjid(3,ipl)
         ntrchoos=ngrapid*ngslow
         do ich=1,ntrchoos
            itrchoos(ich)=icorner+ich-1
         enddo
      endif
c
c   Replace exp. height values with specified vertical coordinate, making
c   sure the quantity and units are consistent with set1 and set2.
c   Also, convert x and y values (which are in coarse dom. grid points)
c   to grid points in the current domain.
c
      do itr=1,ntraj
         ichanged(itr)=0
      enddo
      vcmax=-9e9
      vcmin=9e9
      do ich=1,ntrchoos
         itr=itrchoos(ich)
         if (ichanged(itr).eq.1) goto 49
         do itm=iplbt,iplet
            if (stortr(itm,itr,1).ne.rmsg) then
               stortr(itm,itr,3)=finterp(eh,vc3d,1,miy,mjx,mkzh,
     &            stortr(itm,itr,1),stortr(itm,itr,2),
     &            stortr(itm,itr,3),refrat,yicorn,xjcorn,rmsg,iup)
               if (stortr(itm,itr,3).eq.rmsg) then
                  stortr(itm,itr,1)=rmsg
                  stortr(itm,itr,2)=rmsg
               else
                  stortr(itm,itr,1)=
     &               1.+refrat*(stortr(itm,itr,1)-yicorn)
                  stortr(itm,itr,2)=
     &               1.+refrat*(stortr(itm,itr,2)-xjcorn)
                  if (cvcor(ipl).eq.'z'.or.cvcor(ipl).eq.'f') then ! z in km
                     stortr(itm,itr,3)=
     &                  -.001*sclht*log(stortr(itm,itr,3)) 
                  elseif (cvcor(ipl).eq.'l') then ! ln of pressure
                     stortr(itm,itr,3)=log(stortr(itm,itr,3)) 
                  elseif (cvcor(ipl).eq.'x') then ! pressure**gamma
                     stortr(itm,itr,3)=(stortr(itm,itr,3))**gamma
                  endif
                  vcmax=max(vcmax,stortr(itm,itr,3))
                  vcmin=min(vcmin,stortr(itm,itr,3))
               endif
            endif
         enddo
         ichanged(itr)=1
 49      continue
      enddo
c
c   Convert rcrag (which is in coarse dom. grid points) to
c   grid points in the current domain.
c
      crax=1.+refrat*(rcrag(2,ipl)-xjcorn)
      cray=1.+refrat*(rcrag(1,ipl)-yicorn)
c
c   Plot the storm center position with an "L".
c
      if (itjns(ipl).gt.0.and.cnohl(ipl).ne.' ') then
         call gsplci(ilchl(ipl))
         call gstxci(ilchl(ipl))
         xveca=xposnow-crax
         yveca=yposnow-cray
         dposnow=(xveca*xdist+yveca*ydist-
     &      (yveca*xdist-xveca*ydist)*tanb)*xsecleni
         vposnow=set1+(.015+.5*rtshl(ipl))*ftouv
         call plchhq(dposnow,vposnow,'L',rtshl(ipl),0.,0.)
         call gsplci(1)
         call gstxci(1)
      endif
c
c   Trajectory ribbons or arrows.
c
      if (cfeld(1,ipl).eq.'ribbon    '.or.
     &    cfeld(1,ipl).eq.'arrow     ') then
c
      if (cfeld(1,ipl).eq.'ribbon    ') then
         arrowtilt=45.0*rpd
      else
         arrowtilt=60.0*rpd
      endif
c
      do ich=1,ntrchoos
c
c      Create two 1-d arrays containing the trajectory horizontal
c      and vertical coordinates.
c
         itr=itrchoos(ich)
         do itm=iplbt,iplet
            if (stortr(itm,itr,1).ne.rmsg) then
               ifirst=itm
               goto 123
            endif
         enddo
         ifirst=iplet+1
 123     continue
         do itm=iplet,iplbt,-1
            if (stortr(itm,itr,1).ne.rmsg) then
               ilast=itm
               goto 125
            endif
         enddo
         ilast=iplbt-1
 125     continue
c
         if (iplbt.gt.iplet.or.ifirst.gt.ilast.or.
     &       (iplbt.lt.iplet.and.ifirst.eq.ilast)) goto 189
c
         ishow=0
         do itm=ifirst,ilast
            tmx=stortr(itm,itr,2)-stormpx(itm)
            tmy=stortr(itm,itr,1)-stormpy(itm)
            xveca=tmx-crax
            yveca=tmy-cray
            tmd(itm)=(xveca*xdist+yveca*ydist-
     &         (yveca*xdist-xveca*ydist)*tanb)*xsecleni
            tmv(itm)=stortr(itm,itr,3)
            if (ishow.eq.0.and.tmv(itm).ge.3.) then
c               print*,'traj ht, dist=',tmv(itm),tmd(itm)*dskm
               ishow=1
            endif
c            xtimetraj=trbegtime+(itm-1)*dttraj/3600.
c            write(69,'(f10.4,f10.3,f10.1)')
c     &         xtimetraj,tmd(itm)*dskm,ftalt(itm)
         enddo
c
         if (iplbt.eq.iplet.and.ifirst.eq.ilast) goto 169
c
c      Create arrays for "side curves".
c
         do itm=ifirst,ilast
c
c         Determine angles for left and right sides of arrowhead.
c
            angla=angle( tmd(min(itm,ilast-1)), tmv(min(itm,ilast-1)) ,
     +                   tmd(min(itm+1,ilast)), tmv(min(itm+1,ilast)) ,
     &                   ftoud,ftouv)
            anglb=angle( tmd(max(itm,ifirst+1)),tmv(max(itm,ifirst+1)),
     +                   tmd(max(itm-1,ifirst)),tmv(max(itm-1,ifirst)),
     &                   ftoud,ftouv)
            anglavg=(angla+anglb)/2.
            if (angla.le.anglb) then
               anglleft=anglavg + arrowtilt
            else
               anglleft=anglavg + 3.14159 + arrowtilt
            endif
            anglrit=anglleft + 3.14159 - 2.*arrowtilt
c
c         Determine size of arrowhead
c
            hafar=rtjar(2,ipl)
c
c         Determine arrowhead endpoints.
c
            tld(itm)=tmd(itm)+hafar*ftoud*cos(anglleft)
            tlv(itm)=tmv(itm)+hafar*ftouv*sin(anglleft)
            trd(itm)=tmd(itm)+hafar*ftoud*cos(anglrit)
            trv(itm)=tmv(itm)+hafar*ftouv*sin(anglrit)
         enddo
c
c      Set line width, color, and dash pattern.
c
         call setusv('LW',lwidth)
         call gsplci(icolr(ipl))
         call gstxci(icolr(ipl))
         call dashdb(ndot)
c
c      Plot the trajectory.
c
         if (cfeld(1,ipl).eq.'ribbon    ') then
            call curved(tld(ifirst),tlv(ifirst),ilast-ifirst+1)
            call curved(trd(ifirst),trv(ifirst),ilast-ifirst+1)
         elseif (cfeld(1,ipl).eq.'arrow     ') then
c           call curved(tmd(ifirst),tmv(ifirst),ilast-ifirst+1)
c            write(0,*) 'ifirst = ',ifirst,' ilast = ',ilast
c            write(0,*) 'loop end = ',ilast-ifirst+1
c            write(0,*) 'xtime = ',xtime,' inow = ',inow
c            write(0,*) 'trbegtime = ',trbegtime,' dttraj = ',dttraj
             call frstd(tmd(ifirst),tmv(ifirst))
             do n = ifirst, ilast
               call sflush
	       if (icolr(ipl) .eq. 103) then   !if it is set wheat4, then color by time
               if (float(n) / float(ntrajtime) .lt. .2) then
                 call gsplci(37)   ! tail of the back trajectory
                 call gstxci(37)
               else if (float(n) / float(ntrajtime) .lt. .4) then
                 call gsplci(24)
                 call gstxci(24)
               else if (float(n) / float(ntrajtime) .lt. .6) then
                 call gsplci(30)
                 call gstxci(30)
               else if (float(n) / float(ntrajtime) .lt. .8) then
                 call gsplci(13)
                 call gstxci(13)
               else   ! head of the back trajectory
                 call gsplci(1)
                 call gstxci(1)
               endif
               endif
               CALL VECTD (tmd(n),tmv(n))
             enddo
             CALL LASTD
         endif
c
c      Plot the arrow heads.
c
         if (rtjti(ipl).ge.0.) then
            do itm=ifirst,ilast
	     call sflush
	     if (icolr(ipl) .eq. 103) then
             if (float(itm) / float(ntrajtime) .lt. .2) then
               call gsplci(37)   ! tail of the back trajectory
               call gstxci(37)
             else if (float(itm) / float(ntrajtime) .lt. .4) then
               call gsplci(24)
               call gstxci(24)
             else if (float(itm) / float(ntrajtime) .lt. .6) then
               call gsplci(30)
               call gstxci(30)
             else if (float(itm) / float(ntrajtime) .lt. .8) then
               call gsplci(13)
               call gstxci(13)
             else   ! head of the back trajectory
               call gsplci(1)
               call gstxci(1)
             endif
	     endif
               if (itm.eq.inow) then
                  call line(tmd(itm),tmv(itm),tld(itm),tlv(itm))
                  call line(tmd(itm),tmv(itm),trd(itm),trv(itm))
                  call line(tld(itm),tlv(itm),trd(itm),trv(itm))
               elseif (itm.eq.ifirst.or.itm.eq.ilast.or.
     &                 mod(itm-izero,iait).eq.0)then
                  call line(tmd(itm),tmv(itm),tld(itm),tlv(itm))
                  call line(tmd(itm),tmv(itm),trd(itm),trv(itm))
               endif
            enddo
         else
            call line(tmd(ilast),tmv(ilast),tld(ilast),tlv(ilast))
            call line(tmd(ilast),tmv(ilast),trd(ilast),trv(ilast))
         endif
c
c      re-set line width, color, and dash pattern for labels
c
 169     call setusv('LW',1000) ! default (thinnest) line width
         call gsplci(ilcll(ipl))
         call gstxci(ilcll(ipl))
         call dashdb(65535)  ! solid
c
c      Plot trajectory labels.
c
         if (.not.lnolb(ipl)) then
            write(trlab,'(i4)')itr
            if (itr.lt.10) then
               ilabstart=4
            elseif (itr.lt.100) then
               ilabstart=3
            elseif (itr.lt.1000) then
               ilabstart=2
            else
               ilabstart=1
            endif
            if (cfeld(1,ipl).eq.'ribbon    ') then
               tailadd=1.1
            else
               tailadd=1.6
            endif
            if (ifirst.lt.ilast) then  ! normal state of affairs
               rleno2=.5*(4.-ilabstart)
               anglh=angle( tmd(ilast-1), tmv(ilast-1) ,
     +                      tmd(ilast)  , tmv(ilast),
     &                      ftoud,ftouv)
               anglt=angle( tmd(ifirst+1), tmv(ifirst+1)      ,
     +                      tmd(ifirst)  , tmv(ifirst),
     &                      ftoud,ftouv)
               hlabd=tmd(ilast)+ftoud*(.5+rleno2)*rtslb(ipl)*
     &            cos(anglh)
               hlabv=tmv(ilast)+ftouv*(.5+rleno2)*rtslb(ipl)*
     &            sin(anglh)
               tlabd=tmd(ifirst)+ftoud*(tailadd+rleno2)*rtslb(ipl)*
     &            cos(anglt)
               tlabv=tmv(ifirst)+ftouv*(tailadd+rleno2)*rtslb(ipl)*
     &            sin(anglt)
               call plchhq(hlabd,hlabv,trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
               call plchhq(tlabd,tlabv,trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
            else ! only one time specified,
c                  so just show label, on traj. pos.
               call plchhq(tmd(ilast),tmv(ilast),trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
            endif
         endif
 189     continue
      enddo
c
c   Trajectory swarms
c
      elseif (cfeld(1,ipl).eq.'swarm     '.or.
     &        cfeld(1,ipl).eq.'gridswarm ') then
c
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
      call dashdb(ndot)
c
c   Plot loop.
c
      inow=nint((xtimenow-trbegtime)*3600./dttraj) + 1
      if (cfeld(1,ipl).eq.'swarm     ') then
         do ich=1,ntrchoos-1
            itr=itrchoos(ich)
            itrn=itrchoos(ich+1)
            if ( stortr(inow,itr,1)   .ne. rmsg .and.
     &         stortr(inow,itrn,1) .ne. rmsg ) then
               x1 = stortr(inow,itr,2)-stormpx(inow)
               y1 = stortr(inow,itr,1)-stormpy(inow)
               xveca=x1-crax
               yveca=y1-cray
               d1 = (xveca*xdist+yveca*ydist-
     &            (yveca*xdist-xveca*ydist)*tanb)*xsecleni
               v1 = stortr(inow,itr,3)
               x2 = stortr(inow,itrn,2)-stormpx(inow)
               y2 = stortr(inow,itrn,1)-stormpy(inow)
               xveca=x2-crax
               yveca=y2-cray
               d2 = (xveca*xdist+yveca*ydist-
     &            (yveca*xdist-xveca*ydist)*tanb)*xsecleni
               v2 = stortr(inow,itrn,3)
               call lined(d1,v1,d2,v2)
            endif
         enddo
      elseif (cfeld(1,ipl).eq.'gridswarm ') then
         do j=1,ngslow
         do i=1,ngrapid
            itr=itrchoos((j-1)*ngrapid+i)
            if (j.lt.ngslow) then
               itrn=itrchoos(j*ngrapid+i)
               if ( stortr(inow,itr,1)   .ne. rmsg .and.
     &            stortr(inow,itrn,1) .ne. rmsg ) then
                  x1 = stortr(inow,itr,2)-stormpx(inow)
                  y1 = stortr(inow,itr,1)-stormpy(inow)
                  xveca=x1-crax
                  yveca=y1-cray
                  d1 = (xveca*xdist+yveca*ydist-
     &               (yveca*xdist-xveca*ydist)*tanb)*xsecleni
                  v1 = stortr(inow,itr,3)
                  x2 = stortr(inow,itrn,2)-stormpx(inow)
                  y2 = stortr(inow,itrn,1)-stormpy(inow)
                  xveca=x2-crax
                  yveca=y2-cray
                  d2 = (xveca*xdist+yveca*ydist-
     &               (yveca*xdist-xveca*ydist)*tanb)*xsecleni
                  v2 = stortr(inow,itrn,3)
                  call lined(d1,v1,d2,v2)
               endif
            endif
            if (i.lt.ngrapid) then
               itrn=itrchoos((j-1)*ngrapid+i+1)
               if ( stortr(inow,itr,1)   .ne. rmsg .and.
     &            stortr(inow,itrn,1) .ne. rmsg ) then
                  x1 = stortr(inow,itr,2)-stormpx(inow)
                  y1 = stortr(inow,itr,1)-stormpy(inow)
                  xveca=x1-crax
                  yveca=y1-cray
                  d1 = (xveca*xdist+yveca*ydist-
     &               (yveca*xdist-xveca*ydist)*tanb)*xsecleni
                  v1 = stortr(inow,itr,3)
                  x2 = stortr(inow,itrn,2)-stormpx(inow)
                  y2 = stortr(inow,itrn,1)-stormpy(inow)
                  xveca=x2-crax
                  yveca=y2-cray
                  d2 = (xveca*xdist+yveca*ydist-
     &               (yveca*xdist-xveca*ydist)*tanb)*xsecleni
                  v2 = stortr(inow,itrn,3)
                  call lined(d1,v1,d2,v2)
               endif
            endif
         enddo
         enddo
      endif
c
      endif
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      call dashdb(65535)
c
      return
      end
