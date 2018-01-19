c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine htrajdraw(ilinw,vc3d,eh,ixwin,iywin,
     &         idash,cnohl,lnolb,cvcor,icolr,icomg,
     &         ilcll,ilchl,rtslb,rtshl,
     &         icong,ilwng,idang,lnmsg,rvwin,rtjsp,itjns,itjid,itjni,
     &         rtjar,cfeld,rtjst,rtjen,rtjti,rstrm,rtim,ctim,dttraj,
     &         ntraj,ntrajtime,vtickinc,isense,xtime,maxpl,
     &         miy,mjx,mkzh,ipl,irota)
c
      dimension vc3d(miy,mjx,mkzh),eh(miy,mjx,mkzh),
     &   ixwin(2,maxpl),
     &   iywin(2,maxpl),idash(maxpl),ilinw(maxpl),
     &   icong(maxpl),ilwng(maxpl),idang(maxpl),
     &   icolr(maxpl),icomg(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),rvwin(2,maxpl),rtjsp(3,50,maxpl),
     &   itjns(maxpl),itjid(30,maxpl),itjni(maxpl),
     &   rtjar(2,maxpl),rtjst(maxpl),rtjen(maxpl),
     &   rtjti(maxpl),rstrm(2,maxpl),irota(maxpl)
      logical lnolb(maxpl),lnmsg(maxpl)
      character cvcor(maxpl)*1,cfeld(3,maxpl)*10,cnohl(maxpl)*1
c
      parameter (maxtraj=7000,maxtrajtime=200)
c
      dimension stortr(maxtrajtime,maxtraj,3),stormpx(maxtrajtime),
     &   stormpy(maxtrajtime),tmx(maxtrajtime),tmy(maxtrajtime),
     &   tmv(maxtrajtime),itrchoos(maxtraj),tlx(maxtrajtime),
     &   tly(maxtrajtime),trx(maxtrajtime),try(maxtrajtime),
     &   ichanged(maxtraj)
      parameter (nra=4,nim=nra,nst=nra+nim,nnd=nra+2*nim)
      dimension xra(nra),yra(nra),dst(nst),ind(nnd)
      character trlab*4, leglablow*16,leglabhigh*16
c
      dimension circlex(37),circley(37)
c
      include 'comconst'
c
c     write(0,*) 'begin htraj'
      if (ntraj.gt.maxtraj) then
         write(iup,*)'In htrajdraw, ntraj,maxtraj=',ntraj,maxtraj
         write(iup,*)'Increase maxtraj parameter to exceed ntraj.'
         stop
      elseif (ntrajtime.gt.maxtrajtime) then
         write(iup,*)'In htrajdraw, ntrajtime,maxtrajtime=',
     &      ntrajtime,maxtrajtime
         write(iup,*)'Increase maxtrajtime parameter to exceed',
     &      ' ntrajtime.'
         stop
      endif
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
      xintervs=ixwin(2,ipl)-ixwin(1,ipl)
      yintervs=iywin(2,ipl)-iywin(1,ipl)
      if (irota(ipl).ne.90.and.irota(ipl).ne.-90) then
         aspect=yintervs/xintervs
      else
         aspect=xintervs/yintervs
      endif
      faspect=(ftmax-fbmin)/(frmax-flmin)
      if (aspect.lt.faspect) then
         fl=flmin
         fr=frmax
         fextra=.5*((ftmax-fbmin)-aspect*(frmax-flmin))
         fb=fbmin+fextra
         ft=ftmax-fextra
      else
         fb=fbmin
         ft=ftmax
         fextra=.5*((frmax-flmin)-1./aspect*(ftmax-fbmin))
         fl=flmin+fextra
         fr=frmax-fextra
      endif
c
c   Set user coordinate limits for unrotated view.
c
      ul=float(ixwin(1,ipl))
      ur=float(ixwin(2,ipl))
      ub=float(iywin(1,ipl))
      ut=float(iywin(2,ipl))
c
c   Transpose ul,ur,ub,ut if view is rotated +/- 90 degrees.
c   Trajectory positions will be adjusted for rotation further down.
c
      if (irota(ipl).eq.90.or.irota(ipl).eq.-90) then
         usv=ul
         ul=ub
         ub=usv
         usv=ur
         ur=ut
         ut=usv
      endif
c
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
      ftou=(ur-ul)/(fr-fl)
c
c   Convert plspecs to usable values for line width and dash pattern.
c
      lwidth=ilinw(ipl)*1000
      call getdash(idash(ipl),ndot)
      if (cfeld(1,ipl).eq.'circle    ') then
         lwidthng=ilwng(ipl)*1000
         call getdash(idang(ipl),ndotng)
      endif
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
c     write(0,*) 'ntraj = ',ntraj,' itm1 = ',itm1
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
c
         xposnow=stormpx( nint((xtime-trbegtime)/dttraj*3600.)+1 )
         yposnow=stormpy( nint((xtime-trbegtime)/dttraj*3600.)+1 )
         do itm=1,ntrajtime
            stormpx(itm)=stormpx(itm)-xposnow
            stormpy(itm)=stormpy(itm)-yposnow
         enddo
      else
         do itm=1,ntrajtime
            xtimetraj=trbegtime+(itm-1)*dttraj/3600.
            stormpx(itm)=rstrm(2,ipl)*(xtimetraj-xtime)*3600./ds
            stormpy(itm)=rstrm(1,ipl)*(xtimetraj-xtime)*3600./ds
         enddo
      endif
c
c   Convert headhr, tailhr arinthr to timesteps.
c
      if (rtjst(ipl).eq.rmsg) then
         plbegtime=trbegtime
      else
         plbegtime=max(rtjst(ipl),trbegtime)
      endif
      if (rtjen(ipl).eq.rmsg) then
         plendtime=trendtime
      else
         plendtime=min(rtjen(ipl),trendtime)
      endif
      iplbt=nint((plbegtime-trbegtime)*3600./dttraj) + 1
      iplet=nint((plendtime-trbegtime)*3600./dttraj) + 1
      iait=nint(rtjti(ipl)*3600./dttraj)
      izero=nint((-trbegtime)*3600./dttraj) + 1
      inow=nint((xtime-trbegtime)*3600./dttraj) + 1
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
         if (itjni(ipl).ne.3) then
            write(iup,*)'In htrajdraw, you specified feld=gridswarm,'
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
c   sure the quantity and units are consistent with vtickinc.
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
c            write(iup,*)'ich,itr,itm=',ich,itr,itm
c            write(iup,*)'stortr1,2,3=',stortr(itm,itr,1),
c    &         stortr(itm,itr,2),stortr(itm,itr,3)
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
      if (rvwin(1,ipl).eq.rmsg) then
         if (isense.eq.1) then
            vclow=vcmin
         else
            vclow=vcmax
         endif
      else
         vclow=rvwin(1,ipl)
      endif
      if (rvwin(2,ipl).eq.rmsg) then
         if (isense.eq.1) then
            vchigh=vcmax
         else
            vchigh=vcmin
         endif
      else
         vchigh=rvwin(2,ipl)
      endif
      vclow=nint(vclow/vtickinc)*vtickinc
      vchigh=nint(vchigh/vtickinc)*vtickinc
      if (vclow.eq.vchigh) vchigh=vchigh+isense*vtickinc
c
c   Plot the storm center position with an "L".
c
      if (itjns(ipl).gt.0.and.cnohl(ipl).ne.' ') then
         call gsplci(ilchl(ipl))
         call gstxci(ilchl(ipl))
         call plchhq(xposnow,yposnow,'L',rtshl(ipl),0.,0.)
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
c      Calculate the storm-relative trajectory positions, and rotate
c      if called for
c
         do itm=ifirst,ilast
            tmx(itm)=stortr(itm,itr,2)-stormpx(itm)
            tmy(itm)=stortr(itm,itr,1)-stormpy(itm)
            tmv(itm)=stortr(itm,itr,3)
            if (irota(ipl).eq.90) then
               tempo=tmx(itm)
               tmx(itm)=iywin(2,ipl)-(tmy(itm)-iywin(1,ipl))
               tmy(itm)=tempo
            elseif (irota(ipl).eq.-90) then
               tempo=tmx(itm)
               tmx(itm)=tmy(itm)
               tmy(itm)=ixwin(2,ipl)-(tempo-ixwin(1,ipl))
            elseif (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
               tmx(itm)=ixwin(2,ipl)-(tmx(itm)-ixwin(1,ipl))
               tmy(itm)=iywin(2,ipl)-(tmy(itm)-iywin(1,ipl))
            endif
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
            angla=angle( tmx(min(itm,ilast-1)), tmy(min(itm,ilast-1)) ,
     +                   tmx(min(itm+1,ilast)), tmy(min(itm+1,ilast)) ,
     &                   ftou,ftou)
            anglb=angle( tmx(max(itm,ifirst+1)),tmy(max(itm,ifirst+1)),
     +                   tmx(max(itm-1,ifirst)),tmy(max(itm-1,ifirst)),
     &                   ftou,ftou)
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
            fac=(tmv(itm)-vclow)/(vchigh-vclow)
            hafar=fac*rtjar(2,ipl)+(1.-fac)*rtjar(1,ipl)
c
c         Determine arrowhead endpoints.
c
            tlx(itm)=tmx(itm)+hafar*ftou*cos(anglleft)
            tly(itm)=tmy(itm)+hafar*ftou*sin(anglleft)
            trx(itm)=tmx(itm)+hafar*ftou*cos(anglrit)
            try(itm)=tmy(itm)+hafar*ftou*sin(anglrit)
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
            call curved(tlx(ifirst),tly(ifirst),ilast-ifirst+1)
            call curved(trx(ifirst),try(ifirst),ilast-ifirst+1)
         elseif (cfeld(1,ipl).eq.'arrow     ') then
c           call curved(tmx(ifirst),tmy(ifirst),ilast-ifirst+1)
c            write(0,*) 'ifirst = ',ifirst,' ilast = ',ilast
c            write(0,*) 'loop end = ',ilast-ifirst+1
c            write(0,*) 'xtime = ',xtime,' inow = ',inow
c            write(0,*) 'trbegtime = ',trbegtime,' dttraj = ',dttraj
c            write(0,*) 'icolr(ipl) = ',icolr(ipl)
             call frstd(tmx(ifirst),tmy(ifirst))
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
               CALL VECTD (tmx(n),tmy(n))
             enddo
             CALL LASTD

         endif
c
c      Plot the arrow heads.
c
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
               call line(tmx(itm),tmy(itm),tlx(itm),tly(itm))
               call line(tmx(itm),tmy(itm),trx(itm),try(itm))
               call line(tlx(itm),tly(itm),trx(itm),try(itm))
            elseif (itm.eq.ifirst.or.itm.eq.ilast.or.
     &              mod(itm-izero,iait).eq.0)then
               call line(tmx(itm),tmy(itm),tlx(itm),tly(itm))
               call line(tmx(itm),tmy(itm),trx(itm),try(itm))
            endif
         enddo
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
               anglh=angle( tmx(ilast-1), tmy(ilast-1) ,
     +                      tmx(ilast)  , tmy(ilast),
     &                      ftou,ftou)
               anglt=angle( tmx(ifirst+1)     , tmy(ifirst+1)      ,
     +                      tmx(ifirst)       , tmy(ifirst),
     &                      ftou,ftou)
               hlabx=tmx(ilast)+ftou*(.5+rleno2)*rtslb(ipl)*
     &            cos(anglh)
               hlaby=tmy(ilast)+ftou*(.5+rleno2)*rtslb(ipl)*
     &            sin(anglh)
               tlabx=tmx(ifirst)+ftou*(tailadd+rleno2)*rtslb(ipl)*
     &            cos(anglt)
               tlaby=tmy(ifirst)+ftou*(tailadd+rleno2)*rtslb(ipl)*
     &            sin(anglt)
               call plchhq(hlabx,hlaby,trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
               call plchhq(tlabx,tlaby,trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
            else
               call plchhq(tmx(ilast),tmy(ilast),trlab(ilabstart:),
     &            rtslb(ipl),0.,0.)
            endif
         endif
 189     continue
      enddo
c        
c   Draw scale legend.  First figure out labels for legend.
c
      if (lnmsg(ipl)) goto 832
      leglablow=' '
      leglabhigh=' '
      if (cvcor(ipl).eq.'s') then   ! model vertical level index
         write(leglablow,'(f6.3)') vclow
         write(leglabhigh,'(f6.3)') vchigh
      elseif (cvcor(ipl).eq.'z'.or.cvcor(ipl).eq.'f') then    ! height (geop.)
         write(leglablow,'(f5.2,a3)') vclow,' km'
         write(leglabhigh,'(f5.2,a3)') vchigh,' km'
      elseif (cvcor(ipl).eq.'p'.or.
     &        cvcor(ipl).eq.'l'.or.
     &        cvcor(ipl).eq.'x') then   ! pressure
         write(leglablow,'(i4,a4)') nint(vclow),' hPa'
         write(leglabhigh,'(i4,a4)') nint(vchigh),' hPa'
      elseif (cvcor(ipl).eq.'t') then   ! theta
         write(leglablow,'(i3,a3)') nint(vclow),' K '
         write(leglabhigh,'(i3,a3)') nint(vchigh),' K '
      elseif (cvcor(ipl).eq.'m') then   ! temperature
         write(leglablow,'(i3,a5)') nint(vclow),' dg C'
         write(leglabhigh,'(i3,a5)') nint(vchigh),' dg C'
      elseif (cvcor(ipl).eq.'e') then   ! theta_e
         write(leglablow,'(i3,a3)') nint(vclow),' K '
         write(leglabhigh,'(i3,a3)') nint(vchigh),' K '
      elseif (cvcor(ipl).eq.'q') then   ! PV
         write(leglablow,'(f4.1,a4)') vclow,' PVU'
         write(leglabhigh,'(f4.1,a4)') vchigh,' PVU'
      endif
      ifchlow=13
      ifchhigh=13
      ilchlow=0
      ilchhigh=0
      do i=1,12
         if (leglablow(i:i).ne.' ') ifchlow=min(ifchlow,i)
         if (leglabhigh(i:i).ne.' ') ifchhigh=min(ifchhigh,i)
         j=13-i
         if (leglablow(j:j).ne.' ') ilchlow=max(ilchlow,j)
         if (leglabhigh(j:j).ne.' ') ilchhigh=max(ilchhigh,j)
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call pcseti('TE',1)
      call plchhq(.5,.5,leglablow(ifchlow:ilchlow),.008,360.,-1.)
      call pcgetr('DR',flenleglablow)
      call plchhq(.5,.5,leglabhigh(ifchhigh:ilchhigh),.008,360.,-1.)
      call pcgetr('DR',flenleglabhigh)
      call pcseti('TE',0)
c
c   Now draw legend.  Points pd-pa-pg are the upper arrowhead,
c   points pe-pb-ph are the middle arrowhead, and points pf-pc-pi
c   are the lower arrowhead.  First, set line width, color,
c   and dash pattern.
c
      call setusv('LW',lwidth)
      call gsplci(icomg(ipl))
      call gstxci(icomg(ipl))
      call dashdb(ndot)
c
c   Fill a box with the default background color and draw a perimeter
c   with the trajectory colr specified by icolr.
c
      xright=fr-.015
      xleft=xright-max(2.*rtjar(1,ipl)+.025+flenleglablow,
     &                 2.*rtjar(2,ipl)+.025+flenleglabhigh)-.03
      ybot=fb+.015
      heightbox=.08
      ytop=ybot+heightbox
      xra(1)=xleft
      yra(1)=ybot
      xra(2)=xleft
      yra(2)=ytop
      xra(3)=xright
      yra(3)=ytop
      xra(4)=xright
      yra(4)=ybot
      call gqfais(ier,ifais)
      call gsfais(1)
      call sfsgfa(xra,yra,nra,dst,nst,ind,nnd,0)
      call gsfais(0)
      call sfsgfa(xra,yra,nra,dst,nst,ind,nnd,icomg(ipl))
      call gsfais(ifais)
c
c   Figure out the points
c
      pax=xleft+.015+rtjar(2,ipl)
      pay=ybot+.9*heightbox
      pdx=pax-rtjar(2,ipl)*cos(arrowtilt)
      pdy=pay-rtjar(2,ipl)*sin(arrowtilt)
      pgx=pax+rtjar(2,ipl)*cos(arrowtilt)
      pgy=pdy
      pfy=ybot+.15*heightbox
      piy=pfy
      pcy=pfy+rtjar(1,ipl)*sin(arrowtilt)
      pcx=pax
      pfx=pcx-rtjar(1,ipl)*cos(arrowtilt)
      pix=pcx+rtjar(1,ipl)*cos(arrowtilt)
      pex=(pdx+pfx)/2.
      pey=(pdy+pfy)/2.
      pbx=pax
      pby=(pay+pcy)/2.
      phx=(pgx+pix)/2.
      phy=pey
      call line(pdx,pdy,pax,pay)
      call line(pax,pay,pgx,pgy)
      call line(pex,pey,pbx,pby)
      call line(pbx,pby,phx,phy)
      call line(pfx,pfy,pcx,pcy)
      call line(pcx,pcy,pix,piy)
      if (cfeld(1,ipl).eq.'ribbon    ') then
         call lined(pdx,pdy,pfx,pfy)
         call lined(pgx,pgy,pix,piy)
      elseif (cfeld(1,ipl).eq.'arrow     ') then
         call lined(pax,pay,pcx,pcy)
      endif
c
c   Write labels for legend
c
      call setusv('LW',1000)  ! default (thinnest) line width
      call dashdb(65535)  ! solid
      call plchhq(pgx+.025,pgy,leglabhigh(ifchhigh:ilchhigh),
     &   .008,0.,-1.)
      call plchhq(pix+.025,piy,leglablow(ifchlow:ilchlow),
     &      .008,0.,-1.)
 832  continue
c
c   Trajectory circles showing net ascent.
c
      elseif (cfeld(1,ipl).eq.'circle    ') then
c
c   First determine avcnetmax
c
      avcnetmax=0.0
      do ich=1,ntrchoos
         itr=itrchoos(ich)
         do itm=iplbt,iplet
            if (stortr(itm,itr,1).ne.rmsg) then
               ifirst=itm
               goto 223
            endif
         enddo
         ifirst=iplet+1
 223     continue
         do itm=iplet,iplbt,-1
            if (stortr(itm,itr,1).ne.rmsg) then
               ilast=itm
               goto 225
            endif
         enddo
         ilast=iplbt-1
 225     continue
         if (ifirst.lt.ilast) then
            vcnet=stortr(ilast,itr,3)-stortr(ifirst,itr,3)
            avcnetmax=max(avcnetmax,abs(vcnet))
         endif
      enddo
      if (rvwin(1,ipl).ne.rmsg) then
         vcnetref=abs(rvwin(1,ipl))
      else
         vcnetref=avcnetmax
      endif
c
c   Now make circles
c
      do ich=1,ntrchoos
c
         itr=itrchoos(ich)
         do itm=iplbt,iplet
            if (stortr(itm,itr,1).ne.rmsg) then
               ifirst=itm
               goto 227
            endif
         enddo
         ifirst=iplet+1
 227     continue
         do itm=iplet,iplbt,-1
            if (stortr(itm,itr,1).ne.rmsg) then
               ilast=itm
               goto 229
            endif
         enddo
         ilast=iplbt-1
 229     continue
         if (ifirst.ge.ilast) goto 289
c
c      Calculate the storm-relative trajectory positions, and rotate
c      if called for
c
         do itm=ifirst,ilast
            tmx(itm)=stortr(itm,itr,2)-stormpx(itm)
            tmy(itm)=stortr(itm,itr,1)-stormpy(itm)
            tmv(itm)=stortr(itm,itr,3)
            if (irota(ipl).eq.90) then
               tempo=tmx(itm)
               tmx(itm)=iywin(2,ipl)-(tmy(itm)-iywin(1,ipl))
               tmy(itm)=tempo
            elseif (irota(ipl).eq.-90) then
               tempo=tmx(itm)
               tmx(itm)=tmy(itm)
               tmy(itm)=ixwin(2,ipl)-(tempo-ixwin(1,ipl))
            elseif (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
               tmx(itm)=ixwin(2,ipl)-(tmx(itm)-ixwin(1,ipl))
               tmy(itm)=iywin(2,ipl)-(tmy(itm)-iywin(1,ipl))
            endif
         enddo
c
c      Create array for circle radius
c
         vcnet=tmv(ilast)-tmv(ifirst)
         radius=vcnet/vcnetref*rtjar(2,ipl)
         if (radius.eq.0.0) radius=.01*rtjar(2,ipl)
c
c      Set line width, color, and dash pattern.
c
         if (radius.ge.0.0) then
            call setusv('LW',lwidth)
            call gsplci(icolr(ipl))
            call gstxci(icolr(ipl))
            call dashdb(ndot)
         elseif (radius.lt.0.0) then
            call setusv('LW',lwidthng)
            call gsplci(icong(ipl))
            call gstxci(icong(ipl))
            call dashdb(ndotng)
         endif
c
c      Draw the circle.
c
         if (tmx(inow).ne.rmsg) then
            do iperpt=0,36
               circlex(iperpt+1)=tmx(inow)+
     &            ftou*radius*cos(iperpt*.174533)
               circley(iperpt+1)=tmy(inow)+
     &            ftou*radius*sin(iperpt*.174533)
            enddo
            call curved(circlex,circley,37)
         endif
c
 289     continue
      enddo
c        
c   Draw scale legend.  First figure out labels for legend.
c
      if (lnmsg(ipl)) goto 932
      leglabhigh=' '
      if (cvcor(ipl).eq.'s') then   ! model vertical level index
         write(leglabhigh,'(f6.3)') avcnetmax
      elseif (cvcor(ipl).eq.'z'.or.cvcor(ipl).eq.'f') then    ! height (geop.)
         write(leglabhigh,'(f5.2,a3)') avcnetmax,' km'
      elseif (cvcor(ipl).eq.'p'.or.
     &        cvcor(ipl).eq.'l'.or.
     &        cvcor(ipl).eq.'x') then   ! pressure
         write(leglabhigh,'(i4,a4)') nint(avcnetmax),' hPa'
      elseif (cvcor(ipl).eq.'t') then   ! theta
         write(leglabhigh,'(i3,a3)') nint(avcnetmax),' K '
      elseif (cvcor(ipl).eq.'m') then   ! temperature
         write(leglabhigh,'(i3,a5)') nint(avcnetmax),' dg C'
      elseif (cvcor(ipl).eq.'e') then   ! theta_e
         write(leglabhigh,'(i3,a3)') nint(avcnetmax),' K '
      elseif (cvcor(ipl).eq.'q') then   ! PV
         write(leglabhigh,'(f4.1,a4)') avcnetmax,' PVU'
      endif
      ifchhigh=13
      ilchhigh=0
      do i=1,12
         if (leglabhigh(i:i).ne.' ') ifchhigh=min(ifchhigh,i)
         j=13-i
         if (leglabhigh(j:j).ne.' ') ilchhigh=max(ilchhigh,j)
      enddo
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call pcseti('TE',1)
      call plchhq(.5,.5,leglabhigh(ifchhigh:ilchhigh),.008,360.,-1.)
      call pcgetr('DR',flenleglabhigh)
      call pcseti('TE',0)
c
c   Now draw legend.  First, set line width, color,
c   and dash pattern.
c
      call setusv('LW',lwidth)
      call gsplci(icomg(ipl))
      call gstxci(icomg(ipl))
      call dashdb(ndot)
c
c   Fill a box with the default background color and draw a perimeter
c   with the trajectory colr specified by icolr.
c
      xright=fr-.015
      xleft=xright-(2.*rtjar(2,ipl)+.025+flenleglabhigh)-.03
      ybot=fb+.015
      heightbox=2.*rtjar(2,ipl)+.03
      ytop=ybot+heightbox
      xra(1)=xleft
      yra(1)=ybot
      xra(2)=xleft
      yra(2)=ytop
      xra(3)=xright
      yra(3)=ytop
      xra(4)=xright
      yra(4)=ybot
      call gqfais(ier,ifais)
      call gsfais(1)
      call sfsgfa(xra,yra,nra,dst,nst,ind,nnd,0)
      call gsfais(0)
      call sfsgfa(xra,yra,nra,dst,nst,ind,nnd,icolr(ipl))
      call gsfais(ifais)
c
c   Draw the legend circle
c
      xcenter=xleft+.015+rtjar(2,ipl)
      ycenter=.5*(ybot+ytop)
      radius=avcnetmax/vcnetref*rtjar(2,ipl)
      do iperpt=0,36
         circlex(iperpt+1)=xcenter+
     &      radius*cos(iperpt*.174533)
         circley(iperpt+1)=ycenter+
     &      radius*sin(iperpt*.174533)
      enddo
      call curved(circlex,circley,37)
c
c   Write label for legend
c
      call setusv('LW',1000)  ! default (thinnest) line width
      call dashdb(65535)  ! solid
      call plchhq(xcenter+radius+.025,ycenter,
     &   leglabhigh(ifchhigh:ilchhigh),.008,0.,-1.)
 932  continue
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
      do itm=iplbt,iplet,iait
c
c      inow=nint((xtime-trbegtime)*3600./dttraj) + 1
      if (cfeld(1,ipl).eq.'swarm     ') then
         do ich=1,ntrchoos-1
            itr=itrchoos(ich)
            itrn=itrchoos(ich+1)
            if ( stortr(itm,itr,1)   .ne. rmsg .and.
     &         stortr(itm,itrn,1) .ne. rmsg ) then
               x1 = stortr(itm,itr,2)-stormpx(itm)
               y1 = stortr(itm,itr,1)-stormpy(itm)
               x2 = stortr(itm,itrn,2)-stormpx(itm)
               y2 = stortr(itm,itrn,1)-stormpy(itm)
               call lined(x1,y1,x2,y2)
            endif
         enddo
      elseif (cfeld(1,ipl).eq.'gridswarm ') then
         do j=1,ngslow
         do i=1,ngrapid
            itr=itrchoos((j-1)*ngrapid+i)
            if (j.lt.ngslow) then
               itrn=itrchoos(j*ngrapid+i)
               if ( stortr(itm,itr,1)   .ne. rmsg .and.
     &            stortr(itm,itrn,1) .ne. rmsg ) then
                  x1 = stortr(itm,itr,2)-stormpx(itm)
                  y1 = stortr(itm,itr,1)-stormpy(itm)
                  x2 = stortr(itm,itrn,2)-stormpx(itm)
                  y2 = stortr(itm,itrn,1)-stormpy(itm)
                  call lined(x1,y1,x2,y2)
               endif
            endif
            if (i.lt.ngrapid) then
               itrn=itrchoos((j-1)*ngrapid+i+1)
               if ( stortr(itm,itr,1)   .ne. rmsg .and.
     &               stortr(itm,itrn,1) .ne. rmsg ) then
                  x1 = stortr(itm,itr,2)-stormpx(itm)
                  y1 = stortr(itm,itr,1)-stormpy(itm)
                  x2 = stortr(itm,itrn,2)-stormpx(itm)
                  y2 = stortr(itm,itrn,1)-stormpy(itm)
                  call lined(x1,y1,x2,y2)
               endif
            endif
         enddo
         enddo
      endif
c
      enddo
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
