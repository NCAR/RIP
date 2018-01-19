c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine qgomg(prs,omg,tmk,qvp,ght,bvfsqd,bvfsqm,
     &   rhi,cor,xmap,ter,numpas,ivar,ivar2,iqvecforc,
     &   itopobc,iekmnbc,imo,rhithresh,ihrip,rhrip,chrip,
     &   vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &   nxt,ncxc,maxtavl,miy,mjx,mkzh)
c
      parameter (errmin=.0005)   ! .0005 hPa/s. or .5 microbar/s
c      parameter (errmin=.000001)   ! .000001 hPa/s. or .001 microbar/s
      parameter (itmax=100, alpha=1.8, mkp=30)
      parameter (frnew=.1, frold=1.-frnew)
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*(*)
c
      dimension ter(miy,mjx),
     &   prs(miy,mjx,mkzh),omg(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),
     &   bvfsqd(miy,mjx,mkzh),bvfsqm(miy,mjx,mkzh),rhi(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),cor(miy,mjx),xmap(miy,mjx)
      character temp*4
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24,varname*10
c
      dimension omgqg(miy,mjx,mkp),
     &   ghtqg(miy,mjx,mkp),
     &   alphaqg(miy,mjx,mkp),prsqg(mkp),
     &   divq(miy,mjx,mkp),ipttp(miy,mjx,mkp),
     &   fac(miy,mjx,mkp),rhs(miy,mjx,mkp),
     &   stbldqg(miy,mjx,mkp),stblmqg(miy,mjx,mkp),
     &   stblqg(miy,mjx,mkp),
     &   rhiqg(miy,mjx,mkp),terqg(miy,mjx),work(miy,mjx)
c
      include 'comconst'
c
      if (nproj.eq.4) then
         write(iup,*) 'Routine qgomg will not work exactly right'
         write(iup,*) 'for the NMM SRCE map projection because'
         write(iup,*) 'it is not conformal and I was too lazy'
         write(iup,*) 'to account for different map factors in the'
         write(iup,*) 'x and y directions in this routine.'
         write(iup,*) 'However, errors will be quite small.'
      endif
c
      dzpbl=1000.   ! Assume 1000-m-thick PBL.
      perpfac=.20  ! PBL-avg. comp. of actual wind perp. to geos. wind
      parlfac=.70  ! PBL-avg. comp. of actual wind parallel to geos. wind
      betarot=atan2(perpfac,parlfac)
      speedreduc=sqrt(perpfac*perpfac+parlfac*parlfac)
c
c   Try to get pressure levels, as well as omg, ght, alpha, stbld,
c   stblm, and rhi on pressure levels.
c
c   prsqg is embedded in a 2D array.  temporarily use k=1 slab of divq
c   to hold prsqg data
c
      write(temp,'(i4)')1000+numpas
c
      varname='prsqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,2,0,divq,istatprsqg)
      if (istatprsqg.eq.1) then
         j=1
         i=0
         do k=1,mkp
            i=i+1
            if (i.gt.miy-1) then
               j=j+1
               i=1
            endif
            prsqg(k)=divq(i,j,1)
         enddo
         dp=prsqg(2)-prsqg(1)
         dpi=1./dp
      endif
c
      varname='omgqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,3,0,omgqg,istatomgqg)
      varname='ghtqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,3,0,ghtqg,istatghtqg)
      varname='alphaqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,3,0,alphaqg,istatalphaqg)
      varname='stbldqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,3,0,stbldqg,istatstbldqg)
      if (imo.eq.1) then
         varname='stblmqg'//temp(2:4)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &        nxt, ncxc,varname,
     &        miy,mjx,mkp,maxtavl,3,0,stblmqg,istatstblmqg)
         varname='rhiqg'//temp(2:4)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &        nxt, ncxc,varname,
     &        miy,mjx,mkp,maxtavl,3,0,rhiqg,istatrhiqg)
      else
         istatstblmqg=1
         istatrhiqg=1
      endif
      varname='terqg'//temp(2:4)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,
     &     nxt, ncxc,varname, miy,
     &     mjx,mkp,maxtavl,2,0,terqg,istatterqg)
c
      if (istatomgqg.eq.1.and.istatghtqg.eq.1.and.
     &    istatalphaqg.eq.1.and.istatstbldqg.eq.1.and.
     &    (imo.eq.0.or.(istatstblmqg.eq.1.and.
     &    istatrhiqg.eq.1))) then
         do kp=1,mkp
         do j = 1, mjx-1
         do i = 1, miy-1
            if (omgqg(i,j,kp).ne.rmsg) then
               ipttp(i,j,kp)=0  ! point is within domain
            else
               ipttp(i,j,kp)=-1  ! point is below ground (missing data)
            endif
         enddo
         enddo
         enddo
         goto 473
      endif
c
c   If any files did not exist, calculate all of the variables.
c   First, calculate pressure levels.
c
      pbotqg=-9e9
      ptopqg=-9e9
      do j=1,mjx-1
      do i=1,miy-1
         pbotqg=max(pbotqg,prs(i,j,mkzh))
         ptopqg=max(ptopqg,prs(i,j,1))
      enddo
      enddo
c
c   Add 1 hPa to ptopqg, to guarantee that no pres. dom. points are
c   above model domain
c
      ptopqg=ptopqg+1.
c
      hscale=rgas*celkel/grav
      hscali=1./hscale
      dp=(pbotqg-ptopqg)/(mkp-1.)
      dpi=1./dp
c
c   Define pressure levels.
c
      do kp=1,mkp
         prsqg(kp)=ptopqg+(kp-1.)*dp
      enddo
c
c   Other variables
c
      do kp=1,mkp
      do j = 1, mjx-1
      do i = 1, miy-1
         if (prsqg(kp).gt.prs(i,j,mkzh)) then
            ipttp(i,j,kp)=-1  ! point is below ground (missing data)
            if (istatomgqg.eq.-1) omgqg(i,j,kp)=rmsg
            ttqg=rmsg
            if (istatghtqg.eq.-1) ghtqg(i,j,kp)=rmsg
            if (istatalphaqg.eq.-1) alphaqg(i,j,kp)=rmsg
            if (istatstbldqg.eq.-1) stbldqg(i,j,kp)=rmsg
            if (istatstblmqg.eq.-1) stblmqg(i,j,kp)=rmsg
            if (istatrhiqg.eq.-1) rhiqg(i,j,kp)=rmsg
         elseif (prsqg(kp).lt.prs(i,j,1)) then
            write(iup,*)'At point i,j,kp=',i,j,kp,
     &         '  pressure is too low.'
            write(iup,*)'prsqg,prs=',prsqg(kp),prs(i,j,1)
            stop
         else
            ipttp(i,j,kp)=0  ! point is within domain
            do ks=1,mkzh-1
               if (prsqg(kp).le.prs(i,j,ks+1).and.
     &             prsqg(kp).ge.prs(i,j,ks)) then
                  fac1=prsqg(kp)-prs(i,j,ks)
                  fac2=prs(i,j,ks+1)-prsqg(kp)
                  denom=prs(i,j,ks+1)-prs(i,j,ks)
                  if (istatomgqg.eq.-1) omgqg(i,j,kp)=
     &               (fac1*omg(i,j,ks+1)+fac2*omg(i,j,ks))/
     &               denom*.001   ! microbar/sec to hPa/sec
                  ttqg=(fac1*tmk(i,j,ks+1)+fac2*tmk(i,j,ks))/
     &               denom
                  qvpqg=(fac1*qvp(i,j,ks+1)+fac2*qvp(i,j,ks))/
     &               denom
                  tvqg=virtual(ttqg,qvpqg)
                  prat=(fac1*exp(-ght(i,j,ks+1)*hscali)+
     &                  fac2*exp(-ght(i,j,ks)*hscali)     )/denom
                  if (istatghtqg.eq.-1) ghtqg(i,j,kp)=-hscale*log(prat)
                  if (istatalphaqg.eq.-1) alphaqg(i,j,kp)=
     &               rgas*tvqg/prsqg(kp)
                  if (istatstbldqg.eq.-1) stbldqg(i,j,kp)=
     &               (alphaqg(i,j,kp)/grav)**2*(fac1*bvfsqd(i,j,ks+1)+
     &               fac2*bvfsqd(i,j,ks))/denom
                  if (istatstblmqg.eq.-1) stblmqg(i,j,kp)=
     &               (alphaqg(i,j,kp)/grav)**2*(fac1*bvfsqm(i,j,ks+1)+
     &               fac2*bvfsqm(i,j,ks))/denom
                  if (istatrhiqg.eq.-1) rhiqg(i,j,kp)=
     &               (fac1*rhi(i,j,ks+1)+fac2*rhi(i,j,ks))/denom
                  goto 30
               endif
            enddo
 30         continue
         endif
cc
cc      Set omgqg to rmsg in PBL, as a marker of which points are in
cc      the PBL.
cc
c         pmaxallowed=prs(i,j,mkzh)*(1.-dzpbl*hscali)+.5*dp
c         if (prsqg(kp).gt.pmaxallowed) omgqg(i,j,kp)=rmsg
c
c       Note: The above has been commented out because ground-level
c       will be used as the bottom boundary.  The Ekman BC will not
c       be completely accurate because, although it is calculated
c       assuming a PBL depth of dzpbl, it is applied at ground
c       level.  In setting the lower boundary to ground level, the
c       Ekman inaccuracy was accepted in exchange for the benefit of
c       allowing the Q-vector and terrain forcing to act all the way
c       down to the ground.
c
      enddo
      enddo
         if (numpas.gt.0) then
            write(iup,*)'Smoothing QG input fields.  kp,numpas=',kp,
     &         numpas
            if (istatghtqg.eq.-1) then
               call smooth(ghtqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new ghtqg field)'
            endif
            if (istatalphaqg.eq.-1) then
               call smooth(alphaqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new alphaqg field)'
            endif
            if (istatstbldqg.eq.-1) then
               call smooth(stbldqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new stbldqg field)'
            endif
            if (istatstblmqg.eq.-1) then
               call smooth(stblmqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new stblmqg field)'
            endif
            if (istatrhiqg.eq.-1) then
               call smooth(rhiqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new rhiqg field)'
            endif
            if (istatomgqg.eq.-1) then
               call smooth(omgqg(1,1,kp),work,numpas,miy,miy-1,mjx-1)
               if (kp.eq.1) write(iup,*)'   (new omgqg field)'
            endif
         endif
      enddo
c
c  If we just calculated the pressure-level fields, write them out.
c  Embed prsqg in a 2D array.  temporarily use k=1 slab of divq
c  to hold prsqg data, to write it out.
c
      if (istatprsqg.eq.-1) then
         j=1
         i=0
         do k=1,mkp
            i=i+1
            if (i.gt.miy-1) then
               j=j+1
               i=1
            endif
            divq(i,j,1)=prsqg(k)
         enddo
         vardesc='Pressure levels'
         plchun='hPa'
         varname='prsqg     '
         call writefile(divq,varname,numpas,
     &      2,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatomgqg.eq.-1) then
         vardesc='Omega on pressure levels'
         plchun='hPa/s'
         varname='omgqg     '
         call writefile(omgqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatghtqg.eq.-1) then
         vardesc='Geo. Height on pressure levels'
         plchun='m'
         varname='ghtqg     '
         call writefile(ghtqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatalphaqg.eq.-1) then
         vardesc='Alpha (spec. volume) on pressure levels'
         plchun='J/kg/hPa'
         varname='alphaqg   '
         call writefile(alphaqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatstbldqg.eq.-1) then
         vardesc='Dry stability (Holton''s sigma) on pressure levels'
         plchun='(m/s/hPa)**2'
         varname='stbldqg  '
         call writefile(stbldqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatstblmqg.eq.-1) then
         vardesc='Moist stability (Holton''s sigma) on pressure levels'
         plchun='(m/s/hPa)**2'
         varname='stblmqg  '
         call writefile(stblmqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
      if (istatrhiqg.eq.-1) then
         vardesc='RH w.r.t. ice'
         plchun='%'
         varname='rhiqg     '
         call writefile(rhiqg,varname,numpas,
     &      3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
c
 473  continue
c
      if (istatterqg.eq.-1) then
         do j=1,mjx
         do i=1,miy
            terqg(i,j)=ter(i,j)
         enddo
         enddo
         if (numpas.gt.0) then
            write(iup,*)'Smoothing terrain.'
            call smooth(terqg,work,numpas,miy,miy-1,mjx-1)
         endif
         vardesc='Smoothed terrain height'
         plchun='m'
         varname='terqg     '
         call writefile(terqg,varname,numpas,
     &      2,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
     &      cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
      endif
c
c   Identify interior points.
c
      do j = 3, mjx-3
      do i = 3, miy-3
      do k=2,mkp-1
         nmsg=0
         do jj=j-2,j+2
            idi=2-abs(jj-j)
         do ii=i-idi,i+idi
            if (ipttp(ii,jj,k).eq.-1) nmsg=1
         enddo
         enddo
         if (ipttp(i,j,k+1).eq.-1.or.ipttp(i,j,k-1).eq.-1)
     &      nmsg=1
         if (nmsg.eq.0) ipttp(i,j,k)=1   ! interior point
      enddo
      enddo
      enddo
c
c   Identify boundary points.
c
      do j = 1, mjx-1
      do i = 1, miy-1
      do k=1,mkp
         if (ipttp(i,j,k).eq.0) then
            if (k.eq.1) then
               ipttp(i,j,k)=3   ! top boundary point
            elseif (j.le.2.or.j.ge.mjx-2.or.i.le.2.or.i.ge.miy-2) then
               ipttp(i,j,k)=4   ! lateral boundary point
            else
               ipttp(i,j,k)=2   ! bottom boundary point
            endif
         endif
         rhs(i,j,k)=rmsg
         fac(i,j,k)=rmsg
      enddo
      enddo
      enddo
c
c   Make stability no less than .0005.  This value is
c   equivalent to d(theta)/dp of about 0.5K per 300 hPa.
c
      do k=1,mkp
      do j = 1, mjx-1
      do i = 1, miy-1
         if (ipttp(i,j,k).ge.1) then
            stbldqg(i,j,k)=max(stbldqg(i,j,k),.0005)
            stblmqg(i,j,k)=max(stblmqg(i,j,k),.0005)
         endif
      enddo
      enddo
      enddo
c
c   Get fbar
c
      fbar=0.
      do j = 1, mjx-1
      do i = 1, miy-1
         fbar=fbar+cor(i,j)
      enddo
      enddo
      fbar=fbar/((miy-1)*(mjx-1))
c
c   Create effective stability
c
      numpstb=(numpas/100)*100+mod(numpas,100)/4
c      numpstb=0
      do k=1,mkp
      do j = 1, mjx-1
      do i = 1, miy-1
         if (ipttp(i,j,k).ge.1) then
            if (imo.eq.1.and.omgqg(i,j,k).lt.0.0.and.
     &          rhiqg(i,j,k).gt.rhithresh) then
               stblqg(i,j,k)=stblmqg(i,j,k)
            else
               stblqg(i,j,k)=stbldqg(i,j,k)
            endif
         else
            stblqg(i,j,k)=rmsg
         endif
      enddo
      enddo
         if (imo.eq.1) then
            call smooth(stblqg(1,1,k),work,numpstb,miy,miy-1,mjx-1)
         endif
      enddo
c
c   Create right-hand-side term (div-Q follows Reed's notes)
c
      do j = 1, mjx-1
      do i = 1, miy-1
         dsri=xmap(i,j)/ds
         p5dsri=.5*dsri
         dsri2=dsri*dsri
         p25dsri2=.25*dsri2
         p5dsri3=p5dsri*dsri2
         dsr2=1./dsri2
      do k=1,mkp
         if (ipttp(i,j,k).eq.1) then
            dadx=(alphaqg(i,j+1,k)-alphaqg(i,j-1,k))*p5dsri
            dady=(alphaqg(i+1,j,k)-alphaqg(i-1,j,k))*p5dsri
            d2zdx2=(ghtqg(i,j+1,k)-2.*ghtqg(i,j,k)+
     &         ghtqg(i,j-1,k))*dsri2
            d2zdy2=(ghtqg(i+1,j,k)-2.*ghtqg(i,j,k)+
     &         ghtqg(i-1,j,k))*dsri2
            d2zdxdy=(ghtqg(i+1,j+1,k)+ghtqg(i-1,j-1,k)-
     &         ghtqg(i+1,j-1,k)-ghtqg(i-1,j+1,k))*p25dsri2
            if (ivar.eq.0) then  ! full omega, u-bar/s
               omgqg(i,j,k)=omgqg(i,j,k)*1000.   ! hPa/sec to microbar/sec
            elseif (ivar.eq.1) then  ! Q-vec-x, 10**-6 m/s**3/hPa
               omgqg(i,j,k)=1.e6*grav/fbar*(d2zdxdy*dadx-d2zdx2*dady)
            elseif (ivar.eq.2) then  ! Q-vec-y, 10**-6 m/s**3/hPa
               omgqg(i,j,k)=1.e6*grav/fbar*(d2zdy2*dadx-d2zdxdy*dady)
c               if (j.eq.92.and.i.eq.79) then
c                  print*,'k,prsqg,grav,fbar,d2zdy2,d2zdxdy,dadx,dady='
c                  print*,k,prsqg(k),grav,fbar,d2zdy2,d2zdxdy,dadx,dady
c               endif
            else
               d2adx2=(alphaqg(i,j+1,k)-2.*alphaqg(i,j,k)+
     &            alphaqg(i,j-1,k))*dsri2
               d2ady2=(alphaqg(i+1,j,k)-2.*alphaqg(i,j,k)+
     &            alphaqg(i-1,j,k))*dsri2
               d2adxdy=(alphaqg(i+1,j+1,k)+alphaqg(i-1,j-1,k)-
     &            alphaqg(i+1,j-1,k)-alphaqg(i-1,j+1,k))*p25dsri2
               d3zdx2dy=(ghtqg(i+1,j-1,k)-2.*ghtqg(i+1,j,k)+
     &            ghtqg(i+1,j+1,k)-ghtqg(i-1,j-1,k)+2.*ghtqg(i-1,j,k)-
     &            ghtqg(i-1,j+1,k))*p5dsri3
               d3zdxdy2=(ghtqg(i-1,j+1,k)-2.*ghtqg(i,j+1,k)+
     &            ghtqg(i+1,j+1,k)-ghtqg(i-1,j-1,k)+2.*ghtqg(i,j-1,k)-
     &            ghtqg(i+1,j-1,k))*p5dsri3
               d3zdx3=(2.*ghtqg(i,j-1,k)-2.*ghtqg(i,j+1,k)+
     &            ghtqg(i,j+2,k)-ghtqg(i,j-2,k))*p5dsri3
               d3zdy3=(2.*ghtqg(i-1,j,k)-2.*ghtqg(i+1,j,k)+
     &            ghtqg(i+2,j,k)-ghtqg(i-2,j,k))*p5dsri3
               divq(i,j,k)=grav/fbar*(
     &            d2zdxdy*(d2adx2-d2ady2)+d2adxdy*(d2zdy2-d2zdx2)+
     &            dadx*(d3zdx2dy+d3zdy3)-dady*(d3zdx3+d3zdxdy2))
               if (ivar.eq.3) then
                 if (ivar2.eq.2) then
                  omgqg(i,j,k)=stblqg(i,j,k)
                 elseif (ivar2.eq.3) then
                  omgqg(i,j,k)=-2.*dsr2/stblqg(i,j,k)*divq(i,j,k)
                 elseif (ivar2.eq.4) then
                  omgqg(i,j,k)=fbar*fbar*dsr2/(stblqg(i,j,k)*dp*dp)
                 else
                  omgqg(i,j,k)=1.e12*divq(i,j,k)  ! 10**-12 (s**3 hPa)**-1
                 endif
               elseif (ivar.eq.4) then
                  if (iqvecforc.eq.1) then
                     rhs(i,j,k)=-2.*dsr2/stblqg(i,j,k)*divq(i,j,k)
                  else
                     rhs(i,j,k)=0.0
                  endif
                  fac(i,j,k)=fbar*fbar*dsr2/(stblqg(i,j,k)*dp*dp)
                  if (ipttp(i,j,k-1).eq.3) then ! top bound. above
                     omgqg(i,j,k-1)=0.
                  endif
                  if (ipttp(i,j,k+1).eq.2) then ! bottom bound. below
                     dzdx=(ghtqg(i,j+1,k)-ghtqg(i,j-1,k))*p5dsri
                     dzdy=(ghtqg(i+1,j,k)-ghtqg(i-1,j,k))*p5dsri
                     dztdx=(terqg(i,j+1)-terqg(i,j-1))*p5dsri
                     dztdy=(terqg(i+1,j)-terqg(i-1,j))*p5dsri
c
c                  Set the lower boundary condition (topographic
c                     and Ekman)
c
                     omgqg(i,j,k+1)=0.0
                     if (itopobc.ge.1) then
                        dzdx=(ghtqg(i,j+1,k)-ghtqg(i,j-1,k))*p5dsri
                        dzdy=(ghtqg(i+1,j,k)-ghtqg(i-1,j,k))*p5dsri
                        dztdx=(terqg(i,j+1)-terqg(i,j-1))*p5dsri
                        dztdy=(terqg(i+1,j)-terqg(i-1,j))*p5dsri
                        if (itopobc.eq.1) then
                           omgqg(i,j,k+1)=omgqg(i,j,k+1)-
     &                        grav*grav/(fbar*alphaqg(i,j,k))*
     &                        (-dzdy*dztdx+dzdx*dztdy)
                        elseif (itopobc.eq.2) then
                           zgraddirr=atan2(dzdy,dzdx)+betarot
                           zgradmagr=speedreduc*
     &                        sqrt(dzdx*dzdx+dzdy*dzdy)
                           dzdxr=zgradmagr*cos(zgraddirr)
                           dzdyr=zgradmagr*sin(zgraddirr)
                           omgqg(i,j,k+1)=omgqg(i,j,k+1)-
     &                        grav*grav/(fbar*alphaqg(i,j,k))*
     &                        (-dzdyr*dztdx+dzdxr*dztdy)
                        endif
                     endif
                     if (iekmnbc.eq.1) omgqg(i,j,k+1)=omgqg(i,j,k+1)-
     &                  grav*grav/(fbar*alphaqg(i,j,k))*
     &                     perpfac*dzpbl*(d2zdx2+d2zdy2)
                     do kb=k+2,mkp
                        if (ipttp(i,j,kb).eq.2) then ! bnd. pts. below
                           omgqg(i,j,kb)=omgqg(i,j,k+1)
                        endif
                     enddo
                  endif
               endif
            endif
         else
            if (ivar.le.3) omgqg(i,j,k)=rmsg
         endif
      enddo
      enddo
      enddo
c
c   Do the over-relaxation
c
      if (ivar.eq.4) then
c
c   Over-relaxation iteration loop
c
      iter=0
   40 resmax=0.
      iter=iter+1
c
c   For moist case, recompute stability (and fac and rhs) after every
c   8th iteration.
c
      if (imo.eq.1.and.iter.gt.1.and.mod(iter-1,8).eq.0) then
         write(iup,*)'   Recomputing stability at iter=',iter
c
c      First put new stability in the array "fac", then mix it
c      with old stability and put in stblqg.
c
         do k=1,mkp
         do j=1,mjx-1
         do i=1,miy-1
            if (ipttp(i,j,k).ge.1) then
               if (omgqg(i,j,k).lt.0.0.and.
     &             rhiqg(i,j,k).gt.rhithresh) then
                  fac(i,j,k)=stblmqg(i,j,k)
               else
                  fac(i,j,k)=stbldqg(i,j,k)
               endif
            else
               fac(i,j,k)=rmsg
            endif
         enddo
         enddo
            call smooth(fac(1,1,k),work,numpstb,miy,miy-1,mjx-1)
            do j=1,mjx-1
            do i=1,miy-1
               stblqg(i,j,k)=fac(i,j,k)
c               if (ipttp(i,j,k).ge.1.and.
c     &             fac(i,j,k).ne.stblqg(i,j,k)) then
c                  stblqg(i,j,k)=frnew*fac(i,j,k)+frold*stblqg(i,j,k)
c               endif
            enddo
            enddo
         enddo
         do j=1,mjx-1
         do i=1,miy-1
            dsri=xmap(i,j)/ds
            dsri2=dsri*dsri
            dsr2=1./dsri2
         do k=1,mkp
            if (ipttp(i,j,k).eq.1) then
               if (iqvecforc.eq.1) then
                  rhs(i,j,k)=-2.*dsr2/stblqg(i,j,k)*divq(i,j,k)
               endif
               fac(i,j,k)=fbar*fbar*dsr2/(stblqg(i,j,k)*dp*dp)
            endif
         enddo
         enddo
         enddo
      endif
c
c   Do the over-relaxation
c
      do k=2,mkp-1
      do j=3,mjx-3
      do i=3,miy-3
         if (ipttp(i,j,k).eq.1) then
            aterm=4.+2.*fac(i,j,k)
c            res=omgqg(i-1,j,k)+omgqg(i+1,j,k)+omgqg(i,j-1,k)+
c     &         omgqg(i,j+1,k)+fac(i,j,k)*(omgqg(i,j,k+1)+
c     &         omgqg(i,j,k-1))-aterm*omgqg(i,j,k)-rhs(i,j,k)
            res=1./stblqg(i,j,k)*(stblqg(i-1,j,k)*omgqg(i-1,j,k)+
     &         stblqg(i+1,j,k)*omgqg(i+1,j,k)+stblqg(i,j-1,k)*
     &         omgqg(i,j-1,k)+stblqg(i,j+1,k)*omgqg(i,j+1,k))+
     &         fac(i,j,k)*(omgqg(i,j,k+1)+
     &         omgqg(i,j,k-1))-aterm*omgqg(i,j,k)-rhs(i,j,k)
c            if (i.gt.4.and.i.lt.miy-4.and.j.gt.4.and.j.lt.mjx-4.and.
c     &          k.gt.3) then
               resmax=max(resmax,abs(res))
c            endif
            omgqg(i,j,k)=omgqg(i,j,k)+alpha/aterm*res
         endif
      enddo
      enddo
      enddo
      write(iup,'(a,i4,f13.5)')'iter,resmax=',iter,resmax
      if (resmax.gt.errmin.and.iter.le.itmax) goto 40
      if (resmax.gt.errmin.and.iter.gt.itmax) then
         write(iup,*)'   In QGOMG: Didn''t converge in ',itmax,
     &      ' iterations.'
      else
         write(iup,*)'   In QGOMG: Converged in ',iter,' iterations.'
      endif
c
      endif
c
c      vardesc='qgomg field  '
c      write(vardesc(13:13),'(i1)')ivar
c      plchun='?'
c      write(varname,'(a5,i1)') 'omgqg',ivar
c      call writefile(omgqg,varname,numpas,
c     &   3,1,vardesc,plchun,ihrip,rhrip,chrip,casename,iendc,
c     &   cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,miy,mjx,mkp)
c
c   Interpolate back to original model levels.
c
      do j = 1, mjx-1
      do i = 1, miy-1
c
      kvmax=-10000
      kvmin=10000
      do kp=1,mkp
         if ((ivar.eq.4.and.ipttp(i,j,kp).gt.0).or.
     &       (ivar.lt.4.and.ipttp(i,j,kp).eq.1)) then
            kvmax=max(kvmax,kp)
            kvmin=min(kvmin,kp)
         endif
      enddo
      if (kvmax-kvmin.lt.1) then
         do ks=1,mkzh
            omg(i,j,ks)=rmsg
         enddo
      else
         do ks=1,mkzh
            if (prs(i,j,ks).le.prsqg(kvmax).and.
     &          prs(i,j,ks).ge.prsqg(kvmin)) then
               do kp=kvmin,kvmax-1
                  if (prs(i,j,ks).le.prsqg(kp+1).and.
     &                prs(i,j,ks).ge.prsqg(kp)) then
                     omg(i,j,ks)=((prs(i,j,ks)-prsqg(kp))*
     &                   omgqg(i,j,kp+1)+(prsqg(kp+1)-prs(i,j,ks))*
     &                   omgqg(i,j,kp))*dpi
                     goto 50
                  endif
               enddo
 50            continue
            elseif (prs(i,j,ks).gt.prsqg(kvmax)) then
               omg(i,j,ks)=omgqg(i,j,kvmax)
            elseif (prs(i,j,ks).lt.prsqg(kvmin)) then
               omg(i,j,ks)=omgqg(i,j,kvmin)
            endif
            if (ivar.eq.4) omg(i,j,ks)=omg(i,j,ks)*1000. ! hPa/s to ubar/s
         enddo
      endif
c
      enddo
      enddo
c
      return
      end
