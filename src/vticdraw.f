c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vticdraw(ilinw,icolr,xseclen,raxtd,raxld,cvcor,
     &         rtslb,nscrs,raxlv,raxtv,lnogd,vcground,
     &         rlata,rlona,rlatb,rlonb,lableft,labright,
     &         vv1,vv2,set1,set2,vtickinc,mabpl,maxpl,ipl)
c
      dimension ilinw(maxpl),icolr(maxpl),raxtd(maxpl),raxld(maxpl),
     &   raxlv(maxpl),raxtv(maxpl),vcground(mabpl),
     &   rtslb(maxpl)
      logical lnogd(maxpl)
      character cvcor(maxpl)*1
c
      dimension arrl(1000),rect(4),fillx(1000),filly(1000),
     &   fillsc1(1010),fillsc2(1020)
      character axlab*5,axtit*32, left*2, right*2
      character(len=40), dimension(maxpl) :: lableft
      character(len=40), dimension(maxpl) :: labright
c
      include 'comconst'
c
c   Set line width and color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
      tickl=.005  ! length of small tick in frac. coords.
      ydist=(.6*rtslb(ipl)+.01)
      ydist2=(3.0*rtslb(ipl)+.01)
      xdist=.01
      xdist2=(5.4*rtslb(ipl)+.01)
c
c   Draw perimeter of x-section.
c
      fb=fbmin
      ft=ftmax
      fl=flmin
      fr=frmax
c
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call line(fl,fb,fl,ft)
      call line(fl,ft,fr,ft)
      call line(fr,ft,fr,fb)
      call line(fr,fb,fl,fb)
c
c   Draw distance tick marks and label distance axis.
c
      rntickinc=dskm*xseclen
      if (raxld(ipl).eq.rmsg) then
         if (rntickinc.lt.25.) then
            iaxld=1
         elseif (rntickinc.lt.250.) then
            iaxld=10
         elseif (rntickinc.lt.2500.) then
            iaxld=100
         else
            iaxld=1000
         endif
      else
         iaxld=nint(raxld(ipl))
      endif
      if (raxtd(ipl).eq.rmsg) then
         if (rntickinc.lt.25.) then
            iaxtd=1
         elseif (rntickinc.lt.250.) then
            iaxtd=1
         elseif (rntickinc.lt.2500.) then
            iaxtd=10
         else
            iaxtd=100
         endif
      else
         iaxtd=nint(raxtd(ipl))
      endif
      if (iaxtd.le.0.and.iaxld.le.0) goto 600
      call setusv('LW',1000)
      call plchhq(.5,fb-ydist2,'Distance (km)',1.2*rtslb(ipl),0.,0.)
      call setusv('LW',lwidth)
      clat = (rlata+rlatb)*.5
      clon = (rlona+rlonb)*.5
      xlat = (rlata - clat)
      xlon = (rlona - clon)
      ang = atan2(-xlon,-xlat) * 180./3.141592653
      if ( ang .lt. 0.) ang = ang + 360.
      if ( ang .ge. 22.5 .and. ang .lt. 67.5) then
	left = 'SW'
	right = 'NE'
      else if ( ang .ge. 67.5 .and. ang .lt. 112.5) then
	left = 'W '
	right = 'E '
      else if ( ang .ge. 112.5 .and. ang .lt. 157.5) then
	left = 'NW'
	right = 'SE'
      else if ( ang .ge. 157.5 .and. ang .lt. 202.5) then
	left = 'N '
	right = 'S '
      else if ( ang .ge. 202.5 .and. ang .lt. 247.5) then
	left = 'NE'
	right = 'SW'
      else if ( ang .ge. 247.5 .and. ang .lt. 292.5) then
	left = 'E '
	right = 'W '
      else if ( ang .ge. 292.5 .and. ang .lt. 337.5) then
	left = 'SE'
	right = 'NW'
      else if ( ang .ge. 337.5 .or. ang .le. 22.5 ) then
	left = 'S '
	right = 'N '
      else 
	left = '  '
	right = '  '
      endif
      call setusv('LW',1000)
c Add directional labels or user-supplied labels to the bottom of the cross-section.
c The label's edges can be adjusted by changing the final argument to the plchhq call.
c By default, the directions are centered, while the user-supplied labels are left-justified
c (left label) or right-justified (right label). Centering them might be OK if they're short.
      if (lableft(ipl).eq. 'default') then
         call plchhq(fl,fb-ydist2,left,1.2*rtslb(ipl),0.,-1.)
      else
         call plchhq(fl,fb-ydist2,trim(lableft(ipl)),
     &               1.2*rtslb(ipl),0.,-1.)
      endif
      if (labright(ipl).eq. 'default') then
         call plchhq(fr,fb-ydist2,right,1.2*rtslb(ipl),0.,0.)
      else
         call plchhq(fr,fb-ydist2,trim(labright(ipl)),
     &               1.2*rtslb(ipl),0.,1.)
      endif
      call setusv('LW',lwidth)
      ntick=int(rntickinc)+1
      gint=(fr-fl)/rntickinc
      if (iaxld.gt.0) then
         do j=1,ntick
            if(mod(j-1,iaxld).eq.0) then
               irr=(j-1)
               rrloc=fl+real(j-1)*gint
               write(axlab,'(i5)')irr
               iaxs=5-int(alog10(float(irr)+.5))
               call setusv('LW',1000)
               call plchhq(rrloc,fb-ydist,axlab(iaxs:),
     &            rtslb(ipl),0.,0.)
               call setusv('LW',lwidth)
               call line(rrloc,fb,rrloc,fb+2.*tickl)
               call line(rrloc,ft,rrloc,ft-2.*tickl)
            endif
         enddo
      endif
      if (iaxtd.gt.0) then
         do j=1,ntick
            if(mod(j-1,iaxtd).eq.0) then
               rrloc=fl+real(j-1)*gint
               call line(rrloc,fb,rrloc,fb+tickl)
               call line(rrloc,ft,rrloc,ft-tickl)
            endif
         enddo
      endif
 600  continue
c
c   Do vertical axis
c
      call set(fl,fr,fb,ft,fl,fr,set1,set2,1)
      call gqclip(ierr,iclp,rect)
      call gsclip(0)
c
      if (raxlv(ipl).eq.rmsg) then
         iaxlv=10
      else
         iaxlv=nint(raxlv(ipl)/vtickinc)
      endif
      if (raxtv(ipl).eq.rmsg) then
         iaxtv=1
      else
         iaxtv=nint(raxtv(ipl)/vtickinc)
      endif
      if (iaxtv.le.0.and.iaxlv.le.0) goto 650
c
c   Write vertical axis title
c
      if (cvcor(ipl).eq.'s') then ! these are in vertical level index
         axtit='Vertical Level Index'
         nch=20
      elseif (cvcor(ipl).eq.'p') then ! these are in hPa
         axtit='Pressure (hPa)'
         nch=14
      elseif (cvcor(ipl).eq.'l') then ! b
         axtit='Pressure (hPa), log scale'
         nch=25
      elseif (cvcor(ipl).eq.'x') then ! these are in hPa
         axtit='Pressure (hPa), Exner scale'
         nch=27
      elseif (cvcor(ipl).eq.'z') then ! these are in km
         axtit='Height (km)'
         nch=11
      elseif (cvcor(ipl).eq.'f') then ! these are in km
         axtit='Height AFL (km)'
         nch=11
      elseif (cvcor(ipl).eq.'t') then ! these are in K
         axtit='Theta (K)'
         nch=9
      elseif (cvcor(ipl).eq.'m') then ! these are in deg. C
         axtit='Temperature (dg C)'
         nch=18
      elseif (cvcor(ipl).eq.'e') then ! these are in K
         axtit='Theta_e (K)'
         nch=11
      elseif (cvcor(ipl).eq.'q') then ! these are in PVU
         axtit='PV (PVU)'
         nch=8
      endif
      call setusv('LW',1000)
      call plchhq(fl-xdist2,.5*(set1+set2),axtit(1:nch),
     &   1.2*rtslb(ipl),90.,0.)
      call setusv('LW',lwidth)
c
c   Make tick marks
c
      htvc=vv2-vv1
      ntick=nint(abs(htvc)/vtickinc)+1
      if (htvc.eq.0.) then
         write(iup,*)'Zero depth of cross section.'
         stop
      endif
      idir=nint(htvc/abs(htvc))
      if (iaxlv.gt.0) then
         do i=1,ntick
            rr=vv1+idir*(i-1)*vtickinc
            irrinc=nint(rr/vtickinc)
            if (cvcor(ipl).eq.'l') then
               vtick=alog(rr)
            elseif (cvcor(ipl).eq.'x') then
               vtick=(rr)**gamma
            else
               vtick=rr
            endif
            if (mod(irrinc,iaxlv).eq.0) then
               if (cvcor(ipl).eq.'z'.or.cvcor(ipl).eq.'f') then
                  write(axlab,'(f4.1)') rr
               else
                  write(axlab,'(i4)') nint(rr)
               endif
               call setusv('LW',1000)
               call plchhq(fl-xdist,vtick,axlab,rtslb(ipl),0.,1.)
               call setusv('LW',lwidth)
               call line(fl,vtick,fl+2.*tickl,vtick)
               call line(fr,vtick,fr-2.*tickl,vtick)
            endif
         enddo
      endif
      if (iaxtv.gt.0) then
         do i=1,ntick
            rr=vv1+idir*(i-1)*vtickinc
            irrinc=nint(rr/vtickinc)
            if (cvcor(ipl).eq.'l') then
               vtick=alog(rr)
            elseif (cvcor(ipl).eq.'x') then
               vtick=(rr)**gamma
            else
               vtick=rr
            endif
            if (mod(irrinc,iaxtv).eq.0) then
               call line(fl,vtick,fl+tickl,vtick)
               call line(fr,vtick,fr-tickl,vtick)
            endif
         enddo
      endif
      call gsclip(iclp)
 650  continue
c
c   Draw surface curve.
c
      call set(fl,fr,fb,ft,1.,float(nscrs),set1,set2,1)
      do 700 l=1,nscrs
         arrl(l)=float(l)
  700 continue
      call gqlwsc (ierr, oldw)
      call gslwsc (2.0)
      call curve(arrl,vcground,nscrs)
      call gslwsc (oldw)
c
c   Fill in the underground area.
c
      if (iplevdata.lt.4.and..not.lnogd(ipl)) then
c
      rleftfill=cfux(fl+3.*tickl)
      rrightfill=cfux(fr-3.*tickl)
      rbotfill=cfuy(fb+3.*tickl)
      do l=1,nscrs-1
         if (rleftfill.ge.arrl(l).and.rleftfill.le.arrl(l+1)) then
            l1=l
            l2=l+1
         endif
         if (rrightfill.ge.arrl(l).and.rrightfill.le.arrl(l+1)) then
            l3=l
            l4=l+1
         endif
      enddo
c
c   Start at upper right, work around clockwise
c
      fillx(1)=rleftfill
      filly(1)=(rleftfill-float(l1))*vcground(l2)+
     &         (float(l2)-rleftfill)*vcground(l1)
      do l=l2,l3
         ll=l-l2+2
         fillx(ll)=float(l)
         filly(ll)=vcground(l)
      enddo
      fillx(l3-l2+3)=rrightfill
      filly(l3-l2+3)=(rrightfill-float(l3))*vcground(l4)+
     &         (float(l4)-rrightfill)*vcground(l3)
      fillx(l3-l2+4)=rrightfill
      filly(l3-l2+4)=rbotfill
      fillx(l3-l2+5)=rleftfill
      filly(l3-l2+5)=rbotfill
      npts=l3-l2+5
c      do l=1,npts
c         print*,'l,x,y=',l,fillx(l),filly(l)
c      enddo
      call sfseti('TY',1)
      call sfseti('DO',0)
      call sfsetr('AN',50.)
      call sfsetr('SP',.005)
      call sfsgfa(fillx,filly,npts,fillsc1,1010,fillsc2,1020,
     &   icolr(ipl))
c
c   Must set 'TY' back to 0 to prevent some subsequent color fills from
c   using hatch pattern instead of solid.
c
      call sfseti('TY',0)
c
      endif   ! done filling underground area
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
c
      return
      end



