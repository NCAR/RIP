c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hcondraw(xtime,ilinw,vc3d,tmk,qvp,
     &         prs,ght,ter,sfp,sfpsm,ixwin,iywin,
     &         ismth,rcint,rcbeg,rcend,rcval,incvl,lmult,larng,
     &         idash,rlevl,rlavl,cnohl,lnolb,lnobr,lnozr,incon,
     &         bottextfloor,cfeld,cvcor,work,icdwk,unwk,ilev,lpslb,
     &         icolr,icoll,ilcll,ilchl,rtslb,rtshl,imjsk,icomg,
     &         lnmsg,icong,iconl,icozr,idimn,lhide,lgrad,lhadv,
     &         ilwll,ilwng,ilwnl,ilwzr,idall,
     &         idang,idanl,idazr,ilcnl,ilczr,lcord,
     &         ilcbr,ipwlb,iorlb,ipwhl,ipwbr,ifclb,ifcnl,
     &         ifczr,ifchl,ilclo,ifclo,ccmth,rwdbr,ihvbr,rcfad,ncfadbin,
     &         idotser,tseryi,tserxj,tserdat,ntsers,ntsert,ntserv,
     &         icosq,rcosq,incsq,fred,fgreen,fblue,nco,icomax,pslab1,
     &         pslabt,ipslab,iam,xcs,ycs,niam,ncs,idiffflag,
     &         maxtserv,maxtsers,maxtsert,maxcosq,mabpl,morpl,
     &         maxlev,maxpl,maxcon,miy,mjx,mkzh,ipl,irota,iwkidcgm,
     &         noplots)
c
      parameter (iwrklen=100000)
c
      dimension vc3d(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ter(miy,mjx),sfp(miy,mjx),sfpsm(miy,mjx),ixwin(2,maxpl),
     &   iywin(2,maxpl),ismth(maxpl),
     &   rcint(maxpl),rcbeg(maxpl),rcend(maxpl),
     &   rcval(maxcon,maxpl),incvl(maxpl),
     &   idash(maxpl),rlevl(maxlev,maxpl),incsq(maxpl),
     &   rlavl(maxlev,maxpl),ilinw(maxpl),incon(maxpl),icomg(maxpl),
     &   icolr(maxpl),icoll(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),imjsk(maxpl),idimn(maxpl),
     &   icosq(maxcosq,maxpl),rcosq(maxcosq,maxpl),
     &   icong(maxpl),iconl(maxpl),icozr(maxpl),ilwll(maxpl),
     &   ilwng(maxpl),ilwnl(maxpl),ilwzr(maxpl),idall(maxpl),
     &   idang(maxpl),idanl(maxpl),idazr(maxpl),
     &   ilcnl(maxpl),ilczr(maxpl),
     &   ilcbr(maxpl),ipwlb(maxpl),iorlb(maxpl),
     &   ipwhl(maxpl),ipwbr(maxpl),ifclb(maxpl),ifcnl(maxpl),
     &   ifczr(maxpl),ifchl(maxpl),
     &   ilclo(maxpl),ifclo(maxpl),rwdbr(maxpl),ihvbr(maxpl),
     &   rcfad(3,maxpl),work(miy,mjx,mkzh),icdwk(maxpl),
     &   tserdat(maxtsert,maxtserv,maxtsers),tseryi(maxtsers),
     &   tserxj(maxtsers),fred(0:255),fgreen(0:255),fblue(0:255),
     &   pslab1(mabpl,morpl),pslabt(mabpl,morpl),ipslab(mabpl,morpl),
     &   irota(maxpl)
      dimension ncfadcount(ncfadbin)
      dimension iam(niam),xcs(ncs),ycs(ncs)
      logical lnolb(maxpl),lnobr(maxpl),lnozr(maxpl),lpslb(maxpl),
     &   lmult(maxpl),larng(maxpl),lnmsg(maxpl),lhide(maxpl),
     &   lgrad(maxpl),lhadv(maxpl),lcord(maxpl),lhidel
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1,ccmth(maxpl)*4,
     &   unwk(maxpl)*24,cnohl(maxpl)*1
c
      dimension redcosq(500),greencosq(500),bluecosq(500)
      dimension rwrk(iwrklen),iwrk(iwrklen),rect(4),
     &   valcon(maxcon),majcon(maxcon),lbbarcon(maxcon),
     &   iaia(500),igia(500),lfin(maxcon+1)
      character messg*126,llbs(maxcon+2)*24
      external drawcl,cpcolr,mask
c
      parameter (nsp=6)
      dimension spx(nsp),spy(nsp),spt(nsp)
c      data (spt(i),i=1,nsp) / 24.,25.,27.,30.,33.,36. /
c      data (spx(i),i=1,nsp) / 85.,87.,90.,95.,98.,99. /
c      data (spy(i),i=1,nsp) / 56.,57.,58.,69.,70.,73. /
      data (spt(i),i=1,nsp) / 0.,1.,3.,6.,9.,12. /
      data (spx(i),i=1,nsp) / 37.,37.,39.,41.,43.,46. /
      data (spy(i),i=1,nsp) / 42.,42.,44.,47.,51.,56. /
c
      include 'comconst'
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
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
      if (icdwk(ipl).eq.0) then
         ul=1.
         ur=1.+xintervs
         ub=1.
         ut=1.+yintervs
         niy=nint(ut)
         njx=nint(ur)
      else
         ul=.5
         ur=.5+xintervs
         ub=.5
         ut=.5+yintervs
         niy=int(ut)
         njx=int(ur)
      endif
      if (noplots.eq.0) then
c
c   Transpose ul,ur,ub,ut if view is rotated +/- 90 degrees.
c   niy, njx, and pslab arrays will not be adjusted for rotation
c   until a later point.
c
      if (irota(ipl).eq.90.or.irota(ipl).eq.-90) then
         usv=ul
         ul=ub
         ub=usv
         usv=ur
         ur=ut
         ut=usv
      endif
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
      endif    ! noplots
c
c   Put appropriate data into horizontal slab.
c
      if ((cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).eq.
     &   rlavl(ilev,ipl)).or.idimn(ipl).eq.2) then
         do 90 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 80 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work(ii,jj,nint(rlevl(ilev,ipl)))
   80       continue
   90    continue
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).ge.0) then
         call fillarray(pslab1,mabpl*morpl,0.)
         lev1=min(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         lev2=max(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         do j=1,njx
            jj=j+ixwin(1,ipl)-1
         do i=1,niy
            ii=i+iywin(1,ipl)-1
            prstot=prs(ii,jj,lev2)-prs(ii,jj,lev1)
            do k=lev1,lev2-1
               pslab1(j,i)=pslab1(j,i)+.5*(work(ii,jj,k)+
     &            work(ii,jj,k+1))*(prs(ii,jj,k+1)-prs(ii,jj,k))/prstot
            enddo
         enddo
         enddo
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).lt.0) then
         lev1=nint(rlevl(ilev,ipl))
         lev2=nint(-rlavl(ilev,ipl))
         do j=1,njx
            jj=j+ixwin(1,ipl)-1
            do i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work(ii,jj,lev1)-work(ii,jj,lev2)
            enddo
         enddo
      else
         if (lgrad(ipl).or.lhadv(ipl)) then
            lhidel=.true.
         else
            lhidel=lhide(ipl)
         endif
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,lhidel,idiffflag,cfeld(1,ipl),work,
     &      pslab1,mabpl,morpl,njx,niy,miy,mjx,mkzh)
      endif
c      PSMAX=-500.
c      PSMIN=500.
c      DO 150 J=1,NJX
c      DO 150 I=1,NIY
c         PSMAX=MAX(PSMAX,PSLAB1(J,I))
c         PSMIN=MIN(PSMIN,PSLAB1(J,I))
c  150 CONTINUE
c      WRITE(IUP,*)'MAX, MIN TEMP.=',PSMAX,PSMIN
c
c      print*,'pslab1(90,86)=',pslab1(90,86)
c
c   Smooth field if called for.
c
      if (noplots .eq. 0) 
     &   call smooth(pslab1,pslabt,ismth(ipl),mabpl,njx,niy)
c
c   Print slab to file if asked for, in column-major order
c
      if (lpslb(ipl)) then
         write(59)((pslab1(j,i),i=1,niy),j=1,njx)
      endif
c
c   Do CFAD counting and print out
c
      if (ncfadbin.gt.0) then
         do ib=1,ncfadbin
            ncfadcount(ib)=0
         enddo
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg) then
               ibin=nint((pslab1(j,i)-rcfad(1,ipl))/rcfad(2,ipl))+1
               if (ibin.ge.1.and.ibin.le.ncfadbin) then
                  ncfadcount(ibin)=ncfadcount(ibin)+1
               endif
            endif
         enddo
         enddo
         do ib=1,ncfadbin
            write(58,*)rlevl(ilev,ipl),
     &                 rcfad(1,ipl)+(ib-1)*rcfad(2,ipl),
     &                 ncfadcount(ib)
         enddo
      endif
c
c   Calculate time series data if called for.
c
      do 300 istn=1,idotser*ntsers
         posy=tseryi(istn)-iywin(1,ipl)+1.-.5*icdwk(ipl)
         posx=tserxj(istn)-ixwin(1,ipl)+1.-.5*icdwk(ipl)
         if (posx.le.1..or.posx.ge.float(njx).or.
     &       posy.le.1..or.posy.ge.float(niy)) goto 290
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         if (pslab1(jl,it).eq.rmsg.or.
     &       pslab1(jr,it).eq.rmsg.or.
     &       pslab1(jl,ib).eq.rmsg.or.
     &       pslab1(jr,ib).eq.rmsg) goto 290
         tserdat(ntsert,ntserv,istn)=
     +      (1.-ratlr)*(   ratbt)*pslab1(jl,it)+
     +      (   ratlr)*(   ratbt)*pslab1(jr,it)+
     +      (1.-ratlr)*(1.-ratbt)*pslab1(jl,ib)+
     +      (   ratlr)*(1.-ratbt)*pslab1(jr,ib)
         goto 300
  290    tserdat(ntsert,ntserv,istn)=rmsg
  300 continue
cc
cc   Do low center thing
cc
c      if (lnobr(ipl)) then
c         do isp=1,nsp
c            if (abs(spt(isp)-xtime).le..15) then
c               ispg=isp
c               goto 543
c            endif
c         enddo
c         stop 333
c 543     continue
c         dist=800./dskm
c         pouter=0.
c         do ideg=0,350,10
c            xring=spx(ispg)+dist*cos(rpd*ideg)
c            yring=spy(ispg)+dist*sin(rpd*ideg)
c            ixring=nint(min(max(xring,1.),mjx-1.))
c            iyring=nint(min(max(yring,1.),miy-1.))
c            pouter=pouter+pslab1(ixring,iyring)
c         enddo
c         pouter=pouter/36.
cc         write(iup,*) 'x,y=',nint(spx(ispg)),nint(spy(ispg))
c         pcentral=pslab1(nint(spx(ispg)),nint(spy(ispg)))
c         write(iup,*) 'time,pcentral,pouter,pdiff=',
c     &      spt(ispg),pcentral,pouter,pcentral-pouter
c      endif
c     if (rwdbr(ipl).ge.100.) then  ! one-time code - ignore this stuff
c        gradient=rcosq(1,ipl)  !  in hPa/km
c        direction=rcosq(2,ipl)  ! degrees
c        valinterior=rcosq(3,ipl)
c        xinterior=41.
c        yinterior=47.
cc mod this too...
c        do j=1,mjx
c        do i=1,miy
c           gval=valinterior+gradient*dskm*(
c    &         (j-xinterior)*cos(rpd*direction)+
c    &         (i-yinterior)*sin(rpd*direction))
c           if (rwdbr(ipl).eq.100.) then
c              pslab1(j,i)=gval
c           else
c              pslab1(j,i)=pslab1(j,i)+gval
c           endif
c        enddo
c        enddo
c     endif
c
c   Determine number of contours and contour values.
c
      if (noplots.eq.0) then
      call getconvals(rcbeg(ipl),rcend(ipl),rcint(ipl),
     &   incon(ipl),lmult(ipl),rcval(1,ipl),incvl(ipl),
     &   imjsk(ipl),pslab1,mabpl,njx,niy,rmsg,
     &   maxcon,valcon,majcon,cintuse,numcon,iup)
c
      call cpseti('SET',0)   ! we'll use our own set call
      call cpseti('CLS',0)   ! we'll use our own cont. values
c
c   Some settings for hi/lo and contour labels
c
      call cpseti('NSD',-4)
      call cpseti('NOF',7)
      call cpsetr ('PC6',.6)
c
c   Filling:
c
      if ((ccmth(ipl).eq.'fill'.or.ccmth(ipl).eq.'cell'.or.
     &     ccmth(ipl).eq.'both'.or.ccmth(ipl).eq.'ceco').and.
     &    numcon.gt.0) then
c
c   Assign the red, green, and blue fractions to the color sequence
c
      do i=1,incsq(ipl)
         do j=0,nco
            if (icosq(i,ipl).eq.j) then
               redcosq(i)=fred(j)
               greencosq(i)=fgreen(j)
               bluecosq(i)=fblue(j)
               goto 320
            endif
         enddo
         redcosq(i)=1.
         greencosq(i)=1.
         bluecosq(i)=1.
 320     continue
      enddo
c
c   Colors are assigned to ranges of values bounded by adjacent
c   contour values, as well as the semi-infinite ranges below
c   the lowest and above the highest contours.  Therefore, numcon+1
c   colors need to be defined. In order to choose a color from the
c   color sequence, each range must be assigned a specific value.
c   That value will be set equal to the average of the two
c   surrounding contour values.  The value for the range
c   <minus infinity to valcon(1)> will be assigned a value of the
c   lowest contour minus one half the difference between the lowest
c   and second lowest contours.  A similar definition applies to
c   the value for the range <valcon(numcon) to plus infinity>.
c   Based on the specific value assigned to each range, a color will
c   be defined based on a linear interpolation from the color sequence.
c
      nconarea=numcon+1
      mm=incsq(ipl)
c
      if (numcon.gt.1) then
         valmin=valcon(1)-.5*(valcon(2)-valcon(1))
         valmax=valcon(numcon)+.5*(valcon(numcon)-valcon(numcon-1))
      else
c
c      If there is only one contour value generated, then
c      set valmin to the average of all grid values that are less than
c      the single contour value, and valmax to the average of all
c      grid values that are greater than the single contour value, 
c
         numgt=0
         avggt=0.
         numlt=0
         avglt=0.
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg.and.
     &          pslab1(j,i).gt.valcon(1)) then
               numgt=numgt+1
               avggt=avggt+pslab1(j,i)
            elseif (pslab1(j,i).ne.rmsg.and.
     &              pslab1(j,i).lt.valcon(1)) then
               numlt=numlt+1
               avglt=avglt+pslab1(j,i)
            endif
         enddo
         enddo
         if (numgt.gt.0) then
            valmax=avggt/numgt
	 else
            valmax=valcon(1)
         endif
         if (numlt.gt.0) then
            valmin=avglt/numlt
	 else
            valmin=valcon(1)
         endif
      endif
c
      do i=1,numcon+1
         red=-1.
c
c      Note, in the following code, if icosc<0, that means the user
c      chose "transparent" as the color, which means the value of
c      "red" should remain -1.
c
         if (.not.lcord(ipl)) then
            if (i.eq.1) then
               val=valmin
            elseif (i.eq.numcon+1) then
               val=valmax
            else
               val=.5*(valcon(i-1)+valcon(i))
            endif
            if (larng(ipl)) val=100.*(val-valmin)/(valmax-valmin)
         else
            val=float(i)
            if (larng(ipl)) val=1.+99.*(val-1.)/float(numcon)
         endif
         if (val.le.rcosq(1,ipl)) then
            if (icosq(1,ipl).gt.0) then
               red=redcosq(1)
               green=greencosq(1)
               blue=bluecosq(1)
            endif
         elseif (val.ge.rcosq(mm,ipl)) then
            if (icosq(mm,ipl).gt.0) then
               red=redcosq(mm)
               green=greencosq(mm)
               blue=bluecosq(mm)
            endif
         else
            do isq=1,incsq(ipl)-1
               if (val.ge.rcosq(isq,ipl).and.
     &             val.le.rcosq(isq+1,ipl)) then
                  if (icosq(isq,ipl).gt.0.and.
     &                icosq(isq+1,ipl).gt.0) then
                     fac=(val-rcosq(isq,ipl))/
     &                  (rcosq(isq+1,ipl)-rcosq(isq,ipl))
                     red=fac*redcosq(isq+1)+(1.-fac)*redcosq(isq)
                     green=fac*greencosq(isq+1)+(1.-fac)*greencosq(isq)
                     blue=fac*bluecosq(isq+1)+(1.-fac)*bluecosq(isq)
                  endif
                  goto 330
               endif
            enddo
 330        continue
         endif
c
c      Assign color index to icoindcp
c
         if (red.lt.0.) then
            icoindcp(i)=-1 ! transparent
         else
c
c         First check if color already exists
c
            do ico=2,icomax
               if (red.eq.fred(ico).and.green.eq.fgreen(ico).and.
     &             blue.eq.fblue(ico)) then
                  icoindcp(i)=ico
                  goto 357
               endif
            enddo
c
c         If not, set the new color and assign its index to icoincdcp
c
            icomax=icomax+1
            if (icomax.eq.256) then
               write(iup,*) 'hcondraw is trying to define color number'
               write(iup,*) '256, but 255 is the limit.  Try increasing'
               write(iup,*) 'the contour interval, or reduce the number'
               write(iup,*) 'of colors defined in the color table.'
               stop
            endif
            call gscr(iwkidcgm,icomax,red,green,blue)
            fred(icomax)=red
            fgreen(icomax)=green
            fblue(icomax)=blue
            icoindcp(i)=icomax
         endif
 357     continue
      enddo
c
c  Increase number of vertical strips for efficiency:
c
      call cpseti('NVS',4)
c
      call cpseti('NCL',numcon)
      do i=1,numcon
         call cpseti('PAI',i)
         call cpsetr('CLV',valcon(i))
         call cpseti('CLU',1)
         call cpseti('AIB',i)
         call cpseti('AIA',i+1)
      enddo
c
c   Rotate pslab1 array, and switch njx,niy, if needed.
c
      call rotpslab(pslab1,pslabt,mabpl,morpl,njx,niy,irota(ipl))
c
c   Initialize contouring
c
      call cprect(pslab1,mabpl,njx,niy,rwrk,iwrklen,iwrk,iwrklen)
c
c   Do filling
c
      if (ccmth(ipl).eq.'fill'.or.ccmth(ipl).eq.'both') then
         call arinam(iam,niam)   ! initialize the area map
         call cpclam(pslab1,rwrk,iwrk,iam)  ! put cont. lines in area map
         call arscam(iam,xcs,ycs,ncs,iaia,igia,500,cpcolr)
      elseif (ccmth(ipl).eq.'cell'.or.ccmth(ipl).eq.'ceco') then
c 
c      Set cell array values and map it to user coordinates
c 
         call cpcica (pslab1,rwrk,iwrk,ipslab,mabpl,njx,niy,
     &      cufx(.5001),cufy(.5001),cufx(float(njx)+.4999),
     &      cufy(float(niy)+.4999))
c
c     What ends up in the ipslab array after the call to cpcica is an
c     integer for each cell array element (which, in the case of RIP's
c     usage, refers to a model grid box).  It turns out that this
c     integer is equivalent to the element of icoindcp that should be
c     used as the color index for that cell.  However, routine gca
c     (called below) requires an array of the color indices themselves,
c     not of indices of icoindcp.  Therefore, the following conversion
c     must be performed:
c
         do j=1,njx
         do i=1,niy
            ipslab(j,i)=icoindcp(ipslab(j,i))
         enddo
         enddo
c
c      Then call GCA to fill the cell array in the plot.
c
         call gca (.5,.5,float(njx)+.5,float(niy)+.5,mabpl,morpl,1,1,
     1        njx,niy,ipslab)
         call sflush
c
      endif
c
c   Make label bar
c
      if (.not.lnobr(ipl)) then
c
         call lbseti('CBL',ilcbr(ipl))
         call lbseti('CLB',ilcbr(ipl))
         call lbsetr('WBL',float(ipwbr(ipl)))
         call lbsetr('WLB',float(ipwbr(ipl)))
c
c      Decide whether to orient the bar horizontally or vertically
c
c      Determine anticipated available space at bottom and on right.
c      This needn't be exact, it's just for roughly determining whether a
c      horizontal or vertical label bar will more efficiently use
c      the available space.  For bottom, assume .02 is needed for
c      tick labels (typically present) and perhaps 3 additional fields
c      will be plotted.  For right, assume .045 is needed for lat/lon
c      labels (typically present).
c
         avlbottom=fb-bottextfloor-.02-3.*.0152
         avlright=1.-fr-.045
c
         if (rwdbr(ipl).eq.-1.) then
            reqbottom=.036
            reqright=.07
         else
            reqbottom=rwdbr(ipl)
            reqright=rwdbr(ipl)
         endif
c
         ihov=ihvbr(ipl)
         if (ihov.eq.-1) then
            if (avlbottom-reqbottom.gt.avlright-reqright) then
               ihov=0
            else
               ihov=1
            endif
         endif
         if (ihov.eq.0) then
            xleb=fl
            xreb=fr
            ybeb=bottextfloor
            yteb=ybeb+reqbottom
            wsfb=1.
            hsfb=.45
         else
            xleb=1.-reqright
            xreb=1.
            ybeb=fb
            yteb=ft
            wsfb=.231
            hsfb=1.
         endif
         iftp=0
         lbab=1
c
c      Get label bar labels and colors
c
         lstep=numcon/20+1  ! don't overcrowd label. skip some values
c
c      Use getconvals to decide which contours (color transitions) should
c      be labeled on the label bar.
c
         call getconvals(rcbeg(ipl),rcend(ipl),rcint(ipl),
     &      incon(ipl),lmult(ipl),rcval(1,ipl),incvl(ipl),
     &      lstep-1,pslab1,mabpl,njx,niy,rmsg,
     &      maxcon,valcon,lbbarcon,cintuse,numcon,iup)
c
         nlbs=1
         nbox=1
         if (icoindcp(1).ge.0) then
            lfin(1)=icoindcp(1)
         else
            lfin(1)=0
         endif
         llbs(1)=' '
         do i = 1,numcon
            nlbs=nlbs+1
            nbox=nbox+1
            if ((lbbarcon(i).eq.-2.or.lbbarcon(i).eq.0.or.
     &           lbbarcon(i).eq.2).and.nlbs.le.numcon+2-lstep) then
               call cpseti ('PAI',i)          ! set internal array index
               call cpsetr ('ZDV',valcon(i))  ! give value to a converter
               call cpgetc ('ZDV',llbs(nlbs)) ! get a string back
            else
               llbs(nlbs)=' '
            endif
c            icii=min(i+1+lstep/2,nconarea)
c            if (icoindcp(icii).ge.0) then
c               lfin(nbox)=icoindcp(icii)
            if (icoindcp(nbox).ge.0) then
               lfin(nbox)=icoindcp(nbox)
            else
               lfin(nbox)=0
            endif
         enddo
         nlbs=nlbs+1
         llbs(nlbs)=unwk(ipl)
c
         call lblbar(ihov,xleb,xreb,ybeb,yteb,nbox,wsfb,hsfb,lfin,
     &      iftp,llbs,nlbs,lbab)
         if (ihov.eq.0) bottextfloor=bottextfloor+1.2*reqbottom
c
      endif
c
      endif
c
c   Contours:
c
      if ((ccmth(ipl).eq.'cont'.or.ccmth(ipl).eq.'both'.or.
     &     ccmth(ipl).eq.'ceco').and.numcon.gt.0) then
c
      call cpsetc('ILT',' ') ! no CONPACK informational label
      call cpseti('LBC',-1)  ! use current fill color for lab. boxes
c
c   Set a scale factor which keeps all contour and hi/lo labels the same size
c   regardless of the width of the viewport.
c
      call cpsetr('CWM',1./(fr-fl))
c
c   Set up high/low stuff
c
      if (cnohl(ipl).ne.'B'.and.cnohl(ipl).ne.'b') then
c         call cpsetc('HIT','H~PRL~0~PRU~ $ZDV$')
c         call cpsetc('LOT','L~PRL~0~PRU~ $ZDV$')
         if (cnohl(ipl).eq.'H'.or.cnohl(ipl).eq.'h') then
            call cpsetc('HIT',' ')
         else
            call cpsetc('HIT','H~V-1QH-50~ $ZDV$')
         endif
         if (cnohl(ipl).eq.'L'.or.cnohl(ipl).eq.'l') then
            call cpsetc('LOT',' ')
         else
            call cpsetc('LOT','L~V-1QH-50~ $ZDV$')
         endif
c
c      Following is Jim Bresch's options for nohl, which aloow for nonstandard
c      labeling of extrema.
c
         if (cnohl(ipl).eq.'Z'.or.cnohl(ipl).eq.'z') then
c           Highs and Lows
            call cpsetc('HIT','H')   ! show "H" only, no value
            call cpsetc('LOT','L')   ! show "L" only, no value
         else if (cnohl(ipl).eq.'T'.or.cnohl(ipl).eq.'t') then
c           Warm and Cold extrema
            call cpsetc('HIT','W')   ! show "W" only, no value
            call cpsetc('LOT','C')   ! show "C" only, no value
         else if (cnohl(ipl).eq.'A'.or.cnohl(ipl).eq.'a') then
c           Cyclonic and anticyclonic extrema
            call cpsetc('HIT','A')   ! show "A" only, no value
            call cpsetc('LOT','C')   ! show "C" only, no value
         endif
c
         call cpseti('HIC',ilchl(ipl))
         call cpseti('LOC',ilclo(ipl))
         call cpsetr('HLS',rtshl(ipl))
         if (ipwhl(ipl).ne.0) then
            call cpsetr('HLL',float(ipwhl(ipl)))
         else
            call cpsetr('HLL',1.)
         endif
         if (ifchl(ipl).ne.999999) then
            icpfchl=ifchl(ipl)
            if (ifclo(ipl).ne.999999) then
               icpfclo=ifclo(ipl)
            else
               icpfclo=ifchl(ipl)
            endif
         endif
         if ((ipwhl(ipl).eq.0.and.ifchl(ipl).eq.999999)
     &       .or. ipwhl(ipl) .eq. 1000) then
            ihlb=0
         elseif (ipwhl(ipl).ne.0.and.ifchl(ipl).eq.999999) then
            ihlb=1
         elseif (ipwhl(ipl).eq.0.and.ifchl(ipl).ne.999999) then
            ihlb=2
         else
            ihlb=3
         endif
         call cpseti('HLB',ihlb)
      else
         call cpsetc('HLT',' ')
      endif
c
c   Rotate pslab1 array, and switch njx,niy, if needed.
c
      call rotpslab(pslab1,pslabt,mabpl,morpl,njx,niy,irota(ipl))
c
c   Initialize contouring
c
      call cprect(pslab1,mabpl,njx,niy,rwrk,iwrklen,iwrk,iwrklen)
c
c   Set up contour line stuff
c
      if (iorlb(ipl).eq.1) then
         call cpseti('LLP',3)   ! penalty scheme for line labels
         call cpseti('LLO',0)   ! all labels at angle=0 (hor. orient.)
      elseif (iorlb(ipl).eq.2) then
         call cpseti('LLP',3)   ! penalty scheme for line labels
         call cpseti('LLO',1)   ! all labels oriented along contour
      elseif (iorlb(ipl).eq.3) then
         call cpseti('LLP',1)   ! labels drawn with DASHPAT,
c                                 like in old CONREC
      else
         write(iup,*)'For ipl=',ipl,', invalid value for orlb.'
         write(iup,*)'Should be 1, 2, or 3.'
         stop
      endif
      call cpsetr('LLS',rtslb(ipl))
      call cpseti('NCL',numcon)
      if (.not.lnolb(ipl)) then
         if (ipwlb(ipl).ne.0) then
            call cpsetr('LLL',float(ipwlb(ipl)))
         else
            call cpsetr('LLL',1.)
         endif
         if (ifclb(ipl).ne.999999) then
            icpfclb=ifclb(ipl)
            if (ifcnl(ipl).ne.999999) then
               icpfcnl=ifcnl(ipl)
            else
               icpfcnl=ifclb(ipl)
            endif
            if (ifczr(ipl).ne.999999) then
               icpfczr=ifczr(ipl)
            else
               icpfczr=ifclb(ipl)
            endif
         endif
         if ((ipwlb(ipl).eq.0.and.ifclb(ipl).eq.999999)
     &       .or. ipwhl(ipl) .eq. 1000) then
            illb=0
         elseif (ipwlb(ipl).ne.0.and.ifclb(ipl).eq.999999) then
            illb=1
         elseif (ipwlb(ipl).eq.0.and.ifclb(ipl).ne.999999) then
            illb=2
         else
            illb=3
         endif
         call cpseti('LLB',illb)
      endif
c
      do i=1,numcon
         call cpseti('PAI',i)
         call cpsetr('CLV',valcon(i))
         mconcp(i)=majcon(i)
         isetllc=-1
         if (majcon(i).eq.1) then
            isetclu=1
            isetclc=icolr(ipl)
            setcll=float(ilinw(ipl))
            call getdash(idash(ipl),isetcld)
         elseif (majcon(i).eq.2) then
            if (.not.lnolb(ipl)) then
               isetclu=3
               isetllc=ilcll(ipl)
            else
               isetclu=1
            endif
            isetclc=icoll(ipl)
            setcll=float(ilwll(ipl))
            call getdash(idall(ipl),isetcld)
         elseif (majcon(i).eq.-1) then
            isetclu=1
            isetclc=icong(ipl)
            setcll=float(ilwng(ipl))
            call getdash(idang(ipl),isetcld)
         elseif (majcon(i).eq.-2) then
            if (.not.lnolb(ipl)) then
               isetclu=3
               isetllc=ilcnl(ipl)
            else
               isetclu=1
            endif
            isetclc=iconl(ipl)
            setcll=float(ilwnl(ipl))
            call getdash(idanl(ipl),isetcld)
         elseif (majcon(i).eq.0) then
            if (.not.lnozr(ipl)) then
               if (.not.lnolb(ipl)) then
                  isetclu=3
                  isetllc=ilczr(ipl)
               else
                  isetclu=1
               endif
               isetclc=icozr(ipl)
               setcll=float(ilwzr(ipl))
               call getdash(idazr(ipl),isetcld)
            else
               isetclu=0
            endif
         endif
c
         if (isetclc.lt.0) then  ! "transparent" color was chosen
            isetclu=0            ! Don't draw contour at all.
         endif
c
         call cpseti('CLU',isetclu)
         call cpseti('LLC',isetllc)
         call cpseti('CLC',isetclc)
         call cpsetr('CLL',setcll)
         call cpseti('CLD',isetcld)
      enddo
c
      call arinam(iam,niam)   ! initialize the area map
      call cplbam(pslab1,rwrk,iwrk,iam)  ! put label boxes in area map
      if(ipwhl(ipl) .eq. 1000) then
        call cpcldm(pslab1,rwrk,iwrk,iam,mask)   ! draw contour lines
      else
        call cpcldm(pslab1,rwrk,iwrk,iam,drawcl)   ! draw contour lines
      endif
      call cplbdr(pslab1,rwrk,iwrk)   ! make labels
c
c   Write message at bottom
c
      if (.not.lnmsg(ipl)) then
c
      messg ='                                           LOW=123'//
     &       '456789012  HIGH=123456789012  INTERVAL= 1234567890'//
     &       '12  SCALE=123456789012'
c
c     messg ='CONTOURS:  UNITS=123456789012345678901234  LOW=123'//
c    &       '456789012  HIGH=123456789012  INTERVAL= 1234567890'//
c    &       '12  SCALE=123456789012'
C             12345678901234567890123456789012345678901234567890
C                      1         2         3         4         5
      ilch=lennonblank(unwk(ipl))
      ibwk=1+(24-ilch)
      messg(ibwk:ibwk+16)='CONTOURS:  UNITS='
      write(messg(ibwk+17:ibwk+16+ilch),'(a)') unwk(ipl)(1:ilch)
      write(messg(48:59),'(g12.5)') valcon(1)
      write(messg(67:78),'(g12.5)') valcon(numcon)
      if (.not.lmult(ipl)) then
         write(messg(91:102),'(g12.5)') cintuse
      else
         write(messg(90:102),'(a1,g12.5)') 'X',cintuse
      endif
      call cpgetr('SFU',ash)
      write(messg(111:122),'(g12.5)') ash
      nchar=122
      if (ash .eq. 1.) nchar = 102
      call gqclip (ierr,iclp,rect)
      call gsclip (0)
      call gstxci(icomg(ipl))
      call gsplci(icomg(ipl))
      chsize=.008
      ypos=bottextfloor+.5*chsize
      call pcgeti ('QU',ntextqq)
      call pcseti ('QU',0)
      call plchhq (cfux(.5),cfuy(ypos),messg(ibwk:nchar),chsize,0.,0.)
      call pcseti ('QU',ntextqq)
      call gsclip (iclp)
      bottextfloor=bottextfloor+1.9*chsize
c
      endif
c
      endif
      endif  ! noplots
c
      return
      end

      subroutine mask (xwrk, ywrk, nwrk, iarea, igrp, ngrps)
      integer iarea(ngrps), igrp(ngrps)
      real xwrk(nwrk), ywrk(nwrk)

      if (nwrk .lt. 2) return
      call curved (xwrk, ywrk, nwrk)
      return
      end
