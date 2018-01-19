c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hstrdraw(ilinw,iintv,workspc,vc3d,tmk,qvp,
     &   prs,ght,ter,sfp,sfpsm,
     &   ixwin,iywin,ismth,icolr,lhide,
     &   rlevl,rlavl,cfeld,cvcor,idimn,work1,work2,icdwk,ilev,
     &   idiffflag,pslab1,pslab2,pslabt,mabpl,morpl,maxlev,maxpl,
     &   miy,mjx,mkzh,ipl,irota)
c
      dimension workspc(miy,mjx,mkzh),vc3d(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),
     &   ter(miy,mjx),sfp(miy,mjx),sfpsm(miy,mjx),
     &   prs(miy,mjx,mkzh),icolr(maxpl),idimn(maxpl),
     &   ixwin(2,maxpl),iywin(2,maxpl),ismth(maxpl),ilinw(maxpl),
     &   iintv(maxpl),irota(maxpl),
     &   rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   pslab1(mabpl,morpl),pslab2(mabpl,morpl),pslabt(mabpl,morpl)
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1
      logical lhide(maxpl)
c
      include 'comconst'
      common /str03/  inita , initb , arowl , iterp , iterc , igflg
     +             ,  imsg , uvmsg , icyc , displ , dispc , cstop
      common / stpar /
     +                iud1       ,ivd1       ,ipd1       ,
     +                ixd1       ,ixdm       ,iyd1       ,iydn       ,
     +                ixm1       ,iym1       ,ixm2       ,iym2       ,
     +                iwkd       ,iwku       ,iset       ,ierr       ,
     +                ixin       ,iyin       ,imsk       ,icpm       ,
     +                nlvl       ,ipai       ,ictv       ,wdlv       ,
     +                uvmn       ,uvmx       ,pmin       ,pmax       ,
     +                ithn       ,iplr       ,isst       ,
     +                iclr(64)           ,tvlu(64)
c
      imsg=1            ! needed for strmln common block
      uvmsg=rmsg        ! needed for strmln common block
c
c   Set color and common block variable that controls line width,
c      and common block variable that controls density.
c
      call gsplci(icolr(ipl))
      wdlv=float(ilinw(ipl))
      inita=iintv(ipl)
      initb=inita
c  make the arrow size a function of the number of x gridpoints
      arowl = mjx / 50.
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
c
c   Put appropriate data into horizontal slab.
c
      if ((cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).eq.
     &   rlavl(ilev,ipl)).or.idimn(ipl).eq.2) then
         do 90 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 90 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work1(ii,jj,nint(rlevl(ilev,ipl)))
               pslab2(j,i)=work2(ii,jj,nint(rlevl(ilev,ipl)))
   90    continue
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).ge.0) then
         call fillarray(pslab1,mabpl*morpl,0.)
         call fillarray(pslab2,mabpl*morpl,0.)
         lev1=min(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         lev2=max(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         do j=1,njx
            jj=j+ixwin(1,ipl)-1
         do i=1,niy
            ii=i+iywin(1,ipl)-1
            prstot=prs(ii,jj,lev2)-prs(ii,jj,lev1)
            do k=lev1,lev2-1
               pslab1(j,i)=pslab1(j,i)+.5*(work1(ii,jj,k)+
     &            work1(ii,jj,k+1))*(prs(ii,jj,k+1)-prs(ii,jj,k))/prstot
               pslab2(j,i)=pslab1(j,i)+.5*(work2(ii,jj,k)+
     &            work2(ii,jj,k+1))*(prs(ii,jj,k+1)-prs(ii,jj,k))/prstot
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
               pslab1(j,i)=work1(ii,jj,lev1)-work1(ii,jj,lev2)
               pslab2(j,i)=work2(ii,jj,lev1)-work2(ii,jj,lev2)
            enddo
         enddo
      else
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,lhidel,idiffflag,cfeld(1,ipl),work1,
     &      pslab1,mabpl,morpl,njx,niy,miy,mjx,mkzh)
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,lhidel,idiffflag,cfeld(2,ipl),work2,
     &      pslab2,mabpl,morpl,njx,niy,miy,mjx,mkzh)
      endif
c
c   Smooth data if necessary
c
      call smooth(pslab1,pslabt,ismth(ipl),mabpl,njx,niy)
      call smooth(pslab2,pslabt,ismth(ipl),mabpl,njx,niy)
c
c   Rotate pslab1 and pslab2 arrays, and switch njx,niy, if needed.
c
      njx_sv2=njx
      niy_sv2=niy
      call rotpslab(pslab1,pslabt,mabpl,morpl,njx,niy,irota(ipl))
      call rotpslab(pslab2,pslabt,mabpl,morpl,njx_sv2,niy_sv2,
     &   irota(ipl))
c
c     Not only do the data have to be rearranged within each array, but
c     the wind components also have to be rotated.
c
      if (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg.and.pslab2(j,i).ne.rmsg) then
               pslab1(j,i)=-pslab1(j,i)
               pslab2(j,i)=-pslab2(j,i)
            endif
         enddo
         enddo
      elseif (irota(ipl).eq.90) then
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg.and.pslab2(j,i).ne.rmsg) then
               tempo=pslab1(j,i)
               pslab1(j,i)=-pslab2(j,i)
               pslab2(j,i)= tempo
            endif
         enddo
         enddo
      elseif (irota(ipl).eq.-90) then
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg.and.pslab2(j,i).ne.rmsg) then
               tempo=pslab1(j,i)
               pslab1(j,i)= pslab2(j,i)
               pslab2(j,i)=-tempo
            endif
         enddo
         enddo
      endif
c
c   Call strmln
c
      call strmln(pslab1,pslab2,workspc,mabpl,njx,niy,1,istrerr)
      call setusv('LW',1000)
      call gsplci(1)
c
      return
      end
