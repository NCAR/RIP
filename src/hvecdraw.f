c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hvecdraw(ilinw,vc3d,tmk,qvp,
     &   prs,ght,ter,sfp,sfpsm,
     &   icolr,ixwin,iywin,ismth,rvcmx,cfulb,unwk,lhide,icomg,
     &   iintv,rlevl,rlavl,cfeld,cvcor,idimn,idiffflag,
     &   work1,work2,icdwk,ilev,
     &   lnmsg,bottextfloor,pslab1,pslab2,pslabt,ipslab,ipslabt,
     &   mabpl,morpl,maxlev,maxpl,miy,mjx,mkzh,ipl,irota)
c
      dimension vc3d(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ter(miy,mjx),sfp(miy,mjx),sfpsm(miy,mjx),
     &   ixwin(2,maxpl),idimn(maxpl),icomg(maxpl),
     &   iywin(2,maxpl),ismth(maxpl),ilinw(maxpl),rvcmx(maxpl),
     &   iintv(maxpl),rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   icolr(maxpl),pslab1(mabpl,morpl),pslab2(mabpl,morpl),
     &   pslabt(mabpl,morpl),ipslab(mabpl,morpl),ipslabt(mabpl,morpl),
     &   irota(maxpl)
      logical lnmsg(maxpl),lhide(maxpl)
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1,cfulb(maxpl)*5,
     &   unwk(maxpl)*24
c
      dimension vecskip(2)
      character string*48
c
      include 'comconst'
c
      vecskip(1)=rmsg
      vecskip(2)=rmsg
c
c   Set line width
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
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
     &      prs,ght,ter,sfp,sfpsm,lhide(ipl),idiffflag,cfeld(1,ipl),
     &      work1,pslab1,mabpl,morpl,njx,niy,miy,mjx,mkzh)
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,lhide(ipl),idiffflag,cfeld(2,ipl),
     &      work2,pslab2,mabpl,morpl,njx,niy,miy,mjx,mkzh)
      endif
c
c   Smooth data if necessary
c
      call smooth(pslab1,pslabt,ismth(ipl),mabpl,njx,niy)
      call smooth(pslab2,pslabt,ismth(ipl),mabpl,njx,niy)
c
c   Put in special values where we don't want vectors
c
      vmagmax=0.
      intp=iintv(ipl)
      if (intp.gt.0) then  ! regular staggering
         intph=0
      else   ! "diamond pattern" staggering
         intp=-intp
         intph=intp/2
      endif
      do j=1,njx
      do i=1,niy
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       pslab1(j,i).ne.rmsg.and.pslab2(j,i).ne.rmsg) then
            vmag=sqrt(pslab1(j,i)**2+pslab2(j,i)**2)
            vmagmax=max(vmagmax,vmag)
         else
            pslab1(j,i)=rmsg
            pslab2(j,i)=rmsg
         endif
      enddo
      enddo
c
c   If vmagmax=0, then why the hell are we doing this?
c
      if (vmagmax.eq.0.) goto 200
c
      call getusv('XF',ixpau)
      ixpau=2**ixpau
      gskip=float(abs(iintv(ipl)))
      if (iintv(ipl).lt.0) gskip=gskip*.7071
      rpaubetv=(fr-fl)/(ur-ul)*ixpau*gskip
      if (rvcmx(ipl).gt.0) then
         rpaupervmag=rpaubetv/rvcmx(ipl)
      elseif (rvcmx(ipl).eq.0) then
         rpaupervmag=rpaubetv/vmagmax
      endif
      if (rvcmx(ipl).ge.0) then
         rpaumax=rpaupervmag*vmagmax
         npaumax=int(rpaumax)+1
         vmagpmax=vmagmax*npaumax/rpaumax
      endif
c
c   Rotate pslab1 and pslab2 arrays, and switch njx,niy, if needed.
c
      njx_sv2=njx
      njx_svi=njx
      niy_sv2=niy
      niy_svi=niy
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
c   Call velvct
c
      if (rvcmx(ipl).ge.0.) then   ! vectors
c
      imxvpl=0
      call velvctmts(pslab1,mabpl,pslab2,mabpl,njx,niy,
     &      0.,vmagpmax,1,npaumax,4,vecskip,imxvpl,icolr(ipl))
c
      if (.not.lnmsg(ipl)) then
         call setusv('LW',1000)
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
         call gsplci(icomg(ipl))
         call gstxci(icomg(ipl))
         chsize=.008
         ypos=bottextfloor+.5*chsize
         xleft=.3
         string=' '
         write(string,'(a16,f5.1,1x,a24)')
     &      'MAXIMUM VECTOR: ',vmagmax,unwk(ipl)
         call pcseti('TE',1)
         call pcgeti ('QU',ntextqq)
         call pcseti ('QU',0)
         call plchhq(xleft,ypos,string,chsize,0.,-1.)
         call pcseti ('QU',ntextqq)
         call pcgetr('DR',distr)
         call pcseti('TE',0)
         vecmaxfcor=vmagmax*rpaupervmag/ixpau
         xstart=xleft+distr+.04
         xend=xstart+vecmaxfcor
         call setusv('LW',lwidth)
         call line(xstart,ypos,xend,ypos)
         dxarrow=.24*vecmaxfcor*cos(.45)
         dyarrow=.24*vecmaxfcor*sin(.45)
         call line(xend,ypos,xend-dxarrow,ypos+dyarrow)
         call line(xend,ypos,xend-dxarrow,ypos-dyarrow)
         call gsplci(1)
         call gstxci(1)
         bottextfloor=bottextfloor+1.9*chsize
      endif
c
      else                  ! wind barbs
c
c   First make hemisphere indicator array, put in ipslab
c
      do j=1,njx_svi
         rjx=xjcorn+(.5*icdwk(ipl)+j+ixwin(1,ipl)-2.)/refrat
      do i=1,niy_svi
         riy=yicorn+(.5*icdwk(ipl)+i+iywin(1,ipl)-2.)/refrat
         call maptform(riy,rjx,rlat,rlon,1)
         if (rlat.ge.0.) then
            ipslab(j,i)=1
         else
            ipslab(j,i)=-1
         endif
      enddo
      enddo
c
c   Rotate ipslab array
c
      call rotipslab(ipslab,ipslabt,mabpl,morpl,njx_svi,niy_svi,
     &   irota(ipl))
c
c   Convert to appropriate units so that a full barb represents the
c   desired magnitude.  The barb routine (called by velbrb) always
c   assumes that a full barb = 10 units.  Hence, if cfulb=10mps, the
c   conversion factor is 1.  If cfulb=5mps, the conversion factor is
c   2.  If cfulb=10kts, the conversion factor is 1.94.
c
      string=' '
      if (index(cfulb(ipl),'10mps').ne.0) then
         barbfac=1.
         write(string,'(a41)')
     &      'BARB VECTORS:  FULL BARB = 10 m s~S~-1~N~'
         nch=41
      elseif (index(cfulb(ipl),'5mps').ne.0) then
         barbfac=2.
         write(string,'(a40)')
     &      'BARB VECTORS:  FULL BARB = 5 m s~S~-1~N~'
         nch=40
      else
         barbfac=rktpmps
         write(string,'(a33)')
     &      'BARB VECTORS:  FULL BARB = 10 kts'
         nch=33
      endif
c
      do 140 j=1,njx
      do 140 i=1,niy
         if (pslab1(j,i).ne.vecskip(1))
     &      pslab1(j,i)=barbfac*pslab1(j,i)
         if (pslab2(j,i).ne.vecskip(2))
     &      pslab2(j,i)=barbfac*pslab2(j,i)
  140 continue
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
      call velbrb(pslab1,mabpl,pslab2,mabpl,ipslab,mabpl,njx,niy,
     &    .4*gskip,4,vecskip)
      if (.not.lnmsg(ipl)) then
         call setusv('LW',1000)
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
         chsize=.008
         ypos=bottextfloor+.5*chsize
         call pcgeti ('QU',ntextqq)
         call pcseti ('QU',0)
         call gsplci(icomg(ipl))
         call gstxci(icomg(ipl))
         call plchhq(.5,ypos,string(1:nch),chsize,0.,0.)
         call pcseti ('QU',ntextqq)
         bottextfloor=bottextfloor+1.9*chsize
      endif
      call gsplci(1)
      call gstxci(1)
c
      endif
c
  200 continue
      call setusv('LW',1000)
c
      return
      end
