c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vvecdraw(ilinw,vervl,rcrag,rcrbg,ismth,
     &         icolr,rvcmx,rvvms,pavprof,ivvnx,ixavg,work1,work2,
     &         icdwk,icomg,lnmsg,nscrs,set1,set2,xdist,ydist,xseclen,
     &         cvcor,cfeld,vv1,vv2,bottextfloor,pslab1,pslab2,pslabt,
     &         mabpl,morpl,maxpl,miy,mjx,mkzh,ipl)
c
      dimension vervl(miy,mjx,mkzh),rcrag(2,maxpl),rcrbg(2,maxpl),
     &   ismth(maxpl),rvcmx(maxpl),rvvms(maxpl),ilinw(maxpl),
     &   ivvnx(maxpl),ixavg(maxpl),icolr(maxpl),icomg(maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   pslab1(mabpl,morpl),pslab2(mabpl,morpl),pslabt(mabpl,morpl),
     &   pavprof(1000)
      logical lnmsg(maxpl)
      character cvcor(maxpl)*1,cfeld(3,maxpl)*10
c
      dimension vecskip(2),idolev(100)
      character string*96
c
      include 'comconst'
      include 'comvctran'
c
      vecskip(1)=rmsg
      vecskip(2)=rmsg
c
c   Set line width
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
c
      fb=fbmin
      ft=ftmax
      fl=flmin
      fr=frmax
c
c   Make appropriate mutiplying factors for work array
c
      fac1 = xdist/xseclen
      fac2 = ydist/xseclen
c
c   Make set call for velvct.
c
      call set(fl,fr,fb,ft,1.,float(nscrs),set1,set2,1)
      cfac=(ft-fb)*xseclen*ds/((fr-fl)*(vv2-vv1))
      if (cvcor(ipl).eq.'z') cfac=cfac*.00001
      if (cvcor(ipl).eq.'p') cfac=cfac*.001
c
c   If field is Sawyer-Eliassen streamfunction, data is already
c   arranged into cross-section array.
c
      if (cfeld(2,ipl)(1:2).eq.'se'.or.cfeld(2,ipl)(1:2).eq.'sm') then
         do k=1,mkzh
         do ls=1,nscrs
            n1d=(k-1)*nscrs+ls
            kk=1+(n1d-1)/(miy*mjx)
            n1dleft=n1d-(kk-1)*miy*mjx
            jj=1+(n1dleft-1)/miy
            ii=n1dleft-(jj-1)*miy
            pslab1(ls,k)=-work2(ii,jj,kk)
            pslab2(ls,k)=cfac*vervl(ii,jj,kk)
         enddo
         enddo
         goto 431
      endif
c
c   Set up x-sec averaging parameters
c
      cosangle=xdist/xseclen
      sinangle=ydist/xseclen
      rnavg=1./(2.*ixavg(ipl)+1.)
      do 50 k=1,mkzh
      do 50 ls=1,nscrs
         pslab1(ls,k)=0.
         pslab2(ls,k)=0.
   50 continue
c
c   Interpolate gridded data to x-section.
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
      do 103 islab=-ixavg(ipl),ixavg(ipl)
         xj1t=caxgn+islab*sinangle
         xj2t=cbxgn+islab*sinangle
         yi1t=caygn-islab*cosangle
         yi2t=cbygn-islab*cosangle
         if (xj1t.le.1.5.or.xj1t.ge.mjx-.5.or.
     &       xj2t.le.1.5.or.xj2t.ge.mjx-.5.or.
     &       yi1t.le.1.5.or.yi1t.ge.miy-.5.or.
     &       yi2t.le.1.5.or.yi2t.ge.miy-.5) then
            write(iup,*)'Cross sec. endpoints must be between 1.5'
            write(iup,*)'and (miy-.5) or (mjx-.5).'
            stop
         endif
c
      do 100 k=1,mkzh
      do 100 ls=1,nscrs
         posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t) -
     &        icdwk(ipl)*.5
         posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t) -
     &        icdwk(ipl)*.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         if (work1(it,jl,k).eq.rmsg.or.
     &       work1(it,jr,k).eq.rmsg.or.
     &       work1(ib,jl,k).eq.rmsg.or.
     &       work1(ib,jr,k).eq.rmsg.or.
     &       work2(it,jl,k).eq.rmsg.or.
     &       work2(it,jr,k).eq.rmsg.or.
     &       work2(ib,jl,k).eq.rmsg.or.
     &       work2(ib,jr,k).eq.rmsg.or.
     &       pslab1(ls,k).eq.rmsg) then
            pslab1(ls,k)=rmsg
         else
            wk1=fac1*work1(it,jl,k)+
     &         fac2*work2(it,jl,k)
            wk2=fac1*work1(it,jr,k)+
     &         fac2*work2(it,jr,k)
            wk3=fac1*work1(ib,jl,k)+
     &         fac2*work2(ib,jl,k)
            wk4=fac1*work1(ib,jr,k)+
     &         fac2*work2(ib,jr,k)
            pslab1(ls,k)=pslab1(ls,k)+(
     &                   (1.-ratlr)*(   ratbt)*wk1+
     &                   (   ratlr)*(   ratbt)*wk2+
     &                   (1.-ratlr)*(1.-ratbt)*wk3+
     &                   (   ratlr)*(1.-ratbt)*wk4 )*rnavg
         endif
         posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t) - .5
         posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t) - .5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         pslab2(ls,k)=pslab2(ls,k)+cfac*(
     &                   (1.-ratlr)*(   ratbt)*vervl(it,jl,k)+
     &                   (   ratlr)*(   ratbt)*vervl(it,jr,k)+
     &                   (1.-ratlr)*(1.-ratbt)*vervl(ib,jl,k)+
     &                   (   ratlr)*(1.-ratbt)*vervl(ib,jr,k)  )*rnavg
  100 continue
c
  103 continue
c
 431  continue
c
c   Smooth data if necessary
c
      call smooth(pslab1,pslabt,ismth(ipl),mabpl,nscrs,mkzh)
      call smooth(pslab2,pslabt,ismth(ipl),mabpl,nscrs,mkzh)
c
c   Put in special values where we don't want vectors
c
      dprsmin=rvvms(ipl)
      idolev(mkzh)=1
      dprs=0.
      do 125 k=mkzh-1,1,-1
         dprs=dprs+pavprof(k+1)-pavprof(k)
         if (dprs.gt.dprsmin) then
            dprs=0.
            idolev(k)=1
         else
            idolev(k)=0
         endif
  125 continue
      interval=max(1,nint(float(nscrs)/ivvnx(ipl)))
      vmagmax=0.
      hormax=0.
      vermax=0.
      setmin=min(set1,set2)
      setmax=max(set1,set2)
      do 130 k=1,mkzh
      do 130 ls=1,nscrs
         if (idolev(k).ne.1.or.mod(ls-1,interval).ne.0.or.
     &       vc2d(ls,k).lt.setmin.or.vc2d(ls,k).gt.setmax) then
            pslab1(ls,k)=vecskip(1)
            pslab2(ls,k)=vecskip(2)
         elseif (pslab1(ls,k).ne.rmsg.and.pslab2(ls,k).ne.rmsg) then
            vmag=sqrt(pslab1(ls,k)**2+pslab2(ls,k)**2)
            vmagmax=max(vmagmax,vmag)
            hormax=max(hormax,abs(pslab1(ls,k)))
            vermax=max(vermax,abs(pslab2(ls,k)))
         endif
  130 continue
      vermaxu=abs(vermax/cfac)
c
c   If vmagmax=0, then why the hell are we doing this?
c
      if (vmagmax.eq.0.) goto 200
c
      call getusv('XF',ixpau)
      ixpau=2**ixpau
      gskip=float(interval)
      rpaubetv=(fr-fl)/(nscrs-1.)*ixpau*gskip
      if (rvcmx(ipl).gt.0.) then
         rpaupervmag=rpaubetv/rvcmx(ipl)
      else
         rpaupervmag=rpaubetv/vmagmax
      endif
      rpaumax=rpaupervmag*vmagmax
      npaumax=int(rpaumax)+1
      vmagpmax=vmagmax*npaumax/rpaumax
c
c   Call velvct
c
      imxvpl=0
      ivcs=1
      call velvctmts(pslab1,mabpl,pslab2,mabpl,nscrs,mkzh,
     &            0.,vmagpmax,1,npaumax,4,vecskip,imxvpl,icolr(ipl))
      ivcs=0
c
  200 continue
      call setusv('LW',1000)
c
c   Print max vector components at bottom of frame.
c
      if (.not.lnmsg(ipl)) then
         call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
         call gsplci(icomg(ipl))
         call gstxci(icomg(ipl))
         chsize=.008
         ypos=bottextfloor+.5*chsize
         string=' '
         if     (cvcor(ipl).eq.'s') then
            write(string,'(a16,f5.1,a27,f8.7,a17)')
     &         'MAXIMUM VECTOR: ',hormax,' m s~S~-1~N~ (HORIZ)       ',
     &         vermaxu,' s~S~-1~N~ (VERT)'
            nch=73
         elseif (cvcor(ipl).eq.'p') then
            write(string,'(a16,f5.1,a27,f6.1,a21)')
     &         'MAXIMUM VECTOR: ',hormax,' m s~S~-1~N~ (HORIZ)       ',
     &         vermaxu,' dPa s~S~-1~N~ (VERT)'
            nch=82
         elseif (cvcor(ipl).eq.'z') then
            write(string,'(a16,f5.1,a27,f6.1,a20)')
     &         'MAXIMUM VECTOR: ',hormax,' m s~S~-1~N~ (HORIZ)       ',
     &         vermaxu,' cm s~S~-1~N~ (VERT)'
            nch=74
         endif
         call pcgeti ('QU',ntextqq)
         call pcseti ('QU',0)
         call plchhq(.5,ypos,string(1:nch),chsize,0.,0.)
         call pcseti ('QU',ntextqq)
         vecmaxfcor=hormax*rpaupervmag/ixpau
         xstart=.98-.5*ypos
         xend=xstart-vecmaxfcor
         call setusv('LW',lwidth)
         call line(xstart,ypos,xend,ypos)
         dxarrow=.24*vecmaxfcor*cos(.45)
         dyarrow=.24*vecmaxfcor*sin(.45)
         call line(xend,ypos,xend+dxarrow,ypos+dyarrow)
         call line(xend,ypos,xend+dxarrow,ypos-dyarrow)
         vecmaxfcor=vermax*rpaupervmag/ixpau
         xpos=xstart
         ystart=ypos
         yend=ystart+vecmaxfcor
         call line(xpos,ystart,xpos,yend)
         dxarrow=.24*vecmaxfcor*sin(.45)
         dyarrow=.24*vecmaxfcor*cos(.45)
         call line(xpos,yend,xpos-dxarrow,yend-dyarrow)
         call line(xpos,yend,xpos+dxarrow,yend-dyarrow)
         call gsplci(1)
         call gstxci(1)
         call gsplci(1)
         call gstxci(1)
         bottextfloor=bottextfloor+1.9*chsize
      endif
c
      call setusv('LW',1000)
c
      return
      end
