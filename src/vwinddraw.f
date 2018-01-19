c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vwinddraw(ilinw,rcrag,rcrbg,ismth,unwk,icomg,
     &         icolr,rvcmx,cfulb,rvvms,pavprof,ivvnx,ixavg,
     &         work1,work2,icdwk,lnmsg,nscrs,set1,set2,xdist,ydist,
     &         xseclen,cvcor,vv1,vv2,bottextfloor,unorth,vnorth,
     &         pslab1,pslab2,pslabt,ipslab,mabpl,morpl,maxpl,
     &         miy,mjx,mkzh,ipl)
c
      dimension rcrag(2,maxpl),rcrbg(2,maxpl),
     &   ismth(maxpl),rvcmx(maxpl),rvvms(maxpl),ilinw(maxpl),
     &   ivvnx(maxpl),ixavg(maxpl),icolr(maxpl),icomg(maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   pslab1(mabpl,morpl),pslab2(mabpl,morpl),
     &   pslabt(mabpl,morpl),ipslab(mabpl,morpl),
     &   unorth(miy,mjx),vnorth(miy,mjx),pavprof(1000)
      logical lnmsg(maxpl)
      character cvcor(maxpl)*1,cfulb(maxpl)*5,unwk(maxpl)*24
c
      dimension vecskip(2),idolev(100)
      character string*48
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
c   Make set call for velvct.
c
      call set(fl,fr,fb,ft,1.,float(nscrs),set1,set2,1)
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
            pslab2(ls,k)=rmsg
         else
            vele1 = vnorth(it,jl)*work1(it,jl,k)-
     &              unorth(it,jl)*work2(it,jl,k)
            vele2 = vnorth(it,jr)*work1(it,jr,k)-
     &              unorth(it,jr)*work2(it,jr,k)
            vele3 = vnorth(ib,jl)*work1(ib,jl,k)-
     &              unorth(ib,jl)*work2(ib,jl,k)
            vele4 = vnorth(ib,jr)*work1(ib,jr,k)-
     &              unorth(ib,jr)*work2(ib,jr,k)
            veln1 = unorth(it,jl)*work1(it,jl,k)+
     &              vnorth(it,jl)*work2(it,jl,k)
            veln2 = unorth(it,jr)*work1(it,jr,k)+
     &              vnorth(it,jr)*work2(it,jr,k)
            veln3 = unorth(ib,jl)*work1(ib,jl,k)+
     &              vnorth(ib,jl)*work2(ib,jl,k)
            veln4 = unorth(ib,jr)*work1(ib,jr,k)+
     &              vnorth(ib,jr)*work2(ib,jr,k)
            pslab1(ls,k)=pslab1(ls,k)+ (
     &          (1.-ratlr)*(   ratbt)*vele1+
     +          (   ratlr)*(   ratbt)*vele2+
     +          (1.-ratlr)*(1.-ratbt)*vele3+
     +          (   ratlr)*(1.-ratbt)*vele4 )*rnavg
            pslab2(ls,k)=pslab2(ls,k)+ (
     &          (1.-ratlr)*(   ratbt)*veln1+
     +          (   ratlr)*(   ratbt)*veln2+
     +          (1.-ratlr)*(1.-ratbt)*veln3+
     +          (   ratlr)*(1.-ratbt)*veln4 )*rnavg
         endif
  100 continue
c
  103 continue
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
         endif
  130 continue
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
      elseif (rvcmx(ipl).eq.0) then
         rpaupervmag=rpaubetv/vmagmax
      endif
      if (rvcmx(ipl).ge.0) then
         rpaumax=rpaupervmag*vmagmax
         npaumax=int(rpaumax)+1
         vmagpmax=vmagmax*npaumax/rpaumax
      endif
c
c   Call velvct
c
      if (rvcmx(ipl).ge.0) then   ! vectors
         imxvpl=0
         ivcs=1
         call velvctmts(pslab1,mabpl,pslab2,mabpl,nscrs,mkzh,
     &         0.,vmagpmax,1,npaumax,4,vecskip,imxvpl,icolr(ipl))
         ivcs=0
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
     &         'MAXIMUM VECTOR: ',vmagmax,unwk(ipl)
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
      do ls=1,nscrs
         rjx=xjcorn+(caxgn+(ls-1.)/(nscrs-1.)*xdist-1.)/refrat
         riy=yicorn+(caygn+(ls-1.)/(nscrs-1.)*ydist-1.)/refrat
         call maptform(riy,rjx,rlat,rlon,1)
         do k=1,mkzh
            if (rlat.ge.0.) then
               ipslab(ls,k)=1
            else
               ipslab(ls,k)=-1
            endif
         enddo
      enddo
c
c     Convert to appropriate units so that a full barb represents the
c     desired magnitude.  The barb routine (called by velbrb) always
c     assumes that a full barb = 10 units.  Hence, if cfulb=10mps, the
c     conversion factor is 1.  If cfulb=5mps, the conversion factor is
c     2.  If cfulb=10kts, the conversion factor is 1.94.
c
         string=' '
         if (index(cfulb(ipl),'10mps').ne.0) then
            barbfac=1.
            write(string,'(a41)')
     &         'BARB VECTORS:  FULL BARB = 10 m s~S~-1~N~'
            nch=41
         elseif (index(cfulb(ipl),'5mps').ne.0) then
            barbfac=2.
            write(string,'(a40)')
     &         'BARB VECTORS:  FULL BARB = 5 m s~S~-1~N~'
            nch=40
         else
            barbfac=rktpmps
            write(string,'(a33)')
     &         'BARB VECTORS:  FULL BARB = 10 kts'
            nch=33
         endif
c
         do 140 j=1,nscrs
         do 140 i=1,mkzh
            if (pslab1(j,i).ne.vecskip(1))
     &         pslab1(j,i)=barbfac*pslab1(j,i)
            if (pslab2(j,i).ne.vecskip(2))
     &         pslab2(j,i)=barbfac*pslab2(j,i)
  140    continue
         call gsplci(icolr(ipl))
         call gstxci(icolr(ipl))
         ivcs=1
         call velbrb(pslab1,mabpl,pslab2,mabpl,ipslab,mabpl,nscrs,mkzh,
     &            .4*gskip,4,vecskip)
         ivcs=0
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
      endif
c
  200 continue
      call setusv('LW',1000)
c
      return
      end
