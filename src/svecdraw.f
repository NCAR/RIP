c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine svecdraw(ilinw,prs,rslcg,icolr,icdwk,
     &   work1,work2,ipl,rslcgprv,unorth,vnorth,rvvms,pavprof,prssou,
     &   ipolar,
     &   cfulb,lhodo,ins,flminsou,frmaxsou,fbminsou,ftmaxsou,icomax,
     &   maxpl,miy,mjx,mkzh,iwkidcgm)
c
      dimension prs(miy,mjx,mkzh),
     &   rslcg(2,maxpl),ilinw(maxpl),icolr(maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),rslcgprv(2),
     &   unorth(miy,mjx),vnorth(miy,mjx),rvvms(maxpl),
     &   prssou(mkzh),icdwk(maxpl),pavprof(1000)
      logical lhodo(maxpl)
      character cfulb(maxpl)*5
c
      dimension velnorth(200),veleast(200),rect(4)
      character string*16
c
      include 'comconst'
c
c   Convert to appropriate units so that a full barb represents
c      the desired magnitude.   The barb routine always assumes
c      that a full barb = 10 units.  Hence, if cfulb=10mps, the
c      conversion factor is 1.  If cfulb=5mps, the conversion factor
c      is 2.  If cfulb=10kts, the conversion factor is 1.94.
c
      if (index(cfulb(ipl),'10mps').ne.0) then
         barbfac=1.
         string='10 m s~S~-1~N~  '
         nch=14
      elseif (index(cfulb(ipl),'5mps').ne.0) then
         barbfac=2.
         string='5 m s~S~-1~N~   '
         nch=13
      else
         barbfac=rktpmps
         string='10 kts          '
         nch=6
      endif
c
c   Interpolate data to sounding location
c
      sxgn=1.+(rslcg(2,ipl)-xjcorn)*refrat
      sygn=1.+(rslcg(1,ipl)-yicorn)*refrat
      if (sxgn.le..5.or.sxgn.ge.mjx-.5.or.
     &    sygn.le..5.or.sygn.ge.miy-.5) then
         write(iup,*)'I don''t do soundings outside the'
         write(iup,*)'cross-point domain.'
         stop
      endif
c
c   Make prssou if not already made for this location
c
      if (rslcg(1,ipl).ne.rslcgprv(1).or.
     &    rslcg(2,ipl).ne.rslcgprv(2)) then
         posx=sxgn-.5
         posy=sygn-.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         do 5 k=1,mkzh
            prssou(k)= (    (1.-ratlr)*(   ratbt)*prs(it,jl,k)+
     &                      (   ratlr)*(   ratbt)*prs(it,jr,k)+
     &                      (1.-ratlr)*(1.-ratbt)*prs(ib,jl,k)+
     &                      (   ratlr)*(1.-ratbt)*prs(ib,jr,k) )
    5    continue
      endif
c
c   Make velocity components
c
      posx=sxgn-.5*icdwk(ipl)
      posy=sygn-.5*icdwk(ipl)
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      unorths= (
     &                (1.-ratlr)*(   ratbt)*unorth(it,jl)+
     &                (   ratlr)*(   ratbt)*unorth(it,jr)+
     &                (1.-ratlr)*(1.-ratbt)*unorth(ib,jl)+
     &                (   ratlr)*(1.-ratbt)*unorth(ib,jr) )
      vnorths= (
     &                (1.-ratlr)*(   ratbt)*vnorth(it,jl)+
     &                (   ratlr)*(   ratbt)*vnorth(it,jr)+
     &                (1.-ratlr)*(1.-ratbt)*vnorth(ib,jl)+
     &                (   ratlr)*(1.-ratbt)*vnorth(ib,jr) )
      do 8 k=1,mkzh
         uuus = (
     &                (1.-ratlr)*(   ratbt)*work1(it,jl,k)+
     &                (   ratlr)*(   ratbt)*work1(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*work1(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*work1(ib,jr,k) )
         vvvs = (
     &                (1.-ratlr)*(   ratbt)*work2(it,jl,k)+
     &                (   ratlr)*(   ratbt)*work2(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*work2(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*work2(ib,jr,k) )
         velnorth(k)=(unorths*uuus+vnorths*vvvs)*barbfac
         veleast(k)= (vnorths*uuus-unorths*vvvs)*barbfac
    8 continue
c
c   Plot sounding
c
c   Note: "toppress" can be changed to stretch the vertical scale
c   of the sounding, but labels and upper border of sounding plot are
c   not changed accordingly, so it's not too pretty.  Background is
c   designed for toppress=100 hPa.
c
      toppress=100.
      toplim=132.182-44.061*alog10(toppress)
      space=3.
      spacev=3.*(toplim+.9346)*.02222
      call set(flminsou,frmaxsou,fbminsou,ftmaxsou,
     &   -19.0-space,27.1+space,-.9346217-.2*spacev,toplim+spacev,1)
c    &   -19.0-space,27.1+space,-.9346217-.2*space,44.061+space,1)
      call gqclip (ierr,iclp,rect)
      call gsclip (1)
c
      if (ipolar.eq.0) then
         xm=24.2
      elseif (ipolar.eq.1) then
         xm=26.2
      endif
      call line(xm,-.9346217,xm,44.061)
c
c  Set line width and color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
c
c   Plot wind vectors
c
      yval=132.182-44.061*alog10(prssou(mkzh))
      call barb(xm,yval,veleast(mkzh),velnorth(mkzh),.02,ins)
      dprsmin=rvvms(ipl)
      dprs=0.
      do 75 k=mkzh-1,1,-1
         dprs=dprs+pavprof(k+1)-pavprof(k)
         if (dprs.gt.dprsmin .and. prssou(k) .ge. toppress) then
            yval=132.182-44.061*alog10(prssou(k))
            call barb(xm,yval,veleast(k),velnorth(k),.02,ins)
            dprs=0.
         endif
   75 continue
      if (lhodo(ipl)) call hodograph (veleast, velnorth, prssou, mkzh,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,icomax,barbfac,iwkidcgm)
c
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      chsize=.008
      call pcgeti ('QU',ntextqq)
      call pcseti ('QU',0)
      ypos=.878
      xpos=.9
      call setusv('LW',1000)
      call plchhq(xpos,ypos,'Full barb:',chsize,0.,0.)
      ypos=ypos-.013
      call plchhq(xpos,ypos,string(1:nch),chsize,0.,0.)
      call pcseti ('QU',ntextqq)
c
      call setusv('LW',1000)
      call gsplci(1)
      call gsclip (iclp)
c
      return
      end
