c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine scondraw(ilinw,prs,rslcg,idash,icolr,icdwk,
     &   work,ipl,rslcgprv,prssou,ipolar,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,miy,mjx,mkzh)
c
c   This routine expects a variable in deg. Celsius
c
      dimension prs(miy,mjx,mkzh),
     &   rslcg(2,maxpl),ilinw(maxpl),idash(maxpl),icolr(maxpl),
     &   work(miy,mjx,mkzh),rslcgprv(2),
     &   prssou(mkzh),icdwk(maxpl)
c
      dimension temp(200),rect(4)
      character*10 dashpat
c
      include 'comconst'
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
c   Make temperature array
c
      posx=sxgn-.5*icdwk(ipl)
      posy=sygn-.5*icdwk(ipl)
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      do 8 k=1,mkzh
         if (work(it,jl,k).eq.rmsg.or.work(it,jr,k).eq.rmsg.or.
     &       work(ib,jl,k).eq.rmsg.or.work(ib,jr,k).eq.rmsg) then
            temp(k)=rmsg
         else
            temp(k)= (
     &                (1.-ratlr)*(   ratbt)*work(it,jl,k)+
     &                (   ratlr)*(   ratbt)*work(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*work(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*work(ib,jr,k) )
         endif
c
         if (ipolar.eq.1) then
c
c         Polar mod: Add 30C to temp to account for shifted background
c         produced by modified sticdraw.
c
            temp(k)=temp(k) + 30.
         endif
c
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
c
c  Set line width and color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gqclip (ierr,iclp,rect)
      call gsclip (1)
c
c  Set dash pattern
c
      call getdash(idash(ipl),ndot)
      inew=ndot/abs(ndot)
      ndot=abs(ndot)
      if (inew.lt.0) then
         do 10 id=1,10
            ndiv=2**(id)
            ndmod=mod(ndot,ndiv)/(2**(id-1))
            if (ndmod.eq.1) then
               dashpat(11-id:11-id)='$'
            else
               dashpat(11-id:11-id)=''''
            endif
            ndot=ndot-ndmod
   10    continue
         call dashdc(dashpat,4,1 )
      else
         call dashdb(ndot)
      endif
c
c   Plot the temperature curve
c
      ido=0
      ymax=132.182-44.061*alog10(toppress)
      y=0.
      do 60 k=mkzh,1,-1
         if (temp(k).ne.rmsg) then
            ido=ido+1
            yprev=y
            y=132.182-44.061*alog10(prssou(k))
            temp1=temp(k)
            if (y.gt.ymax) then
               temp1=temp(k+1)+(ymax-yprev)/(y-yprev)*(temp1-temp(k+1))
               y=ymax
            endif
            x=0.54*temp1+0.90692*y
            if (ido.eq.1) call frstd(x,y)
            call vectd(x,y)
            if (y.eq.ymax) goto 62
         endif
   60 continue
 62   continue
      call dashdb(65535)
      call setusv('LW',1000)
      call gsplci(1)
      call gsclip (iclp)
c
      return
      end
