c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine maptick(refrat,yicorn,xjcorn,iywin,ixwin,llint,iside,
     &   rtslb,irota,iup)
c
c   This program labels latitude and longitude lines where they
c   intersect the boundaries of the map.  For each side of the domain,
c   whichever type of line (lat or lon) has more intersections with
c   that boundary is the one that gets labeled.
c   This routine should be called immediately after the map is drawn,
c   since it uses the "set" parameters specified by the map drawing
c   routine.  This routine calls the routine maptform, which converts
c   between x/y coordinates in the model's centered coarse domain and
c   lat/lon.
c
c   Input:
c      refrat: refinement ratio between this nest and the coarse
c         (i.e. centered) domain
c      yicorn,xjcorn: position in coarse (centered) domain of the
c         lower left corner
c      iywin: 2-element array containing min and max limits
c         of the window in the y-dir. (in nest domain grid values).
c      ixwin: 2-element array containing min and max limits
c         of the window in the x-dir. (in nest domain grid values).
c      llint: interval (in degrees) for lat/lon lines.  If zero,
c         no labels are drawn.  Must be an integer.
c      iside: value of 2 means labels are drawn only on top and right;
c         value of -2 means labels are drawn only on bottom and left;
c         any other value means labels are drawn on all 4 sides.
c      rtslb: text size (in fractional coordinates) of labels.
c
      parameter (nres=200)
c   
      dimension iywin(2),ixwin(2),xx(nres),yy(nres),
     &   rlatbot(nres),rlattop(nres),rlatleft(nres),rlatright(nres),
     &   rlonbot(nres),rlontop(nres),rlonleft(nres),rlonright(nres)
c
      logical dolatbot,dolattop,dolatleft,dolatright,
     &        dolonbot,dolontop,dolonleft,dolonright
c
      dimension rect(4)
      character axlab*5
c
      call gqclip (ierr,iclp,rect)
      call gsclip (0)  ! turn off clipping for now
      call getusv ('LW',lwsv)
      call setusv ('LW',1000)
c
c     Determine what the corners of the unrotated map are in terms of
c     coarse domain grid points.
c
      yllc=yicorn+(iywin(1)-1.)/refrat
      xllc=xjcorn+(ixwin(1)-1.)/refrat
      yurc=yicorn+(iywin(2)-1.)/refrat
      xurc=xjcorn+(ixwin(2)-1.)/refrat
c
c     The above set parameters can remain the normal values if rotation
c     is 0 or +/- 180 degrees, but x and y values must be switched with
c     each other if rotation is +/- 90 degrees.
c
      if (irota.eq.90.or.irota.eq.-90) then
         temp=yllc
         yllc=xllc
         xllc=temp
         temp=yurc
         yurc=xurc
         xurc=temp
      endif
c
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
c
c     Note, in getset above, fl,fr,fb,ft apply to rotated view, assuming
c     mapset was called for rotated view.  This is what we want.  When
c     we call set in the next line, the view will be the correct
c     orientation for the rotated map, but the user coordinate system
c     will increase rightward and upward, regardless of rotation angle.
c     This must be carefully accounted for prior to calls to maptform,
c     in order to get the proper lat/lon values for the rotated map.
c
      call set(fl,fr,fb,ft,xllc,xurc,yllc,yurc,ll)
c
      rat=(xurc-xllc)/(fr-fl)
      ydist=(.6*rtslb+.01)*rat
      xdist=.01*rat
c
      ddy=(yurc-yllc)/(nres-1.)
      ddx=(xurc-xllc)/(nres-1.)
c
      dabslatbot=0.
      dabslattop=0.
      dabslatleft=0.
      dabslatright=0.
      dabslonbot=0.
      dabslontop=0.
      dabslonleft=0.
      dabslonright=0.
c
      ladd=(1000/llint)*llint
      laddlonbot=ladd
      laddlontop=ladd
      laddlonleft=ladd
      laddlonright=ladd
c
      iswlonbot=0
      iswlontop=0
      iswlonleft=0
      iswlonright=0
c
      do i=1,nres
         yy(i)=yllc+(i-1.)*ddy
         xx(i)=xllc+(i-1.)*ddx
c
c      Account for rotated view when feeding x/y values to maptform.
c
         if (irota.eq.0) then
            yyleft=yy(i)
            yyright=yyleft
            yybot=yllc
            yytop=yurc
            xxleft=xllc
            xxright=xurc
            xxbot=xx(i)
            xxtop=xxbot
         elseif (irota.eq.90) then
            yyleft=xurc
            yyright=xllc
            yybot=xurc-(i-1.)*ddx
            yytop=yybot
            xxleft=yy(i)
            xxright=xxleft
            xxbot=yllc
            xxtop=yurc
         elseif (irota.eq.-90) then
            yyleft=xllc
            yyright=xurc
            yybot=xx(i)
            yytop=yybot
            xxleft=yurc-(i-1.)*ddy
            xxright=xxleft
            xxbot=yurc
            xxtop=yllc
         elseif (irota.eq.180.or.irota.eq.-180) then
            yyleft=yurc-(i-1.)*ddy
            yyright=yyleft
            yybot=yurc
            yytop=yllc
            xxleft=xurc
            xxright=xllc
            xxbot=xurc-(i-1.)*ddx
            xxtop=xxbot
         endif
         call maptform(yyleft,xxleft,rlatleft(i),rlonleft(i),1)
         call maptform(yyright,xxright,rlatright(i),rlonright(i),1)
         call maptform(yybot,xxbot,rlatbot(i),rlonbot(i),1)
         call maptform(yytop,xxtop,rlattop(i),rlontop(i),1)
         rlatleft(i)=rlatleft(i)+ladd
         rlonleft(i)=rlonleft(i)+laddlonleft
         rlatright(i)=rlatright(i)+ladd
         rlonright(i)=rlonright(i)+laddlonright
         rlatbot(i)=rlatbot(i)+ladd
         rlonbot(i)=rlonbot(i)+laddlonbot
         rlattop(i)=rlattop(i)+ladd
         rlontop(i)=rlontop(i)+laddlontop
         if (i.gt.1) then
            if (rlonleft(i)-rlonleft(i-1).gt.180.) then
               rlonleft(i)=rlonleft(i)-360.
               laddlonleft=laddlonleft-360
               iswlonleft=iswlonleft+1
            elseif (rlonleft(i)-rlonleft(i-1).lt.-180.) then
               rlonleft(i)=rlonleft(i)+360.
               laddlonleft=laddlonleft+360
               iswlonleft=iswlonleft+1
            endif
            if (iswlonleft.gt.1) then
               write(iup,*)'Too many sign switches for lonleft.'
               stop
            endif
            if (rlonright(i)-rlonright(i-1).gt.180.) then
               rlonright(i)=rlonright(i)-360.
               laddlonright=laddlonright-360
               iswlonright=iswlonright+1
            elseif (rlonright(i)-rlonright(i-1).lt.-180.) then
               rlonright(i)=rlonright(i)+360.
               laddlonright=laddlonright+360
               iswlonright=iswlonright+1
            endif
            if (iswlonright.gt.1) then
               write(iup,*)'Too many sign switches for lonright.'
               stop
            endif
            if (rlonbot(i)-rlonbot(i-1).gt.180.) then
               rlonbot(i)=rlonbot(i)-360.
               laddlonbot=laddlonbot-360
               iswlonbot=iswlonbot+1
            elseif (rlonbot(i)-rlonbot(i-1).lt.-180.) then
               rlonbot(i)=rlonbot(i)+360.
               laddlonbot=laddlonbot+360
               iswlonbot=iswlonbot+1
            endif
            if (iswlonbot.gt.1) then
               write(iup,*)'Too many sign switches for lonbot.'
               stop
            endif
            if (rlontop(i)-rlontop(i-1).gt.180.) then
               rlontop(i)=rlontop(i)-360.
               laddlontop=laddlontop-360
               iswlontop=iswlontop+1
            elseif (rlontop(i)-rlontop(i-1).lt.-180.) then
               rlontop(i)=rlontop(i)+360.
               laddlontop=laddlontop+360
               iswlontop=iswlontop+1
            endif
            if (iswlontop.gt.1) then
               write(iup,*)'Too many sign switches for lontop.'
               stop
            endif
            dabslatbot=dabslatbot+abs(rlatbot(i)-rlatbot(i-1))
            dabslattop=dabslattop+abs(rlattop(i)-rlattop(i-1))
            dabslatleft=dabslatleft+abs(rlatleft(i)-rlatleft(i-1))
            dabslatright=dabslatright+abs(rlatright(i)-rlatright(i-1))
            dabslonbot=dabslonbot+abs(rlonbot(i)-rlonbot(i-1))
            dabslontop=dabslontop+abs(rlontop(i)-rlontop(i-1))
            dabslonleft=dabslonleft+abs(rlonleft(i)-rlonleft(i-1))
            dabslonright=dabslonright+abs(rlonright(i)-rlonright(i-1))
         endif
      enddo
c
      if (dabslatbot.gt.dabslonbot) then
         dolatbot=.true.
         dolonbot=.false.
      else
         dolatbot=.false.
         dolonbot=.true.
      endif
c
      if (dabslattop.gt.dabslontop) then
         dolattop=.true.
         dolontop=.false.
      else
         dolattop=.false.
         dolontop=.true.
      endif
c
      if (dabslatleft.gt.dabslonleft) then
         dolatleft=.true.
         dolonleft=.false.
      else
         dolatleft=.false.
         dolonleft=.true.
      endif
c
      if (dabslatright.gt.dabslonright) then
         dolatright=.true.
         dolonright=.false.
      else
         dolatright=.false.
         dolonright=.true.
      endif
c
      if (iside.eq.2) then
         dolatbot=.false.
         dolonbot=.false.
         dolatleft=.false.
         dolonleft=.false.
      elseif (iside.eq.-2) then
         dolattop=.false.
         dolontop=.false.
         dolatright=.false.
         dolonright=.false.
      endif
c
      do i=1,nres-1
c
         if (dolatbot) then
            if (rlatbot(i+1).ge.rlatbot(i)) then
               l1=int(rlatbot(i)/llint)+1
               l2=int(rlatbot(i+1)/llint)
               idr=1
            else
               l1=int(rlatbot(i)/llint)
               l2=int(rlatbot(i+1)/llint)+1
               idr=-1
            endif
            ypos=yllc-ydist
            do l=l1,l2,idr
               lval=l*llint
               xpos=xx(i)+(lval-rlatbot(i))/
     &            (rlatbot(i+1)-rlatbot(i))*ddx
               lval=lval-ladd
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' S'
                  iaxe=3
               elseif (lval.le.-10) then
                  write(axlab(1:4),'(i2,a2)') -lval,' S'
                  iaxe=4
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' N'
                  iaxe=3
               elseif (lval.ge.10) then
                  write(axlab(1:4),'(i2,a2)') lval,' N'
                  iaxe=4
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,0.)
            enddo
         endif
         if (dolattop) then
            if (rlattop(i+1).ge.rlattop(i)) then
               l1=int(rlattop(i)/llint)+1
               l2=int(rlattop(i+1)/llint)
               idr=1
            else
               l1=int(rlattop(i)/llint)
               l2=int(rlattop(i+1)/llint)+1
               idr=-1
            endif
            ypos=yurc+ydist
            do l=l1,l2,idr
               lval=l*llint
               xpos=xx(i)+(lval-rlattop(i))/
     &            (rlattop(i+1)-rlattop(i))*ddx
               lval=lval-ladd
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' S'
                  iaxe=3
               elseif (lval.le.-10) then
                  write(axlab(1:4),'(i2,a2)') -lval,' S'
                  iaxe=4
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' N'
                  iaxe=3
               elseif (lval.ge.10) then
                  write(axlab(1:4),'(i2,a2)') lval,' N'
                  iaxe=4
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,0.)
            enddo
         endif
         if (dolatleft) then
            if (rlatleft(i+1).ge.rlatleft(i)) then
               l1=int(rlatleft(i)/llint)+1
               l2=int(rlatleft(i+1)/llint)
               idr=1
            else
               l1=int(rlatleft(i)/llint)
               l2=int(rlatleft(i+1)/llint)+1
               idr=-1
            endif
            xpos=xllc-xdist
            do l=l1,l2,idr
               lval=l*llint
               ypos=yy(i)+(lval-rlatleft(i))/
     &            (rlatleft(i+1)-rlatleft(i))*ddy
               lval=lval-ladd
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' S'
                  iaxe=3
               elseif (lval.le.-10) then
                  write(axlab(1:4),'(i2,a2)') -lval,' S'
                  iaxe=4
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' N'
                  iaxe=3
               elseif (lval.ge.10) then
                  write(axlab(1:4),'(i2,a2)') lval,' N'
                  iaxe=4
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,1.)
            enddo
         endif
         if (dolatright) then
            if (rlatright(i+1).ge.rlatright(i)) then
               l1=int(rlatright(i)/llint)+1
               l2=int(rlatright(i+1)/llint)
               idr=1
            else
               l1=int(rlatright(i)/llint)
               l2=int(rlatright(i+1)/llint)+1
               idr=-1
            endif
            xpos=xurc+xdist
            do l=l1,l2,idr
               lval=l*llint
               ypos=yy(i)+(lval-rlatright(i))/
     &            (rlatright(i+1)-rlatright(i))*ddx
               lval=lval-ladd
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' S'
                  iaxe=3
               elseif (lval.le.-10) then
                  write(axlab(1:4),'(i2,a2)') -lval,' S'
                  iaxe=4
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' N'
                  iaxe=3
               elseif (lval.ge.10) then
                  write(axlab(1:4),'(i2,a2)') lval,' N'
                  iaxe=4
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,-1.)
            enddo
         endif
         if (dolonbot) then
            if (rlonbot(i+1).ge.rlonbot(i)) then
               l1=int(rlonbot(i)/llint)+1
               l2=int(rlonbot(i+1)/llint)
               idr=1
            else
               l1=int(rlonbot(i)/llint)
               l2=int(rlonbot(i+1)/llint)+1
               idr=-1
            endif
            ypos=yllc-ydist
            do l=l1,l2,idr
               lval=l*llint
               xpos=xx(i)+(lval-rlonbot(i))/
     &            (rlonbot(i+1)-rlonbot(i))*ddx
               lval=mod(lval-ladd+899,360)-179
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.eq.180) then
                  write(axlab(1:3),'(a3)') '180'
                  iaxe=3
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' W'
                  iaxe=3
               elseif (lval.le.-10.and.lval.gt.-100) then
                  write(axlab(1:4),'(i2,a2)') -lval,' W'
                  iaxe=4
               elseif (lval.le.-100) then
                  write(axlab(1:5),'(i3,a2)') -lval,' W'
                  iaxe=5
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' E'
                  iaxe=3
               elseif (lval.ge.10.and.lval.lt.100) then
                  write(axlab(1:4),'(i2,a2)') lval,' E'
                  iaxe=4
               elseif (lval.ge.100) then
                  write(axlab(1:5),'(i3,a2)') lval,' E'
                  iaxe=5
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,0.)
            enddo
         endif
         if (dolontop) then
            if (rlontop(i+1).ge.rlontop(i)) then
               l1=int(rlontop(i)/llint)+1
               l2=int(rlontop(i+1)/llint)
               idr=1
            else
               l1=int(rlontop(i)/llint)
               l2=int(rlontop(i+1)/llint)+1
               idr=-1
            endif
            ypos=yurc+ydist
            do l=l1,l2,idr
               lval=l*llint
               xpos=xx(i)+(lval-rlontop(i))/
     &            (rlontop(i+1)-rlontop(i))*ddx
               lval=mod(lval-ladd+899,360)-179
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.eq.180) then
                  write(axlab(1:3),'(a3)') '180'
                  iaxe=3
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' W'
                  iaxe=3
               elseif (lval.le.-10.and.lval.gt.-100) then
                  write(axlab(1:4),'(i2,a2)') -lval,' W'
                  iaxe=4
               elseif (lval.le.-100) then
                  write(axlab(1:5),'(i3,a2)') -lval,' W'
                  iaxe=5
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' E'
                  iaxe=3
               elseif (lval.ge.10.and.lval.lt.100) then
                  write(axlab(1:4),'(i2,a2)') lval,' E'
                  iaxe=4
               elseif (lval.ge.100) then
                  write(axlab(1:5),'(i3,a2)') lval,' E'
                  iaxe=5
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,0.)
            enddo
         endif
         if (dolonleft) then
            if (rlonleft(i+1).ge.rlonleft(i)) then
               l1=int(rlonleft(i)/llint)+1
               l2=int(rlonleft(i+1)/llint)
               idr=1
            else
               l1=int(rlonleft(i)/llint)
               l2=int(rlonleft(i+1)/llint)+1
               idr=-1
            endif
            xpos=xllc-xdist
            do l=l1,l2,idr
               lval=l*llint
               ypos=yy(i)+(lval-rlonleft(i))/
     &            (rlonleft(i+1)-rlonleft(i))*ddy
               lval=mod(lval-ladd+899,360)-179
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.eq.180) then
                  write(axlab(1:3),'(a3)') '180'
                  iaxe=3
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' W'
                  iaxe=3
               elseif (lval.le.-10.and.lval.gt.-100) then
                  write(axlab(1:4),'(i2,a2)') -lval,' W'
                  iaxe=4
               elseif (lval.le.-100) then
                  write(axlab(1:5),'(i3,a2)') -lval,' W'
                  iaxe=5
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' E'
                  iaxe=3
               elseif (lval.ge.10.and.lval.lt.100) then
                  write(axlab(1:4),'(i2,a2)') lval,' E'
                  iaxe=4
               elseif (lval.ge.100) then
                  write(axlab(1:5),'(i3,a2)') lval,' E'
                  iaxe=5
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,1.)
            enddo
         endif
         if (dolonright) then
            if (rlonright(i+1).ge.rlonright(i)) then
               l1=int(rlonright(i)/llint)+1
               l2=int(rlonright(i+1)/llint)
               idr=1
            else
               l1=int(rlonright(i)/llint)
               l2=int(rlonright(i+1)/llint)+1
               idr=-1
            endif
            xpos=xurc+xdist
            do l=l1,l2,idr
               lval=l*llint
               ypos=yy(i)+(lval-rlonright(i))/
     &            (rlonright(i+1)-rlonright(i))*ddx
               lval=mod(lval-ladd+899,360)-179
               if (lval.eq.0) then
                  write(axlab(1:1),'(a1)') '0'
                  iaxe=1
               elseif (lval.eq.180) then
                  write(axlab(1:3),'(a3)') '180'
                  iaxe=3
               elseif (lval.lt.0.and.lval.gt.-10) then
                  write(axlab(1:3),'(i1,a2)') -lval,' W'
                  iaxe=3
               elseif (lval.le.-10.and.lval.gt.-100) then
                  write(axlab(1:4),'(i2,a2)') -lval,' W'
                  iaxe=4
               elseif (lval.le.-100) then
                  write(axlab(1:5),'(i3,a2)') -lval,' W'
                  iaxe=5
               elseif (lval.gt.0.and.lval.lt.10) then
                  write(axlab(1:3),'(i1,a2)') lval,' E'
                  iaxe=3
               elseif (lval.ge.10.and.lval.lt.100) then
                  write(axlab(1:4),'(i2,a2)') lval,' E'
                  iaxe=4
               elseif (lval.ge.100) then
                  write(axlab(1:5),'(i3,a2)') lval,' E'
                  iaxe=5
               endif
               call plchhq (xpos,ypos,axlab(1:iaxe),rtslb,0.,-1.)
            enddo
         endif
c
      enddo
c
      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)  ! restore set call
      call setusv ('LW',lwsv)
      call gsclip (iclp)  ! turn clipping back on
      return
      end
