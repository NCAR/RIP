c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine sticdraw(icolr,igray,tcolor,mcolor,
     &   ilinw,ngfbuf,icolrgf,ilinwgf,
     &   lhodo,lmand,lsndg,lhodogf,lmandgf,lsndggf,iwhatgf,maxbuf,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,ipl)
c
c  Plots grid for a skewt, log p thermodynamic diagram
c
      dimension icolr(maxpl),ilinw(maxpl),icolrgf(maxbuf),
     &   ilinwgf(maxbuf),iwhatgf(maxbuf)
      logical lhodo(maxpl),lsndg(maxpl),lmand(maxpl)
      logical lhodogf(maxbuf), lsndggf(maxbuf), lmandgf(maxbuf)
      integer tcolor, mcolor
c
      include 'comconst'
c
c  Encoding buffer
c
      character*4 itit
c      character*12 ititl
c
c  Skewt border
c
      dimension xb(7),yb(7)
c
c  Pressure line specs
c
      dimension plvm(20),iplvmtype(20),plvi(20),iplvitype(20),pln(11,2)
c
c  Temperature line specs
c
      dimension tp(15,2)
c
c  Mixing ratio specs
c
      real rat(10)
      character*2 lrat(10)
c
c  Dry/saturated adiabat buffers
c
      dimension sx(162),sy(162),y45(162)
c
      data xb /-19.,27.1,27.1,18.6,18.6,-19.,-19./
      data yb /-.9346217,-.9346217,9.,17.53,44.061,44.061,-.9346217/
      data plvm /100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,
     &   650.,700.,750.,800.,850.,900.,925.,950.,1000./
      data iplvmtype / 1,1,1,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,0,1/
      data plvi /100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,
     &   650.,700.,750.,800.,850.,900.,950.,1000.,1000./
      data iplvitype / 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1/
      data pln /11*-19.,4*18.6,22.83,26.306,5*27.1/
      data tp /8*1050.,855.,625.,459.,337.,247.,181.,132.,730.,
     &   580.,500.,430.,342.,251.,185.,135.,7*100./
      data rat /24.,20.,16.,12.,8.,5.,3.,2.,1.,0.4/
      data lrat /'24','20','16','12',' 8',' 5',' 3',' 2',' 1','.4'/
c
c  Check if skewt background already saved,in a gflash buffers, and,
c      if so, flush that buffer to this frame.
c
      do 5 ib=1,ngfbuf
         if (icolr(ipl).eq.icolrgf(ib).and.ilinw(ipl).eq.
     &       ilinwgf(ib).and.(lhodo(ipl).eqv.lhodogf(ib)).and.
     &       (lmand(ipl).eqv.lmandgf(ib)).and.
     &       (lsndg(ipl).eqv.lsndggf(ib)).and.
     &       iwhatgf(ib).eq.2) then
            call gflas3(ib)
            goto 120
         endif
   5  continue
c
c   If no skewt background saved in gflash buffer, then draw the
c   background to a new gflash buffer, and flush this new buffer to the
c   current frame.
c
      ngfbuf=ngfbuf+1
      if (ngfbuf.gt.maxbuf) then
         write(iup,*)'In sticdraw: too many gflash buffers.'
         stop
      endif
      call gflas1(ngfbuf)
c
c  Set line width, color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
c
c  Draw skewt border
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
c
      call curve(xb,yb,7)
c
c  Draw the pressure lines
c
c      data pln /11*-19.,4*18.6,22.83,26.306,5*27.1/
c  400  -625 500
      call gqplci(ierr,iplcisv)
      call gqtxci(ierr,itxcisv)
      do k=1,20
         if (lmand(ipl)) then
            plv=plvm(k)
            iplvtype=iplvmtype(k)
         else
            plv=plvi(k)
            iplvtype=iplvitype(k)
         endif
         if ((lsndg(ipl) .and. plv.gt.500.) .or. (lhodo(ipl).and.
     &       plv.lt.200)) then
            x1=-5.6
         else
            x1=-19.
         endif
         if (plv.le.400.) then
            x2=18.6
         elseif (plv.ge.625.) then
            x2=27.1
         else
            x2=18.6+alog10(plv/400.)/.19*8.5
         endif
         y1=132.182-44.061*alog10(plv)
         if (iplvtype.eq.1) then
            call gsplci(iplcisv)
            call gstxci(itxcisv)
            call line(x1,y1,x2,y1)
            its=nint(plv)
            write (itit,17) its
            call setusv('LW',1000)
            call plchhq (-19.2,y1,itit,.011,0.,1.)
            call setusv('LW',lwidth)
         else
            call gsplci(igray)
            call gstxci(igray)
            call line(x1,y1,x2,y1)
         endif
      enddo
   17 format(i4)
      call gsplci(iplcisv)
      call gstxci(itxcisv)
c
c  Draw temperature lines
c
      call gsplci(tcolor)   
      if (lsndg(ipl)) then
        ntemps = 5
      else if (lhodo(ipl)) then
        ntemps = 11
      else
        ntemps = 15
      endif
      t=40.
      do 20 i=1,ntemps
         y1=132.182-44.061*alog10(tp(i,1))
         y2=132.182-44.061*alog10(tp(i,2))
         x1=0.54*t+0.90692*y1
         x2=0.54*t+0.90692*y2
         call line(x1,y1,x2,y2)
         its=nint(t)
         if(its.eq.20) go to 19
         if (lhodo(ipl) .and. its .le. -80) go to 19
c         x2=x2+0.4
c         y2=y2+0.441
         x2=x2+0.1
         y2=y2+0.11
         write (itit,17) its
         call setusv('LW',1000)
         call plchhq (x2,y2,itit,.012,47.,-1.)
         call setusv('LW',lwidth)
   19    t=t-10.
   20 continue
c
      if (lhodo(ipl)) then
c
c draw temperature lines around hodograph box
c
        y1 = 132.182-44.061*alog10(337.)
        y2 = 132.182-44.061*alog10(200.)
        x1 = 0.54*(-70.)+0.90692*y1
        x2 = 0.54*(-70.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(158.)
        y2 = 132.182-44.061*alog10(100.)
        x1 = 0.54*(-70.)+0.90692*y1
        x2 = 0.54*(-70.)+0.90692*y2
        call line (x1,y1,x2,y2)
c        x2 = x2 + 0.4
c        y2 = y2 + 0.441
        x2 = x2 + 0.1
        y2 = y2 + 0.11
        its = -70
        write (itit,17) its
        call setusv('LW',1000)
        call plchhq (x2,y2,itit,.012,47.,-1.)
        call setusv('LW',lwidth)
        y1 = 132.182-44.061*alog10(247.)
        y2 = 132.182-44.061*alog10(200.)
        x1 = 0.54*(-80.)+0.90692*y1
        y1 = 132.182-44.061*alog10(115.)
        y2 = 132.182-44.061*alog10(100.)
        x1 = 0.54*(-80.)+0.90692*y1
        x2 = 0.54*(-80.)+0.90692*y2
        call line (x1,y1,x2,y2)
c        x2 = x2 + 0.4
c        y2 = y2 + 0.441
        x2 = x2 + 0.1
        y2 = y2 + 0.11
        its = -80
        write (itit,17) its
        call setusv('LW',1000)
        call plchhq (x2,y2,itit,.012,47.,-1.)
        call setusv('LW',lwidth)
      endif
c
      if (lsndg(ipl)) then
c
c draw temperature lines around sndg analysis box
c
        y1 = 132.182-44.061*alog10(1020.)
        y2 = 132.182-44.061*alog10(tp(6,2))
        x1 = 0.54*(-10.)+0.90692*y1
        x2 = 0.54*(-10.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(748.)
        y2 = 132.182-44.061*alog10(tp(7,2))
        x1 = 0.54*(-20.)+0.90692*y1
        x2 = 0.54*(-20.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(550.)
        y2 = 132.182-44.061*alog10(tp(8,2))
        x1 = 0.54*(-30.)+0.90692*y1
        x2 = 0.54*(-30.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(500.)
        y2 = 132.182-44.061*alog10(tp(9,2))
        x1 = 0.54*(-40.)+0.90692*y1
        x2 = 0.54*(-40.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(500.)
        y2 = 132.182-44.061*alog10(tp(10,2))
        x1 = 0.54*(-50.)+0.90692*y1
        x2 = 0.54*(-50.)+0.90692*y2
        call line (x1,y1,x2,y2)
        y1 = 132.182-44.061*alog10(tp(11,1))
        y2 = 132.182-44.061*alog10(tp(11,2))
        x1 = 0.54*(-60.)+0.90692*y1
        x2 = 0.54*(-60.)+0.90692*y2
        call line (x1,y1,x2,y2)
        if (.not.lhodo(ipl)) then
          do its = -60,-100,-10
            ii = (abs(its + 10)*.1) + 6
            y1=132.182-44.061*alog10(tp(ii,1))
            y2=132.182-44.061*alog10(tp(ii,2))
            x1=0.54*float(its)+0.90692*y1
            x2=0.54*float(its)+0.90692*y2
            call line(x1,y1,x2,y2)
          enddo
        endif
c
        if (lhodo(ipl)) then
          lastt = -60    ! We've already written -70 and -80 above
        else
          lastt = -100
        endif
        do 223 its = -10, lastt, -10
          ii = (abs(its + 10)*.1) + 6
c          y2=132.182-44.061*alog10(tp(ii,2)) + .441
c          x2=0.54*float(its)+0.90692*y2 + .4
          y2=132.182-44.061*alog10(tp(ii,2)) + .11
          x2=0.54*float(its)+0.90692*y2 + .1
          write (itit,17) its
          call setusv('LW',1000)
          call plchhq (x2,y2,itit,.012,47.,-1.)
          call setusv('LW',lwidth)
  223   continue
      endif
      call gsplci(iplcisv)
c
c  Tick marks at 500 mb
c
      y1=13.2627
      y2=13.75
      t=-52.
      do 25 i=1,31
         t=t+2.
         if(amod(t,10.).eq.0.)go to 25
         x1=0.54*t+0.90692*y1
         x2=0.54*t+0.90692*y2
         call line(x1,y1,x2,y2)
   25 continue
c
c  Draw mixing ratio lines
c
      pmin = 550.
      call gsplci(mcolor) 
      call dashdb (13107)     ! pattern = 0011001100110011
c     call dashdb (3855)      ! pattern = 0000111100001111
      y1=132.182-44.061*alog10(1050.)
      y2=132.182-44.061*alog10(pmin)
      if (lsndg(ipl)) then
	ist = 8
      else
	ist = 10
      endif
      do 30 i = 1,ist
         es1=.001*rat(i)*1050./(eps+.001*rat(i))
         elog1=alog(es1/ezero)
         tmr1=(eslcon2*elog1-eslcon1*celkel)/(elog1-eslcon1)
         es2=.001*rat(i)* pmin/(eps+.001*rat(i))
         elog2=alog(es2/ezero)
         tmr2=(eslcon2*elog2-eslcon1*celkel)/(elog2-eslcon1)
         x1=0.54*(tmr1-celkel)+0.90692*y1
         x2=0.54*(tmr2-celkel)+0.90692*y2
         call lined(x1,y1,x2,y2)
         call setusv('LW',1000)
         call plchhq (x2,y2+0.6,lrat(i),.010,0.,0.)
         call setusv('LW',lwidth)
   30 continue
      call gsplci(iplcisv)
c
c  Draw saturated adiabats
c
      call gsplci(mcolor)  
c     call dashdb (31710)     ! pattern = 0111101111011110
      call dashdb (65535)     ! pattern = 1111111111111111
      ts=32.
      do 40 i=1,14
         tk=ts+celkel
         es = ezero * exp( eslcon1*(tk-celkel)/(tk-eslcon2) )
         ws=eps*es/(1000.-es)
         aos=tk*exp((thtecon1/tk-thtecon2)*ws*(1.+thtecon3*ws))
	 k = 0
         p=1060.
         do 35 j=1,86
            p=p-10.
            atsa=tonpsadiabat(aos,p)-celkel
            ypd=132.182-44.061*alog10(p)
            xpd=0.54*atsa+0.90692*ypd
            if (xpd .lt. -19.0) go to 36
	    if (lsndg(ipl) .and. xpd .le. -5.75 .and.
     &         ypd .lt. 13.26) GO TO 36
            k = k + 1
            sy(k)=ypd
            sx(k)=xpd
   35    continue
   36    call curved(sx,sy,k)
         its=nint(ts)
c        if (its.lt.10) then
c           write (ititl,'(i1,a1,f5.1,a1)') its,'(',aos,')'
c           iend=8
c        else
c           write (ititl,'(i2,a1,f5.1,a1)') its,'(',aos,')'
c           iend=9
c        endif
         write (itit,102) its
         call setusv('LW',1000)
         x = sx(k)
         y = sy(k)
c        call plchhq (sx(86),sy(86)+0.6,ititl(1:iend),.008,0.,0.)
         if (its .gt. -10 .and. y .gt. 13.26) 
     &      call plchhq (x+.2,y-0.7,itit(1:2),.010,0.,0.)
         call setusv('LW',lwidth)
         ts=ts-4.0
   40 continue
  102 format(i2)
      call gsplci(iplcisv)
c
c  Draw dry adiabat lines
c
      call gsplci(tcolor) 
      call dashdb (65535)     ! pattern = 1111111111111111
c     call dashdb (21845)     ! pattern = 0101010101010101
cc    call dashd(4444b)
      t=51.
      do 45 i=1,162
         y45(i)=66.67*(5.7625544-alog(t+celkel))
         t=t-1.0
   45 continue
      t=450.
      td=52.
      do 55 i=1,20
         t=t-10.
         k=0
	 kkk = 0
         yd=66.67*(alog(t)-5.7625544)
         do 50 j=1,162
            ypd=y45(j)+yd
            tx=td-float(j)
            if(ypd.gt.44.061) go to 54
            if(ypd.lt.-.9346217) go to 50
            xpd=0.54*tx+0.90692*ypd
            if(xpd.lt.-19.0)go to 54
            if(xpd.gt.27.1)go to 50
            if(xpd.gt.18.6.and.t.gt.350.0)go to 50
            if (lhodo(ipl) .and. xpd .le. -5.75 .and.
     &         ypd .gt. 30.87) GO TO 50
c           if (lsndg(ipl).and. nint(t) .eq. 270 ) then
c           write(6,*) 't = ',nint(t),' xpd = ',xpd,' ypd = ',ypd
c           endif
c For dry adiabats behind the parcel box, we need to split them into
c two pieces for drawing.
c
            if (lsndg(ipl) .and. xpd .le. -5.75 .and.
     &         ypd .lt. 13.26) then
	      if (k .gt. 0 .and. kkk .eq. 0) then
		k=k+1
		sx(k)=xpd
		sy(k)=ypd
		call curved(sx,sy,k)
		kkk = 1
		k = 0
	      endif
	      GO TO 50
            endif
c    if (lsndg(ipl) .and. nint(t) .eq. 270 .or.
c    &          (lsndg(ipl) .and. nint(t) .eq. 280)) then
c    write(6,*) 't = ',t,' j = ',j,' xpd = ',xpd,' ypd = ',ypd
c  endif
            k=k+1
            sx(k)=xpd
            sy(k)=ypd
   50    continue
c
   54    call curved(sx,sy,k)
	 if (k .lt. 4) go to 55
         its=nint(t)
         write (itit,103) its
         x=sx(k-3)
         y=sy(k-3)
         if(x.lt.-15.0) x = -17.95
         if(y.gt.40.0)  y = 42.9
         call setusv('LW',1000)
         if (lhodo(ipl)) then
           if (its .ge. 370) then
             call plchhq (x,42.9,itit(1:3),.010,0.,0.)
           else if (its .le. 290) then
             call plchhq (-17.95,y,itit(1:3),.010,0.,0.)
           endif
         else
           call plchhq (x,y,itit(1:3),.010,0.,0.)
         endif
         call setusv('LW',lwidth)
   55 continue
  103 format(i3)
      call gsplci(iplcisv)
c
      call gflas2
      call gflas3(ngfbuf)
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      icolrgf(ngfbuf)=icolr(ipl)
      ilinwgf(ngfbuf)=ilinw(ipl)
      lhodogf(ngfbuf)=lhodo(ipl)
      lmandgf(ngfbuf)=lmand(ipl)
      lsndggf(ngfbuf)=lsndg(ipl)
      iwhatgf(ngfbuf)=2
c      write(iup,*)'   Made skewt background in buffer ',ngfbuf,
c     &   ' for plot ',ipl
c
  120 continue
      return
      end
