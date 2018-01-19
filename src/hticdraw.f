c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hticdraw(ngfbuf,ilinw,ixwin,iywin,raxlg,raxtg,icolr,
     &         rtslb,ilinwgf,ixwingf,iywingf,raxlggf,raxtggf,icolrgf,
     &         rtslbgf,iwhatgf,maxbuf,maxpl,ipl,irota,irotagf)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),raxtg(maxpl),
     &   icolr(maxpl),raxlg(maxpl),ilinw(maxpl),irotagf(maxbuf),
     &   rtslb(maxpl),ixwingf(2,maxbuf),rtslbgf(maxbuf),
     &   iywingf(2,maxbuf),raxtggf(maxbuf),raxlggf(maxbuf),
     &   ilinwgf(maxbuf),icolrgf(maxbuf),iwhatgf(maxbuf),irota(maxpl)
c
      character axlab*5
c
      include 'comconst'
c
c   Check if this tic background is identical to any of those saved
c      in the GFLASH buffers, and if so, flush that buffer to
c      this frame.
c
      do 10 ib=1,ngfbuf
         if (iwhatgf(ib).ne.3) goto 10
         if ( ixwin(1,ipl).eq.ixwingf(1,ib).and.ixwin(2,ipl).eq.
     &        ixwingf(2,ib).and.iywin(1,ipl).eq.iywingf(1,ib).and.
     &        iywin(2,ipl).eq.iywingf(2,ib).and.raxtg(ipl).eq.
     &        raxtggf(ib).and.raxlg(ipl).eq.raxlggf(ib).and.
     &        ilinw(ipl).eq.ilinwgf(ib).and.icolr(ipl).eq.
     &        icolrgf(ib).and.rtslb(ipl).eq.
     &        rtslbgf(ib).and.irota(ipl).eq.
     &        irotagf(ib)) then
            call gflas3(ib)
            goto 120
         endif
   10 continue
c
c   If identical tic background was not found in GFLASH buffer,
c      draw the tic background to a new GFLASH buffer, flush this
c      new buffer to the current frame, and save the parameters for
c      this new tic background.
c
      ngfbuf=ngfbuf+1
      if (ngfbuf.gt.maxbuf) then
         write(iup,*)'In hticdraw: too many gflash buffers.'
         stop
      endif
      call gflas1(ngfbuf)
c
c   Set line width, color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
c
c   Assign values that depend on rotation of view
c
      if (irota(ipl).eq.0) then
         xintervs=ixwin(2,ipl)-ixwin(1,ipl)
         yintervs=iywin(2,ipl)-iywin(1,ipl)
         j1=ixwin(1,ipl)
         j2=ixwin(2,ipl)
         i1=iywin(1,ipl)
         i2=iywin(2,ipl)
      elseif (irota(ipl).eq.90) then
         xintervs=iywin(2,ipl)-iywin(1,ipl)
         yintervs=ixwin(2,ipl)-ixwin(1,ipl)
         j1=iywin(2,ipl)
         j2=iywin(1,ipl)
         i1=ixwin(1,ipl)
         i2=ixwin(2,ipl)
      elseif (irota(ipl).eq.-90) then
         xintervs=iywin(2,ipl)-iywin(1,ipl)
         yintervs=ixwin(2,ipl)-ixwin(1,ipl)
         j1=iywin(1,ipl)
         j2=iywin(2,ipl)
         i1=ixwin(2,ipl)
         i2=ixwin(1,ipl)
      elseif (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
         xintervs=ixwin(2,ipl)-ixwin(1,ipl)
         yintervs=iywin(2,ipl)-iywin(1,ipl)
         j1=ixwin(2,ipl)
         j2=ixwin(1,ipl)
         i1=iywin(2,ipl)
         i2=iywin(1,ipl)
      endif
      jdir=(j2-j1)/abs(j2-j1)
      idir=(i2-i1)/abs(i2-i1)
c
c   Make set call and draw perimeter.
c
      faspect=(ftmax-fbmin)/(frmax-flmin)
      aspect=yintervs/xintervs
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
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call line(fl,fb,fl,ft)
      call line(fr,fb,fr,ft)
      call line(fl,ft,fr,ft)
      call line(fl,fb,fr,fb)
c
c   Draw tick marks and label axes.
c
      ydist=(.6*rtslb(ipl)+.01)
      xdist=.01
      tickl=.005
      gint=(fr-fl)/xintervs
      if (raxlg(ipl).eq.rmsg) then
         iaxlg=10
      else
         iaxlg=nint(raxlg(ipl))
      endif
      if (raxtg(ipl).eq.rmsg) then
         iaxtg=1
      else
         iaxtg=nint(raxtg(ipl))
      endif
c
      if (iaxlg.ne.0) then
         do j=j1,j2,jdir
            xtick=fl+abs(j-j1)*gint
            if(mod(j,iaxlg).eq.0) then
               write(axlab,910) j
               iaxs=5-int(alog10(float(j)+.5))
               call setusv('LW',1000)
               call plchhq(xtick,fb-ydist,axlab(iaxs:),
     &            rtslb(ipl),0.,0.)
               call setusv('LW',lwidth)
               call line(xtick,fb,xtick,fb+2.*tickl)
               call line(xtick,ft,xtick,ft-2.*tickl)
            endif
         enddo
         do i=i1,i2,idir
            ytick=fb+abs(i-i1)*gint
            if(mod(i,iaxlg).eq.0) then
               write(axlab,910) i
               call setusv('LW',1000)
               call plchhq(fl-xdist,ytick,axlab,rtslb(ipl),0.,1.)
               call setusv('LW',lwidth)
               call line(fl,ytick,fl+2.*tickl,ytick)
               call line(fr,ytick,fr-2.*tickl,ytick)
            endif
         enddo
      endif
      if (iaxtg.ne.0) then
         do j=j1,j2,jdir
            xtick=fl+abs(j-j1)*gint
            if(mod(j,iaxtg).eq.0) then
               call line(xtick,fb,xtick,fb+tickl)
               call line(xtick,ft,xtick,ft-tickl)
            endif
         enddo
         do i=i1,i2,idir
            ytick=fb+abs(i-i1)*gint
            if(mod(i,iaxtg).eq.0) then
               call line(fl,ytick,fl+tickl,ytick)
               call line(fr,ytick,fr-tickl,ytick)
            endif
         enddo
      endif
  910 format(i5)
c
      call gflas2
      call gflas3(ngfbuf)
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      iywingf(1,ngfbuf)=iywin(1,ipl)
      iywingf(2,ngfbuf)=iywin(2,ipl)
      ixwingf(1,ngfbuf)=ixwin(1,ipl)
      ixwingf(2,ngfbuf)=ixwin(2,ipl)
      raxtggf(ngfbuf)=raxtg(ipl)
      raxlggf(ngfbuf)=raxlg(ipl)
      ilinwgf(ngfbuf)=ilinw(ipl)
      icolrgf(ngfbuf)=icolr(ipl)
      rtslbgf(ngfbuf)=rtslb(ipl)
      irotagf(ib)=irota(ipl)
      iwhatgf(ngfbuf)=3
c      write(iup,*)'   Made new tic background number ',ngfbuf,
c     &   ' for plot ',ipl
c
  120 continue
      return
      end
