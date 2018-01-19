c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hsidsdraw(rtslb,ilinw,ixwin,iywin,icolr,rcrag,
     &   icaoid,iwmo,maxpl,ipl,irota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   icolr(maxpl),rtslb(maxpl),ilinw(maxpl),irota(maxpl)
      character icaoid*4, wmo*5
c
c   Routine to plot station IDs on the map background (ptyp=hb)
c
      include 'comconst'
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
c
c   Set user coordinate limits for unrotated view.
c
      ul=float(ixwin(1,ipl))
      ur=float(ixwin(2,ipl))
      ub=float(iywin(1,ipl))
      ut=float(iywin(2,ipl))
c
c   Transpose ul,ur,ub,ut if view is rotated +/- 90 degrees.
c   Station ID position will be adjusted for rotation further down.
c
      if (irota(ipl).eq.90.or.irota(ipl).eq.-90) then
         usv=ul
         ul=ub
         ub=usv
         usv=ur
         ur=ut
         ut=usv
      endif
c
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
c
      call gstxci(icolr(ipl))
      call gsplci(icolr(ipl))
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
c
c   Adjust marker location for rotation.
c
      if (irota(ipl).eq.90) then
         tempo=caxgn
         caxgn=iywin(2,ipl)-(caygn-iywin(1,ipl))
         caygn=tempo
      elseif (irota(ipl).eq.-90) then
         tempo=caxgn
         caxgn=caygn
         caygn=ixwin(2,ipl)-(tempo-ixwin(1,ipl))
      elseif (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
         caxgn=ixwin(2,ipl)-(caxgn-ixwin(1,ipl))
         caygn=iywin(2,ipl)-(caygn-iywin(1,ipl))
      endif
c
c   Draw the text
c
      if (ilinw(ipl) .lt. 0) then   ! draw a circle and offset sid
         size=(ur-ul)/(fr-fl)*rtslb(ipl)
         call ngdots(caxgn,caygn,1,size,icolr(ipl))
         if (icaoid.ne.'----') then
            call plchhq(caxgn+size,caygn,icaoid,rtslb(ipl),0.,-1.)
         else
            write(wmo,'(i5.5)') iwmo
            call plchhq(caxgn+size,caygn,wmo,rtslb(ipl),0.,-1.)
         endif
      else
         if (icaoid.ne.'----') then
            call plchhq(caxgn,caygn,icaoid,rtslb(ipl),0.,0.)
         else
            write(wmo,'(i5.5)') iwmo
            call plchhq(caxgn,caygn,wmo,rtslb(ipl),0.,0.)
         endif
      endif
c
      call gsplci(1)
      call gstxci(1)
      call gspmci(1)
      call gsfaci(1)
c
      return
      end
