c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hboxdraw(ilinw,idash,ixwin,iywin,icolr,rcrag,rcrbg,
     &         maxpl,ipl,irota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   rcrbg(2,maxpl),icolr(maxpl),ilinw(maxpl),idash(maxpl),
     &   irota(maxpl)
c
      include 'comconst'
c
c   Convert plspecs to usable values for line width and dash pattern.
c
      lwidth=ilinw(ipl)*1000
      call getdash(idash(ipl),ndot)
c
c   Set line width, color, and dash pattern.
c
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
      call dashdb(ndot)
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
c   Box position will be adjusted for rotation further down.
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
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
c
c   Adjust box corner locations for rotation.
c
      if (irota(ipl).eq.90) then
         tempo=caxgn
         caxgn=iywin(2,ipl)-(caygn-iywin(1,ipl))
         caygn=tempo
         tempo=cbxgn
         cbxgn=iywin(2,ipl)-(cbygn-iywin(1,ipl))
         cbygn=tempo
      elseif (irota(ipl).eq.-90) then
         tempo=caxgn
         caxgn=caygn
         caygn=ixwin(2,ipl)-(tempo-ixwin(1,ipl))
         tempo=cbxgn
         cbxgn=cbygn
         cbygn=ixwin(2,ipl)-(tempo-ixwin(1,ipl))
      elseif (irota(ipl).eq.180.or.irota(ipl).eq.-180) then
         caxgn=ixwin(2,ipl)-(caxgn-ixwin(1,ipl))
         caygn=iywin(2,ipl)-(caygn-iywin(1,ipl))
         cbxgn=ixwin(2,ipl)-(cbxgn-ixwin(1,ipl))
         cbygn=iywin(2,ipl)-(cbygn-iywin(1,ipl))
      endif
c
c   Draw the box
c
      call lined(caxgn,caygn,caxgn,cbygn)
      call lined(caxgn,cbygn,cbxgn,cbygn)
      call lined(cbxgn,cbygn,cbxgn,caygn)
      call lined(cbxgn,caygn,caxgn,caygn)
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      call dashdb(65535)
c
      return
      end
