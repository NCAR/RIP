c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine barb(xlocu,ylocu,xmag,ymag,arrowlength,ihem)
c
c   Draws a wind barb at point (xlocu,ylocu) in user coordinate
c      system, with magnitude <xmag,ymag>. The length of a 10-unit
c      wind (arrow with one full barb) is given by arrowlength
c      (in fractional coordinates).
c   ihem=1 for N. Hem., =-1 for S. Hem.
c
      xloc=cufx(xlocu)
      yloc=cufy(ylocu)
      call getset(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      call set(fl,fr,fb,ft,fl,fr,fb,ft,1)
      pi=3.141592654
      dstem=arrowlength
      xnorm=dstem/3.8
      radb=xnorm
      dstembb=1.2*xnorm
      dstembf=1.6*xnorm
      dstemtf=.8*xnorm
      dhbarb=1.3*xnorm
      dbarb=2.6*xnorm
      angbarb=ihem*70.*pi/180.
      flagwid=.8*dstembf
      speed=sqrt(xmag**2+ymag**2)
      nspeed=nint(speed/5.)
      xcentr=xloc
      ycentr=yloc
c
      if (speed.lt..5) then ! We'll say "calm" is < 0.5 "units",
c                           ! where a full barb = 10 "units"
c
c      Circle
c
         xf=xcentr+radb
         yf=ycentr
         do 900 ip=1,32
            ang=2.*pi*ip/32.
            xt=xcentr+radb*cos(ang)
            yt=ycentr+radb*sin(ang)
            call line(xf,yf,xt,yt)
            xf=xt
            yf=yt
900      continue
      else
         angle=atan2(ymag,xmag)+pi
         nflag=nspeed/10
         nbarb=(nspeed-nflag*10)/2
         nhbarb=(nspeed-nflag*10-nbarb*2)
         xf=xcentr
         yf=ycentr
         xt=xf+dstem*cos(angle)
         yt=yf+dstem*sin(angle)
         call line(xf,yf,xt,yt)
         xf=xt
         yf=yt
         xt=xf+dhbarb*cos(angle-angbarb)
         yt=yf+dhbarb*sin(angle-angbarb)
         if (nhbarb.eq.1) call line(xf,yf,xt,yt)
         if (nspeed.eq.1) then
            xt=xf+dstembb*cos(angle)
            yt=yf+dstembb*sin(angle)
            call line(xf,yf,xt,yt)
         endif
         do 904 ib=1,nbarb
            xt=xf+dstembb*cos(angle)
            yt=yf+dstembb*sin(angle)
            call line(xf,yf,xt,yt)
            xf=xt
            yf=yt
            xt=xf+dbarb*cos(angle-angbarb)
            yt=yf+dbarb*sin(angle-angbarb)
            call line(xf,yf,xt,yt)
 904     continue
         if (nflag.gt.0) then
            xt=xf+dstemtf*cos(angle)
            yt=yf+dstemtf*sin(angle)
            call line(xf,yf,xt,yt)
            xf=xt
            yf=yt
         endif
         do 908 ib=1,nflag
            xt=xf+dstembf*cos(angle)
            yt=yf+dstembf*sin(angle)
            call line(xf,yf,xt,yt)
            xf=xt
            yf=yt
            xf2=xf-flagwid*cos(angle)
            yf2=yf-flagwid*sin(angle)
            xt=xf+dbarb*cos(angle-angbarb)-.5*flagwid*cos(angle)
            yt=yf+dbarb*sin(angle-angbarb)-.5*flagwid*sin(angle)
            call line(xf,yf,xt,yt)
            call line(xf2,yf2,xt,yt)
 908     continue
c
      endif
c
      call set(fl,fr,fb,ft,ul,ur,ub,ut,ll)
      return
      end
