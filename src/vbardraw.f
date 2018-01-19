c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vbardraw(ilinw,idash,icolr,set1,set2,xdist,ydist,
     &   rcrag,rslcg,xseclen,rtslb,ctitl,igdir,maxpl,ipl)
c
      dimension idash(maxpl),rcrag(2,maxpl),rslcg(2,maxpl),
     &   ilinw(maxpl),icolr(maxpl),rtslb(maxpl),igdir(maxpl)
      character ctitl(maxpl)*82
c
      dimension rect(4)
c
      include 'comconst'
c
      fb=fbmin
      ft=ftmax
      fl=flmin
      fr=frmax
c
c   Make set call.
c
      call set(fl,fr,fb,ft,0.,xseclen,set1,set2,1)
      ftoud=xseclen/(fr-fl)
      ftouv=(set2-set1)/(ft-fb)
      xsecleni=1./xseclen
      if (igdir(ipl).ne.362) then
         tanb=tan(float(igdir(ipl))*rpd)
      else
         tanb=0.
      endif
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
c   Convert rcrag and rslcg (which are in coarse dom. grid points) to
c   grid points in the current domain.
c
      crax=1.+refrat*(rcrag(2,ipl)-xjcorn)
      cray=1.+refrat*(rcrag(1,ipl)-yicorn)
      slcx=1.+refrat*(rslcg(2,ipl)-xjcorn)
      slcy=1.+refrat*(rslcg(1,ipl)-yicorn)
c
c   Draw the vertical bar
c
      xveca=slcx-crax
      yveca=slcy-cray
      dpos=(xveca*xdist+yveca*ydist-(yveca*xdist-xveca*ydist)*tanb)*
     &   xsecleni
      call lined(dpos,set1,dpos,set2)
c
      if (ctitl(ipl)(1:8).ne.'auto    ') then
         call gqclip(ierr,iclp,rect)
         call gsclip(0)
         iendct=lennonblank(ctitl(ipl))
         if (iendct.lt.1) goto 30
         ypos=set2+.015*ftouv
         call plchhq(dpos,ypos,ctitl(ipl)(1:iendct),
     &                 rtslb(ipl),0.,0.)
         call gsclip(iclp)
      endif
 30   continue
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      call dashdb(65535)
c
      return
      end
