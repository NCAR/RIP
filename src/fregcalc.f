c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine fregcalc(clg,vis,meth,freg,miy,mjx)
c
c   This routine determines the flight regulation category based
c   on the given ceiling (meth=c), visibility (v), or ceiling AND
c   visibility (meth=b for both).
c      VFR:  1
c      MVFR: 2
c      IFR:  3
c      LIFR: 4
c
      dimension clg(miy,mjx),vis(miy,mjx),freg(miy,mjx)
      character meth*1
c
c   Transitions between flight categories:
c     Ceiling:     500 ft, 1000 ft, 3000 ft
c     Visibility: 1 mile, 3 miles, 5 miles
c     Values given are in meters for ceiling and km for visibility.
c
      include 'comconst'
c
      data clg1,clg2,clg3 / 152.4, 304.8, 914.4 /
      data vis1,vis2,vis3 / 1.609, 4.828, 8.046 /
c
      iclg=0
      ivis=0
      if (meth.eq.'c'.or. meth.eq.'b') then ! ceiling or both
         iclg=1
      elseif (meth.eq.'v'.or.meth.eq.'b') then ! visibility or both
         ivis=1
      endif
      do 405 j=1,mjx-1
      do 405 i=1,miy-1
         if ((clg(i,j).eq.rmsg.and.iclg.eq.1).or.
     &           (vis(i,j).eq.rmsg.and.ivis.eq.1)) then
            freg(i,j)=rmsg
         elseif ((clg(i,j).lt.clg1.and.iclg.eq.1).or.
     &           (vis(i,j).lt.vis1.and.ivis.eq.1)) then
            freg(i,j)=4.
         elseif ((clg(i,j).lt.clg2.and.iclg.eq.1).or.
     &           (vis(i,j).lt.vis2.and.ivis.eq.1)) then
            freg(i,j)=3.
         elseif ((clg(i,j).le.clg3.and.iclg.eq.1).or.
     &           (vis(i,j).le.vis3.and.ivis.eq.1)) then
            freg(i,j)=2.
         else
            freg(i,j)=1.
         endif
  405 continue
      return
      end
