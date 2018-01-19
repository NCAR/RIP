c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine tdpcalc(qvp,prs,tdp,miy,mjx,mkzh)
c
c   This routine calculates dewpoint temperature (in Celsius)
c
      dimension qvp(miy,mjx,mkzh),prs(miy,mjx,mkzh),tdp(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
         q=max(qvp(i,j,k),1.e-15)
         p=prs(i,j,k)
         e=q*p/(eps+q)
         elog=alog(e/ezero)
         tdp(i,j,k)=(eslcon2*elog-eslcon1*celkel)/(elog-eslcon1)-
     &      celkel
 1000 continue
      return
      end
