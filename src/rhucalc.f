c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine rhucalc(qvp,tmk,prs,icerel,rhu,miy,mjx,mkzh)
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),rhu(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 200 k=1,mkzh
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
         q = qvp(i,j,k)
         t = tmk(i,j,k)
         e = q*prs(i,j,k)/(eps+q)
c
c   es should be wrt ice if user has specified icerel and temperature
c   is below freezing.
c
         if (icerel.eq.1.and.t.lt.celkel) then
            es = ezero * exp( esicon1-esicon2/t )
         else
            es = ezero * exp( eslcon1*(t-celkel)/(t-eslcon2) )
         endif
         rhu(i,j,k)=100.*(e*(prs(i,j,k)-es))/(es*(prs(i,j,k)-e))
  200 continue
      return
      end
