c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine tfpcalc(qv,prs,tdp,miy,mjx,mkzh)
c
c   This routine calculates frostpoint temperature (in Celsius)
c
      dimension qv(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),tdp(miy,mjx,mkzh)
c
      include 'comconst'
c
      do k = 1, mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
        q=max(qv(i,j,k),1.e-15)
        p=prs(i,j,k)
        e=q*p/(eps+q)
c ice
        elog = alog(e/6.132)
        tdp(i,j,k) = (6133.-.61*elog)/(22.452-elog) - celkel
c water
c       elog = alog(e/6.133)
c       tdp(i,j,k) = (4780.8-32.19*elog)/(17.502-elog) - celkel
      enddo
      enddo
      enddo
      return
      end
