c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pfcalc(prs,sfp,pf,miy,mjx,mkzh)
c
c     Historically, this routine calculated the pressure at full sigma
c     levels when RIP was specifically designed for MM4/MM5 output.
c     With the new generalized RIP (Feb '02), this routine is still
c     intended to calculate a set of pressure levels that bound the
c     layers represented by the vertical grid points, although no such
c     layer boundaries are assumed to be defined.  The routine simply
c     uses the midpoint between the pressures of the vertical grid
c     points as the bounding levels.  The array only contains mkzh
c     levels, so the pressure of the top of the uppermost layer is
c     actually excluded.  The kth value of pf is the lower bounding
c     pressure for the layer represented by kth data level.  At the
c     lower bounding level of the lowest model layer, it uses the
c     surface pressure, unless the data set is pressure-level data, in
c     which case it assumes the lower bounding pressure level is as far
c     below the lowest vertical level as the upper bounding pressure
c     level is above.
c
      dimension prs(miy,mjx,mkzh),sfp(miy,mjx),pf(miy,mjx,mkzh)
c
      include 'comconst'
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         if (k.eq.mkzh) then
            if (iplevdata.ge.4) then ! terrain-following data
               pf(i,j,k)=sfp(i,j)
            else ! pressure-level data
               pf(i,j,k)=.5*(3.*prs(i,j,k)-prs(i,j,k-1))
            endif
         else
            pf(i,j,k)=.5*(prs(i,j,k+1)+prs(i,j,k))
         endif
      enddo
      enddo
      enddo
c
      return
      end
