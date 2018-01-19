c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine velbrb(u,lu,v,lv,ihem,lihem,m,n,rlength,ispv,spv)
c
c   Everything is as in velvct except:
c      There is no hi and lo specification.
c      NSET is hardwired to 1
c      RLENGTH is the length, in grid spaces, of a one-full-barb
c      wind arrow.  If the grid is irregular, a grid length is taken
c      as an x-direction grid length in the center of the grid.
c
      dimension u(lu,n),v(lv,n),ihem(lihem,n),spv(2)
      include 'comconst'
c
c   Determine vector length in fractional coordinates.
c
      xhalf=.5*m
      yhalf=.5*n
      gridfrac=cufx(fx(xhalf+.5,yhalf))-cufx(fx(xhalf-.5,yhalf))
      do 1000 iy=1,n
      do 1000 ix=1,m
c
      if (ispv.eq.1.and.u(ix,iy).eq.spv(1)) goto 999
      if (ispv.eq.2.and.v(ix,iy).eq.spv(2)) goto 999
      if (ispv.eq.3.and.(u(ix,iy).eq.spv(1).or.
     &    v(ix,iy).eq.spv(2))) goto 999
      if (ispv.eq.4.and.u(ix,iy).eq.spv(1).and.
     &    v(ix,iy).eq.spv(2)) goto 999
c
      rix=fx(float(ix),float(iy))
      riy=fy(float(ix),float(iy))
      call barb(rix,riy,u(ix,iy),v(ix,iy),rlength*gridfrac,ihem(ix,iy))
c
 999  continue
c
 1000 continue
      return
      end
