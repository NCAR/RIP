c--------------------------------------------------------------------
      subroutine bshear (iw, ib, u, v, ght, ter, wk, miy, mjx, mkzh)
c calculate bulk shear 
c iw is the km level of the upper layer. iw > 0 is u, iw < 0 is v
c ib is the lower level (0, 500, 1000 m, etc.)
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh)
     & ,ter(miy,mjx), wk(miy,mjx,1)
      do 15 j = 1, mjx-1
	do 15 i = 1, miy-1
	  sdh = 0.
	  su = 0.
	  sv = 0.
c find the indices nearest .5 and 6 km AGL
	  k5 = 0
	  k6 = 0
	  top = abs(iw) * 1000.
	  bot = ib
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. top) .and. 
     &             (k6 .eq.0)) then
	      k6 = k
	      go to 8
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. bot) .and. 
     &           (k5 .eq. 0)) k5 = k
    6     continue
    8     continue
	  if (k6 .le. 1) k6=2
c calculate a weighted mean wind AGL over the upper layer
          do k = k5, k6, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
            su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
            sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c     if ( i .eq. 59 .and. j .eq. 88 ) then
c         write(6,*) 'sdh = ',sdh,' su = ',su,' sv = ',sv
c         write(6,*) 'ter = ',ter(i,j),' k5 = ',k5,' ht = ',
c    & ght(i,j,k5),' k6 = ',k6,' ht = ',ght(i,j,k6)
c     endif
	  ua = su / sdh
	  va = sv / sdh
	  sdh = 0.
	  su = 0.
	  sv = 0.
c calculate a weighted mean wind AGL over the lower layer
	  do k = mkzh, k5, -1
            dh = ght(i,j,k-1) - ght(i,j,k)
            sdh = sdh + dh
            su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
            sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
	  u0 = su / sdh
	  v0 = sv / sdh
c     if ( i .eq. 59 .and. j .eq. 88 ) then
c       write(6,*) 'i = ',i,' j = ',j,' k5 = ',k5,' k6 = ',k6
c       write(6,*) 'sfc u = ',u0,' sfc v = ',v0
c       write(6,*) 'bottom level u = ',u(i,j,mkzh),' bottom level v = ',
c    &    v(i,j,mkzh)
c       write(6,*) '6km u = ',ua,' 6km v = ',va
c      endif
	  if (iw .gt. 0) then
	    wk(i,j,1) = ua - u0
	  else
	    wk(i,j,1) = va - v0
	  endif
   15 continue
      end
