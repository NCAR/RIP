c--------------------------------------------------------------------
      subroutine vgp (u, v, ght, ter, cape, g, miy, mjx, mkzh)
c Routine to compute vorticity-generation potential
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh)
     & ,ter(miy,mjx), cape(miy,mjx), g(miy,mjx)
      do j = 1, mjx-1
	do i = 1, miy-1
c find the index nearest 3 km AGL
	  k3 = 0
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. 3000.)) then
	      k3 = k
	      go to 8
	    endif
    6     continue
    8     continue
c calculate the 0-3 km hodograph length
	  totshr = 0.
	  depth=ght(i,j,k3)-ght(i,j,mkzh)
	  do k = mkzh-1, k3, -1
	    su = abs(u(i,j,k+1) - u(i,j,k))
	    sv = abs(v(i,j,k+1) - v(i,j,k))
	    totshr = sqrt(su*su + sv*sv) + totshr
	  enddo
	  totshr = totshr/depth
	  g(i,j) = sqrt(max(cape(i,j),0.)) * totshr
        enddo
      enddo
      end
