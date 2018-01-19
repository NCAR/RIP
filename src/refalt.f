c--------------------------------------------------------------------
      subroutine reflalt (z, p, alt, miy, mjx, mkzh)
c Routine to compute highest altitude of 10 dBZ reflectivity
      real z(miy,mjx,mkzh), p(miy,mjx,mkzh)
     & , alt(miy,mjx)

      do j = 1, mjx-1
	do i = 1, miy-1
	  alt(i,j) = 9999.
	  do 8 k = 1, mkzh
	    if (z(i,j,k) .gt. 10.) then
	      alt(i,j) = p(i,j,k)
	      go to 10
	    endif
    8     continue  ! trop_loop
   10     continue
        enddo
       enddo
       end
