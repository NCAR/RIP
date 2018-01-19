c--------------------------------------------------------------------
      subroutine trop (t, z, p, trp, miy, mjx, mkzh)
c Routine to compute tropopause. Based on code from Cameron Homeyer and WMO definition
      real t(miy,mjx,mkzh), z(miy,mjx,mkzh), p(miy,mjx,mkzh)
     & , trp(miy,mjx)
      integer dtdztest(mkzh), ztest(mkzh)

      do j = 1, mjx-1
	do i = 1, miy-1
	  trp(i,j) = 0.
c    trop_loop
	  do 8 k = mkzh, 2, -1
            if (p(i,j,k) .le. 500.) then
	      dtdz = (t(i,j,k+1) - t(i,j,k))/(z(i,j,k) - z(i,j,k+1))    ! Compute lapse rate (-dT/dz)
	    else
	      dtdz = 999.9                              ! Set lapse rate for p > 500 hPa
	    endif
	    if (dtdz .le. 0.002) then
              do m = 1, mkzh
	        dtdztest = 0                            ! Initialize lapse rate test array
		ztest = 0                               ! Initialize altitude test array
	      enddo
	      do k2 = k-1, 1, -1
		! Compute lapse rate at levels above current candidate
	        dtdz = (t(i,j,k2) - t(i,j,k))/(z(i,j,k) - z(i,j,k2))       
		if ((dtdz .le. 0.002) .and. 
     &                ((z(i,j,k) - z(i,j,k2)) .le. 2000.)) THEN
		  ! If lapse rate <= 2 K/km and z <= trop + 2 km, set pass flag
		  dtdztest(k2) = 1                        
		endif
		IF ((z(i,j,k) - z(i,j,k2)) .le. 2000.0) THEN
		  ! If z <= trop + 2 km, set pass flag
		  ztest(k2) = 1                           
		endif
	      enddo   !k2 loop

	      IF (SUM(dtdztest) .ne. SUM(ztest)) THEN
	    ! If the number of lapse rate passes not equal to number of altitude passes, go on
	       go to 8
	      ELSE
	       ! If qualified tropopause, set altitude index and return value
	        ktrp = k  
	        go to 10
	      ENDIF
	    ENDIF
    8     continue  ! trop_loop
   10     continue
          trp(i,j) = p(i,j,ktrp)
        enddo
       enddo
       end
