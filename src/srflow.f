c--------------------------------------------------------------------
      subroutine srflo (u, v, ght, ter, ilow, g, miy, mjx, mkzh)
c Routine to compute storm-relative flow
c iuv is 0 for u, 1 for v;  ilow - 0 for low, 1 for high
c bsp = storm speed (m/s) (est to be 75 % of mean wind)
c bdr = storm direction (deg) (est. to be 30 degrees to the right of mean wind)
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh)
     & ,ter(miy,mjx), g(miy,mjx)
      parameter (pi=3.14159265, dtr=pi/180., dpr=180./pi)
c     rewind(58)
c     read(58,*) iii,jjj
      do 15 j = 1, mjx-1
	do 15 i = 1, miy-1
	  sdh = 0.
	  su = 0.
	  sv = 0.
c find the indices nearest 2 and 9 km AGL
	  k3 = 0
	  k10 = 0
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. 9000.) .and. 
     &             (k10 .eq. 0)) then
	      k10 = k
	      go to 8
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 2000.) .and. 
     &           (k3 .eq. 0)) k3 = k
    6     continue
    8     continue
	  if (k10 .eq. 0) k10=2
c calculate a 2-9 km AGL mean wind for storm motion
	  do k = k3, k10, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
	    su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
	    sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sdh = ',sdh,' su = ',su,' sv = ',sv
c           write(6,*) 'ter = ',ter(i,j),' k3 = ',k3,' ht = ',
c    &       ght(i,j,k3),' k10 = ',k10,' ht = ',ght(i,j,k10)
c         endif
	  ua = su / sdh
	  va = sv / sdh
	  asp = sqrt(ua*ua + va*va)
	  if (ua .eq. 0. .and. va .eq. 0.) then
	    adr = 0.
	  else
	    adr = dpr * (pi + atan2(ua,va))
	  endif
	  bsp = 0.75 * asp
	  bdr = adr + 30.
	  if (bdr .gt. 360.) bdr = bdr-360.
	  cu = -bsp * sin(bdr*dtr)
	  cv = -bsp * cos(bdr*dtr)
          sum = 0.
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'i = ',i,' j = ',j,' k3 = ',k3,
c    & ' k10 = ',k10
c          write(6,*) 'asp = ',asp,' adr = ',adr
c          write(6,*) 'bsp = ',bsp,' bdr = ',bdr
c         endif
	  sdh = 0.
	  su = 0.
	  sv = 0.
	  if (ilow .eq. 0) then
	    do k = mkzh-1, k3, -1
	      dh = ght(i,j,k) - ght(i,j,k+1)
	      sdh = sdh + dh
	      su = su + 0.5 * dh * 
     &                ((u(i,j,k)-cu)+(u(i,j,k+1)-cu))
	      sv = sv + 0.5 * dh * 
     &                ((v(i,j,k)-cv)+(v(i,j,k+1)-cv))
	    enddo
	    su = su / sdh
	    sv = sv / sdh
	    g(i,j) = sqrt(su*su + sv*sv)
	  else
	    do k = k3, k10, -1
	      dh = ght(i,j,k) - ght(i,j,k+1)
	      sdh = sdh + dh
	      su = su + 0.5 * dh * 
     &                ((u(i,j,k)-cu)+(u(i,j,k+1)-cu))
	      sv = sv + 0.5 * dh * 
     &                ((v(i,j,k)-cv)+(v(i,j,k+1)-cv))
	    enddo
	    su = su / sdh
	    sv = sv / sdh
	    g(i,j) = sqrt(su*su + sv*sv)
	  endif
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sreh = ',-sum
c         endif
   15 continue
      end
