c--------------------------------------------------------------------
      subroutine srflo2 (u, v, ght, ter, ilow, iuv, g, miy, mjx, mkzh)
c Routine to compute storm-relative flow for estimating supercell
c type. 
c iuv - 0 for 0, 1 for v;  ilow - 0 for storm motion, 1 for 
c storm-relative flow
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
          if (ilow .eq. 0) then
	    if (iuv .eq. 0) then
	      g(i,j) = cu
	    else
	      g(i,j) = cv
	    endif
	    go to 15
          endif
c         sum = 0.
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'i = ',i,' j = ',j,' k3 = ',k3,
c    & ' k10 = ',k10
c          write(6,*) 'asp = ',asp,' adr = ',adr
c          write(6,*) 'bsp = ',bsp,' bdr = ',bdr
c         endif
	  sdh = 0.
	  su = 0.
	  sv = 0.
	  do k = k10, k10-1, -1
	    dh = ght(i,j,k) - ght(i,j,k+1)
	    sdh = sdh + dh
	    su = su + 0.5 * dh * 
     &              ((u(i,j,k)-cu)+(u(i,j,k+1)-cu))
	    sv = sv + 0.5 * dh * 
     &              ((v(i,j,k)-cv)+(v(i,j,k+1)-cv))
	  enddo
	  su = su / sdh
	  sv = sv / sdh
	  g(i,j) = sqrt(su*su + sv*sv)
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sreh = ',-sum
c         endif
   15 continue
      end
c--------------------------------------------------------------------
      subroutine srflo3 (u, v, ght, ter, ilow, iuv, g, miy, mjx, mkzh)
c uses Rasmussen and Straka's scheme
c Routine to compute storm-relative flow for estimating supercell
c type. 
c iuv - 0 for 0, 1 for v;  ilow - 0 for storm motion, 1 for 
c storm-relative flow
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
c find the indices nearest 1 and 4 km AGL
	  k1 = 0
	  k4 = 0
	  k9 = 0
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. 9000.) .and. 
     &             (k9 .eq. 0)) then
	      k9 = k
	      go to 8
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 4000.) .and. 
     &           (k4 .eq. 0)) then
	      k4 = k
	      go to 6
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 1000.) .and. 
     &           (k1 .eq. 0)) k1 = k
    6     continue
    8     continue
	  if (k9 .eq. 0) k9=2
	  if (k4 .eq. 0) k4=2
c calculate a 0-1 km AGL mean wind for storm motion
	  do k = mkzh, k1, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
	    su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
	    sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sdh = ',sdh,' su = ',su,' sv = ',sv
c           write(6,*) 'ter = ',ter(i,j),' k1 = ',k1,' ht = ',
c    &       ght(i,j,k1),' k4 = ',k4,' ht = ',ght(i,j,k4)
c         endif
	  ua = su / sdh
	  va = sv / sdh
c Now subtract out the 4km wind and take 60% of it
          ub = .6 * (u(i,j,k4) - ua)
	  vb = .6 * (v(i,j,k4) - va)
c         asp = sqrt(ub*ub + vb*vb)
	  if (ub .eq. 0. .and. vb .eq. 0.) then
	    adr = 0.
	  else
	    adr = dpr * (pi + atan2(ub,vb))
	  endif
c Now find the point that's 8 m/s orthogonal to the BL-4km shear
	  bsp = 8.
	  bdr = adr + 90.
	  if (bdr .gt. 360.) bdr = bdr-360.
	  cu = -bsp * sin(bdr*dtr)
	  cv = -bsp * cos(bdr*dtr)
C Here's our motion vector
	  cu = ua + ub + cu
	  cv = va + vb + cv
          if (ilow .eq. 0) then
	    if (iuv .eq. 0) then
	      g(i,j) = cu
	    else
	      g(i,j) = cv
	    endif
	    go to 15
          endif
c         sum = 0.
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'i = ',i,' j = ',j,' k1 = ',k1,
c    & ' k9 = ',k9
c          write(6,*) 'asp = ',asp,' adr = ',adr
c          write(6,*) 'bsp = ',bsp,' bdr = ',bdr
c         endif
	  sdh = 0.
	  su = 0.
	  sv = 0.
	  do k = k9, k9-1, -1
	    dh = ght(i,j,k) - ght(i,j,k+1)
	    sdh = sdh + dh
	    su = su + 0.5 * dh * 
     &              ((u(i,j,k)-cu)+(u(i,j,k+1)-cu))
	    sv = sv + 0.5 * dh * 
     &              ((v(i,j,k)-cv)+(v(i,j,k+1)-cv))
	  enddo
	  su = su / sdh
	  sv = sv / sdh
	  g(i,j) = sqrt(su*su + sv*sv)
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sreh = ',-sum
c         endif
   15 continue
      end
c--------------------------------------------------------------------
      subroutine srflo4 (u, v, ght, ter, ilow, iuv, g, miy, mjx, mkzh)
c uses modified Bunkers scheme
c Routine to compute storm-relative flow for estimating supercell
c type. 
c iuv - 0 for u, 1 for v;  ilow - 0 for storm motion, 1 for 
c storm-relative flow
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh)
     & ,ter(miy,mjx), g(miy,mjx)
      parameter (pi=3.14159265, dtr=pi/180., dpr=180./pi)
c     rewind(58)
c     read(58,*) iii,jjj
c     write(6,*) 'miy = ',miy,' mjx = ',mjx
      do 15 j = 1, mjx-1
	do 15 i = 1, miy-1
c find the indices nearest 1 and 4 km AGL
c  write(6,*) 'i = ',i,' j = ',j
	  khalf = 0
	  k55 = 0
	  k6 = 0
	  k9 = 0
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. 9000.) .and. 
     &             (k9 .eq. 0)) then
	      k9 = k
	      go to 8
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 6000.) .and. 
     &             (k6 .eq. 0)) then
	      k6 = k
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 5500.) .and. 
     &             (k55 .eq. 0)) then
	      k55 = k
	      go to 6
	    endif
	    if (((ght(i,j,k) - ter(i,j)) .gt. 500.) .and. 
     &           (khalf .eq. 0)) khalf = k
    6     continue
    8     continue
	  if (k9 .eq. 0) k9=2
	  if (k6 .eq. 0) k6=2
	  if (k55 .eq. 0) k55=2
	  if (khalf .eq. 0) khalf=2
c         write(6,*) 'k9 = ',k9,' k6 = ',k6
c         write(6,*) 'k55 = ',k55,' khalf = ',khalf
c calculate a 0-.5 km AGL mean wind for storm motion
	  sdh = 0.
	  su = 0.
	  sv = 0.
	  do k = mkzh, khalf, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
	    su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
	    sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'sdh = ',sdh,' su = ',su,' sv = ',sv
c           write(6,*) 'ter = ',ter(i,j),' k1 = ',k1,' ht = ',
c    &       ght(i,j,k1),' k4 = ',k4,' ht = ',ght(i,j,k4)
c         endif
c  write(6,*) 'khalf sdh = ',sdh
	  u1 = su / sdh
	  v1 = sv / sdh
c         write(6,*) 'u1 = ',u1,' v1 = ',v1
c calc a 5.5 - 6 km wind
	  sdh = 0.
	  su = 0.
	  sv = 0.
	  do k = k55, k6, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
	    su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
	    sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c  write(6,*) 'k55 sdh = ',sdh
	  u55 = su / sdh
	  v55 = sv / sdh
          ush = u55 - u1
	  vsh = v55 - v1
c  write(6,*) 'u55 = ',u55,' u1 = ',u1,' ush = ',ush
c calc a 0 - 6 km wind
	  sdh = 0.
	  su = 0.
	  sv = 0.
c  write(6,*) 'k6 = ',k6,' mkzh = ',mkzh
	  do k = mkzh, k6, -1
	    dh = ght(i,j,k-1) - ght(i,j,k)
	    sdh = sdh + dh
	    su = su + 0.5*dh*(u(i,j,k-1)+u(i,j,k))
	    sv = sv + 0.5*dh*(v(i,j,k-1)+v(i,j,k))
	  enddo
c  write(6,*) 'sdh = ',sdh
	  umean = su / sdh
	  vmean = sv / sdh
C Here's our motion vector
	  denom = sqrt(ush*ush + vsh*vsh)
c  write(6,*) 'denom = ',denom
	  if (denom .ne. 0.) then
	    cu = umean + (7.5/denom) * vsh
	    cv = vmean - (7.5/denom) * ush
	  else
            cu = umean
            cv = vmean
	  endif
c  write(6,*) 'cu = ',cu,' cv = ',cv
          if (ilow .eq. 0) then
	    if (iuv .eq. 0) then
	      g(i,j) = cu
	    else
	      g(i,j) = cv
	    endif
	    go to 15
          endif
c         sum = 0.
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'i = ',i,' j = ',j,' k1 = ',k1,
c    & ' k9 = ',k9
c          write(6,*) 'asp = ',asp,' adr = ',adr
c          write(6,*) 'bsp = ',bsp,' bdr = ',bdr
c         endif
	  sdh = 0.
	  su = 0.
	  sv = 0.
c  write(6,*) 'starting k9 loop'
	  do k = k9, k9-1, -1
	    dh = ght(i,j,k) - ght(i,j,k+1)
	    sdh = sdh + dh
	    su = su + 0.5 * dh * 
     &              ((u(i,j,k)-cu)+(u(i,j,k+1)-cu))
	    sv = sv + 0.5 * dh * 
     &              ((v(i,j,k)-cv)+(v(i,j,k+1)-cv))
	  enddo
c  write(6,*) 'k9 sdh = ',sdh
	  su = su / sdh
	  sv = sv / sdh
	  g(i,j) = sqrt(su*su + sv*sv)
   15 continue
      end
