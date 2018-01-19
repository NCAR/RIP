c--------------------------------------------------------------------
      subroutine brnshr (u, v, ght, prsf, ter, wk, miy, mjx, mkzh)
c
c   Bulk Richardson Number Shear, as found in Stensrud et al. (1997)
c
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh),
     &   prsf(miy, mjx,mkzh), ter(miy,mjx), wk(miy,mjx)
c
      dimension ucross(500),vcross(500)
c
      include 'comconst'
c
      do j = 1, mjx-1
      do i = 1, miy-1
         sdh = 0.
         su = 0.
         sv = 0.
c
c      Find the first level below 500 m AGL and the first level
c      below 6000 m AGL
c
         k500 = 0
         k6000 = 0
         do k = 1,mkzh
            if (((ght(i,j,k)-ter(i,j)).lt.6000.).and.(k6000.eq.0)) then
               k6000 = k
            endif
            if (((ght(i,j,k)-ter(i,j)).lt.500.).and.(k500.eq.0)) then
               k500=k
            endif
         enddo
         if (k6000.eq.0) then
            write(iup,*)'In brnshr, couldn''t get the first'
            write(iup,*)' vertical model level below 6000 m AGL.'
            stop
         endif
         if (k500.eq.0) then
            write(iup,*)'In brnshr, couldn''t get the first'
            write(iup,*)' vertical model level below 500 m AGL.'
            stop
         endif
c
c   Calculate u and v at cross points
c
         do k = 1,mkzh
            ucross(k)=.25*(u(i,j,k)+u(i+1,j,k)+
     &                  u(i,j+1,k)+u(i+1,j+1,k))
            vcross(k)=.25*(v(i,j,k)+v(i+1,j,k)+
     &                  v(i,j+1,k)+v(i+1,j+1,k))
         enddo
c
c      Calculate a 0-6 km AGL mean wind
c
         do k = mkzh, k6000, -1
            dh = prsf(i,j,k) - prsf(i,j,k-1)
            sdh = sdh + dh
            su = su + dh*ucross(k)
            sv = sv + dh*vcross(k)
         enddo
         ua = su / sdh
         va = sv / sdh
c
c      Calculate a 0-500 m AGL mean wind
c
         sdh = 0.
         su = 0.
         sv = 0.
         do k = mkzh, k500, -1
            dh = prsf(i,j,k) - prsf(i,j,k-1)
            sdh = sdh + dh
            su = su + dh*ucross(k)
            sv = sv + dh*vcross(k)
         enddo
         u0 = su / sdh
         v0 = sv / sdh
c
c   Bulk Richardson Number Shear:
c
         wk(i,j) = 0.5 * ( (ua - u0)**2 + (va - v0)**2 )
c
      enddo
      enddo
c
      return
      end
