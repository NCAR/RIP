c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine omgcalc(qvp,tmk,www,prs,omg,miy,mjx,mkzh)
c
c   Calculate approximate omega, based on vertical velocity w (dz/dt).
c   It is approximate because it cannot take into account the vertical
c   motion of pressure surfaces.
c
      include 'comconst'
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   www(miy,mjx,mkzh),prs(miy,mjx,mkzh),omg(miy,mjx,mkzh)
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         omg(i,j,k)=-grav*prs(i,j,k)/
     &      (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))*www(i,j,k)*
     &      10. ! www is in cm/sec, want it in dPa/sec
      enddo
      enddo
      enddo
c
      return
      end
