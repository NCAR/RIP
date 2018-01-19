c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pvocalc(xmap,uuu,vvv,cor,theta,prs,pv,
     &   miy,mjx,mkzh)
c
      dimension xmap(miy,mjx),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),cor(miy,mjx),theta(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),pv(miy,mjx,mkzh)
c
      include 'comconst'
c
      icn=0
      do 200 k=1,mkzh
c
      kp1=min(k+1,mkzh)
      km1=max(k-1,1)
         do 100 j=1,mjx-1
         jp1=min(j+1,mjx-1)
         jm1=max(j-1,1)
         do 100 i=1,miy-1
         ip1=min(i+1,miy-1)
         im1=max(i-1,1)
            dxtr=ds/xmap(i,j)
            if (nproj.ne.4) then
               dytr=ds/xmap(i,j)
            else
               dytr=ds
            endif
            dudy=.5*(uuu(i+1,j,k)+uuu(i+1,j+1,k)-
     &               uuu(i,j,k)-uuu(i,j+1,k))/dytr
            dvdx=.5*(vvv(i+1,j+1,k)+vvv(i,j+1,k)-
     &               vvv(i+1,j,k)-vvv(i,j,k))/dxtr
            avort=dvdx-dudy+cor(i,j)
            dp=prs(i,j,kp1)-prs(i,j,km1)
            dudp=.25*(uuu(i,j,kp1)+uuu(i+1,j,kp1)+
     &                uuu(i,j+1,kp1)+uuu(i+1,j+1,kp1)-
     &                uuu(i,j,km1)-uuu(i+1,j,km1)-
     &                uuu(i,j+1,km1)-uuu(i+1,j+1,km1))/dp
            dvdp=.25*(vvv(i,j,kp1)+vvv(i+1,j,kp1)+
     &                vvv(i,j+1,kp1)+vvv(i+1,j+1,kp1)-
     &                vvv(i,j,km1)-vvv(i+1,j,km1)-
     &                vvv(i,j+1,km1)-vvv(i+1,j+1,km1))/dp
            dthdp=(theta(i,j,kp1)-theta(i,j,km1))/dp
            dx=dxtr*(jp1-jm1)
            dy=dytr*(ip1-im1)
            dthdx=(theta(i,jp1,k)-theta(i,jm1,k))/dx
            dthdy=(theta(ip1,j,k)-theta(im1,j,k))/dy
            pvort=-grav*(dthdp*avort-dvdp*dthdx+dudp*dthdy)
c
c         Convert to PVU: 1.e-2 to go from "per hPa" to "per Pa" (SI),
c         and 1.e6 to go from SI to PVU
c
            pv(i,j,k)=pvort*1.e4
c
  100    continue
  200 continue
      return
      end
