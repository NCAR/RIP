c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine monotonic(vc3d,prs,idir,delta,icorsw,cor,
     &   miy,mjx,mkzh)
c
c   This routine makes vc3d always increase (decrease) with height
c   above the 300 hPa level, and always decrease (increase) with
c   height below that level [if idir is +1 (-1)], by the minimum
c   increment of delta.
c
      dimension vc3d(miy,mjx,mkzh),prs(miy,mjx,mkzh),cor(miy,mjx)
c
      do j = 1, mjx-1
      do i = 1, miy-1
c
      if (icorsw.eq.1.and.cor(i,j).lt.0.) then
         do k=1,mkzh
            vc3d(i,j,k)=-vc3d(i,j,k)
         enddo
      endif
c
c   First find k index that is at or below (height-wise) the 300 hPa
c      level.
c
      do k = 1, mkzh
         if (prs(i,j,k).ge.300.) then
            k300=k
            goto 40
         endif
      enddo
c
 40   continue
c
c   Now adjust vc3d.
c
      do k = k300-1, 1, -1
         if (idir.eq.1) then
            vc3d(i,j,k)=max(vc3d(i,j,k),vc3d(i,j,k+1)+delta)
         elseif (idir.eq.-1) then
            vc3d(i,j,k)=min(vc3d(i,j,k),vc3d(i,j,k+1)-delta)
         endif
      enddo
c
      do k = k300+1, mkzh
         if (idir.eq.1) then
            vc3d(i,j,k)=min(vc3d(i,j,k),vc3d(i,j,k-1)-delta)
         elseif (idir.eq.-1) then
            vc3d(i,j,k)=max(vc3d(i,j,k),vc3d(i,j,k-1)+delta)
         endif
      enddo
c
      enddo
      enddo
      return
      end
