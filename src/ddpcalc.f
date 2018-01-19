c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ddpcalc(prs,arr,ddp,miy,mjx,mkzh)
c
      dimension prs(miy,mjx,mkzh),arr(miy,mjx,mkzh),
     &   ddp(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 200 k=1,mkzh
         kp1=min(k+1,mkzh)
         km1=max(k-1,1)
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
         ddp(i,j,k)=(arr(i,j,kp1)-arr(i,j,km1))/
     &              (prs(i,j,kp1)-prs(i,j,km1))
  200 continue
      return
      end
