c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine wdircalc(uvel,vvel,unorth,vnorth,rstrm,wdr,
     &   miy,mjx,mkzh)
c
      dimension uvel(miy,mjx,mkzh), vvel(miy,mjx,mkzh),
     &    unorth(miy,mjx), vnorth(miy,mjx),rstrm(2),
     &    wdr(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 200 j=1,mjx
      do 200 i=1,miy
         rnorthdir=atan2(vnorth(i,j),unorth(i,j))
      do 200 k=1,mkzh
         wdir=(rnorthdir-atan2((vvel(i,j,k)-rstrm(1)),
     &      (uvel(i,j,k)-rstrm(2))))/rpd
         wdr(i,j,k)=mod(wdir+900.,360.)
  200 continue
      return
      end
