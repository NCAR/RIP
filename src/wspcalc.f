c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine wspcalc(uvel,vvel,cosa,sina,rstrm,wsp,miy,mjx,mkzh)
c
      dimension uvel(miy,mjx,mkzh), vvel(miy,mjx,mkzh),
     &    wsp(miy,mjx,mkzh),rstrm(2)
c
      include 'comconst'
c
      do 200 j=1,mjx
      do 200 i=1,miy
      do 200 k=1,mkzh
         if (cosa.eq.0..and.sina.eq.0.) then
            wsp(i,j,k)=sqrt((uvel(i,j,k)-rstrm(2))**2+
     &         (vvel(i,j,k)-rstrm(1))**2)
         else
            wsp(i,j,k)=cosa*(uvel(i,j,k)-rstrm(2))+
     &         sina*(vvel(i,j,k)-rstrm(1))
         endif
  200 continue
      return
      end
