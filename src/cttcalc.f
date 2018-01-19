c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cttcalc(prs,pf,tmk,qcw,qci,ctt,
     &   miy,mjx,mkzh)
c
c   This routine calculates cloud top brightness temperatures in Cels.
c   It assumes:
c      1) zenith angle is 0
c      2) brightness temperature is roughly the temperature at
c         unity optical depth into the cloud
c      3) cloud absorption coefficient is constant
c
      dimension prs(miy,mjx,mkzh),pf(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),qcw(miy,mjx,mkzh),ctt(miy,mjx),
     &   qci(miy,mjx,mkzh)
c
      include 'comconst'
c
c   Create cloud-top temperature.
c
      do 190 j=1,mjx-1
      do 190 i=1,miy-1
         opdepthd=0.
         k=0
c
c      Integrate downward from model top, calculating path at full
c      model vertical levels.
c
   20    opdepthu=opdepthd
         k=k+1
         if (k.eq.1) then
            dp=200.*(pf(i,j,1)-prs(i,j,1))  ! should be in Pa
         else
            dp=100.*(pf(i,j,k)-pf(i,j,k-1))  ! should be in Pa
         endif
         if (iice.eq.0) then
            if (tmk(i,j,k).lt.celkel) then
c             Note: abscoefi is m**2/g, qcw is g/kg,
c                   so no convrsion needed
               opdepthd=opdepthu+abscoefi*qcw(i,j,k)*dp/grav
            else
               opdepthd=opdepthu+abscoef*qcw(i,j,k)*dp/grav
            endif
         else
            opdepthd=opdepthd+(abscoef*qcw(i,j,k)+
     &                        abscoefi*qci(i,j,k))*dp/grav
         endif
         if (opdepthd.lt.1..and.k.lt.mkzh) then
            goto 20
         elseif (opdepthd.lt.1..and.k.eq.mkzh) then
            prsctt=prs(i,j,mkzh)
         else
            fac=(1.-opdepthu)/(opdepthd-opdepthu)
            prsctt=pf(i,j,k-1)+fac*(pf(i,j,k)-pf(i,j,k-1))
            prsctt=min(prs(i,j,mkzh),max(prs(i,j,1),prsctt))
         endif
         do 30 k=2,mkzh
            if (prsctt.ge.prs(i,j,k-1).and.prsctt.le.prs(i,j,k)) then
               fac=(prsctt-prs(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
               ctt(i,j)=tmk(i,j,k-1)+
     &            fac*(tmk(i,j,k)-tmk(i,j,k-1))-celkel
               goto 40
            endif
   30    continue
   40    continue
 190  continue
      return
      end




