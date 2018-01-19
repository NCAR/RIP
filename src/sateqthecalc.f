c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine sateqthecalc(tmk,prs,sateth,miy,mjx,mkzh)
c
      dimension tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),sateth(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
            t=tmk(i,j,k)
            p=prs(i,j,k)
            esat=ezero*exp(eslcon1*(t-celkel)/
     &         (t-eslcon2))
            qsat=eps*esat/(p-esat)
            sateth(i,j,k)=t*(1000./p)**(gamma*(1.+gammamd*qsat))*
     &         exp((thtecon1/t-thtecon2)*qsat*(1.+thtecon3*qsat))
 1000 continue
      return
      end
