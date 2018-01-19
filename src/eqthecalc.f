c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine eqthecalc(qvp,tmk,prs,eth,miy,mjx,mkzh)
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),eth(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
            q=max(qvp(i,j,k),1.e-15)
            t=tmk(i,j,k)
            p=prs(i,j,k)
            e=q*p/(eps+q)
            tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
            eth(i,j,k)=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &         exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
 1000 continue
      return
      end
