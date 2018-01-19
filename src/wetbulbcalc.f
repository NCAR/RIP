c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine wetbulbcalc(prs,tmk,qvp,twb,miy,mjx,mkzh)
c
c   This routine calculates the wet bulb temperature (in Celsius)
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   twb(miy,mjx,mkzh)
c
      include 'comconst'
c
      do 80 k=1,mkzh
      do 80 j=1,mjx-1
      do 80 i=1,miy-1
            q=max(qvp(i,j,k),1.e-15)
            t=tmk(i,j,k)
            p=prs(i,j,k)
            e=q*p/(eps+q)
            tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
            eth=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &         exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
            twb(i,j,k)=tonpsadiabat(eth,p)-celkel
  80  continue
c
      return
      end
