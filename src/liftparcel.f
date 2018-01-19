c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine liftparcel(prs,tmk,qvp,ght,tlift,kpar,ivr,
     &   miy,mjx,mkzh)
c
c   This routine calculates the temperature (deg. C) of parcels
c   lifted from a specified vertical level, for all points at or above
c   that vertical level.
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),tlift(miy,mjx,mkzh)
c
      include 'comconst'
c
      do j=1,mjx-1
      do i=1,miy-1
c
      do k=kpar+1,mkzh
         tlift(i,j,k)=rmsg
      enddo
      if (ivr.eq.0) then
        tlift(i,j,kpar)=tmk(i,j,kpar)-celkel
      elseif (ivr.eq.1) then
        tlift(i,j,kpar)=virtual(tmk(i,j,kpar),qvp(i,j,kpar))-celkel
      endif
c
c   Calculate temperature and moisture properties of parcel
c
      qvppari=qvp(i,j,kpar)
      tmkpari=tmk(i,j,kpar)
      prspari=prs(i,j,kpar)
      ghtpari=ght(i,j,kpar)
      gammam=gamma*(1.+gammamd*qvppari)
      cpm=cp*(1.+cpmd*qvppari)
c
      e=max(1.e-20,qvppari*prspari/(eps+qvppari))
      tlcl=tlclc1/(log(tmkpari**tlclc2/e)-tlclc3)+tlclc4
      ethpari=tmkpari*(1000./prspari)**(gamma*(1.+gammamd*qvppari))*
     &   exp((thtecon1/tlcl-thtecon2)*qvppari*
     &   (1.+thtecon3*qvppari))
      zlcl=ghtpari+(tmkpari-tlcl)/(grav/cpm)
c
c   Calculate temperature of lifted parcel at all levels.
c
      do k=kpar-1,1,-1
         if (ght(i,j,k).lt.zlcl) then
            tmklift=tmkpari-grav/cpm*(ght(i,j,k)-ghtpari)
            if (ivr.eq.0) then
               tlift(i,j,k)=tmklift-celkel
            elseif (ivr.eq.1) then
               qvplift=qvppari
               tlift(i,j,k)=virtual(tmklift,qvplift)-celkel
            endif
         else
            tmklift=tonpsadiabat(ethpari,prs(i,j,k))
            if (ivr.eq.0) then
               tlift(i,j,k)=tmklift-celkel
            elseif (ivr.eq.1) then
               eslift=ezero*exp(eslcon1*(tmklift-celkel)/
     &            (tmklift-eslcon2))
               qvplift=eps*eslift/(prs(i,j,k)-eslift)
               tlift(i,j,k)=virtual(tmklift,qvplift)-celkel
            endif
         endif
      enddo
c
      enddo
      enddo
c
      return
      end
