c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine condheat(tmk,qvp,www,irhthtype,ilrtype,imalr,prs,
     &   miy,mjx,mkzh,dthdt,rhthresh)
c
c      *** uses temperature and vertical velocity to compute
c      ***    d/dt (theta) due to condensation from lifting,
c      ***    in kelvins per hour
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   www(miy,mjx,mkzh),dthdt(miy,mjx,mkzh)
c
      include 'comconst'
c
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         e = qvp(i,j,k)*prs(i,j,k)/(eps+qvp(i,j,k))
         xlhc=xlhc0-xlhctd*tmk(i,j,k)
         es=ezero*exp(eslcon1*(tmk(i,j,k)-celkel)/(tmk(i,j,k)-eslcon2))
         if (tmk(i,j,k).lt.celkel) then
            esi=ezero*exp(esicon1-esicon2/tmk(i,j,k))
            xlhd=xlhc+xlhf
         else
            esi=es
            xlhd=xlhc
         endif
c
         if (irhthtype.eq.0) then
            rhu=100.*(e*(prs(i,j,k)-es))/(es*(prs(i,j,k)-e))
         else
            rhu=100.*(e*(prs(i,j,k)-esi))/(esi*(prs(i,j,k)-e))
         endif
c
         if (www(i,j,k).gt.0..and.rhu.gt.rhthresh) then
            if (ilrtype.eq.0) then
               xlatent=xlhc
               ws=eps*es/(prs(i,j,k)-es)
            else
               xlatent=xlhd
               ws=eps*esi/(prs(i,j,k)-esi)
            endif
c
c         The following pseudoadiabatic lapse rate can be found in several
c         text books (Bluestein, Cotton and Anthes, Rogers and Yau, etc.).
c         It is probably not exactly consistent with Bolton's
c         relationships for theta_e and e_s, which are used more consistently
c         throughout the rip code.
c
            rmalr=(grav/cp)*(1.+xlatent*ws/(rgas*tmk(i,j,k)))/
     &         (1.+xlatent*xlatent*eps*ws/
     &            (rgas*cp*tmk(i,j,k)*tmk(i,j,k)))
            gammam=gamma*(1.+gammamd*qvp(i,j,k))
            cpm=cp*(1.+cpmd*qvp(i,j,k))
            rgasm=rgas*(1.+rgasmd*qvp(i,j,k))
            theta=tmk(i,j,k)*(1000./prs(i,j,k))**gammam
            if (imalr.eq.0) then
               dthdt(i,j,k)=www(i,j,k)*.01*(1000./prs(i,j,k))**gammam*
     &            (grav/cpm-rmalr)*3600.  ! K/s to K/hr
            else
               dthdt(i,j,k)=1000.*(1000./prs(i,j,k))**gammam*
     &            (grav/cpm-rmalr)  ! MALR (as d(theta)/dz) in K/km
            endif
         else
            dthdt(i,j,k)=0.
         endif
      enddo
      enddo
      enddo
c
      return
      end


