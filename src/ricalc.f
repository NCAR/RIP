c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ricalc(itype,im,prs,ght,tmk,qvp,qls,
     &   uvel,vvel,rim,miy,mjx,mkzh)
c
c     This routine is similar to bvfricalc, but its intent is to exactly
c     reproduce various forms of the moist Richardson number (Rim) that
c     have been (or might be) used in the moist vertical diffusion
c     calculation in the Hi-res and MRF PBL schemes.  Shear is the
c     magnitude of the full shear vector (whereas in bvfricalc it is the
c     shear in one particular direction), and there is no preliminary
c     smoothing of winds before the shear calculation.  Also, the value
c     of Rim calculated will be at a certain k-index, but it is really
c     valid for the layer between k and k-1.
c
c     Note: qvp is vapor mix. rat., qls is sum of cloud water and ice
c     (both in kg/kg)
c
c     itype specifies various tweaks in how Ri is calculated.
c     im specifies whether you want the dry value everywhere (1),
c        the moist value everywhere (2), or dry/moist depending on
c        conditions (3)
c
      dimension prs(miy,mjx,mkzh), ght(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),qls(miy,mjx,mkzh),
     &   uvel(miy,mjx,mkzh),
     &   vvel(miy,mjx,mkzh),rim(miy,mjx,mkzh)
c
      dimension ux(300),vx(300),thvx(300),govrth(300),qs(300),
     &   thil(300)
c
      include 'comconst'
c
      EP1=0.608
      XLV=2.5E6
      RV=461.5
c
      do j=1,mjx-1
      do i=1,miy-1
c
c     Get winds at cross points, and virtual pot. temp.
c
      do k=mkzh,1,-1
         ux(k)=.25*(uvel(i,j,k)+uvel(i+1,j,k)+
     &              uvel(i,j+1,k)+uvel(i+1,j+1,k))
         vx(k)=.25*(vvel(i,j,k)+vvel(i+1,j,k)+
     &              vvel(i,j+1,k)+vvel(i+1,j+1,k))
         thcon=(1000./prs(i,j,k))**gamma
         thx=tmk(i,j,k)*thcon
         tvcon=(1.+ep1*qvp(i,j,k))
         thvx(k)=thx*tvcon
         if (k.eq.mkzh) thxbot=thx
         if (itype.eq.2) then
            govrth(k)=grav/thx
         else
            govrth(k)=grav/thxbot
         endif
         es = ezero * exp( eslcon1*(tmk(i,j,k)-celkel)/
     &      (tmk(i,j,k)-eslcon2) )
         qs(k) = eps*es/(prs(i,j,k)-es)
c         thil(k)=thx * (1. - (2.5e6*qls(i,j,k)+2.834e6*qfrz)
         thil(k)=thx * (1. - (2.5e6*qls(i,j,k))
     +                        /(cp*amax1(tmk(i,j,k),253.)))
      enddo
c
      do k=2,mkzh
         dza=ght(i,j,k-1)-ght(i,j,k)
         ss=((ux(k-1)-ux(k))*(ux(k-1)-ux(k))+(vx(k-1)-
     +      vx(k))*(vx(k-1)-vx(k)))/(dza*dza)+    1.e-9
         if (itype.eq.9) then
            rim(i,j,k)=grav/thil(k)*(thil(k-1)-thil(k))/(ss*dza)
            goto 40
         endif
         rid=govrth(k)*(thvx(k-1)-thvx(k))/(ss*dza)
         rhu=100.*qvp(i,j,k)/qs(k)
         if ( (im.eq.2).or.
     &        (im.eq.3.and.
     &           ((itype.ne.5.and.qls(i,j,k).gt.0.01e-3).or.
     &            (itype.eq.5.and.rhu.gt.90.))                 )
     &                   ) then
            qmean=0.5*(qvp(i,j,k)+qvp(i,j,k-1))
            tmean=0.5*(tmk(i,j,k)+tmk(i,j,k-1))                             
            alph=xlv*qmean/rgas/tmean
            chi=xlv*xlv*qmean/cp/rv/tmean/tmean
            if (itype.eq.7.or.itype.eq.8) then
               dqsdz=(qs(k-1)-qs(k))/dza
               if (itype.eq.7) then
                  dqtdz=dqsdz+(qls(i,j,k-1)-qls(i,j,k))/dza
               else
                  dqtdz=dqsdz
               endif
               rim(i,j,k)=grav/ss*((1.+alph)/(1.+chi)*
     &            (ss*rid/grav+xlv/cp/tmk(i,j,k)*dqsdz)-dqtdz)
            else
               rim(i,j,k)=(1.+alph)*
     &            (rid-grav*grav/ss/tmean/cp*((chi-alph)/(1.+chi)))
            endif
         else
            rim(i,j,k)=rid
         endif
 40      continue
c
      enddo
c
      enddo
      enddo
c
      return
      end
