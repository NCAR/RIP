c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine contrive2(cor,dmap,xmap,ter,prs,ght,
     &   tmk,qvp,uuu,vvv,www,sfp,sfpsm,miy,mjx,mkzh)
c
c   This routine generates artificial model output fields for testing,
c   etc.
c
      parameter(mkzslab=300)
c
      dimension cor(miy,mjx),dmap(miy,mjx),xmap(miy,mjx),ter(miy,mjx),
     &   prs(miy,mjx,mkzh),sfp(miy,mjx),sfpsm(miy,mjx),
     &   qvp(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),tmk(miy,mjx,mkzh),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),www(miy,mjx,mkzh)
      dimension refhp(miy,mjx),refhu(miy,mjx),refhv(miy,mjx),
     &   refhth(miy,mjx)
      dimension slabth(miy,mkzslab),work(miy,mkzslab)
      character rootname*256
      real m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11
c
      include 'comconst'
c
      rewind(iunewdom)
      do i=1,11
         read(iunewdom,*)
      enddo
      read(iunewdom,*)ztop             ! top of model domain, km
      read(iunewdom,*)zjet             ! ht of jet, km
      read(iunewdom,*)thetajet         ! theta at jet, K
      read(iunewdom,*)dthetadzt        ! d(theta)/dz in trop, K/km
      read(iunewdom,*)dthetadzs        ! d(theta)/dz in strat, K/km
      read(iunewdom,*)dthetadyf        ! d(theta)/dy across front, K/km
      read(iunewdom,*)dthetadx         ! d(theta)/dx (along-fr grad), K/km
      read(iunewdom,*)frontslope       ! (dz/dy)_front
      read(iunewdom,*)frontwidth       ! frontal width, km
      read(iunewdom,*)sttpshgt         ! height of ST tropop, km
      read(iunewdom,*)rgtpsslope       ! dz/dy of rev grad part of tropop
      read(iunewdom,*)strattranslope   ! dz/dy of stratos trans sfcs
      read(iunewdom,*)zstabtran        ! hieght of lower trop stab tran, km
      read(iunewdom,*)dthetadzcold     ! d(th)/dz in low trop cold air, K/km
      read(iunewdom,*)dthetadzwarm     ! d(th)/dz in low trop warm air, K/km
      read(iunewdom,*)iycenter         ! i(y) index of jet max
      read(iunewdom,*)nsmooth          ! number of smoothing passes for theta
      read(iunewdom,*)def              ! deformation, in s**-1
      read(iunewdom,*)refh             ! reference height, km
      read(iunewdom,*)refhpmid         ! pressure at refh at mid-domain, hPa
      read(iunewdom,*)vagfac           ! factor for gues at V_ag (~1-5)
c
c   First set up terrain, xmap, dmap, cor
c
      do j=1,mjx-1
      do i=1,miy-1
c         if (i.ge.33.and.i.le.58) then
c            ter(i,j)=800.*(sin(pi*(i-33.)/25.))**2
c         else
            ter(i,j)=0.
c         endif
c         if (i.ge.32.and.i.le.45.and.j.ge.34.and.j.le.59) then
c            ter(i,j)=max(ter(i,j),800.*(sin(pi*(j-34.)/25.))**2)
c         endif
c         if (i.ge.20.and.i.le.32.and.j.ge.34.and.j.le.59) then
c            ter(i,j)=800.*(sin(pi*(j-34.)/25.))**2*
c     &                    (sin(pi*(i-19.5)/25.))**2
c         endif
         xmap(i,j)=1.
         dmap(i,j)=1.
         cor(i,j)=1.e-4
      enddo
      enddo
      call xtodot(dmap,miy,mjx)
c
c   Next, calculate theta.  First, determine line equations.  Assume
c   y-axis points toward colder air, and x-axis points in direction
c   of main tropopause-level jet.
c
c   1: southern frontal boundary
c
      m1=frontslope
      b1=0.
      dthetadzf=dthetadzt-dthetadyf/frontslope
      if (dthetadzf.ge.dthetadzs) then
         write(iup,*)'dthetadzf,dthetadzs=',dthetadzf,dthetadzs
         write(iup,*)'Increase strat. stability.'
         stop
      endif
c
c   2: northern frontal boundary
c
      m2=frontslope
      b2=-frontslope*frontwidth
c
c   3: stratosphere/front transition
c
      m3=dthetadyf/(dthetadzs-dthetadzf)
      b3=0.
c
c   4: polar tropopause (flat)
c
      m4=0.
      b4=m3*b2/(m3-m2)
c
c   6: tropopause beneath reverse-gradient portion of stratosphere
c
      if (rgtpsslope.ge.0.) then
         write(iup,*)'rgtpsslope must be .lt. 0'
         stop
      endif
      m6=rgtpsslope
      b6=0.
c
c   7: flat subtropical tropopause
c
      m7=0.
      b7=sttpshgt-zjet
c
c   5: northern stratospheric transition
c
      if (strattranslope.le.0.) then
         write(iup,*)'strattranslope must be .gt. 0'
         stop
      endif
      m5=strattranslope
      dthetadyrg=(dthetadzt-dthetadzs)/
     &   (1./m6-1./m5)
      dthetadzrg=dthetadzt-dthetadyrg/m6
      b5=0.
c
c   8: southern stratospheric transition
c
      m8=strattranslope
      b8=b7-m8*(b7-b6)/m6
c
c   9: back of lower tropospheric frontal zone
c
      zstabtranrel=zstabtran-zjet
      m9=dthetadyf/(dthetadzcold-dthetadzf)
      b9=zstabtranrel-m9*((zstabtranrel-b2)/m2)
c
c   10: front of lower tropospheric frontal zone
c
      m10=-dthetadyf/(dthetadzf-dthetadzwarm)
      b10=zstabtranrel-m10*((zstabtranrel-b1)/m1)
c
c   11: top of lower-tropospheric layer
c
      m11=0.
      b11=zstabtran-zjet
c
c   Grid stuff
c
      dy=dskm
      dx=dskm
c
      mkzslabuse=3*mkzh
      if (mkzslabuse.gt.mkzslab) then
         write(iup,*)'Increase mkzslab to at least ',mkzslabuse
         stop
      endif
      ztopslab=ztop+2.
      dz=ztopslab/(mkzslabuse-1)
      zgrel=-zjet
c
c   Calculate theta
c
      do i=1,miy-1
         y=(i-iycenter)*dy
         z1=m1*y+b1
         z2=m2*y+b2
         z3=m3*y+b3
         z4=m4*y+b4
         z5=m5*y+b5
         z6=m6*y+b6
         z7=m7*y+b7
         z8=m8*y+b8
         z9=m9*y+b9
         z10=m10*y+b10
         z11=m11*y+b11
         do k=1,mkzslabuse
            z=zgrel+(k-1)*dz
            y1=(z-b1)/m1
            y2=(z-b2)/m2
c            y4=(z-b4)/m4
            y3=(z-b3)/m3
            y5=(z-b5)/m5
            y6=(z-b6)/m6
c            y7=(z-b7)/m7
            y8=(z-b8)/m8
            y9=(z-b9)/m9
            y10=(z-b10)/m10
c            y11=(z-b10)/m11
            slabth(i,k)=100.  ! so it will be obvious if we missed anything
            if (y.ge.y3.and.y.ge.y5.and.z.ge.z4) then     ! polar strat
               slabth(i,k)=thetajet+z*dthetadzs
            elseif (y.le.y5.and.y.ge.y8.and.z.ge.z6) then ! rev grad strat
               slabth(i,k)=thetajet+z*dthetadzrg+y*dthetadyrg
            elseif (y.le.y6.and.y.le.y1.and.
     &              z.ge.z11.and.z.le.z7) then            ! ST trop
               slabth(i,k)=thetajet+z*dthetadzt
            elseif ((y.le.y2.and.y.ge.y1.and.
     &               z.le.z3.and.z.ge.z11).or.
     &              (y.le.y9.and.y.ge.y10.and.z.le.z11)) then ! trop FZ
               slabth(i,k)=thetajet+z*dthetadzf+y*dthetadyf
            elseif (y.ge.y2.and.z.le.z4.and.z.ge.z11) then ! polar trop
               slabth(i,k)=thetajet+z4*dthetadzs+(z-z4)*dthetadzt
            elseif (y.ge.y9.and.y.ge.y10.and.z.le.z11) then  ! lower trop cold
               slabth(i,k)=thetajet+z4*dthetadzs+(z11-z4)*dthetadzt+
     &                     (z-z11)*dthetadzcold
            elseif (y.le.y10.and.y.le.y9.and.z.le.z11) then  ! lower trop warm
               slabth(i,k)=thetajet+z11*dthetadzt+
     &                     (z-z11)*dthetadzwarm
            elseif (y.ge.y9.and.y.le.y10.and.z.le.z11) then ! trop FZ, rev grad
               slabth(i,k)=thetajet+z4*dthetadzs+(z11-z4)*dthetadzt+
     &                     (z-z11)*dthetadzcold-(y-y10)*dthetadyf
            elseif (y.le.y8.and.z.ge.z7) then        ! ST stratosph
               slabth(i,k)=thetajet+z7*dthetadzt+
     &                     (z-z7)*dthetadzs
            endif
         enddo
      enddo
c
cc       yya=(b4-b2)/m2
cc       yyb=0.
cc       yyc=(b7-b6)/m6
cc       yyd=(zgrel-b2)/m2
cc       yye=(zgrel-b1)/m1
cc       zbottom=m2*min(yyc-yyd-20.,0.)
cc       if (zbottom.le.-1000.) then
cc          write(iup,*)'Warning: zbottom = ',zbottom,' km.'
cc          write(iup,*)'This is very low, and may cause inaccuracies.'
cc       endif
cc       zbottomrel=zbottom-zjet
cc       yyd=(zbottomrel-b2)/m2
cc       yye=(zbottomrel-b1)/m1
cc c
cc c   Calculate theta
cc c
cc       bottomthetacold=slthetacold+zbottom*dthetadzt
cc       do i=1,miy-1
cc          y=(i-iycenter)*dy
cc          if (y.ge.yyd) then
cc             bottomtheta=bottomthetacold
cc          elseif (y.le.yyd.and.y.ge.yye) then
cc             bottomtheta=bottomthetacold+dthetadyf*(y-yyd)
cc          else
cc             bottomtheta=bottomthetacold+dthetadyf*(yye-yyd)
cc          endif
cc          z1=m1*y+b1
cc          z2=m2*y+b2
cc          z4=m4*y+b4
cc          z3=m3*y+b3
cc          z5=m5*y+b5
cc          z6=m6*y+b6
cc          z7=m7*y+b7
cc          z8=m8*y+b8
cc          do k=1,mkzslabuse
cc             z=zgrel+(k-1)*dz
cc             if (y.ge.yya) then
cc                if (z.le.z4) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzt
cc                elseif (z.le.z5) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z4-zbottomrel)*dthetadzt+
cc      &               (z-z4)*dthetadzs
cc                elseif (z.le.z8) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z4-zbottomrel)*dthetadzt+
cc      &               (z5-z4)*dthetadzs+
cc      &               (z-z5)*dthetadzrg
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z4-zbottomrel)*dthetadzt+
cc      &               (z5-z4)*dthetadzs+
cc      &               (z8-z5)*dthetadzrg+
cc      &               (z-z8)*dthetadzs
cc                endif
cc             elseif (y.ge.yyb) then
cc                if (z.le.z2) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzt
cc                elseif (z.le.z3) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z-z2)*dthetadzf
cc                elseif (z.le.z5) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z3-z2)*dthetadzf+
cc      &               (z-z3)*dthetadzs
cc                elseif (z.le.z8) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z3-z2)*dthetadzf+
cc      &               (z5-z3)*dthetadzs+
cc      &               (z-z5)*dthetadzrg
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z3-z2)*dthetadzf+
cc      &               (z5-z3)*dthetadzs+
cc      &               (z8-z5)*dthetadzrg+
cc      &               (z-z8)*dthetadzs
cc                endif
cc             elseif (y.ge.yyc) then
cc                if (z.le.z2) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzt
cc                elseif (z.le.z1) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z-z2)*dthetadzf
cc                elseif (z.le.z6) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z1-z2)*dthetadzf+
cc      &               (z-z1)*dthetadzt
cc                elseif (z.le.z8) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z1-z2)*dthetadzf+
cc      &               (z6-z1)*dthetadzt+
cc      &               (z-z6)*dthetadzrg
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z1-z2)*dthetadzf+
cc      &               (z6-z1)*dthetadzt+
cc      &               (z8-z6)*dthetadzrg+
cc      &               (z-z8)*dthetadzs
cc                endif
cc             elseif (y.ge.yyd) then
cc                if (z.le.z2) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzt
cc                elseif (z.le.z1) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z-z2)*dthetadzf
cc                elseif (z.le.z7) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z1-z2)*dthetadzf+
cc      &               (z-z1)*dthetadzt
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z2-zbottomrel)*dthetadzt+
cc      &               (z1-z2)*dthetadzf+
cc      &               (z7-z1)*dthetadzt+
cc      &               (z-z7)*dthetadzs
cc                endif
cc             elseif (y.ge.yye) then
cc                if (z.le.z1) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzf
cc                elseif (z.le.z7) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z1-zbottomrel)*dthetadzf+
cc      &               (z-z1)*dthetadzt
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z1-zbottomrel)*dthetadzf+
cc      &               (z7-z1)*dthetadzt+
cc      &               (z-z7)*dthetadzs
cc                endif
cc             else
cc                if (z.le.z7) then
cc                   slabth(i,k)=bottomtheta+
cc      &               (z-zbottomrel)*dthetadzt
cc                else
cc                   slabth(i,k)=bottomtheta+
cc      &               (z7-zbottomrel)*dthetadzt+
cc      &               (z-z7)*dthetadzs
cc                endif
cc             endif
cc          enddo
cc       enddo
c
      call smooth(slabth,work,nsmooth,miy,miy-1,mkzslabuse)
c
c   Next define ght.  Based on a sigma-z coordinate, with evenly thick
c   layers, and each model level is at the midpoint of a layer
c
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         ght(i,j,k)=(mkzh-float(k)+0.5)/float(mkzh)*
     &      (1000.*ztop-ter(i,j))+ter(i,j)
      enddo
      enddo
      enddo
c
      rimidx=.5*miy
      rjmidx=.5*mjx
      imidx=nint(rimidx)
      jmidx=nint(rjmidx)
      rimidd=.5*(miy+1)
      rjmidd=.5*(mjx+1)
      imidd=nint(rimidd)
      jmidd=nint(rjmidd)
c
c   Interpolate theta to model levels, and include x-gradient of theta
c
      do k=1,mkzh
      do j=1,mjx-1
         dtheta=(j-jmidx)*dx*dthetadx
      do i=1,miy-1
         zmod=ght(i,j,k)/1000.
         if (zmod.le.0.) then
            tmk(i,j,k)=slabth(i,1)+dtheta
         elseif (zmod.ge.ztop) then
            tmk(i,j,k)=slabth(i,mkzslabuse)+dtheta
         else
            do kslab=1,mkzslabuse-1
               z=(kslab-1)*dz
               z1=z+dz
               if (zmod.ge.z.and.zmod.le.z1) then
                  tmk(i,j,k)=((zmod-z)*slabth(i,kslab+1)+
     &                        (z1-zmod)*slabth(i,kslab))/dz+dtheta
                  goto 43
               endif
            enddo
            write(iup,*)'couldn''t interpolate contrived theta.'
            write(iup,*)'i,j,k=',i,j,k
            stop
 43         continue
         endif
      enddo
      enddo
      enddo
c
c   Calculate theta at reference height
c
      refhold=refh
      refhpmidold=refhpmid
      refh=dz*nint(refh/dz)
      refhpmid=refhpmidold*
     &   exp(-grav*1000.*(refh-refhold)/(rgas*celkel))
c      write(iup,*)'refh,refhpmid adjusted to ',refh,refhpmid
      refhm=1000.*refh
      kslabrefh=1+nint(refh/dz)
      do j=1,mjx-1
         dtheta=(j-jmidx)*dx*dthetadx
         do i=1,miy-1
            refhth(i,j)=slabth(i,kslabrefh)+dtheta
         enddo
      enddo
c
c   Next, make wind at reference height level, which is just
c   a geostrophic deformation pattern.
c
      do j = 1, mjx
         x=(j-rjmidd)*ds
      do i = 1, miy
         y=(i-rimidd)*ds
         refhu(i,j)=.5*def*x
         refhv(i,j)=-.5*def*y
      enddo
      enddo
c
c   Next, make pressure at reference height
c
      do j = 1, mjx-1
      do i = 1, miy-1
         refhp(i,j)=refhpmid ! first estimate of refhp
      enddo
      enddo
      do iter=1,5 ! iterative refinement of refhp
         do j=2,mjx-1
            avgcor=.5*(cor(1,j-1)+cor(1,j))
            avgrefhv=.5*(refhv(1,j)+refhv(2,j))
            avgrefhp=.5*(refhp(1,j-1)+refhp(1,j))
            avgrefht=.5*(refhth(1,j-1)+refhth(1,j))*
     &         (1000./avgrefhp)**gamma
            refhp(1,j)=refhp(1,j-1)*
     &         exp(avgcor*avgrefhv/(rgas*avgrefht)*ds)
         enddo
         do j=1,mjx-1
         do i=2,miy-1
            avgcor=.5*(cor(i-1,j)+cor(i,j))
            avgrefhu=.5*(refhu(i,j)+refhu(i,j+1))
            avgrefhp=.5*(refhp(i-1,j)+refhp(i,j))
            avgrefht=.5*(refhth(i-1,j)+refhth(i,j))*
     &         (1000./avgrefhp)**gamma
            refhp(i,j)=refhp(i-1,j)*
     &         exp(-avgcor*avgrefhu/(rgas*avgrefht)*ds)
         enddo
         enddo
         refhperrormid=refhp(imidx,jmidx)-refhpmid
         do j=1,mjx-1
         do i=1,miy-1
            refhp(i,j)=refhp(i,j)-refhperrormid
         enddo
         enddo
      enddo
c
c   Calculate pressure at all levels
c
      const=grav*gamma*1000.**gamma/rgas
      do j = 1, mjx-1
      do i = 1, miy-1
c
c   First get two sigma levels surrounding reference level
c
      do k=1,mkzh-1
         if (refhm.le.ght(i,j,k).and.refhm.ge.ght(i,j,k+1)) then
            kabove=k
            kbelow=k+1
            goto 344
         endif
      enddo
      write(iup,*)'couldn''t find two sigma levels surrounding'
      write(iup,*)'reference level, i,j=',i,j
      stop
 344  continue
c
c   Get p at levels above
c
      prs(i,j,kabove)=(refhp(i,j)**gamma-2.*const*
     &   (ght(i,j,kabove)-refhm)/
     &   (tmk(i,j,kabove)+refhth(i,j)))**(1./gamma)
      do k=kabove-1,1,-1
         prs(i,j,k)=(prs(i,j,k+1)**gamma-2.*const*
     &      (ght(i,j,k)-ght(i,j,k+1))/
     &      (tmk(i,j,k)+tmk(i,j,k+1)))**(1./gamma)
      enddo
c
c   Get p at levels below
c      
      prs(i,j,kbelow)=(refhp(i,j)**gamma-2.*const*
     &   (ght(i,j,kbelow)-refhm)/
     &   (tmk(i,j,kbelow)+refhth(i,j)))**(1./gamma)
      do k=kbelow+1,mkzh
         prs(i,j,k)=(prs(i,j,k-1)**gamma-2.*const*
     &      (ght(i,j,k)-ght(i,j,k-1))/
     &      (tmk(i,j,k)+tmk(i,j,k-1)))**(1./gamma)
      enddo
      sfp(i,j)=(prs(i,j,mkzh)**gamma-const*
     &      (ter(i,j)-ght(i,j,mkzh))/
     &      tmk(i,j,mkzh))**(1./gamma)
      sfpsm(i,j)=sfp(i,j)
c
      enddo
      enddo
c
c   Convert theta to temperature, set www and qvp to 0.
c
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         tmk(i,j,k)=tmk(i,j,k)*(prs(i,j,k)/1000.)**gamma
         www(i,j,k)=0.
         qvp(i,j,k)=0.
      enddo
      enddo
      enddo
c
c   Calculate geostrophic winds aloft.  Use refhu and refhv for scratch
c   space, since they're not needed anymore.
c
      call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &   refhu,refhv,uuu,1,'y',miy,mjx,mkzh)
      call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &   refhu,refhv,vvv,1,'x',miy,mjx,mkzh)
      fbar=cor(miy/2,mjx/2)
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            tv=virtual(tmk(i,j,k),qvp(i,j,k))
            uuu(i,j,k)=-rgas*tv*uuu(i,j,k)/(prs(i,j,k)*fbar)
            vag=vagfac*(-(prs(i,j,k)-600.)/500.)
            vvv(i,j,k)= rgas*tv*vvv(i,j,k)/(prs(i,j,k)*fbar)+vag
         enddo
         enddo
         call xtodot(uuu(1,1,k),miy,mjx)
         call xtodot(vvv(1,1,k),miy,mjx)
      enddo
c
      return
      end
