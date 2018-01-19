c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine capecalc3d(prs,tmk,qvp,ght,ter,prsf,cape,cin,
     &   miy,mjx,mkzh,i3dflag)
c
c   If i3dflag=1, this routine calculates CAPE and CIN (in m**2/s**2,
c   or J/kg) for every grid point in the entire 3D domain (treating
c   each grid point as a parcel).  If i3dflag=0, then it
c   calculates CAPE and CIN only for the parcel with max theta-e in
c   the column, (i.e. something akin to Colman's MCAPE).  By "parcel",
c   we mean a 500-m deep parcel, with actual temperature and moisture
c   averaged over that depth.
c
c   In the case of i3dflag=0,
c   CAPE and CIN are 2D fields that are placed in the k=mkzh slabs of
c   the cape and cin arrays.  Also, if i3dflag=0, LCL and LFC heights
c   are put in the k=mkzh-1 and k=mkzh-2 slabs of the cin array.
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),ter(miy,mjx),cape(miy,mjx,mkzh),
     &   cin(miy,mjx,mkzh),prsf(miy,mjx,mkzh)
c
      dimension buoy(150),zrel(150),benaccum(150)
c
      include 'comconst'
c
      do j=1,mjx-1
      do i=1,miy-1
         cape(i,j,1)=0.
         cin(i,j,1)=0.
c
      if (i3dflag.eq.1) then
         kpar1=2
         kpar2=mkzh
      else
c
c      Find parcel with max theta-e in lowest 3 km AGL.
c
         ethmax=-1.
         do k=mkzh,1,-1
            if (ght(i,j,k)-ter(i,j).lt.3000.) then
               q=max(qvp(i,j,k),1.e-15)
               t=tmk(i,j,k)
               p=prs(i,j,k)
               e=q*p/(eps+q)
               tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
               eth=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &            exp((thtecon1/tlcl-thtecon2)*q*(1.+thtecon3*q))
               if (eth.gt.ethmax) then
                  klev=k
                  ethmax=eth
               endif
            endif
         enddo
         kpar1=klev
         kpar2=klev
c
c      Establish average properties of that parcel
c         (over depth of approximately davg meters)
c
c         davg=.1
         davg=500.
         pavg=davg*prs(i,j,kpar1)*grav/
     &      (rgas*virtual(tmk(i,j,kpar1),qvp(i,j,kpar1)))
         p2=min(prs(i,j,kpar1)+.5*pavg,prsf(i,j,mkzh))
         p1=p2-pavg
         totthe=0.
         totqvp=0.
         totprs=0.
         do k=mkzh,2,-1
            if (prsf(i,j,k).le.p1) goto 35
            if (prsf(i,j,k-1).ge.p2) goto 34
            p=prs(i,j,k)
            pup=prsf(i,j,k)
            pdn=prsf(i,j,k-1)
            q=max(qvp(i,j,k),1.e-15)
            th=tmk(i,j,k)*(1000./p)**(gamma*(1.+gammamd*q))
            pp1=max(p1,pdn)
            pp2=min(p2,pup)
            if (pp2.gt.pp1) then
               deltap=pp2-pp1
               totqvp=totqvp+q*deltap
               totthe=totthe+th*deltap
               totprs=totprs+deltap
            endif
 34         continue
         enddo
 35      continue
         qvppari=totqvp/totprs
         tmkpari=(totthe/totprs)*(prs(i,j,kpar1)/1000.)**
     &      (gamma*(1.+gammamd*qvp(i,j,kpar1)))
      endif
c
      do kpar=kpar1,kpar2
c
c   Calculate temperature and moisture properties of parcel
c     (Note, qvppari and tmkpari already calculated above for 2D case.)
c
      if (i3dflag.eq.1) then
         qvppari=qvp(i,j,kpar)
         tmkpari=tmk(i,j,kpar)
      endif
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
c   Calculate buoyancy and relative height of lifted parcel at
c   all levels, and store in bottom up arrays.  Add a level at the LCL,
c   and at all points where buoyancy is zero.
c
      kk=0 ! for arrays that go bottom to top
      ilcl=0
      if (ghtpari.ge.zlcl) then
c
c      initial parcel already saturated or supersaturated.
c
         ilcl=2
         klcl=1
      endif
      do k=kpar,1,-1
 33      kk=kk+1                ! for arrays that go bottom to top
         if (ght(i,j,k).lt.zlcl) then ! model level is below LCL
            qvplift=qvppari
            tmklift=tmkpari-grav/cpm*(ght(i,j,k)-ghtpari)
            tvenv=virtual(tmk(i,j,k),qvp(i,j,k))
            tvlift=virtual(tmklift,qvplift)
            ghtlift=ght(i,j,k)
         elseif (ght(i,j,k).ge.zlcl.and.ilcl.eq.0) then
c
c         This model level and previous model level straddle the LCL,
c         so first create a new level in the bottom-up array, at the LCL.
c
            tmklift=tlcl
            qvplift=qvppari
            facden=ght(i,j,k)-ght(i,j,k+1)
            fac1=(zlcl-ght(i,j,k+1))/facden
            fac2=(ght(i,j,k)-zlcl)/facden
            tmkenv=tmk(i,j,k+1)*fac2+tmk(i,j,k)*fac1
            qvpenv=qvp(i,j,k+1)*fac2+qvp(i,j,k)*fac1
            tvenv=virtual(tmkenv,qvpenv)
            tvlift=virtual(tmklift,qvplift)
            ghtlift=zlcl
            ilcl=1
         else
            tmklift=tonpsadiabat(ethpari,prs(i,j,k))
            eslift=ezero*exp(eslcon1*(tmklift-celkel)/
     &         (tmklift-eslcon2))
            qvplift=eps*eslift/(prs(i,j,k)-eslift)
            tvenv=virtual(tmk(i,j,k),qvp(i,j,k))
            tvlift=virtual(tmklift,qvplift)
            ghtlift=ght(i,j,k)
         endif
         buoy(kk)=grav*(tvlift-tvenv)/tvenv  ! buoyancy
         zrel(kk)=ghtlift-ghtpari
         if (kk.gt.1) then
         if (buoy(kk)*buoy(kk-1).lt.0.0) then
c
c         Parcel ascent curve crosses sounding curve, so create a new level
c         in the bottom-up array at the crossing.
c
            kk=kk+1
            buoy(kk)=buoy(kk-1)
            zrel(kk)=zrel(kk-1)
            buoy(kk-1)=0.
            zrel(kk-1)=zrel(kk-2)+
     &         buoy(kk-2)/(buoy(kk-2)-buoy(kk))*(zrel(kk)-zrel(kk-2))
         endif
         endif
         if (ilcl.eq.1) then
            klcl=kk
            ilcl=2
            goto 33
         endif
      enddo
      kmax=kk
      if (kmax.gt.150) then
         write(iup,*)'in capecalc3d: kmax got too big. kmax=',kmax
         stop
      endif
c
c   If no LCL was found, set klcl to kmax.  It is probably not really
c   at kmax, but this will make the rest of the routine behave
c   properly.
c
      if (ilcl.eq.0) klcl=kmax
c
c   Get the accumulated buoyant energy from the parcel's starting
c   point, at all levels up to the top level.
c
      benaccum(1)=0.0
      benamin=9e9
      do k=2,kmax
         dz=zrel(k)-zrel(k-1)
         benaccum(k)=benaccum(k-1)+.5*dz*(buoy(k-1)+buoy(k))
         if (benaccum(k).lt.benamin) then
            benamin=benaccum(k)
         endif
      enddo
c
c     Determine equilibrium level (EL), which we define as the highest
c     level of non-negative buoyancy above the LCL. Note, this may be
c     the top level if the parcel is still buoyant there.
c     
      do k=kmax,klcl,-1
         if (buoy(k).ge.0.) then
            kel=k   ! k of equilibrium level
            goto 50
         endif
      enddo
c
c   If we got through that loop, then there is no non-negative
c   buoyancy above the LCL in the sounding.  In these situations,
c   both CAPE and CIN will be set to -0.1 J/kg.  Also, where CAPE is
c   non-zero, CAPE and CIN will be set to a minimum of +0.1 J/kg, so
c   that the zero contour in either the CIN or CAPE fields will
c   circumscribe regions of non-zero CAPE.
c
      cape(i,j,kpar)=-0.1
      cin(i,j,kpar)=-0.1
      klfc=kmax
c
      goto 102
c
 50   continue
c
c   If there is an equilibrium level, then CAPE is positive.  We'll
c   define the level of free convection (LFC) as the point below the
c   EL, but at or above the LCL, where accumulated buoyant energy is a
c   minimum.  The net positive area (accumulated buoyant energy) from
c   the LFC up to the EL will be defined as the CAPE, and the net
c   negative area (negative of accumulated buoyant energy) from the
c   parcel starting point to the LFC will be defined as the convective
c   inhibition (CIN).
c
c   First get the LFC according to the above definition.
c
      benamin=9e9
      klfc=kmax
      do k=klcl,kel
         if (benaccum(k).lt.benamin) then
            benamin=benaccum(k)
            klfc=k
         endif
      enddo
c
c   Now we can assign values to cape and cin
c
      cape(i,j,kpar)=max(benaccum(kel)-benamin,0.1)
      cin(i,j,kpar)=max(-benamin,0.1)
c
c   CIN is uninteresting when CAPE is small (< 100 J/kg), so set
c   CIN to -.1 in that case.
c
      if ( cape(i,j,kpar) .lt. 100. ) cin(i,j,kpar) = -0.1
 102  continue
c
      enddo
c
      if (i3dflag.eq.0) then
         cape(i,j,mkzh)=cape(i,j,kpar1)
         cin(i,j,mkzh)=cin(i,j,kpar1)
         cin(i,j,mkzh-1)=zrel(klcl)+ghtpari-ter(i,j)   ! meters AGL
         cin(i,j,mkzh-2)=zrel(klfc)+ghtpari-ter(i,j)   ! meters AGL
      endif
c
      enddo
      enddo
c
      return
      end
