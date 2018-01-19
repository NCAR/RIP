c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine saweli(ght,tmk,qvp,prs,vge,vag,uuu,omg,bvfsqd,bvfsqm,
     &   rhi,imo,ilhs,rhithresh,imaxlit,imaxbig,
     &   errmin,alphor,ter,tergx,cor,xmap,rrfst,
     &   psixs,rcrag,rcrbg,ixaverage,smfac,smfstb,
     &   nscrs,nscd,xdist,ydist,xseclen,ptopse,pbotse,
     &   cfeld,mkp,miy,mjx,mkzh)
c
c   Note: errmin of 5 means 5 hPa*m/s, or .01 m/s in a 500-hPa layer,
c      or .01 dPa/s in a distance of 500 km.
c
c   Sawyer-Eliassen equation solver.
c      Note: y-axis is assumed to point from right to left in the
c      cross section.  x-axis is assumed to point into
c      the cross section.
c
c   Folloiwng is a list of all the fields that can be asked for.
c   Most are not described in the Users' Guide.
c
c      alp: pert alpa
c      baj: bclinic term (adj for ellip)
c      bcl: bclinic term
c      con: confluence forcing
c      faj: coef of dpsidy (adj for ellip)
c      fr1: overrlx factor 1
c      fr2: overrlx factor 2
c      fr3: overrlx factor 3
c      fr4: overrlx factor 4
c      fr5: overrlx factor 5
c      fr6: overrlx factor 6
c      fst: coef of dpsidy
c      ghp: geop ht pert
c      omb: balanced omega
c      omf: full omega
c      psi: streamfunction
c      pvo: PV
c      rhi: RH wrt ice
c      rsw: RH switch for dry or moist (.5 or 1.5)
c      saj: stability (adj for ellip)
c      she: shear forcing
c      std: dry stability
c      ste: effective stability
c      stm: moist stability
c      tot: confluence plus shear forcing
c      vab: balanced ageostrophic v
c      vag: ageostrophic v
c      vaj: vorticity (adj for ellip)
c      vge: geostrophic v
c      vor: vorticity
c      vtb: balanced v
c      vvv: full v
c 
      dimension ght(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   vge(miy,mjx,mkzh),vag(miy,mjx,mkzh),
     &   uuu(miy,mjx,mkzh),omg(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),ter(miy,mjx),tergx(miy,mjx),
     &   bvfsqd(miy,mjx,mkzh),bvfsqm(miy,mjx,mkzh),rhi(miy,mjx,mkzh),
     &   cor(miy,mjx),xmap(miy,mjx),rcrag(2),rcrbg(2),rrfst(4)
      dimension ghtxs(nscd,mkzh),vgexs(nscd,mkzh),
     &   vagxs(nscd,mkzh),prsxs(nscd,mkzh),
     &   uuuxs(nscd,mkzh),omgxs(nscd,mkzh),thetaxs(nscd,mkzh),
     &   alphaxs(nscd,mkzh),terxs(nscd),
     &   stabdxs(nscd,mkzh),stabmxs(nscd,mkzh),rhixs(nscd,mkzh),
     &   tergxxs(nscd),corxs(nscd),xmapxs(nscd),psixs(nscd,mkzh),
     &   ipttp(nscd,mkp),ghtse(nscd,mkp),stabse(nscd,mkp),
     &   stabdse(nscd,mkp),stabmse(nscd,mkp),rhise(nscd,mkp),
     &   alphase(nscd,mkp),vgese(nscd,mkp),thetase(nscd,mkp),
     &   vagse(nscd,mkp),uuuse(nscd,mkp),omgse(nscd,mkp),
     &   psi(nscd,mkp),resid(nscd,mkp),psitmp(nscd),
     &   facrel1(nscd,mkp),facrel2(nscd,mkp),facrel3(nscd,mkp),
     &   facrel4(nscd,mkp),facrel5(nscd,mkp),facrel6(nscd,mkp),
     &   vorse(nscd,mkp),bclse(nscd,mkp),fstse(nscd,mkp),
     &   qqq(nscd,mkp),rswse(nscd,mkp),thresh(nscd,mkp),plev(mkp),
     &   smtmp(nscd)
      character cfeld*10,fld*3,fldorig*3,frcterms*5
c
      include 'comconst'
c
      if (nproj.eq.4) then
         write(iup,*) 'Routine saweli will not work exactly right'
         write(iup,*) 'for the NMM SRCE map projection because'
         write(iup,*) 'it is not conformal and I was too lazy'
         write(iup,*) 'to account for different map factors in the'
         write(iup,*) 'x and y directions in this routine.'
         write(iup,*) 'However, errors will be quite small.'
      endif
c
      dzpbl=1000.   ! Assume 1000-m-thick PBL.
      perpfac=.20  ! PBL-avg. comp. of actual wind perp. to geos. wind
      imaxbigt=imaxbig
c
c   Get on/off switches for various forcing terms and BCs
c
      fld=cfeld(3:5)
      ifldtype=2
      if (fld.eq.'vvv'.or.fld.eq.'vge'.or.fld.eq.'vag'.or.
     &    fld.eq.'pvo'.or.fld.eq.'con'.or.fld.eq.'she'.or.
     &    fld.eq.'tot'.or.fld.eq.'std'.or.fld.eq.'stm'.or.
     &    fld.eq.'rhi'.or.fld.eq.'ghp'.or.fld.eq.'alp') ifldtype=1
      if (imo.eq.1.and.imaxbigt.gt.1.and.ifldtype.eq.2) then
         fldorig=fld
         fld='psi'
      endif
      if (cfeld(6:10).eq.'     ') then
         frcterms='11111'
      else
         frcterms=cfeld(6:10)
      endif
      read(frcterms(1:1),'(i1)') iconflu
      read(frcterms(2:2),'(i1)') ishear
      read(frcterms(3:3),'(i1)') itopobc
      read(frcterms(4:4),'(i1)') iekmnbc
      read(frcterms(5:5),'(i1)') isidebc
c
c   Set up x-sec averaging parameters
c
      cosangle=xdist/xseclen
      sinangle=ydist/xseclen
      rnavg=1./(2.*ixaverage+1.)
c
      do k=1,mkzh
      do ls=1,nscrs
         ghtxs(ls,k)=0.
         vgexs(ls,k)=0.
         vagxs(ls,k)=0.
         prsxs(ls,k)=0.
         uuuxs(ls,k)=0.
         omgxs(ls,k)=0.
         thetaxs(ls,k)=0.
         alphaxs(ls,k)=0.
         stabdxs(ls,k)=0.
         stabmxs(ls,k)=0.
         rhixs(ls,k)=0.
         if (k.eq.1) then
            terxs(ls)=0.
            tergxxs(ls)=0.
            corxs(ls)=0.
            xmapxs(ls)=0.
         endif
      enddo
      enddo
c
c   Interpolate gridded data to x-section.
c
      caxgn=1.+(rcrag(2)-xjcorn)*refrat
      caygn=1.+(rcrag(1)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1)-yicorn)*refrat
      do islab=-ixaverage,ixaverage
         xj1t=caxgn+islab*sinangle
         xj2t=cbxgn+islab*sinangle
         yi1t=caygn-islab*cosangle
         yi2t=cbygn-islab*cosangle
         if (xj1t.le.1.5.or.xj1t.ge.mjx-.5.or.
     &       xj2t.le.1.5.or.xj2t.ge.mjx-.5.or.
     &       yi1t.le.1.5.or.yi1t.ge.miy-.5.or.
     &       yi2t.le.1.5.or.yi2t.ge.miy-.5) then
            write(iup,*)'xj1t,xj2t,yi1t,yi2t=',xj1t,xj2t,yi1t,yi2t
            write(iup,*)'Cross sec. endpoints must be between 1.5'
            write(iup,*)'and (miy-.5) or (mjx-.5).'
            stop
         endif
c
      do k=1,mkzh
      do ls=1,nscrs
         posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t)-.5
         posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t)-.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         fac1=(1.-ratlr)*(   ratbt)
         fac2=(   ratlr)*(   ratbt)
         fac3=(1.-ratlr)*(1.-ratbt)
         fac4=(   ratlr)*(1.-ratbt)
         ghtxs(ls,k)=ghtxs(ls,k)+
     &      ( fac1*ght(it,jl,k) + fac2*ght(it,jr,k) +
     &        fac3*ght(ib,jl,k) + fac4*ght(ib,jr,k) )*rnavg
         tmkxstmp=
     &      ( fac1*tmk(it,jl,k) + fac2*tmk(it,jr,k) +
     &        fac3*tmk(ib,jl,k) + fac4*tmk(ib,jr,k) )
         qvpxstmp=
     &      ( fac1*qvp(it,jl,k) + fac2*qvp(it,jr,k) +
     &        fac3*qvp(ib,jl,k) + fac4*qvp(ib,jr,k) )
         tvkxstmp=virtual(tmkxstmp,qvpxstmp)
         vgexs(ls,k)=vgexs(ls,k)+
     &      ( fac1*vge(it,jl,k) + fac2*vge(it,jr,k) +
     &        fac3*vge(ib,jl,k) + fac4*vge(ib,jr,k) )*rnavg
         vagxs(ls,k)=vagxs(ls,k)+
     &      ( fac1*vag(it,jl,k) + fac2*vag(it,jr,k) +
     &        fac3*vag(ib,jl,k) + fac4*vag(ib,jr,k) )*rnavg
         uuuxs(ls,k)=uuuxs(ls,k)+
     &      ( fac1*uuu(it,jl,k) + fac2*uuu(it,jr,k) +
     &        fac3*uuu(ib,jl,k) + fac4*uuu(ib,jr,k) )*rnavg
         omgxs(ls,k)=omgxs(ls,k)+
     &      ( fac1*omg(it,jl,k) + fac2*omg(it,jr,k) +
     &        fac3*omg(ib,jl,k) + fac4*omg(ib,jr,k) )*rnavg
         prsxstmp=
     &      ( fac1*prs(it,jl,k) + fac2*prs(it,jr,k) +
     &        fac3*prs(ib,jl,k) + fac4*prs(ib,jr,k) )
         prsxs(ls,k)=prsxs(ls,k)+prsxstmp*rnavg
         alphaxstmp=rgas*tvkxstmp/prsxstmp
         alphaxs(ls,k)=alphaxs(ls,k)+alphaxstmp*rnavg
         thetaxs(ls,k)=thetaxs(ls,k)+
     &      tvkxstmp*(1000./prsxstmp)**gamma*rnavg  ! this is actually theta_v
         stabdxs(ls,k)=stabdxs(ls,k)+(alphaxstmp/grav)**2*
     &      ( fac1*bvfsqd(it,jl,k) + fac2*bvfsqd(it,jr,k) +
     &        fac3*bvfsqd(ib,jl,k) + fac4*bvfsqd(ib,jr,k) )*rnavg
         if (imo.eq.1) then
            stabmxs(ls,k)=stabmxs(ls,k)+(alphaxstmp/grav)**2*
     &         ( fac1*bvfsqm(it,jl,k) + fac2*bvfsqm(it,jr,k) +
     &           fac3*bvfsqm(ib,jl,k) + fac4*bvfsqm(ib,jr,k) )*rnavg
            rhixs(ls,k)=rhixs(ls,k)+
     &         ( fac1*rhi(it,jl,k) + fac2*rhi(it,jr,k) +
     &           fac3*rhi(ib,jl,k) + fac4*rhi(ib,jr,k) )*rnavg
         endif
         if (k.eq.1) then
            terxs(ls)=terxs(ls)+
     &         ( fac1*ter(it,jl) + fac2*ter(it,jr) +
     &           fac3*ter(ib,jl) + fac4*ter(ib,jr) )*rnavg
            tergxxs(ls)=tergxxs(ls)+
     &         ( fac1*tergx(it,jl) + fac2*tergx(it,jr) +
     &           fac3*tergx(ib,jl) + fac4*tergx(ib,jr) )*rnavg
            corxs(ls)=corxs(ls)+
     &         ( fac1*cor(it,jl) + fac2*cor(it,jr) +
     &           fac3*cor(ib,jl) + fac4*cor(ib,jr) )*rnavg
            xmapxs(ls)=xmapxs(ls)+
     &         ( fac1*xmap(it,jl) + fac2*xmap(it,jr) +
     &           fac3*xmap(ib,jl) + fac4*xmap(ib,jr) )*rnavg
         endif
      enddo
      enddo
      enddo
c
c   Some more constants
c
      dp=(pbotse-ptopse)/(mkp-1.)
      dpi=1./dp
      p5dpi=.5*dpi
      hscale=rgas*celkel/grav
      hscali=1./hscale
c
c   Define pressure levels.
c
      do k=1,mkp
         plev(k)=ptopse+(k-1.)*dp
      enddo
c
c   Calculate pressure-level fields.
c
      do kp=1,mkp
      do ls=1,nscrs
         if (plev(kp).gt.prsxs(ls,mkzh)) then
            ghtse(ls,kp)=rmsg
            stabdse(ls,kp)=rmsg
            stabse(ls,kp)=rmsg
            if (imo.eq.1) then
               stabmse(ls,kp)=rmsg
               rhise(ls,kp)=rmsg
            endif
            thetase(ls,kp)=rmsg
            alphase(ls,kp)=rmsg
            vgese(ls,kp)=rmsg
            vagse(ls,kp)=rmsg
            uuuse(ls,kp)=rmsg
            omgse(ls,kp)=rmsg
         elseif (plev(kp).lt.prsxs(ls,1)) then
            write(iup,*)'At point ls,kp=',ls,kp,
     &         '  pressure is too low.'
            write(iup,*)'plev,prsxs=',plev(kp),prsxs(ls,1)
            stop
         else
            do ks=1,mkzh-1
               if (plev(kp).le.prsxs(ls,ks+1).and.
     &             plev(kp).ge.prsxs(ls,ks)) then
                  fac1=plev(kp)-prsxs(ls,ks)
                  fac2=prsxs(ls,ks+1)-plev(kp)
                  denom=prsxs(ls,ks+1)-prsxs(ls,ks)
                  prat=(fac1*exp(-ghtxs(ls,ks+1)*hscali)+
     &                  fac2*exp(-ghtxs(ls,ks)*hscali)     )/denom
                  ghtse(ls,kp)=-hscale*log(prat)
                  alphase(ls,kp)=(fac1*alphaxs(ls,ks+1)+
     &               fac2*alphaxs(ls,ks))/denom
                  stabdse(ls,kp)=(fac1*stabdxs(ls,ks+1)+
     &               fac2*stabdxs(ls,ks))/denom
                  if (imo.eq.1) then
                     stabmse(ls,kp)=(fac1*stabmxs(ls,ks+1)+
     &                  fac2*stabmxs(ls,ks))/denom
c                     if (plev(kp).ge.400.) stabmse(ls,kp)=.0001
                     rhise(ls,kp)=(fac1*rhixs(ls,ks+1)+
     &                  fac2*rhixs(ls,ks))/denom
                  endif
                  thetase(ls,kp)=(fac1*thetaxs(ls,ks+1)+
     &               fac2*thetaxs(ls,ks))/denom
                  vgese(ls,kp)=(fac1*vgexs(ls,ks+1)+fac2*vgexs(ls,ks))/
     &               denom
                  vagse(ls,kp)=(fac1*vagxs(ls,ks+1)+fac2*vagxs(ls,ks))/
     &               denom
                  uuuse(ls,kp)=(fac1*uuuxs(ls,ks+1)+fac2*uuuxs(ls,ks))/
     &               denom
                  omgse(ls,kp)=(fac1*omgxs(ls,ks+1)+fac2*omgxs(ls,ks))/
     &               denom
                  goto 30
               endif
            enddo
            write(iup,*)'shouldn''t have gotten here.'
            stop
 30         continue
         endif
      enddo
      enddo
c
c   Identify points.
c
      do ls=1,nscrs
      do k=1,mkp
         ipttp(ls,k)=-1
      enddo
      enddo
      do k=2,mkp-1
      do ls=2,nscrs-1
         if (ghtse(ls,k).ne.rmsg.and.
     &       ghtse(ls+1,k).ne.rmsg.and.ghtse(ls-1,k).ne.rmsg.and.
     &       ghtse(ls,k+1).ne.rmsg.and.ghtse(ls,k-1).ne.rmsg.and.
     &       ghtse(ls+1,k+1).ne.rmsg.and.ghtse(ls+1,k-1).ne.rmsg.and.
     &       ghtse(ls-1,k+1).ne.rmsg.and.ghtse(ls-1,k-1).ne.rmsg) then
            ipttp(ls,k)=1 !interior point
         endif
      enddo
      enddo
      do ls=1,nscrs
         ipttp(ls,1)=2 ! top boundary
      enddo
      do k=2,mkp
         if (ghtse(1,k).ne.rmsg) then
            ipttp(1,k)=3 ! side bondary
         endif
         if (ghtse(nscrs,k).ne.rmsg) then
            ipttp(nscrs,k)=3 ! side boundary
         endif
      enddo
      do ls=1,nscrs
      do k=1,mkp
         if (ipttp(ls,k).eq.-1.and.ghtse(ls,k).ne.rmsg) then
            ipttp(ls,k)=4 !bottom boundary
         endif
      enddo
      enddo
c
c   Smooth the fields
c
      if (smfac.gt.1.) then
         call smoothslice2d(ghtse,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(alphase,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(stabdse,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         if (imo.eq.1) then
            call smoothslice2d(stabmse,smtmp,ipttp,1,nscrs,mkp,
     &                       smfac)
            call smoothslice2d(rhise,smtmp,ipttp,1,nscrs,mkp,
     &                       smfac)
         endif
         call smoothslice2d(thetase,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(vgese,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(vagse,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(uuuse,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice2d(omgse,smtmp,ipttp,1,nscrs,mkp,
     &      smfac)
         call smoothslice1d(terxs,smtmp,nscrs,smfac)
         call smoothslice1d(tergxxs,smtmp,nscrs,smfac)
      endif
c
c   Get fbar
c
      fbar=0.
      do ls=1,nscrs
         fbar=fbar+corxs(ls)
      enddo
      fbar=fbar/float(nscrs)
c
c   Calculate Q-forcing
c
      do ls=1,nscrs
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
         dyi=1./dy
         p5dyi=.5*dyi
         dyi2=dyi*dyi
      do k=1,mkp
         psi(ls,k)=rmsg
         if (ipttp(ls,k).eq.1) then
            dvgdy=(vgese(ls+1,k)-vgese(ls-1,k))*p5dyi
            d2zdydp=(ghtse(ls+1,k+1)+ghtse(ls-1,k-1)-
     &               ghtse(ls-1,k+1)-ghtse(ls+1,k-1))*p5dyi*p5dpi
            dvgdp=(vgese(ls,k+1)-vgese(ls,k-1))*p5dpi
            d2zdy2=(ghtse(ls-1,k)-2.*ghtse(ls,k)+ghtse(ls+1,k))*dyi2
            confluterm=-2.*grav/fbar*dvgdy*d2zdydp
            shearterm=2.*grav/fbar*dvgdp*d2zdy2
            if (fld.eq.'con') then
               psi(ls,k)=confluterm
            elseif (fld.eq.'she') then
               psi(ls,k)=shearterm
            elseif (fld.eq.'tot') then
               psi(ls,k)=confluterm+shearterm
            endif
            qqq(ls,k)=iconflu*confluterm+ishear*shearterm
         endif
      enddo
      enddo
c
c   Fill output fields that would not change from here on out.
c
      refslp=rrfst(1)
      refslt=rrfst(2)
      reflaps=rrfst(3)
      refstratt=rrfst(4)
      if (fld.eq.'vvv') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=vgese(ls,k)+vagse(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'vge') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=vgese(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'fgg') then
         do ls=2,nscrs-1
            dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
            dyi=1./dy
            p5dyi=.5*dyi
         do k=1,mkp
            if (ipttp(ls,k).eq.1.or.ipttp(ls,k).eq.2) then
               psi(ls,k)=(vgese(ls+1,k)-vgese(ls-1,k))*p5dyi*
     &            (thetase(ls+1,k)-thetase(ls-1,k))*p5dyi*
     &            *1.e5*3600.   ! K per 100 km per hour
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'vag') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=vagse(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'pvo') then
         do ls=2,nscrs-1
            dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
            dyi=1./dy
            p5dyi=.5*dyi
         do k=2,mkp
            if (ipttp(ls,k).eq.1) then
               vvor=-(uuuse(ls+1,k)-uuuse(ls-1,k))*p5dyi+corxs(ls)
               stabil=(thetase(ls,k+1)-thetase(ls,k-1))*p5dpi
               hvor=-(uuuse(ls,k+1)-uuuse(ls,k-1))*p5dpi
               thgrad=-(thetase(ls+1,k)-thetase(ls-1,k))*p5dyi
               psi(ls,k)=-1e4*grav*(vvor*stabil+hvor*thgrad)
               if (k.eq.2) psi(ls,1)=psi(ls,2)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'ghp') then
         cc1=rgas/grav*(-.5)*reflaps
         cc2=rgas/grav*(reflaps*log(.01*refslp)-refslt)
         cc3=rgas/grav*(refslt-.5*reflaps*log(.01*refslp))*
     &      log(.01*refslp)
         alnpreftpause=(refstratt-refslt)/reflaps+log(.01*refslp)
         ztpause=cc1*alnpreftpause*alnpreftpause+
     &      cc2*alnpreftpause+cc3
         do k=1,mkp
            alnpref=log(plev(k))
            if (alnpref.gt.alnpreftpause) then
               refght=cc1*alnpref*alnpref+cc2*alnpref+cc3
            else
               refght=ztpause+rgas*refstratt/grav*
     &            (alnpreftpause-alnpref)
            endif
         do ls=1,nscrs
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=ghtse(ls,k)-refght
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'alp') then
         do k=1,mkp
            reftmk=refslt+reflaps*log(plev(k)/(.01*refslp))
            reftmk=max(reftmk,refstratt)
            refalpha=rgas*reftmk/plev(k)
         do ls=1,nscrs
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=alphase(ls,k)-refalpha
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'std') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=stabdse(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'stm') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=stabmse(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'rhi') then
         do ls=1,nscrs
         do k=1,mkp
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=rhise(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
      endif
      if (fld.eq.'con'.or.fld.eq.'she'.or.fld.eq.'tot'.or.
     &    fld.eq.'pvo') call fillpsi(psi,ipttp,nscrs,mkp)
      if (ifldtype.eq.1) goto 300
c
      ibigiter=0
      inochg=0
c
c   This is the point we come back to if doing big iterations
c   after moist stability adjustment.
c
 987  ibigiter=ibigiter+1
      if (imaxbigt.gt.1) write(iup,*)'Starting big iteration ',ibigiter
c
c   Make effective stability, vorticity, and baroclinic terms.
c   Fill output fields that don't depend
c   on any further calculations, if this is the final big iteration.
c
      if (ibigiter.gt.1) nchg=0
      do ls=1,nscrs
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
         dyi=1./dy
         p5dyi=.5*dyi
         dyi2=dyi*dyi
      do k=1,mkp
         if (ibigiter.gt.1) rswold=rswse(ls,k)
         rswse(ls,k)=0.5
         if (ipttp(ls,k).ge.1) then
            if (imo.eq.1)then
               if (omgse(ls,k).lt.0.0.and.
     &          rhise(ls,k).gt.rhithresh) then
                  rswse(ls,k)=1.5
                  stabse(ls,k)=stabmse(ls,k)
               else
                  stabse(ls,k)=stabdse(ls,k)
               endif
            else
               stabse(ls,k)=stabdse(ls,k)
            endif
            if (ipttp(ls,k).eq.1) then
               vorse(ls,k)=grav/fbar*(ghtse(ls-1,k)-2.*ghtse(ls,k)+
     &            ghtse(ls+1,k))*dyi2+fbar
               bclse(ls,k)=-grav/fbar*(ghtse(ls+1,k+1)+ghtse(ls-1,k-1)-
     &                  ghtse(ls-1,k+1)-ghtse(ls+1,k-1))*p5dyi*p5dpi
               if (fld.eq.'vor') then
                  psi(ls,k)=vorse(ls,k)
               elseif (fld.eq.'bcl') then
                  psi(ls,k)=bclse(ls,k)
               endif
            endif
            if (fld.eq.'omf') then
               psi(ls,k)=omgse(ls,k)
            elseif (fld.eq.'rsw') then
               psi(ls,k)=rswse(ls,k)
            elseif (fld.eq.'ste') then
               psi(ls,k)=stabse(ls,k)
            endif
         else
            if (fld.eq.'omf'.or.fld.eq.'rsw'.or.fld.eq.'ste') then
               psi(ls,k)=psi(ls,k-1)
            endif
         endif
         if (imo.eq.1.and.ibigiter.gt.1..and.
     &                    rswse(ls,k).ne.rswold) nchg=nchg+1
      enddo
      enddo
c
      if (imo.eq.1.and.ibigiter.gt.1.and.nchg.eq.0) inochg=inochg+1
c
      if (fld.eq.'omf'.or.fld.eq.'rsw'.or.fld.eq.'ste') goto 300
      if (fld.eq.'vor'.or.fld.eq.'bcl') then
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Calculate coefficient of d(psi)/dy
c
      do ls=1,nscrs
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
         dyi=1./dy
         p5dyi=.5*dyi
      do k=1,mkp
         if (ipttp(ls,k).eq.1) then
            kp1=k+1
            if (ipttp(ls,kp1).ne.1) kp1=k
            km1=k-1
            if (ipttp(ls,km1).ne.1) km1=k
            dpk=dp*(kp1-km1)
            if (kp1.eq.km1) then
               write(iup,*)'Can''t do p-deriv. of vorticity.'
               write(iup,*)'ls,k=',ls,k
               stop
            endif
            fstse(ls,k)=1./fbar*(stabse(ls+1,k)-stabse(ls-1,k))*p5dyi+
     &         (bclse(ls,kp1)-bclse(ls,km1))/dpk
            if (fld.eq.'fst') psi(ls,k)=fstse(ls,k)
         endif
      enddo
      enddo
      if (fld.eq.'fst') then
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Adjust stab, vor, and bcl for ellipticity, then recalculate fst.
c   Use thresh to hold the threshhold values of each array.
c
c   First do stab (the effective stability term).  Make it greater
c   than a small positive threshold.  If the min. stability is
c   represented by a min. delta theta (DELTH) in a given depth (DELZ),
c   then the corresponding minimum stab term is given by
c   alpha**2/(theta*grav*DELZ)*DELTH
c
      delth=.5 !K
      delz=3000. !km
      do k=1,mkp
      do ls=1,nscrs
         if (ipttp(ls,k).ge.1) thresh(ls,k)=alphase(ls,k)**2/
     &         (thetase(ls,k)*grav*delz)*delth
      enddo
      enddo
      call adjellip(stabse,thresh,ipttp,1,nscrs,mkp,tot,totdiff)
      if (abs(totdiff/tot).gt..1) then
         write(iup,*)'Lots of adjusting of stability term (stab)'
         write(iup,*)'needed to make entire domain stable.'
         write(iup,*)'Results may be suspect, if you get any at all.'
         write(iup,*)'Total, total diff=',tot,totdiff
      endif
c
c   Smooth stab term with linear weighted smoother
c
      if (smfstb.gt.1.0)
     &   call smoothslice2d(stabse,smtmp,ipttp,1,nscrs,mkp,smfstb)
      if (fld.eq.'saj') then
         do k=1,mkp
         do ls=1,nscrs
            if (ipttp(ls,k).ge.1) then
               psi(ls,k)=stabse(ls,k)
            else
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
         enddo
         goto 300
      endif
c
c   Next, do vor (the vorticity term).  Make it greater than a small
c   positive threshold.
c
      vorthresh=.05e-4
      do k=1,mkp
      do ls=1,nscrs
         if (ipttp(ls,k).eq.1) thresh(ls,k)=vorthresh
      enddo
      enddo
      call adjellip(vorse,thresh,ipttp,0,nscrs,mkp,tot,totdiff)
      if (abs(totdiff/tot).gt..1) then
         write(iup,*)'Lots of adjusting of vorticity term (vor)'
         write(iup,*)'needed to make entire domain have pos. abs. vor.'
         write(iup,*)'Results may be suspect, if you get any at all.'
         write(iup,*)'Total, total diff=',tot,totdiff
      endif
      if (fld.eq.'vaj') then
         do k=1,mkp
         do ls=1,nscrs
            if (ipttp(ls,k).eq.1) psi(ls,k)=vorse(ls,k)
         enddo
         enddo
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Finally, do bcl (the baroclinic term).  First put the term
c   stab/fbar*vor-bcl**2 into bcl.  Make this greater than a small
c   positive threshold.  Then extract the new bcl from that adjusted
c   quantity.  Use uuuse to hold the sign of bcl.
c
      do k=1,mkp
      do ls=1,nscrs
         if (ipttp(ls,k).eq.1) then
c            stabthresh=alphase(ls,k)**2/
c     &         (thetase(ls,k)*grav*delz)*delth
c            thresh(ls,k)=.97*stabthresh*vorthresh/fbar
c            uuuse(ls,k)=sign(1.,bclse(ls,k))
c            bclse(ls,k)=stabse(ls,k)*vorse(ls,k)/fbar-
c     &         bclse(ls,k)*bclse(ls,k)
            thresh(ls,k)=-sqrt(.999*stabse(ls,k)*vorse(ls,k)/fbar)
            uuuse(ls,k)=sign(1.,bclse(ls,k))
            bclse(ls,k)=-abs(bclse(ls,k))
         endif
      enddo
      enddo
      call adjellip(bclse,thresh,ipttp,0,nscrs,mkp,tot,totdiff)
      if (abs(totdiff/tot).gt..1) then
         write(iup,*)'Lots of adjusting of therm. wind term (bcl)'
         write(iup,*)'needed to satisfy ellipticity.'
         write(iup,*)'Results may be suspect, if you get any at all.'
         write(iup,*)'Total, total diff=',tot,totdiff
      endif
      do k=1,mkp
      do ls=1,nscrs
         if (ipttp(ls,k).eq.1) then
c            bclse(ls,k)=sqrt(stabse(ls,k)*vorse(ls,k)/fbar-bclse(ls,k))
            bclse(ls,k)=sign(bclse(ls,k),uuuse(ls,k))
         endif
      enddo
      enddo
      if (fld.eq.'baj') then
         do k=1,mkp
         do ls=1,nscrs
            if (ipttp(ls,k).eq.1) psi(ls,k)=bclse(ls,k)
         enddo
         enddo
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Recalculate coefficient of d(psi)/dy
c
      do ls=1,nscrs
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
         dyi=1./dy
         p5dyi=.5*dyi
      do k=1,mkp
         if (ipttp(ls,k).eq.1) then
            kp1=k+1
            if (ipttp(ls,kp1).ne.1) kp1=k
            km1=k-1
            if (ipttp(ls,km1).ne.1) km1=k
            dpk=dp*(kp1-km1)
            if (kp1.eq.km1) then
               write(iup,*)'Can''t do p-deriv. of vorticity.'
               write(iup,*)'ls,k=',ls,k
               stop
            endif
            fstse(ls,k)=1./fbar*(stabse(ls+1,k)-stabse(ls-1,k))*p5dyi+
     &         (bclse(ls,kp1)-bclse(ls,km1))/dpk
            if (fld.eq.'faj') psi(ls,k)=fstse(ls,k)
         endif
      enddo
      enddo
      if (fld.eq.'faj') then
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Calculate various factors for over-relaxation.
c
      do ls=1,nscrs
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
      do k=1,mkp
         if (ipttp(ls,k).eq.1) then
            aterm=stabse(ls,k)/fbar
            f12term=fstse(ls,k)*dy/(2.*aterm)
            if (ilhs.eq.1) then ! SG
               f45term=vorse(ls,k)*dy*dy/(aterm*dp*dp)
               facrel1(ls,k)=1.-f12term
               facrel2(ls,k)=1.+f12term
               facrel3(ls,k)=2.*bclse(ls,k)*dy/(4.*aterm*dp)
            else
               f45term=fbar*dy*dy/(aterm*dp*dp)
               facrel1(ls,k)=1. ! QG
               facrel2(ls,k)=1.
               facrel3(ls,k)=0.
            endif
            facrel4(ls,k)=f45term
            facrel5(ls,k)=2.*(1.+f45term)
            facrel6(ls,k)=dy*dy*qqq(ls,k)/aterm
            if (fld.eq.'fr1') then
               psi(ls,k)=facrel1(ls,k)
            elseif (fld.eq.'fr2') then
               psi(ls,k)=facrel2(ls,k)
            elseif (fld.eq.'fr3') then
               psi(ls,k)=facrel3(ls,k)
            elseif (fld.eq.'fr4') then
               psi(ls,k)=facrel4(ls,k)
            elseif (fld.eq.'fr5') then
               psi(ls,k)=facrel5(ls,k)
            elseif (fld.eq.'fr6') then
               psi(ls,k)=facrel6(ls,k)
            endif
         endif
      enddo
      enddo
      if (fld(1:2).eq.'fr') then
         call fillpsi(psi,ipttp,nscrs,mkp)
         goto 300
      endif
c
c   Set all psi values to zero initially
c
      do ls=1,nscrs
      do k=1,mkp
         psi(ls,k)=0.
      enddo
      enddo
c
c   Set boundary values of psi.  Start in upper right corner
c   (ls=nscrs,k=1) and work around the circuit counterclockwise
c
c   First the top boundary.  Leave it equal to zero.
c
c   Next the left side bondary
c
      do k=2,mkp
         if (ipttp(1,k).eq.3) then
            vagavg=.5*(vagse(1,k-1)+vagse(1,k))
            psi(1,k)=psi(1,k-1)-isidebc*vagavg*dp
            k1last=k
         endif
      enddo
c
c   Next, the bottom boundary, but put in temporary 1-d array
c
c       Note: Ground-level
c       will be used as the bottom boundary.  The Ekman BC will not
c       be completely accurate because, although it is calculated
c       assuming a PBL depth of dzpbl, it is applied at ground
c       level.  In setting the lower boundary to ground level, the
c       Ekman inaccuracy was accepted in exchange for the benefit of
c       allowing the Q-vector and terrain forcing to act all the way
c       down to the ground.
c
      psitmp(1)=psi(1,k1last)
      do ls=2,nscrs
         do k=mkp,1,-1
            if (ipttp(ls,k).ge.1.and.ipttp(ls-1,k).ge.1) then
               kdo=k
               goto 345
            endif
         enddo
 345     continue
         dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
         dyi=1./dy
         psitmp(ls)=psitmp(ls-1)
         if (itopobc.eq.1) then
            psitmp(ls)=psitmp(ls)+
     &       grav*2./(alphase(ls,kdo)+alphase(ls-1,kdo))*(
     &       grav/fbar*(ghtse(ls,kdo)-ghtse(ls-1,kdo))*dyi*tergxxs(ls)-
     &       vgese(ls,kdo)*(terxs(ls)-terxs(ls-1))*dyi  )*dy
         endif
         if (iekmnbc.eq.1) then
            ls1=min(max(ls-2,1),nscrs-3)
            ls2=ls1+1
            ls3=ls1+2
            ls4=ls1+3
            do k=mkp,1,-1
               if (ipttp(ls1,k).ge.1.and.ipttp(ls2,k).ge.1.and.
     &             ipttp(ls3,k).ge.1.and.ipttp(ls4,k).ge.1) then
                  kdo2=k
                  goto 346
               endif
            enddo
 346        continue
            psitmp(ls)=psitmp(ls)-
     &         grav*2./(alphase(ls2,kdo2)+alphase(ls3,kdo2))*(
     &         perpfac*dzpbl*grav/fbar*(ghtse(ls1,kdo2)+
     &         ghtse(ls4,kdo2)-ghtse(ls2,kdo2)-
     &         ghtse(ls3,kdo2))*.5*dyi*dyi )*dy
         endif
      enddo
c
c   Finally, the right side boundary
c
      knlast=-99
      do k=mkp,1,-1
         if (ipttp(nscrs,k).eq.3.or.k.eq.1) then
            if (knlast.eq.-99) then
               knlast=k
               psi(nscrs,k)=psitmp(nscrs)
            else
               vagavg=.5*(vagse(nscrs,k)+vagse(nscrs,k+1))
               psi(nscrs,k)=psi(nscrs,k+1)+isidebc*vagavg*dp
            endif
         endif
      enddo
c
c   Spread the circuit integral "error" over the entire circuit.
c
      error=psi(nscrs,1)
      derr=error/(2.*(nscrs-1.)+k1last+knlast-2.)
      n=0
      do ls=nscrs-1,1,-1
         n=n+1
         psi(ls,1)=psi(ls,1)-derr*n
      enddo
      do k=2,k1last
         n=n+1
         psi(1,k)=psi(1,k)-derr*n
      enddo
      psitmp(1)=psi(1,k1last)
      do ls=2,nscrs
         n=n+1
         psitmp(ls)=psitmp(ls)-derr*n
      enddo
      psi(nscrs,knlast)=psitmp(nscrs)
      do k=knlast-1,1,-1
         n=n+1
         psi(nscrs,k)=psi(nscrs,k)-derr*n
      enddo
      if (psi(nscrs,1).gt..1) then
         write(iup,*)'circuit integral of d(psi)/ds not corrected.'
         stop
      endif
c
c   Transfer the bottom BC from the temp. array to the psi array.
c
      do ls=1,nscrs
      do k=1,mkp
         if (ipttp(ls,k).eq.4.or.ipttp(ls,k).eq.-1)
     &      psi(ls,k)=psitmp(ls)
      enddo
      enddo
c
 207  continue
c
c   Over-relaxation iteration loop
c
      iter=0
   40 resmax=0.
      resav=0.
      nres=0
      iter=iter+1
c
c   Do the over-relaxation
c
      do k=1,mkp
      do ls=1,nscrs
         if (ipttp(ls,k).eq.1) then
            res=facrel1(ls,k)*psi(ls-1,k)+facrel2(ls,k)*psi(ls+1,k)+
     &         facrel3(ls,k)*(psi(ls+1,k+1)+psi(ls-1,k-1)-
     &         psi(ls+1,k-1)-psi(ls-1,k+1))+facrel4(ls,k)*
     &         (psi(ls,k-1)+psi(ls,k+1))-facrel5(ls,k)*psi(ls,k)-
     &         facrel6(ls,k)
            resid(ls,k)=abs(res)
            resmax=max(resmax,abs(res))
            resav=resav+abs(res)
            nres=nres+1
            psi(ls,k)=psi(ls,k)+alphor*res/facrel5(ls,k)
c            psi(ls,k)=min(max(psi(ls,k),-1500.),1500.)
            psi(ls,k)=min(max(psi(ls,k),-20000.),20000.)
         endif
      enddo
      enddo
c
      resav=resav/nres
      write(iup,'(a,i4,2(1x,f18.6))')'     iter,resmax,resav=',
     &   iter,resmax,resav
      if (resmax.gt.errmin.and.iter.lt.imaxlit) goto 40
      if (resmax.gt.errmin.and.iter.ge.imaxlit) then
         write(iup,*)'   In SAWELI: Didn''t converge in ',imaxlit,
     &      ' iterations.'
      else
         write(iup,*)'   In SAWELI: Converged in ',iter,' iterations.'
      endif
c
c   If moist run, recalculate omega (in omgse) from psi, and go through
c   whole process again, starting at the point where effective
c   stability was calculated.
c
      if (imo.eq.1.and.ibigiter.lt.imaxbigt) then
         do ls=1,nscrs
            dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
            dyi=1./dy
            p5dyi=.5*dyi
         do k=1,mkp
            if (ipttp(ls,k).eq.1) then
               omgse(ls,k)=1000.*(psi(ls+1,k)-psi(ls-1,k))*p5dyi
            endif
         enddo
         enddo
         if (inochg.ge.2) imaxbigt=ibigiter+1
         if (ibigiter.eq.imaxbigt-1) then
            fld=fldorig
         endif
         goto 987
      endif
c
c   Calculate fields that are derived from SE streamfunction,
c   if asked for.  Put field into thresh.
c
      if (fld.eq.'vab'.or.fld.eq.'vtb') then
         do ls=1,nscrs
         do k=2,mkp
            if (ipttp(ls,k).eq.1.or.ipttp(ls,k).eq.3) then
               thresh(ls,k)=-(psi(ls,k+1)-psi(ls,k-1))*p5dpi
               if (fld.eq.'vtb') then
                  thresh(ls,k)=thresh(ls,k)+vgese(ls,k)
               endif
               if (k.eq.2) thresh(ls,1)=thresh(ls,2)
            else
               thresh(ls,k)=thresh(ls,k-1)
            endif
         enddo
         enddo
      elseif (fld.eq.'omb') then
         do ls=2,nscrs-1
            dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
            dyi=1./dy
            p5dyi=.5*dyi
         do k=1,mkp
            if (ipttp(ls,k).eq.1.or.ipttp(ls,k).eq.2) then
               thresh(ls,k)=1000.*(psi(ls+1,k)-psi(ls-1,k))*p5dyi
            else
               thresh(ls,k)=thresh(ls,k-1)
            endif
         enddo
         enddo
         do k=1,mkp
            thresh(1,k)=thresh(2,k)
            thresh(nscrs,k)=thresh(nscrs-1,k)
         enddo
      elseif (fld.eq.'fga') then
         do ls=2,nscrs-1
            dy=-(ds*xseclen)/((nscrs-1.)*xmapxs(ls))
            dyi=1./dy
            p5dyi=.5*dyi
         do k=1,mkp
            if (ipttp(ls,k).eq.1.or.ipttp(ls,k).eq.2) then
               thresh(ls,k)=1000.*(psi(ls+1,k)-psi(ls-1,k))*p5dyi
            else
               thresh(ls,k)=thresh(ls,k-1)
            endif
         enddo
         enddo
         do k=1,mkp
            thresh(1,k)=thresh(2,k)
            thresh(nscrs,k)=thresh(nscrs-1,k)
         enddo
      else
         goto 300
      endif
c
c   Move data from thresh back into psi
c
      do ls=1,nscrs
      do k=1,mkp
         psi(ls,k)=thresh(ls,k)
      enddo
      enddo
c
c   Interpolate back to original model vertical levels.
c
 300  continue
      do ls=1,nscrs
      do ks=1,mkzh
         if (prsxs(ls,ks).le.plev(mkp).and.
     &       prsxs(ls,ks).ge.plev(1)) then
            do kp=1,mkp-1
               if (prsxs(ls,ks).le.plev(kp+1).and.
     &             prsxs(ls,ks).ge.plev(kp)) then
                  if (psi(ls,kp+1).eq.rmsg.or.psi(ls,kp).eq.rmsg) then
                     psixs(ls,ks)=rmsg
                  else
                     psixs(ls,ks)=((prsxs(ls,ks)-plev(kp))*
     &                  psi(ls,kp+1)+(plev(kp+1)-prsxs(ls,ks))*
     &                  psi(ls,kp))*dpi
                  endif
                  goto 50
               endif
            enddo
 50         continue
         elseif (prsxs(ls,ks).gt.plev(mkp)) then
            psixs(ls,ks)=psi(ls,mkp)
         elseif (prsxs(ls,ks).lt.plev(1)) then
            psixs(ls,ks)=psi(ls,1)
         endif
      enddo
      enddo
c
      return
      end
c
      subroutine fillpsi(psi,ipttp,nscrs,mkp)
      dimension psi(nscrs,mkp),ipttp(nscrs,mkp)
      do ls=2,nscrs-1
         psi(ls,1)=psi(ls,2)
      enddo
      psi(1,1)=psi(2,1)
      psi(nscrs,1)=psi(nscrs-1,1)
      do k=2,mkp
         if (ipttp(1,k).eq.3.and.ipttp(2,k).ge.1) then
               psi(1,k)=psi(2,k)
         else
            psi(1,k)=psi(1,k-1)
         endif
         if (ipttp(nscrs,k).eq.3.and.ipttp(nscrs-1,k).ge.1) then
               psi(nscrs,k)=psi(nscrs-1,k)
         else
            psi(nscrs,k)=psi(nscrs,k-1)
         endif
         do ls=2,nscrs-1
            if (ipttp(ls,k).eq.4.or.ipttp(ls,k).eq.-1) then
               psi(ls,k)=psi(ls,k-1)
            endif
         enddo
      enddo
      return
      end
c
      subroutine smoothslice2d(arr,smtmp,ipttp,iall,nscrs,mkp,smfac)
      dimension arr(nscrs,mkp), smtmp(nscrs), ipttp(nscrs,mkp)
      ld=int(smfac)
      if (ld.eq.0) return
      do k=1,mkp
         do ls=1,nscrs
            if ((iall.eq.0.and.ipttp(ls,k).eq.1).or.
     &          (iall.eq.1.and.ipttp(ls,k).ge.1)) then
               totw=0.0
               tot=0.0
               do l=ls-ld,ls+ld
                  if (l.ge.1.and.l.le.nscrs) then
                     if ((iall.eq.0.and.ipttp(l,k).eq.1).or.
     &                   (iall.eq.1.and.ipttp(l,k).ge.1)) then
                        wt=smfac-abs(ls-l)
                        totw=totw+wt
                        tot=tot+wt*arr(l,k)
                     endif
                  endif
               enddo
               smtmp(ls)=tot/totw
            endif
         enddo
         do ls=1,nscrs
            if ((iall.eq.0.and.ipttp(ls,k).eq.1).or.
     &          (iall.eq.1.and.ipttp(ls,k).ge.1)) then
               arr(ls,k)=smtmp(ls)
            endif
         enddo
      enddo
      return
      end
c
      subroutine smoothslice1d(arr,smtmp,nscrs,smfac)
      dimension arr(nscrs), smtmp(nscrs)
      ld=int(smfac)
      if (ld.eq.0) return
      do ls=1,nscrs
         totw=0.0
         tot=0.0
         do l=ls-ld,ls+ld
            if (l.ge.1.and.l.le.nscrs) then
               wt=smfac-abs(ls-l)
               totw=totw+wt
               tot=tot+wt*arr(l)
            endif
         enddo
         smtmp(ls)=tot/totw
      enddo
      do ls=1,nscrs
         arr(ls)=smtmp(ls)
      enddo
      return
      end
