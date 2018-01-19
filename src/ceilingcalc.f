c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ceilingcalc(qvp,qcw,qra,qci,qsn,tmk,prs,pf,
     &   prsnow,meth,clgfac,clg,miy,mjx,mkzh)
c
c     This routine computes cloud ceiling.  It uses the same expressions
c     for beta as in routine viscalc.  However, rather than using a
c     constant beta, the routine integrates upward using variable beta,
c     until an upward beam of light would be extinguished to .02 times
c     its original intensity. NOTE!  qvp is passed in in kg/kg, whereas
c     the hydrometeor mixing ratios are passed in as g/kg.
c
      dimension qcw(miy,mjx,mkzh),qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   qci(miy,mjx,mkzh),qsn(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),pf(miy,mjx,mkzh),clg(miy,mjx),
     &   prsnow(miy,mjx)
      character meth*1
c
      include 'comconst'
c
      tice=celkel-10.
      coeflc=144.7
      coeflp=2.24
      coeffc=327.8
      coeffp=clgfac*10.36
      exponlc=0.8800
      exponlp=0.7500
      exponfc=1.0000
      exponfp=0.7776
      thcon=.02
c
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
c
      exting=1.
      clg(i,j)=0.
      pup=pf(i,j,mkzh)
      do 150 k=mkzh,1,-1
         wmix=qvp(i,j,k)
         if (iice.eq.0) then
            qprc=qra(i,j,k)
            qcld=qcw(i,j,k)
            if (tmk(i,j,k).lt.celkel) then
               qrain=0.
               qsnow=qprc
               qclw=0.
               qclice=qcld
            else
               qrain=qprc
               qsnow=0.
               qclw=qcld
               qclice=0.
            endif
         else
            qprc=qra(i,j,k)+qsn(i,j,k)
            qcld=qcw(i,j,k)+qci(i,j,k)
            qrain=qra(i,j,k)
            qsnow=qsn(i,j,k)
            qclw=qcw(i,j,k)
            qclice=qci(i,j,k)
         endif
         tv=virtual(tmk(i,j,k),qvp(i,j,k))
         rhoair=100.*prs(i,j,k)/(rgas*tv)
         if (meth.eq.'d') then
            if (tmk(i,j,k).lt.celkel) then
               vovermd=(1.+wmix)/rhoair+.001*(qprc+qcld)/rhoice
               conclc = 0.
               conclp = 0.
               concfc = qcld/vovermd
               concfp = qprc/vovermd
            else
               vovermd=(1.+wmix)/rhoair+.001*(qprc+qcld)/rhowat
               conclc = qcld/vovermd
               conclp = qprc/vovermd
               concfc = 0.
               concfp = 0.
            endif
         elseif (meth.eq.'b') then
            if (tmk(i,j,k).lt.tice) then
               vovermd=(1.+wmix)/rhoair+.001*(qprc+qcld)/rhoice
               conclc = 0.
               conclp = 0.
               concfc = qcld/vovermd
               concfp = qprc/vovermd
            elseif (prsnow(i,j).ge.50..or.(prsnow(i,j).lt.50..and.
     &              k.ne.mkzh.and.tmk(i,j,k).lt.celkel)) then
               vovermd=(1.+wmix)/rhoair+.001*(qprc/rhoice+qcld/rhowat)
               conclc = qcld/vovermd
               conclp = 0.
               concfc = 0.
               concfp = qprc/vovermd
            else
               vovermd=(1.+wmix)/rhoair+.001*(qprc+qcld)/rhowat
               conclc = qcld/vovermd
               conclp = qprc/vovermd
               concfc = 0.
               concfp = 0.
            endif
         elseif (meth.eq.'r') then
            vovermd=(1.+wmix)/rhoair+.001*((qclw+qrain)/rhowat+
     &         (qclice+qsnow)/rhoice)
            conclc = qclw/vovermd
            conclp = qrain/vovermd
            concfc = qclice/vovermd
            concfp = qsnow/vovermd
         endif
         beta = coeffc * concfc ** exponfc + coeffp * concfp ** exponfp
     &        + coeflc * conclc ** exponlc + coeflp * conclp ** exponlp
     &        + 1.e-10
         pdn=pup
         if (k.eq.1) then
            pup=2.*prs(i,j,k)-pf(i,j,k)
         else
            pup=pf(i,j,k-1)
         endif
         dz=.001*rgas*tv/grav*log(pdn/pup)   ! m to km
         extingt=exting*exp(-beta*dz)
         if (extingt.le.thcon) then
            clg(i,j)=clg(i,j)-log(thcon/exting)/beta
            goto 160
         endif
         clg(i,j)=clg(i,j)+dz
         exting=extingt
  150 continue
  160 continue
      clg(i,j)=min(clg(i,j)*1000.,9000.)  ! km to m, max of 9,000m
c
  200 continue
      return
      end
