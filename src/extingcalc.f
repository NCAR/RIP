c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine extingcalc(qvp,qcw,qra,qci,qsn,tmk,prs,
     &   prsnow,meth,n,exting,miy,mjx,mkzh)
c
c     This routine computes extinction coefficients (per km) at the
c     lowest model layer due to cloud water (n=1), rain (2),
c     cloud ice (3), or snow (4).  It uses the same expressions for
c     extinction coefficient as in routine viscalc.
c
      dimension qcw(miy,mjx,mkzh),qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   qci(miy,mjx,mkzh),qsn(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),exting(miy,mjx),
     &   prsnow(miy,mjx)
      character meth*1
c
      include 'comconst'
c
      tice=celkel-10.
      coeflc=144.7
      coeflp=2.24
      coeffc=327.8
      coeffp=10.36
      exponlc=0.8800
      exponlp=0.7500
      exponfc=1.0000
      exponfp=0.7776
c
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
         wmix=qvp(i,j,mkzh)
         if (iice.eq.0) then
            qprc=qra(i,j,mkzh)
            qcld=qcw(i,j,mkzh)
            if (tmk(i,j,mkzh).lt.celkel) then
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
            qprc=qra(i,j,mkzh)+qsn(i,j,mkzh)
            qcld=qcw(i,j,mkzh)+qci(i,j,mkzh)
            qrain=qra(i,j,mkzh)
            qsnow=qsn(i,j,mkzh)
            qclw=qcw(i,j,mkzh)
            qclice=qci(i,j,mkzh)
         endif
         tv=virtual(tmk(i,j,mkzh),qvp(i,j,mkzh))
         rhoair=100.*prs(i,j,mkzh)/(rgas*tv)
         if (meth.eq.'d') then
            if (tmk(i,j,mkzh).lt.celkel) then
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
            if (tmk(i,j,mkzh).lt.tice) then
               vovermd=(1.+wmix)/rhoair+.001*(qprc+qcld)/rhoice
               conclc = 0.
               conclp = 0.
               concfc = qcld/vovermd
               concfp = qprc/vovermd
            elseif (prsnow(i,j).ge.50.) then
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
         if (n.eq.1) beta = coeflc * conclc ** exponlc
         if (n.eq.2) beta = coeflp * conclp ** exponlp
         if (n.eq.3) beta = coeffc * concfc ** exponfc
         if (n.eq.4) beta = coeffp * concfp ** exponfp
         exting(i,j)=beta  ! in km**-1
  200 continue
      return
      end
