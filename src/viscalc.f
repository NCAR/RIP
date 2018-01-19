c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine viscalc(qvp,qcw,qra,qci,qsn,tmk,prs,prsnow,
     &   meth,vis,miy,mjx,mkzh)
c
c     This routine computes horizontal visibility at the lowest model
c     layer, from qcw, qra, qci, and qsn.  NOTE!  qvp is passed in in
c     kg/kg, whereas the hydrometeor mixing ratios are passed in as
c     g/kg.  At each grid point, the routine assigns values of total
c     liq. eq. precip, total liq.  equiv. cloud, rain, snow, cloud
c     water, and cloud ice, based on the value of iice.
c
c   If iice=0:
c
c      qprc=qra     qrain=qra and qclw=qcw if T>0C
c      qcld=qcw          =0          =0  if T<0C
c                  qsnow=qsn and qclice=qcw  if T<0C
c                       =0            =0   if T>0C
c   If iice=1:
c
c      qprc=qra+qsn   qrain=qra and qclw=qcw
c      qcld=qcw+qci   qsnow=qsn and qclice=qcw
c   
c   Independent of the above definitions, the scheme can use different
c   assumptions of the state of hydrometeors:
c        meth='d': qprc is all frozen if T<0, liquid if T>0
c        meth='b': Bocchieri scheme used to determine whether qprc
c           is rain or snow. A temperature assumption is used to
c           determine whether qcld is liquid or frozen.
c        meth='r': Uses the four mixing ratios qrain, qsnow, qclw,
c           and qclice
c
c   The routine uses the following
c   expressions for extinction coefficient, beta (in km**-1),
c   with C being the mass concentration (in g/m**3):
c
c      cloud water:  beta = 144.7 * C ** (0.8800)
c      rain water:   beta =  2.24 * C ** (0.7500)
c      cloud ice:    beta = 327.8 * C ** (1.0000)
c      snow:         beta = 10.36 * C ** (0.7776)
c
c   These expressions were obtained from the following sources:
c
c      for cloud water: from Kunkel (1984)
c      for rainwater: from M-P dist'n, with No=8e6 m**-4 and
c         rho_w=1000 kg/m**3
c      for cloud ice: assume randomly oriented plates which follow
c         mass-diameter relationship from Rutledge and Hobbs (1983)
c      for snow: from Stallabrass (1985), assuming beta = -ln(.02)/vis
c
c   The extinction coefficient for each water species present is
c   calculated, and then all applicable betas are summed to yield
c   a single beta. Then the following relationship is used to
c   determine visibility (in km), where epsilon is the threshhold
c   of contrast, usually taken to be .02:
c
c      vis = -ln(epsilon)/beta      [found in Kunkel (1984)]
c
      dimension qcw(miy,mjx,mkzh),qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   qci(miy,mjx,mkzh),qsn(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),vis(miy,mjx),prsnow(miy,mjx)
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
      const1=-log(.02)
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
         beta = coeffc * concfc ** exponfc + coeffp * concfp ** exponfp
     &        + coeflc * conclc ** exponlc + coeflp * conclp ** exponlp
     &        + 1.e-10
         vis(i,j)=min(90.,const1/beta)  ! max of 90km
  200 continue
      return
      end
