c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine precprob(tmk,twb,ght,ter,prcprb,
     &   miy,mjx,mkzh)
c
c   This routine calculates the probability of liquid, freezing and ice
c      precip at a grid point based on sounding characteristics,
c      following Bocchieri (1980)
c
c   level:    liquid: 1;   freezing: 2;   ice: 3
c
      parameter (mkz=87,nuor=12)
c
      dimension tmk(miy,mjx,mkzh),twb(miy,mjx,mkzh),ght(miy,mjx,mkzh),
     &   prcprb(miy,mjx,3),ter(miy,mjx)
c
      dimension tmps(100),twbs(100),ghts(100),ghxs(100)
      dimension zz(mkz),zzh(mkz-1),dz(mkz-1),tmpz(mkz-1),twbz(mkz-1)
      dimension lfi(nuor),cx(nuor),cxr(nuor)
      dimension al(nuor),af(nuor),ai(nuor)
c
      include 'comconst'
c
c      qww(g,e)=3.799*exp(17.269*(e-273.16)/(e-35.86))/g
      data cx/ -1.00,  -6.00, 150.00, 300.00, 750.00, 300.00,
     &        550.00, 1200.0,  50.00,  50.00,  800.0,4000.00/
      data bl,bf,bi/20.19, 63.35, 16.46/
      data al/-17.93,  -5.83, -28.48, -30.27, -16.04, +34.90,
     &         -6.38,  +6.55, +23.83, +20.95,  -5.88,  +5.04/
      data af/ +7.99,  -6.82,  +4.88,  -4.56,  -0.13, +18.72,
     &        -22.35, -13.84, -27.23, -25.82, -15.60, +22.96/
      data ai/ +9.94, +12.65, +23.60, +34.83, +16.17, -53.62,
     &        +28.73,  +7.29,  +3.40,  +4.87, +21.49, -28.01/
c
c   Create zz,dz
c
      do 5 k=1,21    ! Every 10 m from 0 to 200 m
         zz(k)=(k-1)*10.
    5 continue
      do 6 k=22,37    ! Every 50 m from 200 to 1000 m
         zz(k)=200.+(k-21)*50.
    6 continue
      do 7 k=38,87    ! Every 100 m from 1000 to 6000 m
         zz(k)=1000.+(k-37)*100.
    7 continue
      do 8 k=1,mkz-1
         zzh(k)=(zz(k+1)+zz(k))/2.
         dz(k)=zz(k+1)-zz(k)
    8 continue
c
      do 1000 j=1,mjx-1
      do 1000 i=1,miy-1
c
c   Make column variables
c
      do 21 k=1,mkzh
         tmps(k)=tmk(i,j,k)-celkel
         twbs(k)=twb(i,j,k)
         ghts(k)=ght(i,j,k)
         ghxs(k)=exp(-ghts(k)/sclht)
  21  continue
c
c   Calculate low-level lapse rates
c
      dtdzlow=(tmps(mkzh)-tmps(mkzh-1))/(ghts(mkzh-1)-ghts(mkzh))
      dtwbdzlow=(twbs(mkzh)-twbs(mkzh-1))/(ghts(mkzh-1)-ghts(mkzh))
c
c   Interpolate tmp, twb to height levels.
c
      do 50 kz=1,mkz-1
         zlev=ter(i,j)+zzh(kz)
         zlevx=exp(-(zlev)/sclht)
         if (zlevx.gt.ghxs(mkzh)) then
            tmpz(kz)=tmps(mkzh)+dtdzlow*(ghts(mkzh)-zlev)
            twbz(kz)=min(twbs(mkzh)+dtwbdzlow*(ghts(mkzh)-zlev),
     &         tmpz(kz))
            goto 40
         endif
         do 30 ks=1,mkzh-1
            vcp1=ghxs(ks+1)
            vcp0=ghxs(ks)
            if (zlevx.ge.vcp0.and.zlevx.le.vcp1) then
               valp0=tmps(ks)
               valp1=tmps(ks+1)
               tmpz(kz)=(zlevx-vcp0)*(valp1-valp0)/(vcp1-vcp0)+
     &                    valp0
               valp0=twbs(ks)
               valp1=twbs(ks+1)
               twbz(kz)=(zlevx-vcp0)*(valp1-valp0)/(vcp1-vcp0)+
     &                    valp0
               goto 40
            endif
   30    continue
         write(iup,*)'Didn''t find z-level in precprob. kz=',kz
         stop
   40    continue
   50 continue
c
c   Determine layer-average temperatures, warm-layer depth, warm-layer
c   area, cold-wet-bulb-layer depth, and cold-wet-bulb-layer area.
c
      tsrfto10=0.
      t5to25=0.
      dzwarm=0.
      dawarm=0.
      dzcoldw=0.
      dacoldw=0.
      lwarmw=0
      do 100 k=1,mkz-1
        if (zzh(k).lt.1000.) tsrfto10=tsrfto10+tmpz(k)*dz(k)
        if (zzh(k).gt.500..and.zzh(k).lt.2500.)
     &     t5to25=t5to25+tmpz(k)*dz(k)
	if (tmpz(k).gt.0.) then
           dzwarm=dzwarm+dz(k)
           dawarm=dawarm+dz(k)*tmpz(k)
        endif
        if (twbz(k).le.0.) then
           if (lwarmw.eq.0) then
              dzcoldw=dzcoldw+dz(k)
              dacoldw=dacoldw-dz(k)*twbz(k)
           endif
        else
           lwarmw=1
        endif
  100 continue
      tsrfto10=.001*tsrfto10
      t5to25=.0005*t5to25
      zr=0
      if (tmpz(1).lt.0..and.dzwarm.gt.0.) zr=1
      dzcoldw=dzcoldw*lwarmw
      dacoldw=dacoldw*lwarmw
c
c   Assign values to predictors
c
      cxr(1)=tsrfto10
      cxr(2)=t5to25
      cxr(3)=dzwarm
      cxr(4)=dzwarm
      cxr(5)=dawarm
      cxr(6)=zr*dzwarm
      cxr(7)=zr*dzwarm
      cxr(8)=zr*dzwarm
      cxr(9)=zr*dawarm
      cxr(10)=dzcoldw
      cxr(11)=dzcoldw
      cxr(12)=dacoldw
c
c   Determine the lfi based on criteria for each predictor
c
      do 260 li=1,nuor
         lfi(li)=0
         if(cxr(li).le.cx(li)) lfi(li)=1
  260 continue
c
c   Calculate the probabilities from the regression equations
c
      prcprb(i,j,1)=bl
      prcprb(i,j,2)=bf
      prcprb(i,j,3)=bi
      do 270 l=1,nuor
         prcprb(i,j,1)=prcprb(i,j,1)+al(l)*lfi(l)
         prcprb(i,j,2)=prcprb(i,j,2)+af(l)*lfi(l)
         prcprb(i,j,3)=prcprb(i,j,3)+ai(l)*lfi(l)
  270 continue
      prbtot=0.
      do 280 it=1,3
         prcprb(i,j,it)=min(max(prcprb(i,j,it),0.),100.)
         prbtot=prbtot+prcprb(i,j,it)
  280 continue
      prbfac=100./prbtot
      do 290 it=1,3
         prcprb(i,j,it)=prbfac*prcprb(i,j,it)
  290 continue
c
 1000 continue
c
      return
      end
