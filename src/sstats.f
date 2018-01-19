c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine sstats(ilinw,prs,rslcg,icolr,icdwk,
     &   uuu,vvv,tmk,qvp,ipl,rslcgprv,unorth,vnorth,prssou,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,miy,mjx,mkzh)
c
c   Below is a table of all the parameters that show up in the lower
c   right corner of the sounding when this routine is called.
c
c   T is initial temperature of lifted parcel (deg. C)
c   LI is lifted index
c   K is K-index
c   TT is Total Totals Index
c   SWI is Showalter Index
c   PW is Precipitable Water (cm)
c   CAPE is Convective Available Potential Energy (J/kg)
c   CIN is Convective Inhibition (J/kg)
c   Tc is Convective Temperature (deg. C)
c   SREH is Storm-Relative Environmental Helicity (J/kg)
c   CELL is Cell Motion (Direction in degrees / Speed in knots)
c   Td is Dewpoint of lifted parcel (deg. C)
c   LCL is Lifted Condensation Level (hPa)
c   LFC is Level of Free Convection (hPa)
c   EL is Equilibrium Level (hPa)
c   CCL is Convective Condensation Level (hPa)
c   VGP is Vorticity Generation Potential (m/s^2)
c   SWEAT Sweat Index
c   HWBZ Wet-Bulb Zero Height (m)
c   SHEAR 0-6 km Shear (m/s)
c   LAPSE 700-500 hPa Lapse Rate (deg. C/km)
c
      dimension prs(miy,mjx,mkzh),rslcg(2,maxpl),ilinw(maxpl),
     &   icolr(maxpl),unorth(miy,mjx),vnorth(miy,mjx),
     &   uuu(miy,mjx,mkzh),vvv(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),
     &   rslcgprv(2),prssou(mkzh),icdwk(maxpl)
      real t(mkzh), r(mkzh), u(mkzh), v(mkzh), p(mkzh), rect(4)
      integer storm_motion
c
      character string*16
c
      include 'comconst'
c
      parameter (xm = 24.2)
c
c   Interpolate data to sounding location
c
      sxgn=1.+(rslcg(2,ipl)-xjcorn)*refrat
      sygn=1.+(rslcg(1,ipl)-yicorn)*refrat
      if (sxgn.le..5.or.sxgn.ge.mjx-.5.or.
     &    sygn.le..5.or.sygn.ge.miy-.5) then
         write(iup,*)'I don''t do soundings outside the'
         write(iup,*)'cross-point domain.'
	 write(iup,*) 'sxgn = ',sxgn,' sygn = ',sygn
         stop
      endif
c
      call gqclip (ierr,iclp,rect)
      call gsclip (1)
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
c
c make the vertical profiles, invert the columns
c
         posx=sxgn-.5
         posy=sygn-.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
c
         do 5 k=1,mkzh
	    kk = mkzh-k+1
	    p(k) = prssou(kk) * 100.     ! convert to Pascals
                 t(k)= (    (1.-ratlr)*(   ratbt)*tmk(it,jl,kk)+
     &                      (   ratlr)*(   ratbt)*tmk(it,jr,kk)+
     &                      (1.-ratlr)*(1.-ratbt)*tmk(ib,jl,kk)+
     &                      (   ratlr)*(1.-ratbt)*tmk(ib,jr,kk) )
                 r(k)= (    (1.-ratlr)*(   ratbt)*qvp(it,jl,kk)+
     &                      (   ratlr)*(   ratbt)*qvp(it,jr,kk)+
     &                      (1.-ratlr)*(1.-ratbt)*qvp(ib,jl,kk)+
     &                      (   ratlr)*(1.-ratbt)*qvp(ib,jr,kk) )
                 r(k)=max(r(k),.000000001) ! set to min value of 10^-6 g/kg
    5    continue
c
c   Make velocity components
c
      posx=sxgn-.5*icdwk(ipl)
      posy=sygn-.5*icdwk(ipl)
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      unorths= (
     &                (1.-ratlr)*(   ratbt)*unorth(it,jl)+
     &                (   ratlr)*(   ratbt)*unorth(it,jr)+
     &                (1.-ratlr)*(1.-ratbt)*unorth(ib,jl)+
     &                (   ratlr)*(1.-ratbt)*unorth(ib,jr) )
      vnorths= (
     &                (1.-ratlr)*(   ratbt)*vnorth(it,jl)+
     &                (   ratlr)*(   ratbt)*vnorth(it,jr)+
     &                (1.-ratlr)*(1.-ratbt)*vnorth(ib,jl)+
     &                (   ratlr)*(1.-ratbt)*vnorth(ib,jr) )
      do 8 k=1,mkzh
	 kk = mkzh-k+1
         uuus = (
     &                (1.-ratlr)*(   ratbt)*uuu(it,jl,kk)+
     &                (   ratlr)*(   ratbt)*uuu(it,jr,kk)+
     &                (1.-ratlr)*(1.-ratbt)*uuu(ib,jl,kk)+
     &                (   ratlr)*(1.-ratbt)*uuu(ib,jr,kk) )
         vvvs = (
     &                (1.-ratlr)*(   ratbt)*vvv(it,jl,kk)+
     &                (   ratlr)*(   ratbt)*vvv(it,jr,kk)+
     &                (1.-ratlr)*(1.-ratbt)*vvv(ib,jl,kk)+
     &                (   ratlr)*(1.-ratbt)*vvv(ib,jr,kk) )
         v(k)= (unorths*uuus+vnorths*vvvs)
         u(k)= (vnorths*uuus-unorths*vvvs)
    8 continue
c
       call sndg_anal (t,p,r,u,v,flminsou,frmaxsou,fbminsou,ftmaxsou,
     &   mkzh)
c
      call gsclip (iclp)
      return
      end
C+---+-----------------------------------------------------------------+
c
      subroutine sndg_anal(t,p,r,u,v,
     &   flminsou,frmaxsou,fbminsou,ftmaxsou,nlevs)
c
c..Subroutine to analyze stability, positive area, level info from a sounding
c
c..Local variables KLCL, KCCL, KLFC, KEL refer to the vertical index above
c..the levels of LCL, CCL, LFC, and EL respectively.  In other words, the
c..exact location of the LCL is between levels KLCL and KLCL-1 in the data.
c
      parameter (max_levs=100, zero=273.15)
      parameter (Rd=287.05, Cp=1004.0, G=9.8, XMS2KTS=1.94)
      dimension t(*),p(*),r(*),u(*),v(*),the(max_levs)
      character*128 label
      logical is_lfc
c
c..Initialize some stuff
c
      RCP = Rd/Cp
      CPR = Cp/Rd
      area_pos = 0.
      area_pos2 = 0.
      area_neg = 0.
      area_neg2 = 0.
      is_lfc = .false.
      prec_w = 0.
      plfc = 0.
      tlfc = 0.
      pel = 0.
      tel = 0.
      cape = 0.
      sreh = 0.
      vgpi = 0.
      alifted = 0.
      akindex = 0.
      totals = 0.
      sweat = 0.
c
c..First lets find all the vertical indicies of important levels
c..to make things easier later. (get array index corresponding to ??? hPa)
c
      k850 = 1
      k700 = 2
      k500 = 3
      do k=1,nlevs
         if(p(k).ge.85000.0) k850=k
         if(p(k).ge.70000.0) k700=k
         if(p(k).ge.50000.0) k500=k
         if(p(k).ge.30000.0) k300=k
      enddo
      k850 = max(1,k850)
      k700 = max(2,k700)
      k300 = min(k300,nlevs)
c
C+---+------------------------------------c
c..Compute Dewpoint, LCL pres, Theta, Theta-e, Wet Bulb Theta for each level
c
c      print11,' Temp | Pres |  W |  Th  | Th-e |  Td  |  RH |'
c     +        ,' T_LCL|P_LCL | Th-Wet '
 11   format(a,a)
      hwbz = -999.
      do k=1,nlevs

         theta = t(k) * (100000.0/p(k))**RCP
         dew_t = t_dew(p(k),r(k))
         tlcl = t_lcl(t(k),dew_t)
         plcl = 100000.0 * (tlcl/theta)**CPR
c write(6,*) 'k = ',k,' p = ',p(k),' t = ',t(k),' r = ',r(k),
c    &    ' tlcl = ',tlcl
         the(k)=theta_e(p(k),t(k),r(k),tlcl)
         rh = e_sub_s(dew_t) / e_sub_s(t(k))
         thwet = theta_wetb(the(k))
         twet = satlft (thwet, p(k))
	 if (hwbz .lt. 0. .and. twet .le. 0.) then
	   hwbz = std_atmos(p(k))
	   if (k .eq. 1) hwbz = 0.
	 endif
cc         print21, t(k)-zero,p(k)*0.01,r(k)*1000.,theta-zero,the(k)-zero
cc     +        ,dew_t-zero,rh*100.,tlcl-zero,plcl*0.01,thwet-zero
 21      format(f6.1,1x,f6.1,1x,f5.2,1x,f6.2,1x,f6.2,1x
     +        ,f6.2,1x,f5.1,1x,f6.2,1x,f6.1,1x,f6.1)
      enddo
c      print*
c
C+---+------------------------------------c
c..Compute precipitable water (meters) in column
c
      prec_w=prec_water(p(1),r(1),nlevs)
c
C+---+------------------------------------c
c..Compute mixed layer info using lowest 110 hPa data
c
      tot_pres = 0.
      tot_mixr = 0.
      tot_temp = 0.
      sum_delta = 0.
      avg_temp = t(1)
      avg_mixr = r(1)
      avg_pres = p(1)
      do k=2,nlevs
         if(p(k).ge.p(1)-11000.0) then
            delta_p = p(k-1) - p(k)
            sum_delta = sum_delta + delta_p
            tot_pres = tot_pres + delta_p*(p(k)+p(k-1))*0.5
            tot_mixr = tot_mixr + delta_p*(r(k)+r(k-1))*0.5
            tot_temp = tot_temp + delta_p*(t(K)+t(k-1))*0.5
         endif
      enddo
      if (sum_delta .gt. 0.0) then
         avg_temp = tot_temp / sum_delta
         avg_mixr = tot_mixr / sum_delta
         avg_pres = tot_pres / sum_delta
      endif
      avg_th = avg_temp * (100000.0/avg_pres)**RCP
      avg_t_dew = t_dew(avg_pres,avg_mixr)
c
C+---+------------------------------------c
c..Compute Lifted Index (C) using SELS method based on 100hPa thick
c..mixed layer temp + 2 degrees C
c Abandoned the +2 for consistency... JFB
c
c     avg_temp2 = avg_temp + 2.
      avg_temp2 = avg_temp
      avg_th2 = avg_temp2 * (100000.0/avg_pres)**RCP
      tlcl2 = t_lcl(avg_temp2,avg_t_dew)
      plcl2 = 100000.0 * (tlcl2/avg_th2)**CPR
      thelcl2 =  theta_e(plcl2,tlcl2,avg_mixr,tlcl2)
      parcel500 = compT_fr_The(thelcl2,50000.0,iup)

      t500 = t(k500) + ((t(k500+1)-t(k500))/(p(k500+1)-p(k500)))
     &                 * (50000.-p(k500))
      alifted = t500 - parcel500
c
c  compute showalter index
c
       if (k850 .eq. 1) then
	 t850 = t(k850)
	 r850 = r(k850)
       else
         t850 = t(k850) + ((t(k850+1)-t(k850))/(p(k850+1)-p(k850)))
     &                 * (85000.-p(k850))
         r850 = r(k850) + ((r(k850+1)-r(k850))/(p(k850+1)-p(k850)))
     &                 * (85000.-p(k850))
       endif
       th850 = t850 * (100000.0/85000.0)**RCP
       d850 = t_dew(85000.,r850)
       tlcl8 = t_lcl(t850,d850)
       plcl8 = 100000.0 * (tlcl8/th850)**CPR
       thelcl8 = theta_e(plcl8,tlcl8,r850,tlcl8)
       pcl500 = compT_fr_The(thelcl8,50000.0,iup)
       swi = t500 - pcl500
c
C+---+------------------------------------c
c..Find LCL, CCL, EL, and LFC from the mixed layer mixing data
c..These findings are exact.  One assumption hurting the
c..computations: theta-w is assumed to vary linearly between
c..two points on the environment temperature curve.  This is
c..used in y=mx + b fit to interpolate any temperature between
c..the two endpoints.  It is a reasonable assumption with
c..small delta-pressure spacing and in the lower levels but
c..gets more unreasonable as the spacing increases and
c..go higher.  Still, it beats the alternative - just assigning
c..the nearest level info.
c

c
c..First LCL - Lifting Condensation Level and save the
c..theta-e and theta-w from here cause we can check positive
c..and negative areas and LFC, EL based on these.
c
      tlcl = t_lcl(avg_temp,avg_t_dew)
      plcl = 100000.0 * (tlcl/avg_th)**CPR
      do k =1,nlevs
         if(p(k).le.plcl) goto 100
      enddo
 100  continue
      klcl = min(k,nlevs)
      thelcl =  theta_e(plcl,tlcl,avg_mixr,tlcl)
      thwlcl =  theta_wetb(thelcl)
c
c..Second, find CCL - Convective Condensation Level and 
c..thus, convective temperature.  Start at 500hPa and go
c..downward along the temperature curve until it crosses
c..the mixing ratio line that goes thru the LCL.  Do so by
c..finding the first time the mixing ratio of the environment
c..air (assuming saturation) exceeds the mixing ratio of
c..the mixed layer parcel (the same as the one going thru
c..the LCL).
c
      do k=k500,2,-1
         tmixr = r_sub_s(p(k),t(k))
         if(tmixr .ge. avg_mixr) goto 110
      enddo
 110  continue
      kccl = min(k+1,nlevs)
      tmixr = r_sub_s(p(kccl),t(kccl))
      if(abs(tmixr - r_sub_s(p(kccl-1),t(kccl-1))) .gt. 1.e-5 ) then       
         slope = (tmixr-avg_mixr) / (tmixr-r_sub_s(p(kccl-1),t(kccl-1)))
         if (slope .lt. 0.) slope = 0.
      else
         slope = 0.5
      endif
      tccl = t(kccl) - slope*(t(kccl)-t(kccl-1))
      pccl = p(kccl) - slope*(p(kccl)-p(kccl-1))
      thccl = tccl * (100000.0/pccl)**RCP
      conv_temp = min(323.15,thccl * (p(1)/100000.0)**RCP)
c
c..Third, find LFC - Level of Free Convection.  Start at 500hPa
c..and go downward along the temperature curve until it crosses
c..the moist adiabat which passes thru the LCL.  Do so by checking
c..the first time the theta-w of the environment air (using temp(k),
c..pres(k)) is greater than the theta-w of the parcel (moist adiabat
c..going thru LCL).
c
      the_scratch = theta_e(p(k500),t(k500)
     +       ,r_sub_s(p(k500),t(k500)),t(k500))
      thw_scratch = theta_wetb(the_scratch)
      if (thw_scratch .gt. thwlcl) then
         do k = k500-1, klcl, -1
            the_scratch = theta_e(p(k),t(k),r_sub_s(p(k),t(k)),t(k))
            thw_scratch = theta_wetb(the_scratch)
            if(thw_scratch .lt. thwlcl) goto 120
         enddo
 120     continue
      else
         is_lfc = .true.
         k = k500 -1
      endif
 
      do kk = k, klcl, -1
         the_scratch = theta_e(p(kk),t(kk),r_sub_s(p(kk),t(kk)),t(kk))
         thw_scratch = theta_wetb(the_scratch)
         if(thw_scratch .ge. thwlcl+0.7) goto 125
      enddo
 125  continue
      if (kk.le.klcl .and. (.not. is_lfc) ) goto 140
      klfc=kk+1
      scratch1 = theta_e(p(klfc),t(klfc)
     +                  ,r_sub_s(p(klfc),t(klfc)),t(klfc))
      scratch2 = theta_wetb(scratch1)
      if(thw_scratch .ne. scratch2) then
         slope = (scratch2-thwlcl) / (scratch2-thw_scratch)
      else
         slope = 0.5
      endif
      tlfc = t(klfc) - slope*(t(klfc)-t(klfc-1))
      plfc = p(klfc) - slope*(p(klfc)-p(klfc-1))
c
c..Fourth, find EL - Equilibrium Level.  Starting at the LFC, go up
c..along the temperature line until it re-crosses the moist adiabat
c..which passes thru the LCL.  Do so by checking the first time the
c..theta-w of the environment air (using temp(k),pres(k)) is less than
c..the theta-w of the parcel (moist adiabat going thru LCL).
c
      do k=klfc,nlevs
         the_scratch = theta_e(p(k),t(k),r_sub_s(p(k),t(k)),t(k))
         thw_scratch = theta_wetb(the_scratch)
         if( (thw_scratch .ge. thwlcl+1.5) .
     +        and. p(k).lt.plfc) goto 130
      enddo
c couldn't find an EL. Set it to the top level and move on.
      pel = p(nlevs)
      tel = t(nlevs)
      kel = nlevs
      go to 131
 130  continue
      kel=min(k,nlevs)
      scratch1 = theta_e(p(kel-1),t(kel-1)
     +                  ,r_sub_s(p(kel-1),t(kel-1)),t(kel-1))
      scratch2 = theta_wetb(scratch1)
      if(thw_scratch .ne. scratch2) then
         slope = (thw_scratch-thwlcl) / (thw_scratch-scratch2)
      else
         slope = 0.5
      endif
      tel = t(kel) - slope*(t(kel)-t(kel-1))
      pel = p(kel) - slope*(p(kel)-p(kel-1))
      if ( pel .lt. p(nlevs)) then
	pel = p(nlevs)
	tel = t(nlevs)
      endif
  131 continue
c
C+---+------------------------------------c
c..Now, compute CAPE, positive/negative areas from mixed layer parcel
c
      if ( (pel.lt.plfc) .and. (kel.gt.klfc) ) then
        npts=kel-klfc+1
        call compute_area(t(klfc),p(klfc),npts,thelcl,area_pos,area_neg)
        cape = area_pos
c
C+---+------------------------------------c
c..Compute Negative area to overcome (below LFC)
c
        call compute_area(t(1),p(1),klfc,thelcl,area_pos2,area_neg2)
      endif
c
C+---+------------------------------------c
c..Continue to here in the case where there is no LFC
c
 140  continue
c
C+---+------------------------------------c
c..Find storm relative helicity
c
      m_storm = 0
      if ( (plcl.gt.p(k300)) .and. (klcl.lt.k300) ) then
c First, compute a storm motion
	m_storm = storm_motion (p,u,v,nlevs)
        sreh = storm_helicity(p,u,v,m_storm,nlevs)
      endif
c
c+---+------------------------------------c
c  find 0-6km shear
c
      k6km = -999
      do k = 2,nlevs
        if (std_atmos(p(k)) .gt. 6000. .and. k6km .eq. -999) then
          k6km = k
        endif
      enddo
      shear = 1.94 * sqrt((u(k6km) - u(1))**2 + (v(k6km) - v(1))**2)
c
C+---+------------------------------------c
c..Find vorticity generation potential (VGP)
c
      vgpi = vgp_gt(p,u,v,cape,nlevs)
c
C+---+------------------------------------c
c..Compute K-index
c..K Index is defined as T @ 850 - T @ 500 + Td @ 850 - (T-Td)@700
c
      t700 = t(k700) + ((t(k700+1)-t(k700))/(p(k700+1)-p(k700)))
     &                 * (70000.-p(k700))
      r700 = r(k700) + ((r(k700+1)-r(k700))/(p(k700+1)-p(k700)))
     &                 * (70000.-p(k700))
      akindex = t850 - t500 + t_dew(85000.,r850)
     +     - (t700-t_dew(70000.,r700) ) -zero
c
c+---+------------------------------------c
c  find 700-500 lapse rate
c
      t500 = t(k500) + ((t(k500+1)-t(k500))/(p(k500+1)-p(k500)))
     &                 * (50000.-p(k500))
      xlaps = 1000. * (t700-t500) /
     &        (std_atmos(50000.) - std_atmos(70000.))
c
C+---+------------------------------------c
c..Compute Totals Totals Index
c..Totals is defined as T @ 850 + Td @ 850 - 2*T @ 500
c
      Totals = t850 + t_dew(85000.,r850) - 2.*t500
c
C+---+------------------------------------c
c..Compute Sweat Index
c..Sweat is defined as 12*D + 20*(TT-49) + 2*(f8) + f5 + 125*(S+0.2)
c..where: D is Td @ 850  (set to zero if < 0)
c..       TT is Totals Index  (set whole term to zero if TT-49 < 0)
c..       f8 is 850 hPa wind speed (kts)
c..       f5 is 500 hPa wind speed (kts)
c..       S  is (500-850 hPa wind direction)
c..           (set whole term =0 if one of:   wdir850 < 130 or > 250 degrees
c..                                           wdir500 < 210 or > 310 degrees
c..                                           wdir500-wdir850 < 0 degrees
c..                                           speed @ 850 or 500 < 15 kts.f)
c
      Td850 = t_dew(85000.,r850)-zero
      if(Td850 .le. 0.0) then
         first = 0.
      else
         first = 12.0 * Td850
      endif
      if(Totals-49.0 .lt. 0.0) then
         second = 0.
      else
         second = 20.0*(Totals-49.0)
      endif
      spd850 = sqrt( u(k850)*u(k850) + v(k850)*v(k850) ) * XMS2KTS
      dir850 = atan2(-u(k850),-v(k850))*180./3.14159
      if(dir850.lt.0.) dir850 = dir850 + 360.
      spd500 = sqrt( u(k500)*u(k500) + v(k500)*v(k500) ) * XMS2KTS
      dir500 = atan2(-u(k500),-v(k500))*180./3.14159
      if(dir500.lt.0.) dir500 = dir500 + 360.
      third = 2.0*spd850
      fourth = spd500
      if(dir850.lt.90. .or. dir850.gt.250. .or. dir500.lt.210.
     +     .or. dir500.gt.310. .or. (dir500-dir850).lt.0.
     +     .or. spd850.lt.7.717 .or. spd500.lt.7.717 ) then
         fifth = 0.
      else
         fifth = 125. * (sin((dir500-dir850)*3.14159/180.) + 0.2)
      endif

      Sweat = first + second + third + fourth + fifth
c
C+---+------------------------------------c
c..Print analysis info
c
c     print*, '     Precipitable water is ',prec_w*1000. ,' mm'
c     print*, '     Mixed layer:  temp = ',avg_temp-zero,' pres = '
c    +                   ,avg_pres/100.,' mix_ratio = ',avg_mixr*1000.
c     print*, '     This provides an LCL level of: ',plcl/100.,' hPa '
c    +                   ,' at temperature: ', tlcl-zero
c     print*, '     Temp, pressure at CCL  and Conv_Temp are: '
c    +                   ,tccl-zero,pccl*0.01,conv_temp-zero
c     print*, '     Temp, pressure at EL are: ',tel-zero,pel*0.01
c     print*, '     Temp, pressure at LFC are: ',tlfc-zero,plfc*0.01
c     print*, '     Lifted Index is: ',alifted
c     print*, '     K Index is: ',akindex
c     print*, '     Totals Totals Index is: ',Totals
c     print*, '     Sweat Index is: ',Sweat
c     print*, '     Positive/Negative Areas, Cape = '
c    +                   ,area_pos,area_neg,cape
c     print*, '     Negative Area to overcome = ',area_neg2
c     print*, '     Storm Relative Helicity value is ',sreh
c     print*, '     Showalter Index is ',swi
c     print*, '     Wet-bulb zero height is ',hwbz
c     print*
c
C+---+------------------------------------c
C+---+-----------------------------------------------------------------+
c..Put analysis info onto skewt plot

      space=3.
      call set(flminsou,frmaxsou,fbminsou,ftmaxsou,
     &   -19.0-space,27.1+space,-.9346217-.2*space,44.061+space,1)
c
c First, blank out the parcel info area
c
      sln = cufx(-19. + .1 )
c     srn = cufx(-4.5)
      srn = .31
      UTB=13.262683 -.1    ! y coord of the 500 hPa level
      utn = -0.93462171 + .1    ! y coord of the 1050 hPa level
      sbn = cufy(utn)
      stn = cufy(utb)
      CALL SET(SLN,SRN,SBN,STN,-1.,1.,-1.,1.,0)
      call ngdots(0.,0.,1,1000.,0)
      call set(flminsou,frmaxsou,fbminsou,ftmaxsou,
     &   -19.0-space,27.1+space,-.9346217-.2*space,44.061+space,1)
c     CALL SET(.05,.95,.05,.95,-19.0,27.1,-.9346217,44.061,1)
c     call plchhq(-18.5,fyjb(plcl*0.01),'LCL->',10.,0.,-1.)
      fyy = fyjb(plcl*0.01)
      fxx = t(klcl)+.5-273.15
c     write(6,*) 't(klcl) = ',t(klcl)-273.15,' zero = ',zero
      call plchhq(fxjb(fxx,fyy),fyy,'--LCL',10.,0.,-1.)
      fxx = t(klcl)+1.-273.15
      call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,-1.)
      fyy = fyjb(pccl*0.01)
      if ( abs(plcl-pccl) .lt. 2500. ) then
        fxx = tlcl-1.-273.15
        call plchhq(fxjb(fxx,fyy),fyy,'CCL--',10.,0.,1.)
        fxx = tlcl-1.5-273.15
        call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,1.)
      else
        fxx = tlcl+.5-273.15
        call plchhq(fxjb(fxx,fyy),fyy,'--CCL',10.,0.,-1.)
        fxx = tlcl+1.-273.15
	call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,-1.)
      endif
      if (plfc.gt.10000. .and. plfc.lt.100000.) then
	fyy = fyjb(plfc*0.01)
	if ( abs(plcl-plfc) .lt. 2000. ) then
	  fxx = tlfc+4.5-273.15
	  call plchhq(fxjb(fxx,fyy),fyy,'--LFC',10.,0.,-1.)
	  fxx = tlfc+5.-273.15
	  call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,-1.)
	else
	  fxx = tlfc+.5-273.15
	  call plchhq(fxjb(fxx,fyy),fyy,'--LFC',10.,0.,-1.)
	  fxx = tlfc+1.-273.15
	  call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,-1.)
	endif
	fyy = fyjb(pel*0.01)
	fxx = tel+.5-273.15
	call plchhq(fxjb(fxx,fyy),fyy,'--EL',10.,0.,-1.)
	fxx = tel+1.-273.15
	call plchhq(fxjb(fxx,fyy),fyy,'-',10.,0.,-1.)
      endif
      te_500 = compT_fr_The(thelcl,50000.,iup) - 273.15
      fy_500 = fyjb(500.)
      call plchhq(fxjb(te_500,fy_500),fy_500,'M',10.,0.,0.)

c     call set(0.01, 0.99, 0.01, 0.99, -21.0, 29.198, -3.0, 47.0, 1)
c     write(label,'(a,2x,a,4x,a,2x,a,3x,a,2x,a,2x,a,2x,a,3x,a,3x,a
c    +,2x,a,3x,a)') 
c    +'T(F)','Td','LI','SWT','K','TT','Pw(cm)','CAPE','Tc','CELL'
c    +,'SREH','VGP'
c     call plchhq(-19.0, 46.6, label, 9., 0., -1.)

c     write(label,'(i3,2x,i3,1x,f5.1,2x,i3,1x,i3,1x,i3,3x,f4.2
c    +              ,3x,i4,2x,i3,2x,i3.3,a1,i2.2,2x,i3,2x,f4.2)')
c    + nint((t(1)-zero)*9./5. + 32.)
c    + ,nint((t_dew(p(1),r(1))-zero)*9./5. + 32.)
c    + ,alifted
c    + ,nint(Sweat)
c    + ,nint(akindex)
c    + ,nint(Totals)
c    + ,prec_w*100.
c    + ,nint(cape)
c    + ,nint((conv_temp-zero)*9./5. + 32.)
c    + ,int(m_storm/100),'/',m_storm-(int(m_storm/100))*100
c    + ,nint(sreh)
c    + ,vgpi
c     call plchhq(-19.0, 45.4, label, 9., 0., -1.)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c     chsize=.008
      chsize=12.
      call pcgeti ('QU',ntextqq)
      call pcseti ('QU',1)     ! 1=medium quality equiv to plchmq
      call gsplci(1)
      xpos=.12
      ypos=.26
      write(label,'(a)') 'Parcel Info'
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      chsize=10.
      xpos=.07
      ypos=.23
      write(label,'(a,f7.1)') 'T  = ',(t(1)-zero)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f7.1)') 'LI = ',alifted
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i7)') 'K  = ',nint(akindex)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i7)') 'TT = ',nint(totals)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f6.1)') 'SWI = ',swi 
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f7.2)') 'PW = ',prec_w*100.
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i5)') 'CAPE = ',nint(cape)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i6)') 'CIN = ',nint(area_neg2)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f7.1)') 'Tc = ',conv_temp-zero
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i5)') 'SREH = ',nint(sreh)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i3.3,a1,i2.2)') 'CELL = ',
     &   int(m_storm/100),'/',m_storm-(int(m_storm/100))*100
      call plchhq(xpos,ypos,label,chsize,0.,-1.)

      xpos=.20
      ypos=.23
      write(label,'(a,f5.1)') 'Td = ',(t_dew(p(1),r(1))-zero)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.0)') 'LCL = ',plcl*.01
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.0)') 'LFC = ',plfc*.01
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.0)') 'EL  = ',pel*.01
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.0)') 'CCL = ',pccl*.01
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.1)') 'VGP = ',vgpi
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,i4)') 'SWEAT= ',nint(Sweat)
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f5.0)') 'HWBZ= ',hwbz
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f4.0)') 'SHEAR= ',shear
      call plchhq(xpos,ypos,label,chsize,0.,-1.)
      ypos=ypos-.02
      write(label,'(a,f4.1)') 'LAPSE= ',xlaps
      call plchhq(xpos,ypos,label,chsize,0.,-1.)

      call pcseti ('QU',ntextqq)
      call setusv('LW',1000)
      call gsplci(1)

      return
      end

C
      function theta_e(pres_Pa,temp_K,w_non,tlcl_K)
c..
c..         The following code was based on Bolton (1980) eqn #43
c..         and claims to have 0.3 K maximum error within -35 < T < 35 C
c..            pres_Pa = Pressure in Pascals
c..            temp_K  = Temperature in Kelvin
c..            w_non   = mixing ratio (non-dimensional = kg/kg)
c..            tlcl_K  = Temperature at Lifting Condensation Level (K)
c..
      pp = pres_Pa
      tt = temp_K
      rr = w_non + 1.e-8
c     write(6,*) 'rr = ',rr,' w_non = ',w_non
      tlc = tlcl_K
c
      power=(0.2854*(1.0 - (0.28*rr) ))
      xx = tt * (100000.0/pp)**power

      p1 = (3.376/tlc) - 0.00254
      p2 = (rr*1000.0) * (1.0 + 0.81*rr)
c     write(6,*) 'p1 = ',p1,' p2 = ',p2
      
      theta_e = xx * exp(p1*p2)

      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function t_lcl(temp_K,tdew_K)
c..
c..         The following code was based on Bolton (1980) eqn #15
c..         and claims to have 0.1 K maximum error within -35 < T < 35 C
c..            temp_K  = Temperature in Kelvin
c..            tdew_K  = Temperature at Lifting Condensation Level (K)
c..
      tt = temp_K
      tttd= tdew_K
      denom= ( 1.0/(tttd-56.0) ) + (log(tt/tttd)/800.)
      t_lcl = ( 1.0 / denom ) + 56.0
      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function t_dew(pres_Pa,w_non)
c..
c..            pres_Pa = Pressure in Pascals
c..            w_non   = mixing ratio (non-dimensional = kg/kg)
c..
      p = pres_Pa
      RR=w_non+1e-8
      ES=P*RR/(.622+RR)
      ESLN=LOG(ES)
      T_Dew=(35.86*ESLN-4947.2325)/(ESLN-23.6837)
      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function r_sub_s(pres_Pa,temp_K)
c..
c..         this calls function e_sub_s which computes saturation
c..         vapor pressure (Pa) and converts to sat. mixing ratio (kg/kg)
c..             pres_Pa - pressure (pa)
c..             temp_K  - temperature (k)
c..
      es=e_sub_s(temp_K)
      r_sub_s=.622*es/(pres_Pa-es)
c
      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function e_sub_s(temp_K)
c..
c..      compute saturation vapor pressure (Pa) over liquid with
c..      polynomial fit of goff-gratch (1946) formulation. (walko, 1991)
c..
      dimension c(0:8)
      data c/610.5851,44.40316,1.430341,.2641412e-1,.2995057e-3
     +      ,.2031998e-5,.6936113e-8,.2564861e-11,-.3704404e-13/

      x=max(-80.,temp_K-273.16)
      e_sub_s = c(0)+x*(c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)
     +              +x*(c(6)+x*(c(7)+x*c(8))))))))

      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function theta_wetb(thetae_K)
c..
c..              Eqn below was gotten from polynomial fit to data in
c..              Smithsonian Meteorological Tables showing Theta-e
c..              and Theta-w
c..
      real*8 c(0:6), d(0:6)
      data c/-1.00922292e-10, -1.47945344e-8, -1.7303757e-6
     +      ,-0.00012709,      1.15849867e-6, -3.518296861e-9
     +      ,3.5741522e-12/
      data d/0.00000000,   -3.5223513e-10, -5.7250807e-8
     +     ,-5.83975422e-6, 4.72445163e-8, -1.13402845e-10
     +     ,8.729580402e-14/

      x=min(475.0,thetae_K)
           
      if( x .le. 335.5 ) then
         answer = c(0)+x*(c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+
     +            x*c(6) )))))
      else
         answer = d(0)+x*(d(1)+x*(d(2)+x*(d(3)+x*(d(4)+x*(d(5)+
     +            x*d(6) )))))
      endif

      theta_wetb = answer + 273.15

      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function prec_water(pres_Pa, w_non, nlevels)
c..
c..            pres_Pa = Pressure in Pascals
c..            w_non   = mixing ratio (non-dimensional = kg/kg)
c..            returned precipitable water value in meters only below 300hPa
c..
      dimension pres_Pa(*), w_non(*)
      parameter (RHO_W=1000.0, G=9.8)

      sum=0.
      do k=nlevels,2,-1
         if(pres_Pa(k).gt.30000.0) goto 10
      enddo
 10   continue
      do kk=k,2,-1
         sum = sum + ( (w_non(kk)+w_non(kk-1))*0.5)
     +         * abs(pres_Pa(kk)-pres_Pa(kk-1))
      enddo

      prec_water = sum / (G*RHO_W)

      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function compT_fr_The(thelcl_K,pres_Pa,iup)
c..
c..            pres_Pa = Pressure in Pascals
c..            thelcl  = Theta-e at LCL (units in Kelvin)
c..
c..            Temperature (K) is returned given Theta-e at LCL
c..            and a pressure.  This describes a moist-adiabat.
c..            This temperature is the parcel temp at level Pres
c..            along moist adiabat described by theta-e.
c..
      guess= (thelcl_K - 0.5 * ( max(thelcl_K-270., 0.))**1.05)
     +     * (pres_Pa/100000.0)**.2
      epsilon=0.01
      do iter=1,100
         w1 = r_sub_s(pres_Pa,guess)
         w2 = r_sub_s(pres_Pa,guess+1.)
         tenu = theta_e(pres_Pa,guess,w1,guess)
         tenup = theta_e(pres_Pa,guess+1,w2,guess+1.)
         cor = (thelcl_K - tenu) / (tenup - tenu)
         guess = guess + cor
         if( (cor.lt.epsilon) .and. (-cor.lt.epsilon) ) then
            compT_fr_The=guess
            return
         endif
      enddo
      write(iup,*)
     &   'Convergence not reached in compT_fr_The, continuing.'
      thwlcl_K=theta_wetb(thelcl_K)
      compT_fr_The = thwlcl_K*((pres_Pa/100000.0)**0.286)

      return
      end
c+---+-----------------------------------------------------------------+
c
      subroutine compute_area(temp_K,pres_Pa,npts,thelcl
     +     ,area_pos,area_neg)
c..
c..            npts    = number of points in vertical column
c..            pres_Pa = Pressure in Pascals
c..            temp_K  = temperature in Kelvin
c..            thelcl  = Theta-e at LCL (units in Kelvin)
c..
c..            area_pos = positive area (J/Kg)
c..            area_neg = negative area (J/Kg)
c..
      parameter(P0=100000.0, RCP=0.2859, g=9.8, zero=273.15)
      dimension temp_K(*),pres_Pa(*)

      area = 0.0
      area_pos = 0.0
      area_neg = 0.0
c
c..Beginning at the LFC, go up to EL and sum positive and negative energies
c..using moist adiabatic lapse rate of parcel.  To do this find temp of parcel
c..given moist adiabat and pressure, then convert to theta then use eqn
c..for CAPE = g * (delta_theta/envirn_theta) * delta_height
c
      do k=1,npts-1
         p1        = pres_Pa(k)
         p2        = pres_Pa(k+1)
         avg_pres  = (p1+p2)*0.5
         z1        = std_atmos(p1)
         z2        = std_atmos(p2)
         delta_z   = abs(z2 - z1)
         th1       = temp_K(k)*(P0/pres_Pa(k))**RCP
         th2       = temp_K(k+1)*(P0/pres_Pa(k+1))**RCP
         envirn_TH = (th1+th2)*0.5
         parcel_T  = compT_fr_The(thelcl,avg_pres,iup)
         parcel_TH = parcel_T*(P0/avg_pres)**RCP
         delta_TH  = parcel_TH - envirn_TH

         area = g * (delta_TH/envirn_TH) * delta_z

cc         print*, 'envTH, delTH, delz, area: '
cc     +        ,envirn_TH,delta_TH,delta_z,area

         if(area.lt.0.0) then
            area_neg = area_neg + area
         else
            area_pos = area_pos + area
         endif
      enddo

      return
      end
c
c+---+-----------------------------------------------------------------+
c

c      function std_atmos(pres_Pa)
cc..
cc..            pres_Pa = Pressure in Pascals
cc..            standard atmos height in meters is returned for given p
cc..
c      pr = pres_Pa*0.01
c      height = 44307.692 * (1.0 - (pr/1013.25)**0.190)
c
c      std_atmos=height
c
c      return
c      end
c
c+---+-----------------------------------------------------------------+
c

      function std_atmos_p(height)
c..
c..            standard atmos pressure in Pascals is returned for given
c..            height in meters
c..

      std_atmos_p = exp(log(1.0-height/44307.692)/0.19)*101325.0

      return
      end
c
c+---+-----------------------------------------------------------------+
c
      function storm_motion (p,u,v,nlevs)
      dimension u(nlevs), v(nlevs), p(nlevs)
      parameter(pi180 = 3.14159/180.0)
c..Estimate a storm motion from cloud base to 10 km wind info
c     k6km = -999
      do k = 2,nlevs
c       if (std_atmos(p(k)) .gt. 6000. .and. k6km .eq. -999) then
c         k6km = k
c       endif
	if (std_atmos(p(k)) .gt. 10000.) then
	  k10km = k
	  go to 10
	endif
      enddo
      k10km = nlevs - 1
   10 delta_z = 0.
      sum_deltaz = 1.e-4
      su = 0.
      sv = 0.
      do kp = 2, k10km
	delta_z = std_atmos(p(kp)) - std_atmos(p(kp-1))
	sum_deltaz = sum_deltaz + delta_z
	su = su + 0.5*delta_z*(u(kp)+u(kp-1))
	sv = sv + 0.5*delta_z*(v(kp)+v(kp-1))
      enddo
      ua = su/sum_deltaz
      va = sv/sum_deltaz
      asp = sqrt(ua*ua + va*va)
      if (ua.eq.0. .and. va.eq.0.) then
         adr = 0.
      else
         adr = (1./pi180) * (3.14159 + atan2(ua,va))
      endif
c     write(6,*) 'adr = ',adr,' asp = ',asp
c     write(6,*) 'k10km = ',k10km,' k6km = ',k6km
      storm_motion = nint(adr)*100 + asp*1.94
      storm_motion = max(0.,storm_motion)
c     write(6,*) 'storm_motion = ',storm_motion
c..bsp and bdr are guesses for storm motion (75% of mean wind)
c..and 30 deg to the right
c     bsp = 0.75 * asp
c     bdr = adr + 30.
c     if (bdr.ge.360.) bdr = bdr-360.
      return
      end
c
c+---+-----------------------------------------------------------------+
c

      function storm_helicity(p,u,v,m_storm,nlevs)
c..
c..            p = pressure (Pascals) and u,v are wind components
c..            an entire vertical column needs to be provided
c..            m_storm is storm motion integer such as 24022
c..            meaning from 240 degrees at 22 knots.
c..
      dimension u(nlevs), v(nlevs), p(nlevs)
      parameter(pi180 = 3.14159/180.0)

c..Estimate a storm motion from cloud base to 10 km wind info
cc      delta_z = 0.
cc      sum_deltaz = 1.e-4
cc      su = 0.
cc      sv = 0.
cc      do kp = kbot+1, k300
cc         if (u(kp).ne.-99.9 .and. u(kp-1).ne.-99.9 .
cc     +        and. v(kp).ne.-99.9 .and. v(kp-1).ne.-99.9) then
cc            delta_z = std_atmos(p(kp)) - std_atmos(p(kp-1))
cc            sum_deltaz = sum_deltaz + delta_z
cc            su = su + 0.5*delta_z*(u(kp)+u(kp-1))
cc            sv = sv + 0.5*delta_z*(v(kp)+v(kp-1))
cc         endif
cc      enddo
cc      ua = su/sum_deltaz
cc      va = sv/sum_deltaz
cc      asp = sqrt(ua*ua + va*va)
cc      if (ua.eq.0. .and. va.eq.0.) then
cc         adr = 0.
cc      else
cc         adr = (1./pi180) * (3.14159 + atan2(ua,va))
cc      endif
c..bsp and bdr are guesses for storm motion (75% of mean wind)
c..and 30 deg to the right
cc      bsp = 0.75 * asp
cc      bdr = adr + 30.
cc      if (bdr.ge.360.) bdr = bdr-360.

      bdr = int(m_storm/100)
      bsp = (m_storm - (int(m_storm/100))*100)/1.94

c     write(6,*) '     Estimated storm motion ',bdr,' degrees at '
c    +       ,bsp,'m/s'

c..Now estimate storm relative helicity for 0-3 km AGL
      cu = -bsp * sin(bdr*pi180)
      cv = -bsp * cos(bdr*pi180)
      sum = 0.
      zbot = std_atmos(p(1))
      z3km = zbot + 3000.
      do k = 1,nlevs
         if (std_atmos(p(k)).gt.z3km) then
	   k3km = k
	   goto 22
	 endif
      enddo
   22 k3km = min(k3km,nlevs)

c     if ( u(1).ne.-99.9 .and. v(1).ne.-99.9 ) then
         x = ((u(1)-cu)*v(1)) - ((v(1)-cv)*u(1))
         sum = sum + x
c     endif
      do kk = 2,k3km
c        if (u(kk).ne.-99.9 .and. u(kk-1).ne.-99.9 .
c    +        and. v(kk).ne.-99.9 .and. v(kk-1).ne.-99.9) then
         x = ((u(kk)-cu)*(v(kk)-v(kk-1))) - ((v(kk)-cv)*(u(kk)-u(kk-1)))
         sum = sum + x
c        endif
      enddo

      storm_helicity = -sum
      storm_helicity = max(0.0, storm_helicity)

      return
      end
c
c+---+-----------------------------------------------------------------+
c

      function vgp_gt(p,u,v,cape,nlevs)
c..
c..            p = pressure (Pascals) and u,v are wind components.
c..            entire vertical column is passed but only sfc to
c..            3 km is used.
c..
      dimension u(nlevs), v(nlevs), p(nlevs)

c..Find k-index for level nearest 3 km AGL
      sum = 0.
      zbot = std_atmos(p(1))
      z3km = zbot + 3000.
      do k = 1,nlevs
         if (std_atmos(p(k)).gt.z3km) goto 100
      enddo
 100  continue
      k3km = min(k,nlevs)
      if (k3km.eq.nlevs) return

c..Now calculate 0-3 km AGL hodograph length
      depth = std_atmos(p(k3km)) - std_atmos(p(1))
      do kk = 2,k3km
c        if (u(kk).ne.-99.9 .and. u(kk-1).ne.-99.9 .
c    +        and. v(kk).ne.-99.9 .and. v(kk-1).ne.-99.9) then
            su = abs(u(kk)-u(kk-1))
            sv = abs(v(kk)-v(kk-1))
            sum = sum + sqrt(su*su+sv*sv)
c        endif
      enddo
      sum = sum/depth
      vgp_gt = sqrt(cape)*sum

      return
      end
C
      FUNCTION FYJB(P)
      FYJB = 132.182-44.061*ALOG10(P)
      RETURN
      END
C
      FUNCTION FXJB(T,Y)
      FXJB = 0.54*T+0.90692*Y
      RETURN
      END
c-------------------------------------------------------------
      function satlft (thw1, p1)
      parameter (cta = 273.16, akap = 0.28541)
      thw = thw1 - 273.15
      p = p1 * .01
      if (p .ne. 1000.) go to 5
      satlft = thw
      return
    5 continue
      pwrp = (p / 1000.) ** akap
      tone = (thw + cta) * pwrp - cta
      eone = wobf(tone) - wobf(thw)
      rate = 1.
      go to 15
   10 continue
      rate = (ttwo - tone) / (etwo - eone)
      tone = ttwo
      eone = etwo
   15 continue
      ttwo = tone - eone * rate
      pt = (ttwo + cta) / pwrp - cta
      etwo = pt + wobf(ttwo) - wobf(pt) - thw
      dlt = etwo * rate
      if (abs(dlt) .gt. 0.1) go to 10
      satlft = ttwo - dlt
      return
      end
c-------------------------------------------------------------
      function wobf (t)
      x = t - 20.
      if (x .gt. 0.) go to 10
      pol = 1.                     + x * (-8.8416605e-03
     &     + x * ( 1.4714143e-04   + x * (-9.6719890e-07
     &     + x * (-3.2607217e-08   + x * (-3.8598073e-10)))))
      wobf = 15.130 / pol ** 4
      return
   10 continue
      pol = 1.                     + x * ( 3.6182989e-03
     &     + x * (-1.3603273e-05   + x * ( 4.9618922e-07
     &     + x * (-6.1059365e-09   + x * ( 3.9401551e-11
     &     + x * (-1.2588129e-13   + x * ( 1.6688280e-16)))))))
      wobf = 29.930 / pol ** 4 + 0.96 * x - 14.8
      return
      end
