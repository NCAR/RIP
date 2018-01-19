c------------------------------------------------------------
      subroutine lifted_index(tmk, qvp, prs, li,
     &   miy,mjx,mkzh)
c
c   Compute LI based on the method in sstats. 
c   sstats functions are used here.   3/31/15
c   T is initial temperature of lifted parcel (deg. C)
c   LI is lifted index
c   Td is Dewpoint of lifted parcel (deg. C)
c
      dimension prs(miy,mjx,mkzh), tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh)
      real t(mkzh), r(mkzh), p(mkzh), li(miy,mjx)
c
      include 'comconst'
c
c make the vertical profiles, invert the columns
c
         do 10 i=1,miy-1
         do 10 j=1,mjx-1
         do 5 k=1,mkzh
	    kk = mkzh-k+1
	    p(k) = prs(i,j,kk) * 100.     ! convert to Pascals
                 t(k)= tmk(i,j,kk)
                 r(k)= qvp(i,j,kk)
                 r(k)=max(r(k),.000000001) ! set to min value of 10^-6 g/kg
    5    continue
c
       call calc_li (t,p,r,li(i,j),mkzh)
   10  continue
c
      return
      end
c--------------------------------------------------------------------
      subroutine calc_li(t,p,r,alifted,nlevs)
c
      parameter (max_levs=100)
      parameter (Rd=287.05, Cp=1004.0, G=9.8 )
      dimension t(*),p(*),r(*),the(max_levs)
c
c..Initialize some stuff
c
      RCP = Rd/Cp
      CPR = Cp/Rd
      alifted = 0.
c
      k500 = 3
      do k=1,nlevs
         if(p(k).ge.50000.0) k500=k
      enddo
c     write(6,*) 'k500 = ',k500
c
c..Compute Dewpoint, LCL pres, Theta, Theta-e for each level
c
      do k=1,nlevs

         theta = t(k) * (100000.0/p(k))**RCP
         dew_t = t_dew(p(k),r(k))
         tlcl = t_lcl(t(k),dew_t)
         plcl = 100000.0 * (tlcl/theta)**CPR
c write(6,*) 'k = ',k,' p = ',p(k),' t = ',t(k),' r = ',r(k),
c    &    ' tlcl = ',tlcl
         the(k)=theta_e(p(k),t(k),r(k),tlcl)
      enddo
c      print*
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
c     write(6,*) 'tlcl2 = ',tlcl2,' plcl2 = ',plcl2,' avg_temp2 = ',
c    & avg_temp2
      parcel500 = compT_fr_The(thelcl2,50000.0,iup)

      t500 = t(k500) + ((t(k500+1)-t(k500))/(p(k500+1)-p(k500)))
     &                 * (50000.-p(k500))
      alifted = t500 - parcel500
      return
      end
