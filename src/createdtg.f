c                                                                     c
c*********************************************************************c
c                                                                     c
      SUBROUTINE createdtg(mdateb,rhourb,xtime,timezone,
     &   iusdaylightrule,inearesth,dtgi,dtgfz,dtgfl)
c
c  Creates labels for initial time in UTC and forecast time in UTC and
c     local time.  Adapted from Ernie Recker's create_dtg routine.
c
      CHARACTER*22 dtgi,dtgfz,dtgfl
      CHARACTER name(12)*3,dayofweek(7)*3,czone*3
c
      data name /'Jan','Feb','Mar','Apr','May','Jun',   
     &          'Jul','Aug','Sep','Oct','Nov','Dec' /
      data dayofweek /'Sun','Mon','Tue','Wed','Thu','Fri','Sat'/
c
      tiny=.0001  ! add this amount (equiv. to .36 s) so that time
c                 ! will round up to next minute if very close to it
c
c   First, initial time in UTC
c
      md=mdateb
      rh=rhourb+tiny
      if (rh.lt.0.) then
         iadd=-int(-rh)-1
         rh=rh-float(iadd)
      elseif (rh.ge.1.) then
         iadd=int(rh)
         rh=rh-float(iadd)
      else
         iadd=0
      endif
      if (inearesth.eq.1.and.rh.gt..5+2.*tiny) iadd=iadd+1
      call mconvert(md,mh,1,1940)
      mh=mh+iadd
      call mconvert(md,mh,-1,1940)
      iyy = md/1000000
      imm = md/10000 - iyy*100
      idd = md/100 - (iyy*10000 + imm*100)
      ihh = md - (iyy*1000000 + imm*10000 + idd*100)
      imin=int(rh*60.)
      idow=mod(mh/24+1,7)+1
      czone='UTC'
      if (inearesth.eq.0) then
         WRITE (dtgi,501)
     &      ihh,imin,czone,dayofweek(idow),idd,name(imm),iyy
      else
         WRITE (dtgi,502)
     &      ihh,czone,dayofweek(idow),idd,name(imm),iyy
      endif
  501 FORMAT (i2.2,i2.2,1x,a3,1x,a3,1x,i2.2,1x,a3,1x,i2.2)
  502 FORMAT (i2.2,1x,a3,1x,a3,1x,i2.2,1x,a3,1x,i2.2,2x)
c
c   Next, valid time in UTC
c
      md=mdateb
      rh=rhourb+xtime+tiny
      if (rh.lt.0.) then
         iadd=-int(-rh)-1
         rh=rh-float(iadd)
      elseif (rh.ge.1.) then
         iadd=int(rh)
         rh=rh-float(iadd)
      else
         iadd=0
      endif
      if (inearesth.eq.1.and.rh.gt..5+2.*tiny) iadd=iadd+1
      call mconvert(md,mh,1,1940)
      mh=mh+iadd
      call mconvert(md,mh,-1,1940)
      iyy = md/1000000
      imm = md/10000 - iyy*100
      idd = md/100 - (iyy*10000 + imm*100)
      ihh = md - (iyy*1000000 + imm*10000 + idd*100)
      imin=int(rh*60.)
      idow=mod(mh/24+1,7)+1
      czone='UTC'
      if (inearesth.eq.0) then
         WRITE (dtgfz,501)
     &      ihh,imin,czone,dayofweek(idow),idd,name(imm),iyy
      else
         WRITE (dtgfz,502)
     &      ihh,czone,dayofweek(idow),idd,name(imm),iyy
      endif
c
c   Next, valid time in local time
c
      md=mdateb
      rh=rhourb+xtime+timezone+tiny
      if (rh.lt.0.) then
         iadd=-int(-rh)-1
         rh=rh-float(iadd)
      elseif (rh.ge.1.) then
         iadd=int(rh)
         rh=rh-float(iadd)
      else
         iadd=0
      endif
      if (inearesth.eq.1.and.rh.gt..5+2.*tiny) iadd=iadd+1
      call mconvert(md,mh,1,1940)
      mh=mh+iadd
      call mconvert(md,mh,-1,1940)
      iyy = md/1000000
c
c   Daylight saving
c
      idsv=0
      if (iusdaylightrule.eq.0) goto 39
c
c   Which ccyy are we dealing with 
c   Need this to adjust the daylight saving rule
      md1=iyy*1000000+10102  ! 2 am LST, 1 Jan of current year
      call mconvertccyy(md1,1940,iyear)

      if ( iyear < 2007 ) then 
         md1=iyy*1000000+40102  ! 2 am LST, 1 Apr of current year
         call mconvert(md1,mh1,1,1940)
         idow1=mod(mh1/24+1,7)+1
         if (idow1.gt.1) then
            mdadd=8-idow1
            mh1=mh1+24*mdadd ! 2 am LST, first Sunday in April
         endif
         md2=iyy*1000000+103101  ! 1 am LST (2 am LDT), 31 Oct of current year
         call mconvert(md2,mh2,1,1940)
         idow2=mod(mh2/24+1,7)+1
         if (idow2.gt.1) then
            mdsub=idow2-1
            mh2=mh2-24*mdsub ! 1 am LST (2 am LDT), last Sunday in October
         endif
      else
c         md1=iyy*1000000+30102  ! 2 am LST, 1 Mar of current year
         md1=iyy*1000000+30802  ! 2 am LST, 8 Mar of current year
         call mconvert(md1,mh1,1,1940)
         idow1=mod(mh1/24+1,7)+1
         if (idow1.gt.1) then
c            mdadd=15-idow1
            mdadd=8-idow1
            mh1=mh1+24*mdadd ! 2 am LST, second Sunday in March
         endif
         md2=iyy*1000000+110101  ! 1 am LST (2 am LDT), 1 Nov of current year
         call mconvert(md2,mh2,1,1940)
         idow2=mod(mh2/24+1,7)+1
         if (idow2.gt.1) then
            mdadd=8-idow2
            mh2=mh2+24*mdadd ! 1 am LST (2 am LDT), first Sunday in November
         endif
      endif
 
c
      if (mh.ge.mh1.and.mh.lt.mh2) then !daylight savings time
         idsv=1
         mh=mh+1
         call mconvert(md,mh,-1,1940)
         iyy = md/1000000
         imm = md/10000 - iyy*100
         idd = md/100 - (iyy*10000 + imm*100)
         ihh = md - (iyy*1000000 + imm*10000 + idd*100)
      endif
 39   continue
c
      call mconvert(md,mh,-1,1940)
      iyy = md/1000000
      imm = md/10000 - iyy*100
      idd = md/100 - (iyy*10000 + imm*100)
      ihh = md - (iyy*1000000 + imm*10000 + idd*100)
      imin=int(rh*60.)
      idow=mod(mh/24+1,7)+1
      if ( timezone .eq. -8. ) then
         czone = 'PST'
      else if ( timezone .eq. -7. ) then
         czone = 'MST'
      else if ( timezone .eq. -6. ) then
         czone = 'CST'
      else if ( timezone .eq. -5. ) then
         czone = 'EST'
      else if ( timezone .eq. -9. ) then   ! Alaska
         czone = 'AKT'
      else if ( timezone .eq. -10. ) then   ! Hawaii
         czone = 'HST'
      else if ( timezone .eq. 10. ) then    ! Guam
         czone = 'ChT'
      else
         czone = 'LST'
      endif
      if (idsv.eq.1) czone(2:2)='D'
      if (inearesth.eq.0) then
         WRITE (dtgfl,501)
     &      ihh,imin,czone,dayofweek(idow),idd,name(imm),iyy
      else
         WRITE (dtgfl,502)
     &      ihh,czone,dayofweek(idow),idd,name(imm),iyy
      endif
c
      return
      end
