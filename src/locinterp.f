c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine locinterp(string,gridx,gridy,rlat,rlon,iwmo,icaoid,
     &   stelev,locdesc,rip_root)
c
      character string(2)*20,icaoid*4,fnm*256,rip_root*256
     & ,icaoid2*4,locdesc*44
c
c   Stationlist common block
c
      common /sl/ igot_sl,ns_sl,icao_sl,iwmo_sl,slat_sl,slon_sl,loc_sl,
     &            stelev_sl
      dimension iwmo_sl(15000),slat_sl(15000),slon_sl(15000),
     &   stelev_sl(15000)
      character icao_sl(15000)*4,loc_sl(15000)*44
c
      logical numeric
      external numeric
c
      include 'comconst'
c
      if (string(1).eq.'missing             '.and.
     &    string(2).eq.'missing             ') then ! Everyth. mssng.
         write(iup,*)'In locinterp, location specification is'
         write(iup,*)'completely missing.  Stopping.'
         stop
      elseif (string(2).eq.'missing             ') then ! WMO or ICAO
         if (igot_sl.eq.0) then  ! stationlist not yet in memory. Read it.
            iendci=index(rip_root,' ')-1
            fnm=rip_root(1:iendci)//'/stationlist'
            open (unit=iustnlist,file=fnm,form='formatted',status='old')
            read(iustnlist,*)
            read(iustnlist,*)
            read(iustnlist,*)
            i=0
  200       i=i+1
            read(iustnlist,206,end=205) iwmo_sl(i),icao_sl(i),loc_sl(i),
     &         ilatd,ilatm,ilats,ilond,ilonm,ilons,stelev_sl(i)
            slat_sl(i)=sign(abs(float(ilatd))+.016667*float(ilatm)+
     &         .000278*float(ilats),float(ilatd))
            slon_sl(i)=sign(abs(float(ilond))+.016667*float(ilonm)+
     &         .000278*float(ilons),float(ilond))
            goto 200
 205        ns_sl=i-1
            igot_sl=1
            close (iustnlist)
 206        format(i5,1x,a4,1x,a44,i3,1x,i2,1x,i2,1x,i4,1x,i2,1x,i2,
     &             1x,f4.0)
         endif
         if (numeric(string(1)(1:5))) then ! WMO
            read(string(1)(1:7),'(i5)') iwmo
            do i=1,ns_sl
               if (iwmo.eq.iwmo_sl(i)) then
                  call maptform(gridy,gridx,slat_sl(i),slon_sl(i),-1)
                  icaoid=icao_sl(i)
                  rlat = slat_sl(i)
                  rlon = slon_sl(i)
                  stelev = stelev_sl(i)
                  locdesc=loc_sl(i)
                  goto 210
               endif
            enddo
            write(iup,*)'Couldn''t find station with WMO # =',iwmo
            stop 'locinterp'
  210       continue
         else ! ICAO
            read(string(1)(1:4),'(a4)') icaoid
            if ( icaoid(4:4) .eq. ' ' ) then
c  assume that if it's a 3-character id, it's a U.S. station
              icaoid2 = 'K'//icaoid(1:3)
            else
              icaoid2 = icaoid
            endif
            do i=1,ns_sl
               if (icaoid2.eq.icao_sl(i)) then
                  call maptform(gridy,gridx,slat_sl(i),slon_sl(i),-1)
                  iwmo=iwmo_sl(i)
                  rlat = slat_sl(i)
                  rlon = slon_sl(i)
                  stelev = stelev_sl(i)
                  locdesc=loc_sl(i)
                  goto 212
               endif
            enddo
            write(iup,*)'Couldn''t find station with ICAO ID = ',icaoid
            stop 'locinterp'
 212        continue
         endif
      else ! Either x/y or lat/lon specification.
         iwmo=99999
         icaoid='XXXX'
         locdesc=' '
         inlat1=index(string(1),'lat')
         inlon1=index(string(1),'lon')
         inlat2=index(string(2),'lat')
         inlon2=index(string(2),'lon')
         if (inlat1.ne.0.or.inlon1.ne.0) then ! Lat/lon specification.
            if (inlat1.ne.0) then ! Lat is first.
               read(string(1)(1:inlat1-1),fmt=*) rlat
               read(string(2)(1:inlon2-1),fmt=*) rlon
            else ! Lon is first.
               read(string(1)(1:inlon1-1),fmt=*) rlon
               read(string(2)(1:inlat2-1),fmt=*) rlat
            endif
            call maptform(gridy,gridx,rlat,rlon,-1)
         else ! x/y specification
            read(string(1)(1:12),fmt=*) gridxn
            read(string(2)(1:12),fmt=*) gridyn
            gridx=(gridxn-1.)/refrat+xjcorn
            gridy=(gridyn-1.)/refrat+yicorn
            call maptform(gridy,gridx,rlat,rlon,1)
         endif
      endif
c      print*,'gridx,gridy,rlat,rlon=',gridx,gridy,rlat,rlon
      return
      end
