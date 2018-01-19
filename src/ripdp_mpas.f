      program ripdp_mpas
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIPDP_MPAS is a modified version of RIPDP_WRFARW that reads       c
c   data from MPAS (Model for Prediction Across Scales) and generates c
c   new data files in ripdp format for the RIP data analysis and      c
c   visualization program.  RIP stands for Read/Interpolate/Plot and  c
c   RIPDP stands for RIP data preparation.                            c     
c                                                                     c
c   RIPDP_MPAS was written from RIPDP_WRFARW by Priscilla Mooney      c
c   in 2017.  In some places the author and contributors to the code  c
c   are noted with with initial identifiers.                          c
c                                                                     c 
c   Priscilla Mooney: PAM                                             c 
c   Robert Fovell: RGF                                                c  
c                                                                     c
c   The original RIPDP_WRFARW introduction is below.                  c
c                                                                     c
c                                                                     c
c   RIPDP_WRFARW ia a data preparation program that reads in data     c
c   from the WRF ARW modeling system and makes data files             c
c   appropriate for input to the RIP data analysis and visualization  c
c   program.  RIP stands for Read/Interpolate/Plot and RIPDP stands   c
c   for RIP data preparation.                                         c
c                                                                     c
c   FYI: Output unit number is 65.                                    c
c                                                                     c
c   Originally written by Mark Stoelinga, Univ. of Washington         c
c                                                                     c
c   Modified Feb 2003 for new, generalized vertical coordinate        c
c   version of RIP.  See code commented with "ccc".                   c
c                                                                     c
c   Modified March 2003 for WRF model output.                         c
c      [Wei Wang (NCAR) and Mark Stoelinga (UW)]                      c
c                                                                     c
c   Modified November 2006 for WRF geogrid and metgrid data.          c
c      [Cindy Bruyere (NCAR) and Mark Stoelinga (UW)]                 c
c                                                                     c
c   Modified June 2007 (and name changed from ripdp_wrf to            c
c   ripdp_wrfarw) to accomodate the separate version of ripdp for     c
c   the WRF NMM model (ripdp_wrfnmm.f)                                c
c      [Mark Stoelinga (UW)]                                          c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The main program does very little "work".  It simply reads
c   gets some basic information from the metadata part of the
c   data file in order to give information to
c   subroutine process, which does most of the "work" in RIPDP.
c   Some of the arguments to subroutine process are used as
c   dimensions of the 2- and 3-D arrays.  The purpose of this
c   set up is for dynamic memory allocation.  This capability is
c   not standard to Fortran 77, but is allowable on some f77
c   compilers (such as Cray).  It is standard to Fortran 90.
c   Therefore, RIPDP may be compiled with f77 compilers that allow
c   adjustable dimensioning of local (non-argument) arrays in
c   subroutines, or with any f90 compiler.
c
c   Map transform common block
c
      include "commptf"
c
c   netcdf variables
c
      integer ncid, dimid, nf_status
      character nf_att_text*128, start_date*19
c
      character argum(256)*256
c
      include "netcdf.inc"

c
c   Get command line arguments.
c
      nargum=iargc()
      do i=1,nargum
         argum(i)=' '
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:13).eq.'ripdp_mpas_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if ((argum(1)(1:3).eq.'-n '.and.nargum.lt.5).or.
     &    (argum(1)(1:3).ne.'-n '.and.nargum.lt.3)) then
         print*,'Usage:'
         print*,'  ripdp_mpas [-n namelist_file] casename basic|all',
     &          ' data_file_1 data_file_2 data_file_3 ...'
         print*
         print*,'Note: "namelist_file" can be any name, either with'
         print*,'an extension of your choosing or without an extension.'
         print*,'The ".in" extension is not assumed, as it is in rip.'
         print*,'"basic" tells RIPDP to process only the basic'
         print*,'variables it expects to find.  "all" tells RIP to'
         print*,'process all variables it encounters.'
         stop
      endif
c
      if (argum(1)(1:3).eq.'-n ') then
         nnl=2
         ncn=3
         ndd=4
         nsets=nargum-4
         nsetsbeg=5
      else
         nnl=0
         ncn=1
         ndd=2
         nsets=nargum-2
         nsetsbeg=3
      endif
c
c   Open first netcdf file and get model dimensions
c
      nf_status = nf_open (argum(nsetsbeg), nf_nowrite, ncid)
      call handle_err(000.,nf_status)
c
      write(*,'(" Input file: ",A)') trim(argum(nsetsbeg))
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'model_name', nf_att_text)
      call handle_err(001.,nf_status)
      if (index(nf_att_text,'mpas').ne.0) then
         print*,'Data is recognized as MPAS data.'
         print*
         iprog = 8       
      else
         print*,'Data does not seem to be from the MPAS system.'
         print*,'nf_att_text = ',nf_att_text
         stop
      endif
c
      nf_status = nf_inq_dimid (ncid, 'nLat_m', dimid)
      call handle_err(002.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, miy)
      call handle_err(003.,nf_status)
      miy=miy+1
c
      nf_status = nf_inq_dimid (ncid, 'nLon_m', dimid)
      call handle_err(004.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, mjx)
      call handle_err(005.,nf_status)
      mjx=mjx+1
c
      nf_status = nf_inq_dimid (ncid, 'StrLen', dimid)
      call handle_err(002.,nf_status)

      if ( iprog .gt. 1 ) then
        nf_status = nf_inq_dimid (ncid, 'nVertLevels', dimid)
        if (nf_status .ne. nf_noerr) then
          iprog = 1 
        else
           nf_status = nf_inq_dimlen (ncid, dimid, mkzh)
           call handle_err(007.,nf_status)
        endif
      else
        mkzh = 2      ! fake for geogrid files
      endif
c
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_start_time', start_date)
      if (nf_status .eq. nf_noerr) then
         print *, 'config_start_time exists, use it as start_date '
      else 
         nf_status = nf_get_att_text (ncid,nf_global,
     &   'START_DATE', start_date)
         print *, 'Pre-2.1 WRF data, use START_DATE instead'
      endif
      call handle_err(008.,nf_status)
c
      print*,'start_date=',start_date
      read(start_date,'(i4,5(1x,i2))')iyr4,imo,idy,ihr,imn,isc
      if (iyr4.eq.1) then  !idealized (non-real-weather) model output
         iyr4=1940
         imo=1
         idy=1
         ihr=0
         imn=0
         isc=0
      endif
      iyr=mod(iyr4,100)
      mdateb=1000000*iyr+10000*imo+
     &   100*idy+ihr
      rhourb=imn/60.+isc/3600.

      call mconvert(mdateb,mhourb,1,1940)
c
c   Close first netcdf file
c
      nf_status = nf_close (ncid)
      call handle_err(009.,nf_status)
c
c   Run a loop to open all netcdf files, for the sole purpose of
c   getting max number of model output times in any one file.
c
      nwrftimes=0
      do i=1,nsets
         ind=nsetsbeg+i-1
         nf_status = nf_open (argum(ind), nf_nowrite, ncid)
         call handle_err(010.,nf_status)
         nf_status = nf_inq_dimid (ncid, 'Time', dimid)
         call handle_err(011.,nf_status)
         nf_status = nf_inq_dimlen (ncid, dimid, nwrftimes_this_file)
         call handle_err(012.,nf_status)
         nwrftimes=max(nwrftimes,nwrftimes_this_file)
         nf_status = nf_close (ncid)
         call handle_err(013.,nf_status)
      enddo
      if (nwrftimes .le. 0) then
        print *, 'No valid wrf times were found. nwrftimes = ',nwrftimes
        stop
      endif 
c
      call process(miy,mjx,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes,mdateb,mhourb,rhourb,iprog)
      stop
      end
c                                                                     c
c*********************************************************************c
c               r                                                     c
      subroutine process(miy,mjx,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes,mdateb,mhourb,rhourb,iprog)
c
c   This subroutine does most of the "work".
c
c   miy, and mjx are dot-point dimensions, in the x and y directions
c      respectively, of the domain to be analyzed.
c   mkzh is number of model levels in the domain.
c   argum carries the names of the model data files.
c   nsets is the number of files in the model dataset.
c   nsetsbeg is the element of argum that holds the first file name
c      of the model dataset.
c
      character argum(256)*256
c
      include "netcdf.inc"
c
      dimension ter(miy,mjx),xmap(miy,mjx),xlat(miy,mjx),xlon(miy,mjx),
     &   xlat_u(miy,mjx),xlon_u(miy,mjx),
     &   xlat_v(miy,mjx),xlon_v(miy,mjx),
     &   xlat_uv_bgrid(miy,mjx),xlon_uv_bgrid(miy,mjx),
     &   cor(miy,mjx),xlus(miy,mjx),rtc(miy,mjx),rte(miy,mjx),
     &   tgk(miy,mjx),dmap(miy,mjx),tmk(miy,mjx,mkzh),
     &   uuu(miy,mjx,mkzh),ght(miy,mjx,mkzh+1),www(miy,mjx,mkzh+1),
     &   vvv(miy,mjx,mkzh),prs(miy,mjx,mkzh),sfp(miy,mjx),
     &   uue(miy,mjx,mkzh),vve(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),scr3wlev(miy,mjx,mkzh+1),
     &   scr3(miy,mjx,mkzh),scr2(miy,mjx)
      real mu0(miy,mjx)
      character varname*64,fname*256,cxtime*10,
     &   cxtimeavl(256)*10, tmpc*20
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   netcdf variables
c
      parameter (max_wrf_vars=300)
      real nf_uarr3(mkzh,mjx,miy),  nf_varr3(mkzh,mjx,miy),
     &     nf_tarr3(mkzh,mjx-1,miy-1),nf_warr3(mkzh+1,mjx-1,miy-1),
     &     nf_tarr2(mjx-1,miy-1), nf_sarr3(mjx-1,miy-1,24),
     &     nf_tarr2uvb(mjx,miy)
      integer ncid, ndims, nvars, ngatts, unlimdimid, dimid,
     &   nf_status, varid, nf_att_int, iprocvarid(max_wrf_vars)
      integer nf_ustart3(4),nf_vstart3(4),nf_tstart3(4),nf_wstart3(4),
     &   nf_tstart2(3),nf_gstart2(3),nf_tstart2uvb(3)
      integer nf_ucount3(4),nf_vcount3(4),nf_tcount3(4),nf_wcount3(4),
     &   nf_tcount2(3),nf_tcount3s(4), nf_tmp_count(4),nf_gcount2(3),
     &   nf_tcount2uvb(3)
      real nf_att_real
      character nf_att_text*128, nf_varname*64,
     &   wrftimes(nwrftimes)*64
      real nf_znu(mkzh,nwrftimes),nf_znw(mkzh+1,nwrftimes)
      real znu(mkzh),znw(mkzh+1),znfac(mkzh),rdnw(mkzh),rdn(mkzh)
      integer dimid_tm, dimid_we, dimid_sn, dimid_bt, dimid_sls
      integer vardimids(20), nf_att_len
c
c minfo string
c
      character minfostring*128
c
      logical needvapor
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character discard(maxptimes)*64,retain(maxptimes)*64,ptimeunits*1
      namelist/userin/ ptimes,iptimes,ptimeunits,tacc,discard,retain,
     &   iexpandedout,iskpd1
c
      print*, 'Welcome to your friendly RIPDP (V4.6.7) output file !'   !April 2015

c   IF we have pressure level data (metgrid) all 3D fields also contain SFC data
      mkzh_out = mkzh
c
c   Define constants.  Many are taken from Bolton (1980, MWR 108,1046-1053). 
c
      rgas=287.04  !J/K/kg
      rgasmd=.608   ! rgas_moist=rgas*(1.+rgasmd*qvp)
      cp=1004.     ! J/K/kg  Note: not using Bolton's value of 1005.7
      cpmd=.887   ! cp_moist=cp*(1.+cpmd*qvp)
      gamma=rgas/cp
      gammamd=rgasmd-cpmd  ! gamma_moist=gamma*(1.+gammamd*qvp)
      grav=9.81           ! m/s**2
      sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
      eps=0.622
      ezero=6.112  ! hPa
      xlhc0=3.1484e6   ! J/kg
      xlhctd=2370.  !
      xlhf=3.34e5
      rktpmps=1.94
      celkel=273.15
      eslcon1=17.67
      eslcon2=29.65
      esicon1=22.514
      esicon2=6.15e3
      thtecon1=3376. ! K
      thtecon2=2.54
      thtecon3=.81
      tlclc1=2840.
      tlclc2=3.5
      tlclc3=4.805
      tlclc4=55.
      rhoice=917.
      rhowat=1000.
      pi=4.*atan(1.)
      rpd=pi/180.
      abscoef=.145      ! cloud water absorption coefficient in m^2/g
      abscoefi=.272     ! cloud ice absorption coefficient in m^2/g
      ussalr=.0065      ! deg C per m
      rmsg=9.0e+9       ! indicates missing data or specification
c
c   Set a constant limit on small microphysical values
c
      qdelta = 1.e-9
c
c   Define unit numbers
c
      iuinput=7     ! input unit# for namelist, color table,
c                        and plspec table
      ifilecount=1
c
c   Read the namelist values.
c
      iskpd1=0
      do i=1,maxptimes
         ptimes(i)=9e9
         iptimes(i)=99999999
         discard(i)=' '
         retain(i)=' '
      enddo
      ptimeunits='h'
      tacc=1.
      iexpandedout=0
c
      if (nnl.ne.0) then
         open (unit=iuinput,file=argum(nnl),form='formatted',
     &         status='old')
         read (iuinput,userin)
      endif
      iexpanded=0
      if (iexpanded.eq.1) then
         print*,'Input data is on an expanded domain.'
         if (iexpandedout.eq.1) then
           print*,'RIPDP will process the full (expanded) domain.'
         else
           print*,'RIPDP will output the standard (unexpanded) domain.'
         endif
      endif
      tacch=tacc/3600.
c
      ibasic=0
      if (argum(ndd)(1:5).eq.'basic') then
        ibasic=1
      else if (argum(ndd)(1:3).eq.'all') then
        ibasic=0
      else
	print*,' '
        print*,'You must specify "basic" or "all" variables to process:'
        print*,'  ripdp_mpas [-n namelist_file] casename basic|all',
     &          ' data_file_1 data_file_2 data_file_3 ...'
	print*,' '
        print*,'I think you said ',argum(ndd)
	stop
      endif
c
c   Determine number of names in "discard" array.
c
      do i=1,maxptimes
         if (discard(i).eq.'         ') then
            ndiscard=i-1
            goto 492
         endif
      enddo
      ndiscard=maxptimes
 492  continue
      do i=1,maxptimes
         if (retain(i).eq.'         ') then
            nretain=i-1
            goto 493
         endif
      enddo
      nretain=maxptimes
 493  continue
      if (nretain.ne.0.and.ndiscard.ne.0) then
         print*,'You cannot specify both a "discard" list and a'
         print*,'"retain" list.  Specify one or the other.'
         stop
      endif
      if (ibasic.eq.1.and.ndiscard.ne.0) then
         print*,'If you specify "basic" on the command line,'
         print*,'you should only specify a "retain" list (or no list)'
         print*,'in the namelist, but not a "discard" list.'
         stop
      endif
      if (ibasic.eq.0.and.nretain.ne.0) then
         print*,'If you specify "all" on the command line,'
         print*,'you should only specify a "discard" list (or no list)'
         print*,'in the namelist, but not a "retain" list.'
         stop
      endif
c
      iendc=index(argum(ncn),' ')-1
      if (iendc.eq.-1) iendc=256
      iendf1=iendc+12
      nxtavl=0
c
c   If using iptimes, convert the mdates in the iptimes array to
c   xtimes in the ptimes array.  Also, determine nptimes.
c
      if (ptimes(1).lt.0..or.iptimes(1).lt.0.or.(ptimes(1).eq.
     &    9e9.and.iptimes(1).eq.99999999)) then !user wants all times
         nptimes=0
         print*,'Note: RIPDP will process all times encountered.'
      elseif (ptimes(1).ne.9e9.and.iptimes(1).ne.99999999) then
         print*,'Can''t use both ptimes and iptimes.'
         stop
      elseif (iptimes(1).ne.99999999) then
         do i=1,maxptimes
            if (iptimes(i).eq.99999999) then
               nptimes=i-1
               goto 259
            else
               if (iptimes(i).lt.0) then
                  call mconvert(-iptimes(i),mhourp,1,1940)
                  ptimes(i)=-float(mhourp-mhourb)
               elseif (i.ge.2.and.iptimes(i-1).lt.0) then
                  ptimes(i)=float(iptimes(i))
               else
                  call mconvert(iptimes(i),mhourp,1,1940)
                  ptimes(i)=float(mhourp-mhourb)
               endif
            endif
         enddo
         nptimes=maxptimes
 259     continue
      else
         if (ptimeunits.eq.'h') then
            tunitfac=1.
         elseif (ptimeunits.eq.'m') then
            tunitfac=1./60.
         elseif (ptimeunits.eq.'s') then
            tunitfac=1./3600.
         endif
         do i=1,maxptimes
            if (ptimes(i).eq.9e9) then
               nptimes=i-1
               goto 261
            else
               ptimes(i)=tunitfac*ptimes(i)
            endif
         enddo
         nptimes=maxptimes
 261     continue
      endif
c
c   Process time sequences in ptimes array.
c
      ii=0
      itime=0
      ptusemax=-9e9
  400 ii=ii+1
      if (ii.gt.nptimes) goto 420
      if (ptimes(ii).ge.0.) then
         itime=itime+1
         if (itime.gt.maxptimes) then
            print*,'Number of times requested exceeds maxptimes.'
            print*,'Increase maxptimes in ripdp code, recompile,'
            print*,'and run ripdp again.'
            stop
         endif
         ptuse(itime)=ptimes(ii)
         ptusemax=max(ptuse(itime),ptusemax)
      else
         ii=ii+1
         if (ptimes(ii).gt.0.) then
            tstart=ptimes(ii-2)
            tend=-ptimes(ii-1)
            if (tend.eq.tstart) tend=tend+1.e-10
            tinc=ptimes(ii)
            tdist=tend-tstart
            isign=nint(tdist/abs(tdist))
            ntseries=int(abs(tdist)/tinc+.00001) + 1
            do i=2,ntseries
               itime=itime+1
               if (itime.gt.maxptimes) then
                  print*,
     &              'Number of times requested exceeds maxptimes.'
                  print*,
     &              'Increase maxptimes in ripdp code, recompile,'
                  print*,'and run ripdp again.'
                  stop
               endif
               ptuse(itime)=ptuse(itime-1)+isign*tinc
               ptusemax=max(ptuse(itime),ptusemax)
            enddo
         else
            print*,'Error in ptimes sequence specification.'
            stop
         endif
      endif
      goto 400
  420 nptuse=itime
      if (nptuse.eq.0) ptusemax=9e9
c
c=================================================================c
      do iwf=1,nsets         ! File loop
c=================================================================c
c
      ind=nsetsbeg+iwf-1
      nf_status = nf_open (argum(ind), nf_nowrite, ncid)
      call handle_err(014.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'Time', dimid)
      call handle_err(015.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nwrftimes_this_file)
      call handle_err(016.,nf_status)
c
c   Get basic netcdf info
c
      nf_status = nf_inq (ncid, ndims, nvars, ngatts, unlimdimid)
      call handle_err(017.,nf_status)
c
c     Initialize iprocvarid array to 0.  For each varid obtained for a
c     variable that gets processed, set that element of iprocvarid to 1.
c     This will allow ripdp to check all unprocessed variables to see if
c     a data file can be created.
c
      do i=1,max_wrf_vars
         iprocvarid(i)=0
      enddo
c
c   Get some dimension IDs
c
      nf_status = nf_inq_dimid (ncid, 'Time', dimid_tm)
      call handle_err(018.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'nLon_m', dimid_we)
      call handle_err(019.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'nLat_m', dimid_sn)
      call handle_err(020.,nf_status)
      if ( iprog .gt. 1 ) then
        nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid_bt)
        if ( nf_status .ne. 0 ) then
          nf_status = nf_inq_dimid (ncid, 'nVertLevels', dimid_bt)
        endif
        if ( iprog .gt. 2 ) then
          call handle_err(021.,nf_status)
          nf_status = nf_inq_dimid (ncid, 'nSoilLevels', dimid_sls)
          call handle_err(022.,nf_status)
        endif
      endif
c
c   Get array of model output dates/times
c
      nf_status = nf_inq_varid (ncid, 'xtime', varid)
      call handle_err(023.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_text (ncid, varid, wrftimes)
      call handle_err(024.,nf_status)

      if ( iprog .le. 8 ) goto 513    ! don't have this for either geogrid for metgrid data
c
c     Get ptop.
c
      nf_status = nf_inq_varid (ncid, 'P_TOP', varid)
      call handle_err(025.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_real (ncid, varid, nf_znw)
      call handle_err(026.,nf_status)
      ptop=nf_znw(1,1)
c      print*,'ptop=',ptop,' Pa'
c
c     Get the vertical "eta" coordinate values of both the "w" and
c     "mass" levels, and re-order in RIP's top-bottom order.  Note, the
c     value of eta increases from the surface to model top, whereas
c     k-index in the RIP system decreases from the surface to model top.
c
      nf_status = nf_inq_varid (ncid, 'ZNW', varid)
      call handle_err(027.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_real (ncid, varid, nf_znw)
      call handle_err(028.,nf_status)
      do k=1,mkzh+1
         znw(k)=nf_znw(mkzh+1-k+1,1)
      enddo
c
      nf_status = nf_inq_varid (ncid, 'ZNU', varid)
      call handle_err(029.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_real (ncid, varid, nf_znu)
      call handle_err(030.,nf_status)
      do k=1,mkzh
         znu(k)=nf_znu(mkzh-k+1,1)
         znfac(k)=(znw(k)-znu(k))/(znw(k)-znw(k+1))
         dnw=(znw(k) - znw(k+1))
         rdnw(k) = 1./dnw
      enddo
      do k=2,mkzh
         dn=.5*( 1./rdnw(k-1) + 1./rdnw(k))
         rdn(k)=1./dn
      enddo

 513  continue     

c
c   Set values of "count" and "start" arrays for reading netcdf data.
c
      nf_ustart3(1)=1
      nf_ustart3(2)=1
      nf_ustart3(3)=1
      nf_ustart3(4)=1     ! will be reset at each time in time loop
      nf_ucount3(1)=mkzh
      nf_ucount3(2)=mjx
      nf_ucount3(3)=miy
      nf_ucount3(4)=1
      nf_vstart3(1)=1
      nf_vstart3(2)=1
      nf_vstart3(3)=1
      nf_vstart3(4)=1     ! will be reset at each time in time loop
      nf_vcount3(1)=mkzh
      nf_vcount3(2)=mjx
      nf_vcount3(3)=miy
      nf_vcount3(4)=1
      nf_tstart3(1)=1
      nf_tstart3(2)=1
      nf_tstart3(3)=1
      nf_tstart3(4)=1     ! will be reset at each time in time loop
      nf_tcount3(1)=mkzh
      nf_tcount3(2)=mjx-1
      nf_tcount3(3)=miy-1
      nf_tcount3(4)=1
      nf_wstart3(1)=1
      nf_wstart3(2)=1
      nf_wstart3(3)=1
      nf_wstart3(4)=1     ! will be reset at each time in time loop
      nf_wcount3(1)=mkzh+1
      nf_wcount3(2)=mjx-1
      nf_wcount3(3)=miy-1
      nf_wcount3(4)=1
      nf_gstart2(1)=1
      nf_gstart2(2)=1
      nf_gstart2(3)=1
      nf_gcount2(1)=mkzh+1
      nf_gcount2(2)=mjx-1
      nf_gcount2(3)=miy-1
      nf_tstart2(1)=1
      nf_tstart2(2)=1
      nf_tstart2(3)=1     ! will be reset at each time in time loop
      nf_tcount2(1)=mjx-1
      nf_tcount2(2)=miy-1
      nf_tcount2(3)=1
      nf_tstart2uvb(1)=1
      nf_tstart2uvb(2)=1
      nf_tstart2uvb(3)=1     ! will be reset at each time in time loop
      nf_tcount2uvb(1)=mjx
      nf_tcount2uvb(2)=miy
      nf_tcount2uvb(3)=1
c
      if ( iprog .gt. 2  ) then
        nf_status = nf_inq_dimid (ncid, 'nSoilLevels', dimid)
        call handle_err(031.,nf_status)
        nf_status = nf_inq_dimlen (ncid, dimid, nsoil)
        call handle_err(032.,nf_status)
        nf_tcount3s(3)=nsoil
      endif
      nf_tcount3s(1)=mjx-1
      nf_tcount3s(2)=miy-1
      nf_tcount3s(4)=1
c
c   Get information from netcdf global attributes
c
c   First check if "STAND_LON" exists.  It is the longitude (deg.) that is
c   parallel to the y-axis of the grid.  If it does not exist, this is an
c   "old" WRF header (prior to nesting capabilities), and things have to be
c   done differently than with a "new" WRF header (after nesting capabilities
c   were available).
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'STAND_LON', nf_att_real)
      if (nf_status .ne. nf_noerr) then
         iheader=1
      else
         iheader=2
         xlonc=nf_att_real
      endif
c
c   "Coarse domain" dimensions are not defined in WRF output.  However, RIP is
c    going to continue to work in this framework for historical purposes.
c    For WRF output prior to nesting, the coarse domain is the same as the
c    current domain.  For WRF output after nesting is avaialable, a "coarse"
c    domain is a contrived of size 100x100 and grid space 50 km, centered
c    at latitude = TRUELAT1, longitude = STAND_LON (already set above).

      if (iheader.eq.1) then
         nf_status = nf_get_att_int (ncid, nf_global,
     &      'SOUTH-NORTH_GRID_DIMENSION', nf_att_int)
         call handle_err(033.,nf_status)
         miycors=nf_att_int

         nf_status = nf_get_att_int (ncid, nf_global,
     &      'WEST-EAST_GRID_DIMENSION', nf_att_int)
         call handle_err(034.,nf_status)
         mjxcors=nf_att_int
      else
         miycors=100
         mjxcors=100
      endif

c      if (iexpanded.eq.1) then
c         miycors=miy
c         mjxcors=mjx
c         ioffexp=??
c         joffexp=??
c      endif
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DX', nf_att_real)   ! DX is in meters
      call handle_err(035.,nf_status)
      dskm=.001*nf_att_real
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LAT', nf_att_real)
      call handle_err(036.,nf_status)
      if (iheader.eq.1) then
         xlatc=nf_att_real
      elseif (iheader.eq.2) then
         cenlat_thisdom=nf_att_real  ! to be used later to get yicorn,xjcorn
      endif
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LON', nf_att_real)
      call handle_err(037.,nf_status)
      if (iheader.eq.1) then
         xlonc=nf_att_real
      elseif (iheader.eq.2) then
         cenlon_thisdom=nf_att_real  ! to be used later to get yicorn,xjcorn
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'MAP_PROJ', nf_att_int)
      call handle_err(038.,nf_status)
      nproj=nf_att_int
      if (nproj.gt.3.or.nproj.lt.0) then
         print*,'   Map proj. #',nproj,' is not recognized.'
         stop
      endif
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'TRUELAT1', nf_att_real)
      call handle_err(039.,nf_status)
      truelat1=nf_att_real
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'TRUELAT2', nf_att_real)
      call handle_err(040.,nf_status)
      truelat2=nf_att_real
c
      if (iheader.eq.2.and.iwf.eq.1) then
c
c      Only set xlatc to truelat1 if this is the first WRF data set.
c      There is some sort of bug in WRF ARW that changes truelat1 in
c      subsequent model output files if the model output is written
c      to more than one file.
c
         xlatc=truelat1
      endif
c
      if (iheader.eq.1) then
         dskmc=dskm  ! grid spacing (in km) of coarsest domain
      elseif (iheader.eq.2) then
         dskmc=50.  ! km
      endif
      dsc=dskmc*1000.
c
c   Call premaptform...it will be needed to do some grid-to-lat/lon
c   conversions later
c
      call premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,truelat1,truelat2,6)
c
c   Things for this domain
c
      ds=dskm*1000.
      refrat=dsc/ds  ! refinement ratio of this grid w.r.t. coarse grid
c
c     Note, yicorn and xjcorn, the y and x locations of the lower left
c     corner of the current domain with respect to the centered coarse
c     domain (i.e. the coarse domain centered on lat=xlatc,lon=xlonc)
c     are not set here.  They are set at each time (within the time loop),
c     to account for the future possibility of a moving nest.
c
c   Land use data set:
c
      nf_att_text='USGS'
c      nf_status = nf_get_att_text (ncid, nf_global,
c     &   'MMINLU', nf_att_text)
c      call handle_err(041.,nf_status)
      if (nf_att_text(1:5).eq.'USGS ') then
         ilandset=2    ! USGS 24-category land use data set
      elseif (nf_att_text(1:5).eq.'MODIF') then
         ilandset=3    ! MODIS 20-category land use data set
      elseif (nf_att_text(1:5).eq.'') then
         ilandset=1    ! Probably idealized data
         print*,'WARNING: If this is not iealized data, 
     &           then your MMINLU variable in your data is incorrect'
      else
         print*,'Unexpected land use data set specified.'
         print*,'MMINLU=',nf_att_text(1:5)
         STOP
      endif
c
c   Other
c
      iplevdata=iprog   ! Anything larger than 3 means terrain-following data
c
c     Determine iice.  In RIP, iice=1 means that there are separate
c     arrays for frozen hydrometeors (as would be created by a
c     mixed-phase bulk scheme), whereas ice=0 means that frozen and
c     liquid hydrometeors are combined into a single array (e.g. the
c     cloud array contains both cloud water and ice, and the rain array
c     contains both rain and snow, as would be created by a "simple"
c     (non-mixed phase) bulk scheme.  iice will be determined by looking
c     for the presence of a snow mixing ratio array.
c     
      iice=0
      nf_status = nf_inq_varid (ncid, 'qs', varid)
      if (nf_status .eq. nf_noerr) then
         print *, 'Found qs variable, set iice = 1'
         iice=1
      endif
c
c   Write out model info to the Jim Bresch-inspired
c   ".minfo" file (only if this is the first data file)
c
      if (iwf.eq.1) then   ! start of minfo writing
c
      fname=argum(ncn)(1:iendc)//'.minfo'
      open (unit=58,file=fname,form='formatted',
     &      status='unknown')
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'model_name', nf_att_text)
      call handle_err(042.,nf_status)
c
      if ( iprog .le. 2 ) then        ! geogrid or metgrid output
      else                            ! mpas output
        if (index(nf_att_text,' REAL_EM') .ne. 0) then      ! real wrfinput file
        else if (index(nf_att_text,' WRFVAR') .ne. 0) then   ! wrfvar wrfinput file
        else            ! mpas output file
          istart = 1
          iend = len_trim(nf_att_text)
          ! istart=index(nf_att_text,'WRF V')+4
          ! iend=index(nf_att_text(istart:),' ')-1+istart
          if (nf_att_text == "mpas") then
             minfostring='Model Info: MPAS'
          else
             minfostring='Model Info: '//nf_att_text
          endif
          nf_att_text=' '
          nf_status = nf_get_att_text (ncid, nf_global,
     &         'git_version', nf_att_text)
          minfostring = trim(minfostring)//nf_att_text
          ia = len_trim(minfostring)+1
       minfostring(ia:)=' CU:'
       ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_convection_scheme', nf_att_text)
         call handle_err(043.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' No Cu'
      else
         minfostring(ia:)=nf_att_text
      endif
c
      ia = len_trim(minfostring)+1
       minfostring(ia:)=' PBL:'
       ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_pbl_scheme', nf_att_text)
      call handle_err(044.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' None '
      else
         minfostring(ia:)=nf_att_text
      endif
c
      ia = len_trim(minfostring)+1
       minfostring(ia:)=' MP:'
       ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_microp_scheme', nf_att_text)
      call handle_err(045.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' No microph'
      else
         minfostring(ia:)=nf_att_text
      endif
c
      ia = len_trim(minfostring)+1
       minfostring(ia:)=' SF:'
      ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_lsm_scheme', nf_att_text)
      call handle_err(046.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' No SFC'
      else
         minfostring(ia:)=nf_att_text
      endif

      ia = len_trim(minfostring)+1
      if (ds.ge.100000.) then
         write(minfostring(ia:),'(1x,i3,'' km '')') nint(.001*ds)
      elseif (ds.ge.10000.) then
         write(minfostring(ia:),'(1x,i2,'' km '')') nint(.001*ds)
      elseif (ds.ge.1000.) then
         write(minfostring(ia:),'(1x,f3.1,'' km '')') .001*ds
      else
         write(minfostring(ia:),'(1x,i3,'' m '')') nint(ds)
      endif
c
      ia = len_trim(minfostring)+1
      write(minfostring(ia:),'(i3,'' levels '')') mkzh
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'config_dt', nf_att_real)
      call handle_err(047.,nf_status)
      ia = len_trim(minfostring)+1
      write(minfostring(ia:),'(1x,f7.3,'' sec'')') nf_att_real
      write(58,'(a)') minfostring(1:128)
c
c   Second line of minfo file
c
      minfostring = ' '
      ia = 22
       minfostring(ia:)=' LW:'
       ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_radt_lw_scheme', nf_att_text)
      call handle_err(048.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)='LW: none'
      else 
         minfostring(ia:)=nf_att_text
      endif
      ia = len_trim(minfostring)+1
       minfostring(ia:)=' SW:'
      ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_radt_sw_scheme', nf_att_text)
      call handle_err(049.,nf_status)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' SW: none'
      else
         minfostring(ia:)=nf_att_text
      endif
c      ia = 44
c      ib = 56
c      nf_status = nf_get_att_int (ncid, nf_global,
c     &   'DIFF_OPT', nf_att_int)
c      call handle_err(050.,nf_status)
c      if (nf_att_int.eq.0) then
c         minfostring(ia:ib)=' DIFF: none'
c      elseif (nf_att_int.eq.1) then
c         minfostring(ia:ib)=' DIFF: simple'
c      elseif (nf_att_int.eq.2) then
c         minfostring(ia:ib)=' DIFF: full'
c      endif


c      ia = 57
c      ib = 70
c      if ( nf_att_int.eq.0 ) then
c         minfostring(ia:ib)='                '
c      else
c        nf_status = nf_get_att_int (ncid, nf_global,
c     &     'KM_OPT', nf_att_int)
c	call handle_err(051.,nf_status)
c	if (nf_att_int.eq.0) then
c	   minfostring(ia:ib)=' KM: undef'
c	elseif (nf_att_int.eq.1) then
c	   minfostring(ia:ib)=' KM: constant'
c	elseif (nf_att_int.eq.2) then
c	   minfostring(ia:ib)=' KM: 3D TKE'
c	elseif (nf_att_int.eq.3) then
c	   minfostring(ia:ib)=' KM: 3D Smagor'
c	elseif (nf_att_int.eq.4) then
c	   minfostring(ia:ib)=' KM: 2D Smagor'
c	endif
c      endif
c      ia = 71
c      ib = 87
c      nf_status = nf_get_att_int (ncid, nf_global,
c     &   'DAMP_OPT', nf_att_int)
c      if (nf_att_int.eq.0) then
c         minfostring(ia:ib)=' DAMP: none'
c      elseif (nf_att_int.eq.1) then
c         minfostring(ia:ib)=' DAMP: diffusive'
c      elseif (nf_att_int.eq.2) then
c         minfostring(ia:ib)=' DAMP: Rayleigh2'
c      elseif (nf_att_int.eq.3) then
c         minfostring(ia:ib)=' DAMP: Rayleigh3'
c      endif
      ia = len_trim(minfostring)+1
      minfostring(ia:)=' SFLAY:'
      ia = len_trim(minfostring)+1
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'config_sfclayer_scheme', nf_att_text)
      if (nf_att_text.eq.'') then
         minfostring(ia:)=' SFLAY: none'
      else
         minfostring(ia:)=nf_att_text
      endif

      write(58,'(a)') minfostring(1:128)
      endif  ! end of real/wrf/wrfvar block
      endif  ! end of iprog block
 514  continue
      close(58,status='keep')
c
      endif     ! end of minfo writing
c
c   Set up stuff for RIP format data files.
c
      do i=1,32
         ihrip(i)=999999999
         rhrip(i)=9e9
         chrip(i)=' '
         chrip(i+32)=' '
      enddo
      chrip(1)=
     &'map projection (0: None/Ideal, 1: LamCon, 2: PolSter, 3: Merc)'
      ihrip(1)=nproj
      chrip(2)=
     &'number of dot points in the y-direction (coarse domain)'
      ihrip(2)=miycors
      chrip(3)=
     &'number of dot points in the x-direction (coarse domain)'
      ihrip(3)=mjxcors
      chrip(4)=
     &'number of dot points in the y-direction (this domain)'
      ihrip(4)=miy
      chrip(5)=
     &'number of dot points in the x-direction (this domain)'
      ihrip(5)=mjx
      if (iexpanded.eq.1.and.iexpandedout.eq.0) then
         ihrip(2)=miy-2*ioffexp
         ihrip(3)=mjx-2*joffexp
         ihrip(4)=ihrip(2)
         ihrip(5)=ihrip(3)
      endif
      chrip(6)=
     &'number of dimensions of this variable (2 or 3)'
      ihrip(6)=999    ! this is set separately for each variable
      chrip(7)=
     &'grid of this variable (1: cross point dom., 0: dot point dom.)'
      ihrip(7)=999    ! this is set separately for each variable
ccc      chrip(8)=
ccc     &'vertical coordinate (0: hydrostatic sigma, 1: nonhyd. sigma)'
ccc      ihrip(8)=inhyd
ccc      chrip(9)=
ccc     &'number of half sigma levels'
      chrip(9)=
     &'number of vertical levels in the data'
      ihrip(9)=mkzh
      if ( iprog .eq. 2 ) ihrip(9)=mkzh-1
      chrip(10)=
     &'mdateb: YYMMDDHH (truncated hour) of hour-0 for this dataset'
      ihrip(10)=mdateb
      chrip(11)=
     &'mdate: YYMMDDHH (truncated hour) of this time'
      ihrip(11)=99999999    ! this is set separately for each time
      chrip(12)=
     &'ice physics (1: sep. arrays for ice fields, 0: no sep. arrays)'
      ihrip(12)=iice
ccc      chrip(13)=
ccc     &'Program #: 1:TER. 2:DG/RG. 3:RAW. 5:INT. 6:MOD. 11:MOD.(MM5V3)'
      chrip(13)=
     &'ver. coord. type: <or=3: hgt. or prs.; >or=4: terrain-following'
      ihrip(13)=iplevdata
      chrip(14)=
     &'landuse dataset (1: old, 13-cat; 2: USGS, 24-cat; 3: SiB, 16 )'
      ihrip(14)=ilandset
c
      ijmp=32
      chrip(ijmp+1)=
     &'first true latitude (deg.)'
      rhrip(1)=truelat1
      chrip(ijmp+2)=
     &'second true latitude (deg.)'
      rhrip(2)=truelat2
      chrip(ijmp+3)=
     &'central latitude of coarse domain (deg.)'
      rhrip(3)=xlatc
      chrip(ijmp+4)=
     &'central longitude of coarse domain(deg.)'
      rhrip(4)=xlonc
      chrip(ijmp+5)=
     &'grid distance of coarse domain (km)'
      rhrip(5)=dskmc
      chrip(ijmp+6)=
     &'grid distance of this domain (km)'
      rhrip(6)=dskm
      chrip(ijmp+7)=
     &'coarse dom. y-position of lower left corner of this domain'
      rhrip(7)=9e9       ! this is set separately for each time
      chrip(ijmp+8)=
     &'coarse dom. x-position of lower left corner of this domain'
      rhrip(8)=9e9       ! this is set separately for each time
ccc      chrip(ijmp+9)=
ccc     &'pressure level (hPa) of the model top'
ccc      rhrip(9)=ptop
ccc      chrip(ijmp+10)=
ccc     &'reference sea-level pressure (Pa)'
ccc      rhrip(10)=refslp
ccc      chrip(ijmp+11)=
ccc     &'reference sea-level temperature (K)'
ccc      rhrip(11)=refslt
ccc      chrip(ijmp+12)=
ccc     &'reference lapse rate (dT/d(ln(p)), K)'
ccc      rhrip(12)=reflaps
      chrip(ijmp+13)=
     &'rhourb: diff (in h) between exact time and mdate of hour-0'
      rhrip(13)=rhourb
      chrip(ijmp+14)=
     &'rhour: diff (in h) between exact time and mdate of this data'
      rhrip(14)=9e9       ! this is set separately for each time
      chrip(ijmp+15)=
     &'xtime: exact time of this data relative to exact hour-0 (in h)'
      rhrip(15)=9e9       ! this is set separately for each time
ccc      chrip(ijmp+16)=
ccc     &'reference stratospheric constant temperature (K)'
ccc      rhrip(16)=refstratt
c
c   Initialize "max" time level. This feature causes RIPDP to
c      ignore standard data which is out of chronological order.
c
      xtimemax=-1.
c
      iprintonce=0
c
c=================================================================c
      do iwt=1,nwrftimes_this_file     ! Time loop
c=================================================================c
c
c   Note the different time specifications:
c
c   mdateb: This refers to the truncated integer hour of the
c     beginning of the model run for model output, or of the
c     starting time for analysis data sets.  It is an 8-digit
c     integer specified as YYMMDDHH.
c
c   mhourb: This integer variable refers to the same time as
c     mdateb, but instead of the YYMMDDHH format, it is specified
c     as the number of hours since 00 UTC 1 January 1 AD, according
c     to the Gregorian calendar with full leap year specification
c     (leap year every 4 years except on century years not dvisible
c     by 400).
c
c   rhourb: This is a real number specifying the fraction of an
c     hour that the exact start time of the model forecast or
c     analysis data set exceeds the truncated integer hour of the
c     start time (mdateb/mhourb).  Typically, mesoscale models
c     are initialized precisely on the hour, so rhourb is
c     typically 0.00.  However, rhourb could, in
c     principal, be in the range 0.0 < or = rhourb < 1.0.
c
c   xtime: This is a real number referring to this particular
c     data time, and is specified as the exact number of hours
c     since the exact model start time or first analysis time,
c     i.e. since the time mhourb+rhourb.
c
c   mdate: This refers to the truncated integer hour of this
c     particular data time.  It is an 8-digit integer specified as
c     YYMMDDHH.
c
c   mhour: This integer variable refers to the same time as
c     mdate, but instead of the YYMMDDHH format, it is specified as
c     the number of hours since 00 UTC 1 January 1 AD (similar to
c     mhourb).
c
c   rhour: This is a real number specifying the fraction of an
c     hour that the exact time of this particular data time exceeds
c     the truncated integer hour of this particular data time
c     (mdate/mhour). rhour is in the range 0.0 < or = rhour < 1.0.
c
c   Set count variables that control time index in netcdf arrays
c
      nf_ustart3(4)=iwt
      nf_vstart3(4)=iwt
      nf_tstart3(4)=iwt
      nf_wstart3(4)=iwt
      nf_tstart2(3)=iwt
c
c   Create mdate, rhour, xtime.
c
      ! print *, 'iyr4 = ',iyr4
      ! print *, 'imo = ',imo
      ! print *, 'idy = ',idy
      ! print *, 'ihr = ',ihr
      ! print *, 'imn = ',imn
      ! print *, 'isc = ',isc
      ! print *, 'iwt = ',iwt
      ! print *, 'wrftimes = ',wrftimes(iwt)
      read(wrftimes(iwt),'(i4,5(1x,i2))')iyr4,imo,idy,ihr,imn,isc
      ! print *, 'new iyr4 = ',iyr4
      ! print *, 'inew imo = ',imo
      ! print *, 'new idy = ',idy
      ! print *, 'new ihr = ',ihr
      ! print *, 'new imn = ',imn
      ! print *, 'new isc = ',isc
      if (iyr4.eq.1) then  !idealized (non-real-weather) model output
         iyr4=1940
      endif
      iyr=mod(iyr4,100)
      mdate=1000000*iyr+10000*imo+
     &   100*idy+ihr
      ! print *, 'mdate = ',mdate
      rhour=imn/60.+isc/3600.
      ! print *, 'rhour = ',rhour
      call mconvert(mdate,mhour,1,1940)
      xtime=float(mhour-mhourb)+rhour-rhourb
c
c   Get latitude, longitude arrays
c
      nf_status = nf_inq_varid (ncid, 'lat', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLAT_M
        nf_status = nf_inq_varid (ncid, 'LAT_M', varid)
      endif
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(053.,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xlat(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'lon', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLONG_M
        nf_status = nf_inq_varid (ncid, 'LONG_M', varid)
      endif
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(055.,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xlon(i,j)=nf_tarr2(j,i)
      enddo
      enddo

      nf_status = nf_inq_varid (ncid, 'LAT_U', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLAT_U
        nf_status = nf_inq_varid (ncid, 'LAT_U', varid)
      endif
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(053.,nf_status)
      do j=1,mjx
      do i=1,miy-1
         xlat_u(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'LONG_U', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLONG_M
        nf_status = nf_inq_varid (ncid, 'LONG_U', varid)
      endif
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(055.,nf_status)
      do j=1,mjx
      do i=1,miy-1
         xlon_u(i,j)=nf_tarr2(j,i)
      enddo
      enddo

      nf_status = nf_inq_varid (ncid, 'LAT_V', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLAT_V
        nf_status = nf_inq_varid (ncid, 'LAT_V', varid)
      endif
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(053.,nf_status)
      do j=1,mjx-1
      do i=1,miy
         xlat_v(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'LONG_V', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLONG_V
        nf_status = nf_inq_varid (ncid, 'LONG_V', varid)
      endif
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(055.,nf_status)
      do j=1,mjx-1
      do i=1,miy
         xlon_v(i,j)=nf_tarr2(j,i)
      enddo
      enddo


      nf_status = nf_inq_varid (ncid, 'LAT_UV_BGRID', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLAT_V
        nf_status = nf_inq_varid (ncid, 'LAT_UV_BGRID', varid)
      endif
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2uvb,
     &   nf_tcount2uvb, nf_tarr2uvb)
      call handle_err(053.,nf_status)
      do j=1,mjx
      do i=1,miy
         xlat_uv_bgrid(i,j)=nf_tarr2uvb(j,i)
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'LONG_UV_BGRID', varid)
      if ( nf_status .ne. 0 ) then   ! maybe we have XLONG_V
        nf_status = nf_inq_varid (ncid, 'LONG_UV_BGRID', varid)
      endif
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2uvb,
     &   nf_tcount2uvb, nf_tarr2uvb)
      call handle_err(055.,nf_status)
      do j=1,mjx
      do i=1,miy
         xlon_uv_bgrid(i,j)=nf_tarr2uvb(j,i)
      enddo
      enddo
c
c     Determine yicorn, xjcorn, the y and x locations of the lower left
c     corner of the current domain with respect to the centered coarse
c     domain (i.e. the coarse domain centered on lat=xlatc,lon=xlonc).
c
      if (iheader.eq.1.and.nproj.gt.0) then

c     Unfortunately, for WRF output prior to nesting capabilities, there
c     is no information in the netcdf attributes to indicate yicorn and
c     xjcorn. It can only be ascertained from the lat/lon data arrays.
c     First, use given lat/lon at cross grid point (y=1,x=1) to
c     determine where that point is in the centered coarse domain.
c
         call maptform(yicorncross,xjcorncross,xlat(1,1),xlon(1,1),-1)
c
c   Set yicorn, xjcorn
c
         yicorn=yicorncross-.5/refrat
         xjcorn=xjcorncross-.5/refrat
c
      elseif (iheader.eq.1.and.nproj.eq.0) then

c     This is idealized output, so no map background, but just to make
c     things reasonable, set up yicorn and xjcorn so that this domain is
c     exactly in the center of the "coarse" (reference) domain.
c
         yicorn=1.+.5*(float(miycors)-1.)-.5*(float(miy)-1.)/refrat
         xjcorn=1.+.5*(float(mjxcors)-1.)-.5*(float(mjx)-1.)/refrat
c
      elseif (iheader.eq.2) then
c
c     The global attributes gave the lat/lon for the center of the
c     "dot-point" domain.  This can be used to get the "coarse domain"
c     location of the center of the "dot-point" domain, and from this
c     can be determined the "coarse domain" location of the lower-left
c     corner of the dot-point domain.
c
         call maptform(yitemp,xjtemp,cenlat_thisdom,cenlon_thisdom,-1)
         yioffset=.5*(float(miy)-1.)/refrat
         xjoffset=.5*(float(mjx)-1.)/refrat
         yicorn=yitemp-yioffset
         xjcorn=xjtemp-xjoffset
c
      endif
c
c  Set RIP header values that are time-dependent
c
      ihrip(11)=mdate
      rhrip(7)=yicorn
      rhrip(8)=xjcorn
      rhrip(14)=rhour
      rhrip(15)=xtime
c
      secondspast=rhour*3600.
      print*
      print*,'****  Reading model output at'
      print*,'      forecast time=',xtime
      if (secondspast.lt..2) then
         write(6,926) mdate
      else
         write(6,927) mdate,secondspast
      endif
 926  format('        (YYMMDDHH = ',i8.8,')')
 927  format('        (YYMMDDHH = ',i8.8,' plus ',f12.5,' seconds)')
c
c   See if the time we just encountered is beyond the latest time
c   requested.  If so, then jump out of time loop here.
c
      if (xtime.gt.ptusemax+tacch) then
         print*,'   But this time is beyond the latest requested time.'
         print*,'   RIPDP is now stopping.'
         goto 1000
      endif
c
c   See if this is a time that is out of chron. order
c
      iskipit=0
      if (xtime.le.xtimemax) iskipit=1  ! this means out of ch. order
      xtimemax=max(xtimemax,xtime)
c
c   See if this is a time that the user doesn't want.
c   Note: if no ptimes or iptimes were specified (nptuse=0), then
c   it is assumed the user wants all encountered times processed
c   (unless they are out of chronological order).
c
      if (iskipit.eq.0.and.nptuse.gt.0) then
         do i=1,nptuse
            if (abs(xtime-ptuse(i)).le.tacch) goto 40
         enddo
         iskipit=2
 40      continue
      endif
c      
      if (iskipit.gt.0) then
c
         if (iskipit.eq.1) then
            print*,'   But this time is chronologically backward,',
     &        ' so we will skip it.'
         else
            print*,'   But you do not want to process this time,',
     &        ' so we will skip it.'
         endif
         goto 240
c
      endif
c
c   Create as much of the data file name as we know at this point.
c
      cxtime=' '
      write(cxtime,'(f10.5)')xtime
      if (cxtime(1:1).eq.' ') cxtime(1:1)='0'
      if (cxtime(2:2).eq.' ') cxtime(2:2)='0'
      if (cxtime(3:3).eq.' ') cxtime(3:3)='0'
      if (cxtime(4:4).eq.' ') cxtime(4:4)='0'
      fname=argum(ncn)(1:iendc)//'_'//cxtime//'_'
      nxtavl=nxtavl+1
      cxtimeavl(nxtavl)=cxtime
c
c   Obtain variables specifically sought by RIPDP, and write them out.
c
      print*,'Processing basic variables.'
c
c   Get U, transfer to "dot-point" grid, and write out.
c
      if ( iprog .gt. 1 ) then
        nf_status = nf_inq_varid (ncid, 'uReconstructZonal', varid)
        call handle_err(056.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_ustart3,
     &     nf_ucount3, nf_uarr3)
        call handle_err(057.,nf_status)
        do k=1,mkzh
           do j=1,mjx
           do i=1,miy
              uue(i,j,k)=nf_uarr3(mkzh-k+1,j,i)
           enddo
           enddo
        enddo
      else  
        call fillarray(uuu,miy*mjx*mkzh,0.)    ! for geogrid files we still need a placeholder
      endif
      plchun='m s~S~-1~N~'
      if ( iprog .eq. 2 ) then
        scr2 = uuu(:,:,mkzh)
        vardesc='Hor. wind (sfc,x-comp.), m/s'
        call writefile_rdp(scr2,'U10       ',2,0,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get V, transfer to "dot-point" grid, and write out.
c
      if ( iprog .gt. 1 ) then
        nf_status = nf_inq_varid (ncid, 'uReconstructMeridional', varid)
        call handle_err(058.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_vstart3,
     &     nf_vcount3, nf_varr3)
        call handle_err(059.,nf_status)
        do k=1,mkzh
           do j=1,mjx
           do i=1,miy
              vve(i,j,k)=nf_varr3(mkzh-k+1,j,i)
           enddo
           enddo
c           call vtodot(vve(1,1,k),miy,mjx)
        enddo
      else
        call fillarray(vvv,miy*mjx*mkzh,0.)    ! for geogrid files we still need a placeholder
      endif
      plchun='m s~S~-1~N~'
      if ( iprog .eq. 2 ) then
        scr2 = vvv(:,:,mkzh)
        vardesc='Hor. wind (sfc,y-comp.), m/s'
        call writefile_rdp(scr2,'V10       ',2,0,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif

c
c     Need to rotate MPAS winds from earth relative coordinates (uue,vve) 
c     to computational domain coordinates (uuu,vvv) - this is what rip functions expect
c
        do k=1,mkzh
           do j=1,mjx
           do i=1,miy
              call rotate_earth2comp(uue(i,j,k),vve(i,j,k),
     &            xlonc,xlon_uv_bgrid(i,j),xlon_uv_bgrid(i,j),
     &            nproj,uuu(i,j,k),vvv(i,j,k))
           enddo
           enddo
        enddo

      vardesc='Horizontal wind (x-comp.), m/s'
      call writefile_rdp(uuu,'uuu       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)

      vardesc='Horizontal wind (y-comp.), m/s'
      call writefile_rdp(vvv,'vvv       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Get pressure 
c   change to hPa (from Pa), and write out.
c
      print*,' ' 
      nf_status = nf_inq_varid (ncid, 'pressure', varid)
c      nf_status = 333
      if (nf_status .ne. nf_noerr .and. iprog .gt. 2) then
         igotpressure=0
         if (iprintonce.eq.0) then
            print*,'Pressure variable (P) not found.'
            print*,'  This indicates the file is probably surface file '
            iprintonce=1
         endif
      else
c      elseif (iprog .eq. 1)
         igotpressure=1
           nf_status = nf_inq_varid (ncid, 'pressure', varid)
           call handle_err(102.,nf_status)
           iprocvarid(varid)=1
           nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &        nf_tcount3, nf_tarr3)
           call handle_err(103.,nf_status)
           do k=1,mkzh
           do j=1,mjx-1
           do i=1,miy-1
              prs(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*0.01
           enddo
           enddo
           enddo
         vardesc='Pressure, hPa'
         plchun='hPa'
         call writefile_rdp(prs,'prs       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get water vapor mixing ratio (qv), convert it to g/kg,
c   and write out.  Otherwise, fill qvp array with 0s but don't write out.
c
      needvapor = .false.
      nf_status = nf_inq_varid (ncid, 'qv', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qv.'
         call fillarray(qvp,miy*mjx*mkzh,0.)
	 if (iprog .eq. 2) then
           nf_status = nf_inq_varid (ncid, 'RH', varid)
           if (nf_status .ne. nf_noerr) then
             print*,'   Did not find RH either.'
	   else
	     iprocvarid(varid)=1
             nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &         nf_tcount3, nf_tarr3)
             call handle_err(113.,nf_status)
	     do k=1,mkzh
	     do j=1,mjx-1     ! for metgrid, store rh in qvp array, convert it later
	     do i=1,miy-1
		qvp(i,j,k)=nf_tarr3(mkzh-k+1,j,i)  ! %
	     enddo
	     enddo
	     enddo
	     needvapor = .true.
	     vardesc='Relative humidity (w.r.t. water), %'
	     plchun='%'
             call writefile_rdp(qvp,'rhu       ',3,1,vardesc,plchun,
     &        fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &        iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)

	   endif
	 endif
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(063.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            qvp(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Water vapor mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qvp,'qvp       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c     Geopotential perturbation and basic state, as well as vertical
c     velocity, are defined at vertical staggered levels, i.e. the
c     "w" levels, similar to so-called "full sigma" levels in the MM5
c     coordinate system.  For the geopotential, since the WRF vertical
c     coordinate is mass based, it would be best to convert the
c     geopotential to an exponential height variable, then linearly
c     interpolate in WRF vert. coord. ("eta") from "w levels" to "mass
c     levels", and then convert back to geopotential height.
c
c   First read geopotential base state (PHB) into the netcdf scratch array
c   for "w-level" variables, and transfer to RIP scratch array for
c   "w-level" variables.
c
      nf_status = nf_inq_varid (ncid, 'zgrid', varid)
      call handle_err(064.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_gstart2,
     &   nf_gcount2, nf_warr3)
      call handle_err(065.,nf_status)
      do k=1,mkzh+1
      do j=1,mjx-1
      do i=1,miy-1
         scr3wlev(i,j,k)=nf_warr3(mkzh+1-k+1,j,i)
      enddo
      enddo
      enddo
c
c   Now interpolate ght to "mass levels" and write out.  
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         ght(i,j,k)=0.5*(scr3wlev(i,j,k)+scr3wlev(i,j,k+1))
      enddo
      enddo
      enddo
c
c   Get terrain height (ter), and write out.
c   Since MPAS doesn't output terrain use lowest ght level
c   PAM
c       scr2 = ght(:,:,mkzh)
        scr2 = scr3wlev(:,:,mkzh+1)
c   PAM
        vardesc = " "
        plchun = " "
        vardesc = "Terrain height AMSL, m"
        plchun = "m"
        call writefile_rdp(scr2,'ter       ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
 515  continue
      vardesc='Geopotential height, m'
      plchun='m'
      call writefile_rdp(ght,'ght       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Now get vertical velocity (convert from m/s to cm/s)
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 'w', varid)
        call handle_err(068.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &     nf_wcount3, nf_warr3)
        call handle_err(069.,nf_status)
        do k=1,mkzh+1
        do j=1,mjx-1
        do i=1,miy-1
           scr3wlev(i,j,k)=100.*nf_warr3(mkzh+1-k+1,j,i)
        enddo
        enddo
        enddo
c
c   Interpolate vert. vel. to "mass levels" and write out
c
        do j=1,mjx-1
        do i=1,miy-1
        do k=1,mkzh
           www(i,j,k)=0.5*(scr3wlev(i,j,k) + scr3wlev(i,j,k+1))
        enddo
        enddo
        enddo
        vardesc='Vertical velocity, cm/s'
        plchun='cm s~S~-1~N~'
        call writefile_rdp(www,'www       ',3,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If pressure was not previously processed (because this is a WRF
c   input file, which doesn't contain pressure arrays), create it
c   now.
c
      if (igotpressure.eq.0) then
c
c      Need to calculate pressure because this is a WRF input file,
c      which doesn't have base-state (PB) and perturbation (P) arrays.
c      First get MU0.  MU0 = MUB + MU
c
         nf_status = nf_inq_varid (ncid, 'MUB', varid)
	 if (nf_status .ne. nf_noerr) then
c  Old SI files only have MU0
           nf_status = nf_inq_varid (ncid, 'MU0', varid)
           call handle_err(072.,nf_status)
	   iprocvarid(varid)=1
	   nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &        nf_tcount2, nf_tarr2)
	   call handle_err(073.,nf_status)
	   do j=1,mjx-1
	   do i=1,miy-1
	      mu0(i,j)=nf_tarr2(j,i)
	   enddo
	   enddo
	 else
c  we've got MUB, so compute MU0 
c  (in wrf3dvar_input files, mu0 is .ne. mub+mu, so recompute it here).
           iprocvarid(varid)=1
           nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &        nf_tcount2, nf_tarr2)
           call handle_err(173.,nf_status)
           do j=1,mjx-1
           do i=1,miy-1
              mu0(i,j)=nf_tarr2(j,i)
           enddo
           enddo
           nf_status = nf_inq_varid (ncid, 'MU', varid)
           call handle_err(272.,nf_status)
           iprocvarid(varid)=1
           nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &        nf_tcount2, nf_tarr2)
           call handle_err(273.,nf_status)
           do j=1,mjx-1
           do i=1,miy-1
              mu0(i,j)=mu0(i,j) + nf_tarr2(j,i)
           enddo
           enddo
	  endif
c
c      Need base state sea-level pressure, sea-level temperature, and
c      temperature difference from slp to 300 hPa.  In an ideal world,
c      these would be available in the WRF input file, but currently
c      the world is less than ideal so we must guess at them.
c
         nf_status = nf_inq_varid (ncid, 'P00', varid)
         if ( nf_status == 0 ) then
	   nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &        nf_tstart2, nf_tarr2)
           p00 = nf_tarr2(1,1)
           print*,"   Found and will use P00 from file"
         else
           p00 = 100000.   ! Pa
           print*,"   Did not find P00 in file - will use 1000mb"
         endif
c
         nf_status = nf_inq_varid (ncid, 'T00', varid)
         if ( nf_status == 0 ) then
	   nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &        nf_tstart2, nf_tarr2)
           t00 = nf_tarr2(1,1)
           print*,"   Found and will use T00 from file"
         else
           t00 = 290.      ! K
           print*,"   Did not find T00 in file - will use 290K"
         endif
           print*,"   "
         a = 50.         ! K
c
c         rgas_wrf=287.
         rgas_wrf=rgas
c
         do j=1,mjx-1
         do i=1,miy-1
c
c         Calculate base state surface pressure (p_surf)
c
            p_surf = p00*exp(-t00/a+((t00/a)**2-2.*grav*
     &         ter(i,j)/a/rgas_wrf)**0.5)     
c
c         Calculate base state column mass (mub)
c
            mub = p_surf - ptop
c
c         Calculate dry-air surface pressure (pd_surf)
c
            pd_surf = mu0(i,j) + ptop
c
c         Calculate perturbation dry-air column mass (mu_2)
c
            mu_2 = pd_surf-p_surf
c
c         Get pressure perturbation at model top
c
            k = 1  ! (uppermost half-eta level)
            qvf1 = 0.5*(qvp(i,j,k)+qvp(i,j,k))*.001
            qvf2 = 1./(1.+qvf1)
            qvf1 = qvf1*qvf2
            prs(i,j,k) = - 0.5*(mu_2+qvf1*mub)/rdnw(k)/qvf2
c
c         Now get pressure perturbation at levels below
c
            do k=2,mkzh
               qvf1 = 0.5*(qvp(i,j,k)+qvp(i,j,k-1))*.001
               qvf2 = 1./(1.+qvf1)
               qvf1 = qvf1*qvf2
               prs(i,j,k) = prs(i,j,k-1)-(mu_2+qvf1*mub)/qvf2/rdn(k)
            enddo
c
c         Finally compute base state pressure and add to pressure perturbation
c         to get total pressure
c
            do k = 1, mkzh
               pb = znu(k)*(p_surf - ptop) + ptop
               prs(i,j,k) = .01*(prs(i,j,k)+pb)   ! .01 converts from Pa to hPa
            enddo
         enddo
         enddo
c
         vardesc='Pressure, hPa'
         plchun='hPa'
         call writefile_rdp(prs,'prs       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get theta (T), convert it to temperature, and write out.
c
        nf_status = nf_inq_varid (ncid, 'theta', varid)
        call handle_err(074.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &     nf_tcount3, nf_tarr3)
        call handle_err(075.,nf_status)
        do k=1,mkzh
        do j=1,mjx-1
        do i=1,miy-1
           tmk(i,j,k)=(nf_tarr3(mkzh-k+1,j,i))*
     &        (prs(i,j,k)/1000.)**gamma
        enddo
        enddo
        enddo
      vardesc='Temperature, K'
      plchun='K'
      call writefile_rdp(tmk,'tmk       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Calculate surface pressure using altimeter equation.
c
      if ( iprog .gt. 2 ) then
        do j=1,mjx-1
        do i=1,miy-1
           nf_status = nf_inq_varid (ncid, 'surface_pressure', varid)  !! First check to see if we have sfc pressure
           if ( nf_status .ne. 0 ) then   ! Don't have PSFC, create one from lowest soiltemp and PMSL
              tv=virtual(tmk(i,j,mkzh),.001*qvp(i,j,mkzh))
              sfp(i,j)=prs(i,j,mkzh)*(tv/(tv+ussalr*
     &        (ght(i,j,mkzh)-ter(i,j))))**(-grav/(rgas*ussalr))
           else
             nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &       nf_tcount2, nf_tarr2)
             call handle_err(112.,nf_status)
             sfp(i,j) = nf_tarr2(j,i)*0.01
           endif
        enddo
        enddo
      endif
      vardesc='Surface pressure, hPa'
      plchun='hPa'
      call writefile_rdp(sfp,'sfp       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)

c   Get mapscale factor on cross points (MAPFAC_M), and write out.
c
      nf_status = nf_inq_varid (ncid, 'MAPFAC_M', varid)
      if ( nf_status .ne. 0 ) then   ! We probably have V3 data - lets look
            nf_status = nf_inq_varid (ncid, 'MAPFAC_MX', varid)
            call handle_err(076.,nf_status)
          endif
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(077.,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xmap(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Map factor on cross points'
      plchun='none'
      call writefile_rdp(xmap,'xmap      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)

c   Transfer xmap to dot-point grid (dmap) and write out
c
      do j=1,mjx-1
      do i=1,miy-1
         dmap(i,j)=xmap(i,j)
      enddo
      enddo
      call xtodot(dmap,miy,mjx)
      vardesc='Map factor on dot points'
      plchun='none'
      call writefile_rdp(dmap,'dmap      ',2,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Get coriolis parameter on cross points (F), and write out.
c
      nf_status = nf_inq_varid (ncid, 'F', varid)
      call handle_err(078.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(079.,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         cor(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Coriolis parameter, per s'
      plchun='s~S~-1~N~'
      call writefile_rdp(cor,'cor       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Already have latitude (XLAT) and longitude (XLONG).  Write them out.
c
      vardesc='Latitude, degrees'
      plchun='degrees'
      call writefile_rdp(xlat,'xlat      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      vardesc='Longitude, degrees'
      plchun='degrees'
      call writefile_rdp(xlon,'xlon      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
c
c   Get total accumulated cumulus preciitation (rainc), and write out.
c
      nf_status = nf_inq_varid (ncid, 'rainc', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find rainc.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(080.,nf_status)
         do j=1,mjx-1
         do i=1,miy-1
            rtc(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         vardesc='Cumulus precip. since h 0, mm'
         plchun='mm'
         call writefile_rdp(rtc,'rtc       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get total accumulated explicit (grid-resolved) preciitation (rainnc),
c      and write out.
c
      nf_status = nf_inq_varid (ncid, 'rainnc', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find rainnc.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(081.,nf_status)
         do j=1,mjx-1
         do i=1,miy-1
            rte(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         vardesc='Explicit precip. since h 0, mm'
         plchun='mm'
         call writefile_rdp(rte,'rte       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get 2m temperature and write out as T2 
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 't2m', varid)
        call handle_err(082.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &     nf_tcount2, nf_tarr2)
        call handle_err(083.,nf_status)
        do j=1,mjx-1
        do i=1,miy-1
           tgk(i,j)=nf_tarr2(j,i)
        enddo
        enddo
        vardesc='Temperature at 2m, K'
        plchun='K'
        call writefile_rdp(tgk,'T2        ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get 2m humidity and write out as Q2 
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 'q2', varid)
        call handle_err(082.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &     nf_tcount2, nf_tarr2)
        call handle_err(083.,nf_status)
        do j=1,mjx-1
        do i=1,miy-1
           tgk(i,j)=nf_tarr2(j,i)
        enddo
        enddo
        vardesc='Humidity at 2m'
        plchun='kg/kg'
        call writefile_rdp(tgk,'Q2        ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get 10m zonal velocity and write out as U10 
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 'u10', varid)
        call handle_err(082.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &     nf_tcount2, nf_tarr2)
        call handle_err(083.,nf_status)
        do j=1,mjx-1
        do i=1,miy-1
           tgk(i,j)=nf_tarr2(j,i)
        enddo
        enddo
        vardesc='10m zonal wind'
        plchun='m/s'
        call writefile_rdp(tgk,'U10       ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get 10m meridional velocity and write out as V10 
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 'v10', varid)
        call handle_err(082.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &     nf_tcount2, nf_tarr2)
        call handle_err(083.,nf_status)
        do j=1,mjx-1
        do i=1,miy-1
           tgk(i,j)=nf_tarr2(j,i)
        enddo
        enddo
        vardesc='10m meridional wind'
        plchun='m/s'
        call writefile_rdp(tgk,'V10       ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get ground temperature (and sea-surface temperature over water) (TSK),
c      and write out.
c
      if ( iprog .gt. 2 ) then
        nf_status = nf_inq_varid (ncid, 'skintemp', varid)
        call handle_err(082.,nf_status)
        iprocvarid(varid)=1
        nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &     nf_tcount2, nf_tarr2)
        call handle_err(083.,nf_status)
        do j=1,mjx-1
        do i=1,miy-1
           tgk(i,j)=nf_tarr2(j,i)
        enddo
        enddo
        vardesc='Ground/sea-surface temperature, K'
        plchun='K'
        call writefile_rdp(tgk,'tgk       ',2,1,vardesc,plchun,
     &     fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &     iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   Get land use (LU_INDEX), and write out.
c
      varid = -999999
      nf_status = nf_inq_varid (ncid, 'LU_INDEX', varid)
c      call handle_err(084.,nf_status)
      if (nf_status .ne. nf_noerr) then
         varid = -999998
         print*,'   Did not find LU_INDEX.... try ivgtyp'
         nf_status = nf_inq_varid (ncid, 'ivgtyp', varid)
         if (nf_status .ne. nf_noerr) then
            varid = -999997
            print*,'   Did not find ivgtyp'
         else
            print*, '   Found ivgtyp'
         endif
      endif
         
      if ( varid > -999997 ) then
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
         call handle_err(085.,nf_status)
         do j=1,mjx-1
         do i=1,miy-1
            xlus(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         vardesc='Land use category'
         plchun='none'
      call writefile_rdp(xlus,'xlus      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
      print*,'Checking for hydrometeor variables.'
c
c   If available, get cloud water mixing ratio (qc), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'qc', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qc.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(scr3,'qcw       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get cloud water number conc. (QNCLOUD), do not convert,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QNCLOUD', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QNCLOUD.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            !scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)! orig
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i) ! RGF fix
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud water number concentration, /kg'
         plchun='kg~S~-1~N~'
         call writefile_rdp(scr3,'qnc       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get rain water mixing ratio (qr), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'qr', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qr.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(087.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Rain water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(scr3,'qra       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get rain water number conc. (QNRAIN), do not convert,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QNRAIN', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QNRAIN.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(087.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           !scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1) ! orig
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Rain water number concentration, /kg'
         plchun='kg~S~-1~N~'
         call writefile_rdp(scr3,'qnr       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get cloud ice mixing ratio (qi), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'qi', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qi.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(088.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud ice mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(scr3,'qci       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get cloud ice mixing ratio (QNICE), do not convert,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QNICE', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QNICE.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(088.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud ice number concentration, /kg'
         plchun='kg~S~-1~N~'
         call writefile_rdp(scr3,'qni       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get snow mixing ratio (qs), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'qs', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qs.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(089.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Snow mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(scr3,'qsn       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get graupel mixing ratio (qg), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'qg', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find qg.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(090.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Graupel mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(scr3,'qgr       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get water-friendly aerosol number conc. (QNWFA),
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QNWFA', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QNWFA.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Water-friendly aerosol number concentration, /kg'
         plchun='kg~S~-1~N~'
         call writefile_rdp(scr3,'qnwfa     ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c   If available, get ice-friendly aerosol number conc. (QNIFA),
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QNIFA', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QNIFA.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i)
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Ice-friendly aerosol number concentration, /kg'
         plchun='kg~S~-1~N~'
         call writefile_rdp(scr3,'qnifa     ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      endif
c
c     Loop through other variables not specifically sought by RIPDP.
c
      if (ibasic.eq.1.and.nretain.eq.0) goto 230
c
      print*,'Processing other variables in output file.'
c
      do ivar=1,nvars
c
      nf_varname=' '
      nf_status = nf_inq_varname (ncid, ivar, nf_varname)
      call handle_err(091.,nf_status)
c
c     First, jump to end of loop if variable was already processed in some
c     way shape or form.
c
      if (iprocvarid(ivar).eq.1) then
         goto 229
      endif
c
c     Next, jump to end of loop if "all" was specified on the command line
c     but the variable name is in the discard list.
c
      if (ibasic.eq.0) then
         do idis=1,ndiscard
            if (nf_varname.eq.discard(idis)) then
               goto 229
            endif
         enddo
      endif
c
c     Next, jump to end of loop if "basic" was specified on the command line
c     but the variable name is NOT in the retain list.
c
      if (ibasic.eq.1) then
         icount=0
         do iret=1,nretain
            if (nf_varname.eq.retain(iret)) icount=icount+1
         enddo
         if (icount.eq.0) goto 229
      endif
c
c     Next, check for certain combinations of variable dimensions
c     that RIP can make use of (i.e., 3D cross point array, 2D cross
c     point array, or soil-layer cross point array).
c
      nf_status = nf_inq_varndims (ncid, ivar, ndims)
      call handle_err(092.,nf_status)
      nf_status = nf_inq_vardimid (ncid, ivar, vardimids)
      call handle_err(093.,nf_status)
c
      ! RGF
      !if(ndims.ge.3)then
       !print *,' ndims= ',ndims,' dimid_tm= ',dimid_tm,' dimid_bt= ',dimid_bt
       !print *,' dimid_sn= ',dimid_sn,' dimid_we= ',dimid_we 
       !print *,' vardimids(4) ',vardimids(4),' dimid_tm ',dimid_tm
       !print *,' vardimids(3) ',vardimids(3),' dimid_bt ',dimid_bt
       !print *,' vardimids(2) ',vardimids(2),' dimid_sn ',dimid_sn
       !print *,' vardimids(1) ',vardimids(1),' dimid_we ',dimid_we
      !endif
      ! RGF
      if (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(3).eq.dimid_bt.and.vardimids(2).eq.dimid_sn.and.
     &    vardimids(1).eq.dimid_we) then  ! original version (does not work)
         itype=1  ! 3D cross point array
! RGF: Alternative code added.  Array dimension ordering above is incorrect.
! Comparisons altered.
      elseif (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(1).eq.dimid_bt.and.vardimids(3).eq.dimid_sn.and.
     &    vardimids(2).eq.dimid_we) then ! RGF
         itype=1  ! 3D cross point array
      elseif (ndims.eq.3.and.vardimids(3).eq.dimid_tm.and.
     &    vardimids(2).eq.dimid_sn.and.vardimids(1).eq.dimid_we) then
         itype=2  ! 2D cross point array
      elseif (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(3).eq.dimid_sls.and.vardimids(2).eq.dimid_sn.and.
     &    vardimids(1).eq.dimid_we) then
         itype=3  ! soil-layer cross point array
      else
         itype=0  ! none of the above
         if (iprog .eq. 1 .and. nf_varname == 'albedo12m') itype = 4 ! Monthly surface albedo  
         if (iprog .eq. 1 .and. nf_varname == 'greenfrac') itype = 4 ! Monthly green fraction  
         if (iprog .eq. 1 .and. nf_varname == 'isltype') itype = 4 ! Top layer soil type     
         if (iprog .eq. 1 .and. nf_varname == 'soilcat_bot') itype = 4 ! Bottom layer soil type     
      endif
c
c     If variable does not have a combination of dimensions
c     that RIP can make use of, skip it.
c
      if (itype.eq.0) goto 229
c
      varname=nf_varname
      inname=0
      do ic=10,1,-1
         if (inname.eq.0) then
            if (varname(ic:ic).ne.' ') inname=1
         else
            if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
         endif
      enddo
      print *,'------------ IMPORTANT NOTE -------------'
      print *,'MPAS does not have variable information (e.g. units)'
      nf_att_text=' '
c      nf_status = nf_get_att_text (ncid, ivar,
c     &   'units', nf_att_text)
c      call handle_err(094.,nf_status)
      plchun=nf_att_text
c      nf_att_text=' '
c      nf_status = nf_get_att_text (ncid, ivar,
c     &   'description', nf_att_text)
c      call handle_err(095.,nf_status)
c      nf_status = nf_inq_attlen (ncid, ivar,
c     &   'description', nf_att_len)
c      call handle_err(096.,nf_status)
c      if (nf_att_len .GT. LEN(nf_att_text)) then
c         print*, 'NOT POSSIBLE TO CONTINUE'
c         print*, ' MEMORY OVERWITE WILL RESULT BECAUSE'
c         print*, ' nf_att_len is greater than declared'
c         print*, ' size of nf_att_text.  Increase to a'
c         print*, ' minimum value of ', nf_att_len
c         STOP 'ABORT, UNABLE TO CONTINUE'
c      endif

C..Added by G. Thompson to work around problem of potentially
C.. exceeding declared length of vardesc due to a very long
C.. description in the wrfout file, variable description.
      n1 = LEN(plchun)
      iendn1 = 1
      do n = n1, 1, -1
         ich = ichar(plchun(n:n))
	 if ( ich .ne. 32 ) then
            iendn1 = n
	    if ( iendn1 .eq. 1 ) then
         if (.not. ( (ich.ge.65 .and. ich.le.90) .or.     ! Letters A-Z
     &               (ich.ge.97 .and. ich.le.122) .or.    ! Letters a-z
     &               (ich.ge.45 .and. ich.le.58) .or.     ! digits 0-9, also [-./]
     &               (ich.eq.95) .or. ich.eq.32) ) then   ! underscore, blank
c  flag fields have a single character garbage unit - set to blank.
	      plchun = ' ' 
	    endif
	    endif
            goto 44
         endif
      enddo
 44   continue
      if (iendn1.le.0) then
         iendn1 = n1
      endif
      n2 = LEN(nf_att_text)
      iendn2 = 1
      do n = n2, 1, -1
         ich = ichar(nf_att_text(n:n))
	 if ( ich .ne. 32 ) then
c        if (      ( (ich.ge.65 .and. ich.le.90) .or.     ! Letters A-Z
c    &               (ich.ge.97 .and. ich.le.122) .or.    ! Letters a-z
c    &               (ich.ge.45 .and. ich.le.58) .or.     ! digits 0-9, also [-./]
c    &               (ich.eq.95) .or. ich.eq.32) ) then   ! underscore, blank
            iendn2 = n
            goto 45
         endif
      enddo
 45   continue
      if (iendn2.le.0) then
         iendn2 = n2
      endif

C..If we may exceed 64 chars of vardesc, then leave units intact
C.. and shorten text description by amount needed to keep from
C.. memory overwrite.
      vardesc = ' '
      if (iendn1+iendn2 .lt. LEN(vardesc)) then
         iendn3 = nf_att_len
         iendt = iendn1+iendn2 + 2
      else
         iendn3 = MAX(0, nf_att_len - iendn1)
         iendt = LEN(vardesc)
      endif

      if ( iendn1 .eq. 1 .and. plchun(1:1) .eq. ' ' ) then
        vardesc(1:iendt)=nf_att_text(1:iendn3)
      else
        vardesc(1:iendt)=nf_att_text(1:iendn3)//', '//plchun(1:iendn1)
      endif


      icd=1
      if (itype.eq.1) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(097.,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           !scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1) ! orig
            scr3(i,j,k)=nf_tarr3(mkzh-k+1,j,i) ! RGF fix
         enddo
         enddo
         enddo
         call writefile_rdp(scr3,varname,3,icd,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      elseif (itype.eq.2) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(098.,nf_status)
         do j=1,mjx-1
         do i=1,miy-1
            scr2(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
      elseif (itype.eq.3) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3s, nf_tarr3)
         call handle_err(099.,nf_status)
         isp=index(varname,' ')
         if (isp.eq.0.or.isp.eq.10) isp=9
         do k=1,nsoil
            do j=1,mjx-1
            do i=1,miy-1
               scr2(i,j)=nf_tarr3(j,i,k)
            enddo
            enddo
            write(varname(isp:isp+1),'(i2.2)') k
            vardesc=nf_att_text(1:nf_att_len)//
     &         ', layer '//varname(isp:isp+1)//', '//plchun
            call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
         enddo
      elseif (itype.eq.4) then
         nf_tmp_count = nf_tcount3s
         nf_tmp_count(3)=16        ! this will take care of SOILCTOP and SOILCBOT
         if (varname == 'LANDUSEF')  then
             varname = 'LANDUSE'
             if ( ilandset == 2 ) nf_tmp_count(3)=24
             if ( ilandset == 3 ) nf_tmp_count(3)=20
         else if (varname == 'ALBEDO12M') then
             varname = 'ALBEDO'
             nf_tmp_count(3)=12
         else if (varname == 'GREENFRAC') then
             varname = 'GREENFRC'
             nf_tmp_count(3)=12
         end if

         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tmp_count, nf_sarr3)
         call handle_err(101.,nf_status)
         isp=index(varname,' ')
         if (isp.eq.0.or.isp.eq.10) isp=9
         do k=1,nf_tmp_count(3)
            do j=1,mjx-1
            do i=1,miy-1
               scr2(i,j)=nf_sarr3(j,i,k)
            enddo
            enddo
            write(varname(isp:isp+1),'(i2.2)') k
            vardesc=nf_att_text(1:nf_att_len)//
     &         ': month '//varname(isp:isp+1)//', '//plchun
            call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh_out)
         enddo
      endif
c
 229  continue  ! point to jump to if variable is being skipped
c
      enddo     ! end of "other variables (do ivar=1,nvars)" loop
c
 230  continue  ! Destination to jump past "other variables" loop
c
      print*,'Finished reading data for this time.'
c
c   Write out the available xtimes in the ".xtimes" file.  Do this after
c   every time, so that if the user kills the program in mid-run, he'll
c   still have a useful .xtimes file
c
      fname=argum(ncn)(1:iendc)//'.xtimes'
      open (unit=57,file=fname,form='formatted',status='unknown')
      write(57,*) nxtavl
      do i=1,nxtavl
         write(57,'(a10)') cxtimeavl(i)
      enddo
      close (57)
c
 240  continue   ! Destination for skipping this time
c
c=================================================================c
      enddo      ! End of time loop.
c=================================================================c
c
c=================================================================c
      enddo      ! End of file loop.
c=================================================================c
c
 1000 continue   ! destination for jumping out of both loops
c
      print*
      print*,'===================================='
      print*,' We''re outta here like Vladimir !! '
      print*,'===================================='
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine handle_err(rmarker,nf_status)
      include "netcdf.inc"
      integer nf_status
      real rmarker
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'NetCDF error in ripdp_mpas.  Marker = ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop
      endif
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vtodot(slab,maxiy,maxjx)
c
c   This routine converts data that is on the C-grid v-velocity
c   staggered grid to the B-grid velocity staggered grid (known in
c   MM5 lingo as "dot points") 
c
      dimension slab(maxiy,maxjx),rleft(2000)
c
c   Extrapolate out to left and right edges.
c
      do i=1,maxiy
         rleft(i)=(3.*slab(i,1)-slab(i,2))/2.
         slab(i,maxjx)=(3.*slab(i,maxjx-1)-slab(i,maxjx-2))/2.
      enddo
c
c   Interpolate in the interior.
c
      do j=maxjx-1,2,-1
      do i=maxiy,1,-1
         slab(i,j)=.5*(slab(i,j-1)+slab(i,j))
      enddo
      enddo
c
c   Put "rleft" values into slab.
c
      do i=1,maxiy
         slab(i,1)=rleft(i)
      enddo
c
      return
      end
cc                                                                     c
