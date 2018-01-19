      program ripdp_wrfnmm
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIPDP_WRFNMM is a data preparation program that reads in data     c
c   from the WRF NMM modeling system and makes data files             c
c   appropriate for input to the RIP data analysis and visulaization  c
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
c   Modified March 2003 for WRF ARW model output.                     c
c      [Wei Wang (NCAR) and Mark Stoelinga (UW)]                      c
c                                                                     c
c   Modified June 2007 for WRF NMM model output, including treatment   c
c   of NMM's E-grid and the stretched rotated cylindrical             c
c   equidistant (SRCE) map projection.                                c
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
c   set up is for dynamic memory allocation.  This capability,
c   known as adjustable dimensioning of local (non-argument) arrays,
c   is not standard to Fortran 77, but is allowable on some f77
c   compilers.  It is standard to Fortran 90.  As of February 2006,
c   the RIP code assumes it will be compiled with this capability.
c
c   Map transform common block
c
      include "commptf"
c
c   netcdf variables
c
      integer ncid, dimid, nf_status
      character nf_att_text*256, start_date*19
c
      character argum(256)*256,rip_root*256,fname*256
c
      include "netcdf.inc"
c
      include "CMASSI.comm"    ! common block for NMM micro lookup table
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character discard(maxptimes)*16,retain(maxptimes)*16,ptimeunits*1
      namelist/userin/ ptimes,iptimes,ptimeunits,tacc,discard,retain,
     &   iexpandedout,iskpd1,
     &   iinterp,dskmcib,miycorsib,mjxcorsib,nprojib,
     &   xlatcib,xloncib,truelat1ib,truelat2ib,miyib,mjxib,yicornib,
     &   xjcornib,dskmib
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
      if (argum(1)(1:12).eq.'ripdp_wrfnmm') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if ((argum(1)(1:3).eq.'-n '.and.nargum.lt.5).or.
     &    (argum(1)(1:3).ne.'-n '.and.nargum.lt.3)) then
         print*,'Usage:'
         print*,'  ripdp_wrfnmm [-n namelist_file] casename basic|all',
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
      iup = 6
c
c   Get rip_root from environment variable RIP_ROOT
c   if rip_root='/dev/null '
c
      call getenv('RIP_ROOT',rip_root)
      if (rip_root(1:10).eq.'          ') then
         write(iup,*)'Either RIP_ROOT environment variable or rip_root'
         write(iup,*)'is not properly set.'
         stop
      endif
      iendrr=index(rip_root,' ')-1
c
c   Open and read NMM micro lookup table file
c
      fname=rip_root(1:iendrr)//'/eta_micro_lookup.dat'
      open (unit=36,file=fname,form='unformatted',status='old')
      do i=1,3
        read(36)
      enddo
      read(36) massr
      do i=1,5
        read(36)
      enddo
      read(36) massi
      close(36)
c
c   Open first netcdf file and get model dimensions
c
      nf_status = nf_open (argum(nsetsbeg), nf_nowrite, ncid)
      call handle_err(000.,nf_status)
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'TITLE', nf_att_text)
      call handle_err(001.,nf_status)
      if (index(nf_att_text,'OUTPUT FROM WRF').ne.0) then
         print*,'Data is recognized as WRF model output data.'
         nf_status = nf_get_att_text (ncid, nf_global,
     &      'GRIDTYPE', nf_att_text)
         call handle_err(001.,nf_status)
         if (nf_att_text(1:1).eq.'E') then
            print*,'Data is recognized as WRF NMM model output data.'
         else
            stop
         endif
         print*
      else
         stop
      endif
c
c     Some background on dimensions in this program:  The E-grid has
c     mass ("H") and velocity ("V") points that are staggered in a
c     diamond-shaped pattern, as shown below on left side.  Note that an
c     E-grid is essentially a B-grid (as used in MM5), rotated by 45
c     degrees, although that is not relevant to RIPDP/RIP.
c
c                                 |     D   D   D   D   D   D     |
c           H   V   H   V   H     |       C   C   C   C   C       |
c                                 |     D   D   D   D   D   D     |
c           V   H   V   H   V     |       C   C   C   C   C       |
c                                 |     D   D   D   D   D   D     |
c           H   V   H   V   H     |       C   C   C   C   C       |
c                                 |     D   D   D   D   D   D     |
c           V   H   V   H   V     |       C   C   C   C   C       |
c                                 |     D   D   D   D   D   D     |
c           H   V   H   V   H     |       C   C   C   C   C       |
c                                 |     D   D   D   D   D   D     |
c
c     E-grid data are stored in compact arrays that contain only the
c     number of points necessary for either the H or V grid.  One
c     dimension is the number of rows in the E-grid, and the other
c     dimension is the number of H points in and odd-numbered row (or V
c     points in and even-numbered row).  For H (V) arrays, the last
c     element in each even (odd) numbered row is blank.
c
c     RIP was written for B-grid data (as in MM5), so the easiest way to
c     deal with E-grid data is to make it look like B-grid data.  This
c     can be done in two ways.  In the first way, we define a B-grid in
c     which the mass ("cross or C") points collocate with all the H and
c     V points in the E-grid, and the velocity ("dot or D") points are
c     staggered in the usual B-grid way, as shown above on right side.
c     The RIPDP output data files are written using the "compact E-grid"
c     arrays, but after the data are read into RIP, H-point data are
c     transferred directly to overlapping C points, and non-overlapping
c     C points and all D points are interpolated from the H and V point
c     grids.  This is the best way to retain as much of the exact
c     original data as possible, but effectively doubles the number of
c     horizontal grid points in RIP, which can be undesirable.  The
c     second way is to define a completely new B-grid that has no
c     relation to the E-grid points, possibly (or even preferably)
c     including a different map background, but presumably with
c     substantial overlap between the two grids, and a horizontal
c     resolution similar to the effective resolution of the E-grid.  The
c     E-grid data is then bilinearly interpolated to the new B-grid in
c     RIPDP and the new B-grid data is then written out to the RIPDP
c     output data files.  With this method, the fact that the original
c     data was on the E-grid is completely transparent to the RIP
c     plotting program.
c
c     With that background, there are four sets of horizontal dimension
c     parameters that appear in this program:
c
c     miyec,mjxec (ec: "E-grid, compact") are the x and y dimensions of
c     the compact E-grid data arrays.  In above diagram, miyec=5,
c     mjxec=3.
c
c     miyef,mjxef (ef: "E-grid, full") are the x and y dimensions of the
c     full E-grid, including all H and V points.  In general,
c     miyef=miyec and mjxef=(2*mjxec)-1.  In above diagram, miyef=5,
c     mjxef=5.
c
c     miyeb,mjxeb (eb: "B-grid corresponding to E-grid") are the x and y
c     dimensions of the dot-point grid of a B-grid whose cross points
c     are collocated with the H and V points of the E-grid.  In general,
c     miyeb=miyef+1=miyec+1 and mjxeb=mjxef+1=2*mjxec.  In above
c     diagram, miyeb=6, mjxeb=6.  miyeb,mjxeb are not used in this
c     program, but are used in RIP routines that read in the output data
c     from this program if method 1 (described above) is followed.
c
c     miyib,mjxib (ib: "interpolatioin B-grid") are the x and y
c     dimensions of the dot-point grid of a B-grid that the user defines
c     if method 2 (interpolation, described above) is followed.  There
c     is no relationship between miyib,mjxib and the E-grid dimensions,
c     since the user is free to define the interpolation grid in any
c     manner he/she chooses.
c
      nf_status = nf_inq_dimid (ncid, 'south_north', dimid)
      call handle_err(002.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, miyec)
      call handle_err(003.,nf_status)
c
      nf_status = nf_inq_dimid (ncid, 'west_east', dimid)
      call handle_err(004.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, mjxec)
      call handle_err(005.,nf_status)
c
      nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid)
      call handle_err(006.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, mkzh)
      call handle_err(007.,nf_status)
c
c   Get start date string and create RIP parameters for start time.
c
      nf_status = nf_get_att_text (ncid,nf_global,
     &   'START_DATE', start_date)
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
c
c     Get dimensions of interpolation domain, if interpolation will be
c     used.  miyib_dim, mjxib_dim are just used as arguments to pass
c     into subroutine process.  They are set to 1 if no interpolation
c     will be done.
c
      if (nnl.ne.0) then
         iinterp=0  ! temporarily set this, just for this chunk of code
         open (unit=7,file=argum(nnl),form='formatted',
     &         status='old')
         read (7,userin)
         if (iinterp.eq.1) then
            miyib_dim=miyib
            mjxib_dim=mjxib
            miyef_dim=miyec
            mjxef_dim=(2*mjxec)-1
         else
            miyib_dim=1
            mjxib_dim=1
            miyef_dim=1
            mjxef_dim=1
         endif
         close (7)
      else
         miyib_dim=1
         mjxib_dim=1
         miyef_dim=1
         mjxef_dim=1
      endif
c
      call process(miyec,mjxec,miyef_dim,mjxef_dim,
     &   miyib_dim,mjxib_dim,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes,mdateb,mhourb,rhourb)
      stop
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine process(miyec,mjxec,miyef_dim,mjxef_dim,
     &   miyib_dim,mjxib_dim,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes,mdateb,mhourb,rhourb)
c
c   This subroutine does most of the "work".
c
c   miyec and mjxec are the E-grid compact data array dimensions.
c   miyef and mjxef are the E-grid full dimensions (all H + V points).
c   miyib_dim and mjxib_dim are the dimensions of the interpolation
c      B-grid (set equal to 1 if no interpolation will occur).
c   mkzh is number of model levels in the domain.
c   argum carries the names of the model data files.
c   nsets is the number of files in the model dataset.
c   nsetsbeg is the element of argum that holds the first file name
c      of the model dataset.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      character argum(256)*256
c
      include "netcdf.inc"
c
      include "CMASSI.comm"    ! common block for NMM micro lookup table
c
      dimension ter(miyec,mjxec),xmap(miyec,mjxec),xlat(miyec,mjxec),
     &   xlon(miyec,mjxec),
     &   cor(miyec,mjxec),rtc(miyec,mjxec),rte(miyec,mjxec),
     &   tgk(miyec,mjxec),tmk(miyec,mjxec,mkzh),
     &   uuu(miyec,mjxec,mkzh),ght(miyec,mjxec,mkzh),
     &   www(miyec,mjxec,mkzh),
     &   vvv(miyec,mjxec,mkzh),prs(miyec,mjxec,mkzh),sfp(miyec,mjxec),
     &   qvp(miyec,mjxec,mkzh),scr3wlev(miyec,mjxec,mkzh+1),
     &   qcw(miyec,mjxec,mkzh),qra(miyec,mjxec,mkzh),
     &   qci(miyec,mjxec,mkzh),
     &   qnr(miyec,mjxec,mkzh),qni(miyec,mjxec,mkzh),
     &   qsn(miyec,mjxec,mkzh),alph(miyec,mjxec,mkzh+1),
     &   scr3(miyec,mjxec,mkzh),scr2(miyec,mjxec),
     &   xsrcc(miyib_dim,mjxib_dim),ysrcc(miyib_dim,mjxib_dim),
     &   xsrcd(miyib_dim,mjxib_dim),ysrcd(miyib_dim,mjxib_dim),
     &   xlatcrib(miyib_dim,mjxib_dim),xloncrib(miyib_dim,mjxib_dim),
     &   xlatdib(miyib_dim,mjxib_dim),xlondib(miyib_dim,mjxib_dim),
     &   unorthib(miyib_dim,mjxib_dim),vnorthib(miyib_dim,mjxib_dim),
     &   tranx(miyib_dim,mjxib_dim),trany(miyib_dim,mjxib_dim),
     &   arrib(miyib_dim,mjxib_dim,mkzh),
     &   arrib2(miyib_dim,mjxib_dim,mkzh),
     &   arref(miyef_dim,mjxef_dim,mkzh)
c
      character varname*10,fname*256,cxtime*10,
     &   cxtimeavl(256)*10
c
      real nlice,n0r,nlimin,nlimax

      real xlatc_mid, xlonc_mid
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   netcdf variables
c
      real nf_uarr3(mjxec,miyec,mkzh),
     &     nf_tarr3(mjxec,miyec,mkzh),nf_warr3(mjxec,miyec,mkzh+1),
     &     nf_tarr2(mjxec,miyec)
      integer ncid, ndims, nvars, ngatts, unlimdimid, dimid,
     &   nf_status, varid, nf_att_int, iprocvarid(200)
      integer nf_ustart3(4),nf_tstart3(4),nf_wstart3(4),
     &   nf_tstart2(3)
      integer nf_ucount3(4),nf_tcount3(4),nf_wcount3(4),
     &   nf_tcount2(3),nf_tcount3s(4)
      real nf_att_real
      character nf_att_text*256, nf_varname*16,
     &   wrftimes(nwrftimes)*19
      integer dimid_tm, dimid_we, dimid_sn, dimid_bt, dimid_sls
      integer vardimids(20), nf_att_len
      dimension hlat(mjxec,miyec),hlon(mjxec,miyec),vlat(mjxec,miyec),
     &   vlon(mjxec,miyec)
c      dimension xlatetall(miyec,mjxec),xlonetall(miyec,mjxec)
c
c minfo string
c
      character minfostring*88
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character discard(maxptimes)*16,retain(maxptimes)*16,ptimeunits*1
      namelist/userin/ ptimes,iptimes,ptimeunits,tacc,discard,retain,
     &   iexpandedout,iskpd1,
     &   iinterp,dskmcib,miycorsib,mjxcorsib,nprojib,
     &   xlatcib,xloncib,truelat1ib,truelat2ib,miyib,mjxib,yicornib,
     &   xjcornib,dskmib
c
      print*,'Welcome to your friendly RIPDP (V4.7) output file !'    ! January 2017
c
c   Set size of full E-grid domain (H + V points).
c
      miyef=miyec
      mjxef=(2*mjxec)-1
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
c   Set default namelist values.
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
      iinterp=0
      dskmcib=50.
      miycorsib=100
      mjxcorsib=100
      nprojib=1
      xlatcib=45.
      xloncib=-90.
      truelat1ib=30.
      truelat2ib=60.
      miyib=75
      mjxib=75
      yicornib=25.
      xjcornib=25.
      dskmib=25.
c
c   Read the namelist values.
c
      if (nnl.ne.0) then
         open (unit=iuinput,file=argum(nnl),form='formatted',
     &         status='old')
         read (iuinput,userin)
      endif
      iexpandedout=0  ! not used
      iexpanded=0     ! not used
      ioffexp=0       ! not used
      joffexp=0       ! not used
      iskpd1=0        ! not used
      tacch=tacc/3600.
c
      ibasic=0
      if (argum(ndd)(1:5).eq.'basic') ibasic=1
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
      do i=1,200
         iprocvarid(i)=0
      enddo
c
c   Get some dimension IDs
c
      nf_status = nf_inq_dimid (ncid, 'Time', dimid_tm)
      call handle_err(018.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'west_east', dimid_we)
      call handle_err(019.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'south_north', dimid_sn)
      call handle_err(020.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid_bt)
      call handle_err(021.,nf_status)
      nf_status = nf_inq_dimid (ncid, 'soil_layers_stag', dimid_sls)
      call handle_err(022.,nf_status)
c
c   Get array of model output dates/times
c
      nf_status = nf_inq_varid (ncid, 'Times', varid)
      call handle_err(023.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_text (ncid, varid, wrftimes)
      call handle_err(024.,nf_status)
c
c   Set values of "count" and "start" arrays for reading netcdf data.
c
      nf_ustart3(1)=1
      nf_ustart3(2)=1
      nf_ustart3(3)=1
      nf_ustart3(4)=1     ! will be reset at each time in time loop
      nf_ucount3(1)=mjxec
      nf_ucount3(2)=miyec
      nf_ucount3(3)=mkzh
      nf_ucount3(4)=1     ! will be reset at each time in time loop
      nf_tstart3(1)=1
      nf_tstart3(2)=1
      nf_tstart3(3)=1
      nf_tstart3(4)=1     ! will be reset at each time in time loop
      nf_tcount3(1)=mjxec
      nf_tcount3(2)=miyec
      nf_tcount3(3)=mkzh
      nf_tcount3(4)=1     ! will be reset at each time in time loop
      nf_wstart3(1)=1
      nf_wstart3(2)=1
      nf_wstart3(3)=1
      nf_wstart3(4)=1     ! will be reset at each time in time loop
      nf_wcount3(1)=mjxec
      nf_wcount3(2)=miyec
      nf_wcount3(3)=mkzh+1
      nf_wcount3(4)=1     ! will be reset at each time in time loop
      nf_tstart2(1)=1
      nf_tstart2(2)=1
      nf_tstart2(3)=1     ! will be reset at each time in time loop
      nf_tcount2(1)=mjxec
      nf_tcount2(2)=miyec
      nf_tcount2(3)=1     ! will be reset at each time in time loop
c
      nf_status = nf_inq_dimid (ncid, 'soil_layers_stag', dimid)
      call handle_err(031.,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nsoil)
      call handle_err(032.,nf_status)
      nf_tcount3s(1)=mjxec
      nf_tcount3s(2)=miyec
      nf_tcount3s(3)=nsoil
      nf_tcount3s(4)=1     ! will be reset at each time in time loop
c
c   Get information from netcdf global attributes
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LON', nf_att_real)
      if (nf_status .eq. nf_noerr) xlonc=nf_att_real
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LAT', nf_att_real)
      if (nf_status .eq. nf_noerr) xlatc=nf_att_real
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DX', nf_att_real)
      if (nf_status .eq. nf_noerr) dxdeg=nf_att_real
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DY', nf_att_real)
      if (nf_status .eq. nf_noerr) dydeg=nf_att_real
      ds=dydeg*111111.111111
      dskm=ds*.001
c
c     True latitudes are not used for the NMM's SRCE map projection, but
c     "truelat1" will be used to hold the important ratio of
c     dxdeg/dydeg, which I refer to as the "stretch factor" for the SRCE
c     projection (normally slightly > 1.0).
c
      truelat1=dxdeg/dydeg
      truelat2=0.
c
c     "Coarse domain" dimensions are not defined in WRF output.
c     However, RIP is going to continue to work in this framework for
c     historical purposes.  A "coarse" domain is contrived here, which
c     is a B-grid whose dot-point grid is of size 100x100 and
c     grid space ds = 1 deg lat = 111.111111111 km, centered at latitude
c     = CEN_LAT and longitude = CEN_LON.  It is based on the same map
c     background (SRCE) as the NMM "nest" domain whose output is being
c     processed, and also has the same ratio of dxdeg/dydeg as the
c     "nest".
c
      miycors=100
      mjxcors=100
      dskmc=111.111111111  ! km
      dsc=dskmc*1000.
      nproj=4   ! For NMM, always SRCE
      refrat=dsc/ds  ! refinement ratio of grid being processed w.r.t.
                     ! contrived coarse grid
c
c   Set up map transformation stuff
c
      call premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,truelat1,truelat2,6)
c
c     As it turns out (after much frustrating investigation), the values
c     of dxdeg and dydeg in the model output global attributes are not
c     necesarilly consistent with the lat/lon arrays in the model
c     output, at least at the time this version of ripdp_wrfnmm was
c     coded.  This unhealthy situation comes about apparently because
c     the model user can set dx and dy at model startup to values that
c     are different from those that were used to create the initial
c     condition grids in the SI.  While I am not sure if the new dx/dy
c     values actually get used in the model run (which would be insane),
c     I *am* sure that they get transferred to the global attributes of
c     the model output.  Therefore, from RIP's perspective, I will
c     assume that the lat/lon arrays are based on the correct dx/dy,
c     whereas the values of dx/dy in the global attributes are not
c     trustworthy.  The values of dx/dy that *are* consistent with the
c     lat/lon arrays can be calculated with some fancy coding:
c
c     First get lat/lon arrays from first model output time
c
      nf_tstart2(3)=1
      nf_status = nf_inq_varid (ncid, 'GLAT', varid)
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(053.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         xlat(i,j)=nf_tarr2(j,i)/rpd   ! radians to degrees
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'GLON', varid)
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(055.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         xlon(i,j)=nf_tarr2(j,i)/rpd   ! radians to degrees
      enddo
      enddo

c
c..Check here if the grid center point is actually nearby the
c.. (CEN_LAT,CEN_LON) position.  For a NMM grid1, they should be,
c.. but if the model ran with a nest and we are processing this grid,
c.. then we could be far from center and things needs adjusting.
c

      xlatc_mid = xlat(nint(miyec*.5),nint(mjxec*.5))
      xlonc_mid = xlon(nint(miyec*.5),nint(mjxec*.5))
      if (ABS(xlatc_mid-xlatc).gt.0.075 .or.
     &    ABS(xlonc_mid-xlonc).gt.0.075) then
c     print*, ' DEBUG-GT, midpt not co-located: ',
c    &     xlatc, xlonc, xlatc_mid, xlonc_mid
      xlatc = xlatc_mid
      xlonc = xlonc_mid
      call premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,truelat1,truelat2,6)
      endif
c
c   Now figure out dx/dy consistent with lat/lon arrays
c
      iym=miyec/2+1
      iyu=iym+10
      iyd=iym-10
      jxm=mjxec/2+1
      jxr=jxm+5
      jxl=jxm-5
c
      call maptform(riy_l,rjx_l,xlat(iym,jxl),xlon(iym,jxl),-1)
      call maptform(riy_r,rjx_r,xlat(iym,jxr),xlon(iym,jxr),-1)
      call maptform(riy_d,rjx_d,xlat(iyd,jxm),xlon(iyd,jxm),-1)
      call maptform(riy_u,rjx_u,xlat(iyu,jxm),xlon(iyu,jxm),-1)
      dxfac=(rjx_r-rjx_l)*refrat/20.
      dyfac=(riy_u-riy_d)*refrat/20.
      dxdegnew=dxdeg*dxfac
      dydegnew=dydeg*dyfac
      print*,'dx,dy in global attributes =',dxdeg,dydeg
      print*,'dx,dy consistent with lat/lon arrays =',
     &   dxdegnew,dydegnew
      if (ABS(dxfac-1.0).gt.1.E-8.or.ABS(dyfac-1.0).gt.1.E-8) then
         print*,'These are not consistent, so RIP will assume the'
         print*,'value consistent with the lat/lon arrays is correct.'
         print*,' computed dxfac,dyfac = ',dxfac,dyfac
         dxdeg=dxdegnew
         dydeg=dydegnew
c
c      Re-do map transformation stuff
c
         ds=dydeg*111111.111111
         dskm=ds*.001
         truelat1=dxdeg/dydeg
         truelat2=0.
         miycors=100
         mjxcors=100
         dskmc=111.111111111  ! km
         dsc=dskmc*1000.
         nproj=4   ! For NMM, always SRCE
         refrat=dsc/ds
         call premaptform(dskmc,miycors,mjxcors,nproj,
     &      xlatc,xlonc,truelat1,truelat2,6)
      else
         print*,'Close enough.  No adjustment performed.'
         print*, ' dxfac,dyfac = ', dxfac,dyfac
      endif
c
c   Land use data set:
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'MMINLU', nf_att_text)
      call handle_err(041.,nf_status)
      if (nf_att_text(1:5).eq.'USGS ') then
         ilandset=2    ! WRF SI only uses USGS 24-category land use data set
                       !   right now.
      else
         print*,'Unexpected land use data set specified.'
         print*,'MMINLU=',nf_att_text(1:5)
         print*,'Setting ilandset to 2 (USGS 24-cat) and continuing.'
      endif
c
c   Other
c
      iplevdata=8   ! Anything larger than 3 means terrain-following data
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
      nf_status = nf_inq_varid (ncid, 'QSNOW', varid)
      if (nf_status .eq. nf_noerr) then
         print *, 'Found QSNOW variable, set iice = 1'
         iice=1
      endif
c
c     If QSNOW was not found, also check if Ferrier micro arrays are
c     there (e.g., CWM, the condensate water mass) in which case iice
c     would also be 1.
c
          iferrier=0
cjkw      if (iice.eq.0) then
cjkw         nf_status = nf_inq_varid (ncid, 'CWM', varid)
cjkw         if (nf_status .eq. nf_noerr) then
cjkw            print*,'Found CWM (condensate water mass) variable.'
cjkw            print*,'This implies Ferrier microphysics was used, so'
cjkw            print*,'setting iice = 1, iferrier=1.'
cjkw            iice=1
cjkw            iferrier=1
cjkw         endif
cjkw      endif
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
     &   'TITLE', nf_att_text)
      call handle_err(042.,nf_status)
      istart=index(nf_att_text,'WRF V')+4
      iendtmp=min(8,index(nf_att_text(istart:),' ')-1)
      iend=iendtmp+istart-1
      minfostring='Model Info: '//nf_att_text(istart:iend)  ! First 20 chars.
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'CU_PHYSICS', nf_att_int)
      call handle_err(043.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(21:28)=' No Cu'
      elseif (nf_att_int.eq.1) then
         minfostring(21:28)=' KF'
      elseif (nf_att_int.eq.2) then
         minfostring(21:28)=' BMJ'
      elseif (nf_att_int.eq.3) then
         minfostring(21:28)=' G-D Ens'
      elseif (nf_att_int.eq.99) then
         minfostring(21:28)=' Smp A-S'
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'BL_PBL_PHYSICS', nf_att_int)
      if (nf_status .ne. nf_noerr) then
         nf_status = nf_get_att_int (ncid, nf_global,
     &      'SF_PBL_PHYSICS', nf_att_int)
      endif
      call handle_err(044.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(29:37)=' No PBL'
      elseif (nf_att_int.eq.1) then
         minfostring(29:37)=' YSU PBL'
      elseif (nf_att_int.eq.2) then
         minfostring(29:37)=' MYJ PBL'
      elseif (nf_att_int.eq.3) then
         minfostring(29:37)=' GFS PBL'
      elseif (nf_att_int.eq.99) then
         minfostring(29:37)=' MRF PBL'
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'MP_PHYSICS', nf_att_int)
      call handle_err(045.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(38:49)=' No microph'
      elseif (nf_att_int.eq.1) then
         minfostring(38:49)=' Kessler'
      elseif (nf_att_int.eq.2) then
         minfostring(38:49)=' Lin et al'
      elseif (nf_att_int.eq.3) then
         minfostring(38:49)=' WSM 3class'
      elseif (nf_att_int.eq.4) then
         minfostring(38:49)=' WSM 5class'
      elseif (nf_att_int.eq.5) then
         minfostring(38:49)=' Ferrier'
         print*,'This implies Ferrier microphysics was used, so'
         print*,'setting iice = 1, iferrier=1.'
         iice=1
         iferrier=1
      elseif (nf_att_int.eq.6) then
         minfostring(38:49)=' WSM 6class'
      elseif (nf_att_int.eq.8) then
         minfostring(38:49)=' Thompson'
      elseif (nf_att_int.eq.99) then
         minfostring(38:49)=' Zhao-Carr'
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'BL_SURFACE_PHYSICS', nf_att_int)
      if (nf_status .ne. nf_noerr) then
         nf_status = nf_get_att_int (ncid, nf_global,
     &      'SF_SURFACE_PHYSICS', nf_att_int)
      endif
      call handle_err(046.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(50:60)=' No SFC'
      elseif (nf_att_int.eq.1) then
         minfostring(50:60)=' Ther-Diff'
      elseif (nf_att_int.eq.2) then
         minfostring(50:60)=' Noah LSM'
      elseif (nf_att_int.eq.3) then
         minfostring(50:60)=' RUC LSM'
      elseif (nf_att_int.eq.99) then
         minfostring(50:60)=' NMM Noah'
      endif

      if (ds.ge.100000.) then
         write(minfostring(61:68),'(i3,'' km, '')') nint(.001*ds)
      elseif (ds.ge.10000.) then
         write(minfostring(61:68),'(i2,'' km, '')') nint(.001*ds)
      elseif (ds.ge.1000.) then
         write(minfostring(61:68),'(f3.1,'' km, '')') .001*ds
      else
         write(minfostring(61:68),'(i3,'' m, '')') nint(ds)
      endif
c
      write(minfostring(69:81),'(i3,'' levels, '')') mkzh
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DT', nf_att_real)
      call handle_err(047.,nf_status)
      write(minfostring(82:88),'(i3,'' sec'')') nint(nf_att_real)
      write(58,'(a)') minfostring(1:88)
c
c   Second line of minfo file
c
      do i = 1, 88
	minfostring(i:i) = ' '
      enddo
      ia = 22
      ib = 29
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'RA_LW_PHYSICS', nf_att_int)
      call handle_err(048.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(ia:ib)='LW: none'
      elseif (nf_att_int.eq.1) then
         minfostring(ia:ib)='LW: RRTM'
      elseif (nf_att_int.eq.99) then
         minfostring(ia:ib)='LW: GFDL'
      endif
      ia = 30
      ib = 41
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'RA_SW_PHYSICS', nf_att_int)
      call handle_err(049.,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(ia:ib)=' SW: none'
      elseif (nf_att_int.eq.1) then
         minfostring(ia:ib)=' SW: Dudhia'
      elseif (nf_att_int.eq.2) then
         minfostring(ia:ib)=' SW: Goddard'
      elseif (nf_att_int.eq.99) then
         minfostring(ia:ib)=' SW: GFDL'
      endif
      write(58,'(a)') minfostring(1:88)
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
     &'map projection (0: none/ideal, 1: LC, 2: PS, 3: ME, 4: SRCE)'
      ihrip(1)=nproj
      chrip(2)=
     &'number of dot points in the y-direction (coarse domain)'
      ihrip(2)=miycors
      chrip(3)=
     &'number of dot points in the x-direction (coarse domain)'
      ihrip(3)=mjxcors
      chrip(4)=
     &'number of E-grid rows (this domain)'
      ihrip(4)=miyec
      chrip(5)=
     &'number of E-grid H points in an odd row (this domain)'
      ihrip(5)=mjxec
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
      nf_tstart3(4)=iwt
      nf_wstart3(4)=iwt
      nf_tstart2(3)=iwt
c
c   Create mdate, rhour, xtime.
c
      read(wrftimes(iwt),'(i4,5(1x,i2))')iyr4,imo,idy,ihr,imn,isc
      if (iyr4.eq.1) then  !idealized (non-real-weather) model output
         iyr4=1940
      endif
      iyr=mod(iyr4,100)
      mdate=1000000*iyr+10000*imo+
     &   100*idy+ihr
      rhour=imn/60.+isc/3600.
      call mconvert(mdate,mhour,1,1940)
      xtime=float(mhour-mhourb)+rhour-rhourb
c
c   Get latitude, longitude arrays
c
      nf_status = nf_inq_varid (ncid, 'GLAT', varid)
      call handle_err(052.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(053.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         xlat(i,j)=nf_tarr2(j,i)/rpd   ! radians to degrees
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'GLON', varid)
      call handle_err(054.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(055.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         xlon(i,j)=nf_tarr2(j,i)/rpd   ! radians to degrees
      enddo
      enddo
c
c   test to see if the "etall" routine, using parameters from the model
c   output, produces the same lat/lon
c   points as in model output arrays
c
c      call ETALL(mjxec,miyec,xlatc,xlonc,dxdeg,dydeg,
c     & HLAT,HLON,VLAT,VLON)
c      print*,'EGRID CORNERS from ETALL'
c      print*,'hlat(1,1),hlon(1,1):',hlat(1,1),hlon(1,1)
c      print*,'hlat(1,miyec),hlon(1,miyec):',hlat(1,miyec),hlon(1,miyec)
c      print*,'hlat(mjxec,1),hlon(mjxec,1):',hlat(mjxec,1),hlon(mjxec,1)
c      print*,'hlat(mjxec,miyec),hlon(mjxec,miyec):',
c     &   hlat(mjxec,miyec),hlon(mjxec,miyec)
c      do j=1,mjxec
c      do i=1,miyec
c         xlatetall(i,j)=hlat(j,i)
c         xlonetall(i,j)=hlon(j,i)
c      enddo
c      enddo
c
c     The contrived "coarse" domain and the "nest" domain have the same
c     center point, so it is easy to calculate yicorn and xjcorn, the y
c     and x locations of the lower left H-grid point of the E-grid
c     "nest" domain with respect to the "coarse" domain dot-point grid.
c
      yitemp=0.5*(1.+miycors)
      xjtemp=0.5*(1.+mjxcors)
      yioffset=.5*(float(miyef)-1.)/refrat
      xjoffset=.5*(float(mjxef)-1.)/refrat
      yicorn=yitemp-yioffset
      xjcorn=xjtemp-xjoffset
c
c   Set RIP header values that are time-dependent
c
      ihrip(11)=mdate
      rhrip(7)=yicorn
      rhrip(8)=xjcorn
      rhrip(14)=rhour
      rhrip(15)=xtime
c
c   Set up things for horizontal domain interpolation, if asked for.
c
      if (iinterp.eq.1) then
c
c      First, temporarily reset the parameters that affect the map
c      transformation to those of the interpolated domain
c
         call premaptform(dskmcib,miycorsib,mjxcorsib,nprojib,
     &      xlatcib,xloncib,truelat1ib,truelat2ib,6)
c
c      Now get the lat/lon values of all the grid points
c      on the interpolation domain.  Also get unorthib,vnorthib, which
c      is a unit vector that points northward at each dot point
c      in the interpolation domain.
c
         refratib=dskmcib/dskmib
         do j=1,mjxib
            rjxc=xjcornib+(j-0.5)/refratib
            rjxd=xjcornib+(j-1.0)/refratib
         do i=1,miyib
            riyc=yicornib+(i-0.5)/refratib
            riyd=yicornib+(i-1.0)/refratib
            call maptform(riyc,rjxc,xlatcrib(i,j),xloncrib(i,j),1)
            call maptform(riyd,rjxd,xlatdib(i,j),xlondib(i,j),1)
c
            xlatdib2=min(89.9999,xlatdib(i,j)+.1)
            call maptform(riyd2,rjxd2,xlatdib2,xlondib(i,j),-1)
            unn=rjxd2-rjxd
            vnn=riyd2-riyd
            dnorth=sqrt(unn**2+vnn**2)
            unorthib(i,j)=unn/dnorth
            vnorthib(i,j)=vnn/dnorth
         enddo
         enddo
c
c      Reset the map-tran parameters back to the source domain values
c
         call premaptform(dskmc,miycors,mjxcors,nproj,
     &      xlatc,xlonc,truelat1,truelat2,6)
c
c      Now get the x/y values (in source-domain coordinate system)
c      of all the grid points on the interpolation domain.  Also,
c      get a unit vector that points northward in the source domain,
c      but valid at the lat/lon locations of the dot points in the 
c      interpolation domain (unorthsrcib,vnorthsrcib).  
c
         do j=1,mjxib
         do i=1,miyib
            call maptform(riyc,rjxc,xlatcrib(i,j),xloncrib(i,j),-1)
            xsrcc(i,j)=1.+(rjxc-xjcorn)*refrat
            ysrcc(i,j)=1.+(riyc-yicorn)*refrat
            call maptform(riyd,rjxd,xlatdib(i,j),xlondib(i,j),-1)
            xsrcd(i,j)=1.+(rjxd-xjcorn)*refrat
            ysrcd(i,j)=1.+(riyd-yicorn)*refrat
c
            xlatdib2=min(89.9999,xlatdib(i,j)+.1)
            call maptform(riyd2,rjxd2,xlatdib2,xlondib(i,j),-1)
            unn=rjxd2-rjxd
            vnn=riyd2-riyd
            dnorth=sqrt(unn**2+vnn**2)
            unorthibe=unn/dnorth
            vnorthibe=vnn/dnorth
            tranx(i,j)= unorthib(i,j)*unorthibe+vnorthib(i,j)*vnorthibe
            trany(i,j)=-unorthib(i,j)*vnorthibe+vnorthib(i,j)*unorthibe
         enddo
         enddo
c
c      Reset the map-tran parameters back to the interpolation domain values
c
         call premaptform(dskmcib,miycorsib,mjxcorsib,nprojib,
     &      xlatcib,xloncib,truelat1ib,truelat2ib,6)
c
c         print*,'Corners:'
c         print*,'xsrcc(1,1),ysrcc(1,1)=',
c    &       '   ',xsrcc(1,1),ysrcc(1,1)
c         print*,'xsrcc(1,mjxib-1),ysrcc(1,mjxib-1)=',
c    &       '   ',xsrcc(1,mjxib-1),ysrcc(1,mjxib-1)
c         print*,'xsrcc(miyib-1,1),ysrcc(miyib-1,1)=',
c    &       '   ',xsrcc(miyib-1,1),ysrcc(miyib-1,1)
c         print*,'xsrcc(miyib-1,mjxib-1),ysrcc(miyib-1,mjxib-1)=',
c    &       '   ',xsrcc(miyib-1,mjxib-1),ysrcc(miyib-1,mjxib-1)
c         print*,'xsrcd(1,1),ysrcd(1,1)=',
c    &       '   ',xsrcd(1,1),ysrcd(1,1)
c         print*,'xsrcd(1,mjxib),ysrcd(1,mjxib)=',
c    &       '   ',xsrcd(1,mjxib),ysrcd(1,mjxib)
c         print*,'xsrcd(miyib,1),ysrcd(miyib,1)=',
c    &       '   ',xsrcd(miyib,1),ysrcd(miyib,1)
c         print*,'xsrcd(miyib,mjxib),ysrcd(miyib,mjxib)=',
c    &       '   ',xsrcd(miyib,mjxib),ysrcd(miyib,mjxib)
c
c      Change rip header values that need to be changed for the
c      interpolation domain
c
         rhrip(7)=yicornib
         rhrip(8)=xjcornib
         chrip(4)=
     &   'number of dot points in the y-direction (this domain)'
         chrip(5)=
     &   'number of dot points in the x-direction (this domain)'
         ihrip(1)=nprojib
         ihrip(2)=miycorsib
         ihrip(3)=mjxcorsib
         ihrip(4)=miyib
         ihrip(5)=mjxib
         rhrip(1)=truelat1ib
         rhrip(2)=truelat2ib
         rhrip(3)=xlatcib
         rhrip(4)=xloncib
         rhrip(5)=dskmcib
         rhrip(6)=dskmib
c
      endif
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
c   Get U and (maybe) write out.
c
      nf_status = nf_inq_varid (ncid, 'U', varid)
      call handle_err(056.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_ustart3,
     &   nf_ucount3, nf_uarr3)
      call handle_err(057.,nf_status)
      do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            uuu(i,j,k)=nf_uarr3(j,i,mkzh-k+1)
         enddo
         enddo
      enddo
      vardesc='Horizontal wind (x-comp.), m/s'
      plchun='m s~S~-1~N~'
      icd=2
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(uuu,'uuu       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(uuu,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcd,ysrcd,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         iwarn=0
         do j=1,mjxib
         do i=1,miyib
            if (arrib(i,j,1).eq.rmsg) then
               iwarn=1
            endif
         enddo
         enddo
         if (iwarn.eq.1) then
            write(iup,*)
            write(iup,*)'Warning: The interpolation domain you have'
            write(iup,*)'defined contains points outside of the NMM'
            write(iup,*)'E-grid.  These points will be assigned the'
            write(iup,*)'"missing data" value.  Some parts of RIP'
            write(iup,*)'recognize the "missing data" value and deal'
            write(iup,*)'with it appropriately, but some parts do not,'
            write(iup,*)'yielding unpredictable results.'
            write(iup,*)'To avoid this problem, redefine the interp-'
            write(iup,*)'olation domain to a smaller area that is'
            write(iup,*)'completely enclosed by the E-grid.  This'
            write(iup,*)'must be done through trial and error.'
            write(iup,*)
         endif
c
c         Will write out data after vvv is interpolated and entire vector
c         is transformed (below).
c
      endif
c
c   Get V and write out.
c
      nf_status = nf_inq_varid (ncid, 'V', varid)
      call handle_err(058.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_ustart3,
     &   nf_ucount3, nf_uarr3)
      call handle_err(059.,nf_status)
      do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            vvv(i,j,k)=nf_uarr3(j,i,mkzh-k+1)
         enddo
         enddo
      enddo
      vardesc='Horizontal wind (y-comp.), m/s'
      plchun='m s~S~-1~N~'
      icd=2
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(vvv,'vvv       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(vvv,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcd,ysrcd,arrib2,miyib,mjxib,mkzh,icd,ndim,rmsg)
c
c      Transform the u/v values from the E-grid reference frame
c      to the interpolation grid reference frame
c
         do k=1,mkzh
         do j=1,mjxib
         do i=1,miyib
            uuue=arrib(i,j,k)
            vvve=arrib2(i,j,k)
            arrib(i,j,k) =tranx(i,j)*uuue-trany(i,j)*vvve
            arrib2(i,j,k)=tranx(i,j)*vvve+trany(i,j)*uuue
         enddo
         enddo
         enddo
c
c      Write out both transformed components.
c
         call writefile_rdp(arrib,'uuu       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         call writefile_rdp(arrib2,'vvv       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Get temperature and write out.
c
      nf_status = nf_inq_varid (ncid, 'T', varid)
      call handle_err(074.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &   nf_tcount3, nf_tarr3)
      call handle_err(075.,nf_status)
      do k=1,mkzh
      do j=1,mjxec
      do i=1,miyec
         tmk(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
      enddo
      enddo
      enddo
      vardesc='Temperature, K'
      plchun='K'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(tmk,'tmk       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(tmk,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'tmk       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c     If available, get mixing ratio (Q) in kg/kg, convert to g/kg,
c     and write out.  Otherwise, fill qvp array with 0s but don't write
c     out.
c
      nf_status = nf_inq_varid (ncid, 'Q', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find mixing ratio.'
         call fillarray(qvp,miyec*mjxec*mkzh,0.)
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(063.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            qvp(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000. ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Water vapor mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qvp,'qvp       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qvp,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qvp       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   Get terrain height (FIS), and write out.
c
      nf_status = nf_inq_varid (ncid, 'FIS', varid)
      call handle_err(070.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(071.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         ter(i,j)=nf_tarr2(j,i)/grav
      enddo
      enddo
      plchun='m'
      vardesc='Terrain height AMSL, m'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(ter,'ter       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(ter,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'ter       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Get pressure at full levels.  Save mean of log pressure.
c   
c
      nf_status = nf_inq_varid (ncid, 'PINT', varid)
      call handle_err(068.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &   nf_wcount3, nf_warr3)
      call handle_err(069.,nf_status)
      do k=1,mkzh+1
      do j=1,mjxec
      do i=1,miyec
         alph(i,j,k)=alog(nf_warr3(j,i,mkzh+1-k+1))
      enddo
      enddo
      enddo
c
c   Convert to pressure at half levels and write out
c
      do k=1,mkzh
      do j=1,mjxec
      do i=1,miyec
         prs(i,j,k)=exp(.5*(alph(i,j,k)+alph(i,j,k+1)))*.01 ! Pa to hPa
      enddo
      enddo
      enddo
      vardesc='Pressure, hPa'
      plchun='hPa'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(prs,'prs       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(prs,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'prs       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Integrate hydrostatic eq. to get heights at full levels
c
      do j=1,mjxec
      do i=1,miyec
         scr3wlev(i,j,mkzh+1)=ter(i,j)
      enddo
      enddo
      do k=mkzh,1,-1
      do j=1,mjxec
      do i=1,miyec
         tv=virtual(tmk(i,j,k),.001*qvp(i,j,k))
         scr3wlev(i,j,k)=scr3wlev(i,j,k+1)+rgas*tv/grav*
     &      (alph(i,j,k+1)-alph(i,j,k))
      enddo
      enddo
      enddo
c
c   Convert to heights at half levels and write out
c
      do k=1,mkzh
      do j=1,mjxec
      do i=1,miyec
         ght(i,j,k)=0.5*(scr3wlev(i,j,k)+scr3wlev(i,j,k+1))
      enddo
      enddo
      enddo
      vardesc='Geopotential height, m'
      plchun='m'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(ght,'ght       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(ght,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'ght       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Now get vertical velocity (convert from m/s to cm/s)
c
      nf_status = nf_inq_varid (ncid, 'W', varid)
      call handle_err(068.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &   nf_wcount3, nf_warr3)
      call handle_err(069.,nf_status)
      do k=1,mkzh+1
      do j=1,mjxec
      do i=1,miyec
         scr3wlev(i,j,k)=100.*nf_warr3(j,i,mkzh+1-k+1)
      enddo
      enddo
      enddo
c
c   Convert to v.v. at half levels and write out
c
      do k=1,mkzh
      do j=1,mjxec
      do i=1,miyec
         www(i,j,k)=0.5*(scr3wlev(i,j,k)+scr3wlev(i,j,k+1))
      enddo
      enddo
      enddo
      vardesc='Vertical velocity, cm/s'
      plchun='cm s~S~-1~N~'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(www,'www       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(www,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'www       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Calculate surface pressure using altimeter equation.
c
      do j=1,mjxec
      do i=1,miyec
         tv=virtual(tmk(i,j,mkzh),.001*qvp(i,j,mkzh))
         sfp(i,j)=prs(i,j,mkzh)*(tv/(tv+ussalr*
     *      (ght(i,j,mkzh)-ter(i,j))))**(-grav/(rgas*ussalr))
      enddo
      enddo
      vardesc='Surface pressure, hPa'
      plchun='hPa'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(sfp,'sfp       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(sfp,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'sfp       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Compute mapscale factor and write out.
c
      if (iinterp.eq.0) then
         do jec=1,mjxec
            jef=2*jec-1
            rjx=xjcorn+(float(jef)-1.)/refrat
         do iec=1,miyec
            ief=iec
            riy=yicorn+(float(ief)-1.)/refrat
            xmap(iec,jec)=xmapcalc(riy,rjx)
         enddo
         enddo
         vardesc='Map factor on cross points'
         plchun='none'
         icd=3
         ndim=2
         call writefile_rdp(xmap,'xmap      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         do j=1,mjxib-1
            rjxib=xjcornib+(float(j)-.5)/refratib
         do i=1,miyib-1
            riyib=yicornib+(float(i)-.5)/refratib
            arrib(i,j,1)=xmapcalc(riyib,rjxib)
         enddo
         enddo
         vardesc='Map factor on cross points'
         plchun='none'
         icd=1
         ndim=2
         call writefile_rdp(arrib,'xmap      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
c
c      Also make and write out map factor on dot points
c
         call xtodot(arrib(1,1,1),miyib,mjxib)
         vardesc='Map factor on dot points'
         plchun='none'
         icd=0
         ndim=2
         call writefile_rdp(arrib,'dmap      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Compute coriolis parameter and write out.
c
      twoomega=4.*pi/(24.*3600)
      do j=1,mjxec
      do i=1,miyec
         cor(i,j)=twoomega*sin(rpd*xlat(i,j))
      enddo
      enddo
      vardesc='Coriolis parameter, per s'
      plchun='s~S~-1~N~'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(cor,'cor       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(cor,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'cor       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c   Already have latitude (GLAT) and longitude (GLON).  Write them out.
c
      vardesc='Latitude, degrees'
      plchun='degrees'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(xlat,'xlat      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(xlat,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'xlat      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
      vardesc='Longitude, degrees'
      plchun='degrees'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(xlon,'xlon      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(xlon,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'xlon      ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c      vardesc='Latitude (etall), degrees'
c      plchun='degrees'
c      icd=3
c      ndim=2
c      if (iinterp.eq.0) then
c         call writefile_rdp(xlatetall,'xlatetall ',ndim,icd,vardesc,
c     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
c      else
c         call hinterp(xlatetall,miyec,mjxec,arref,miyef,mjxef,
c     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
c         call writefile_rdp(arrib,'xlatetall ',ndim,icd,vardesc,
c     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
c      endif
c      vardesc='Longitude (etall), degrees'
c      plchun='degrees'
c      icd=3
c      ndim=2
c      if (iinterp.eq.0) then
c         call writefile_rdp(xlonetall,'xlonetall ',ndim,icd,vardesc,
c     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
c      else
c         call hinterp(xlonetall,miyec,mjxec,arref,miyef,mjxef,
c     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
c         call writefile_rdp(arrib,'xlonetall ',ndim,icd,vardesc,
c     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
c      endif
c
c   Get total accumulated cumulus preciitation (CUPREC), and write out.
c
      nf_status = nf_inq_varid (ncid, 'CUPREC', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find CUPREC.'
         do j=1,mjxec
         do i=1,miyec
            rtc(i,j)=0.
         enddo
         enddo
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(080.,nf_status)
         do j=1,mjxec
         do i=1,miyec
            rtc(i,j)=nf_tarr2(j,i)*1000.  ! (convert to mm)
         enddo
         enddo
         vardesc='Cumulus precip. since h 0, mm'
         plchun='mm'
         icd=3
         ndim=2
         if (iinterp.eq.0) then
            call writefile_rdp(rtc,'rtc       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(rtc,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'rtc       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   Get total accumulated explicit (grid-resolved) + cumulus precipitation
c   (ACPREC), subtract cumulus to yield total acc. expl., and write out.
c
      nf_status = nf_inq_varid (ncid, 'ACPREC', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find ACPREC.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(081.,nf_status)
         do j=1,mjxec
         do i=1,miyec
            rte(i,j)=nf_tarr2(j,i)*1000.-rtc(i,j)  ! (convert to mm)
         enddo
         enddo
         vardesc='Explicit precip. since h 0, mm'
         plchun='mm'
         icd=3
         ndim=2
         if (iinterp.eq.0) then
            call writefile_rdp(rte,'rte       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(rte,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'rte       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   Get ground potential temperature (and sea-surface potential temperature
c   over water) (THS), convert to temperature, and write out.
c
      nf_status = nf_inq_varid (ncid, 'THS', varid)
      call handle_err(082.,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(083.,nf_status)
      do j=1,mjxec
      do i=1,miyec
         gammam=gamma*(1.+gammamd*.001*qvp(i,j,mkzh))
         tgk(i,j)=(nf_tarr2(j,i))*
     &      (sfp(i,j)/1000.)**gammam
      enddo
      enddo
      vardesc='Ground/sea-surface temperature, K'
      plchun='K'
      icd=3
      ndim=2
      if (iinterp.eq.0) then
         call writefile_rdp(tgk,'tgk       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(tgk,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'tgk       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
c     If Ferrier micro (iferrier=1), then need to do special things to
c     create usual micro arrays from Ferrier output.  This follows
c     CALMICT routine in WRF-POST.  Otherwise, do the same thing that is
c     done in ripdp_wrfarw.
c
      if (iferrier.eq.1) then  !   ******** start of "iferrier" if-block
c
c     Obtain Ferrier micro output arrays.
c
      nf_status = nf_inq_varid (ncid, 'CWM', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find CWM.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
c
c         Put CWM (condensate water mass) in "qcw" for now, and leave in kg/kg
c
            qcw(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
      endif
c
      nf_status = nf_inq_varid (ncid, 'F_ICE', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find F_ICE.'
         stop
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
c
c         Put F_ICE (fraction of CWM that is cloud ice) into "qci" for now.
c
            qci(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
      endif
c
      nf_status = nf_inq_varid (ncid, 'F_RAIN', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find F_RAIN.'
         stop
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
c
c         Put F_RAIN (fraction of CWM that is rain) into "qra" for now.
c
            qra(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
      endif
c
      nf_status = nf_inq_varid (ncid, 'F_RIMEF', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find F_RIMEF.'
         stop
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
c
c         Put F_RIMEF ("rime factor") into "qsn" for now.
c
            qsn(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
      endif
c
c   Code below is adapted from WRFPOST routine CALMICT for determining
c   hydrometeor mixing ratios and reflectivities from total condensate
c   and mass fractions.
c
C   INPUT ARGUMENT LIST:
C     P1D     - PRESSURE (PA)
C     T1D     - TEMPERATURE (K)
C     Q1D     - SPECIFIC HUMIDITY (KG/KG)
C     C1D     - TOTAL CONDENSATE (CWM, KG/KG)
C     FI1D    - F_ice (fraction of condensate in form of ice)
C     FR1D    - F_rain (fraction of liquid water in form of rain)
C     FS1D    - F_RimeF ("Rime Factor", ratio of total ice growth 
C                       to deposition growth)
C     CUREFL  - Radar reflectivity contribution from convection (mm**6/m**3)
C
C   OUTPUT ARGUMENT LIST:
C     QW1   - CLOUD WATER MIXING RATIO (KG/KG)
C     QI1   - CLOUD ICE MIXING RATIO (KG/KG)
C     QR1   - RAIN MIXING RATIO (KG/KG)
C     QS1   - "SNOW" (precipitation ice) MIXING RATIO (KG/KG)
C     DBZ1  - Equivalent radar reflectivity factor in dBZ; i.e., 10*LOG10(Z)
c
      tfrz=celkel  
      t_ice=-10.    ! liq./ice cutoff temp. in Ferrier micro (deg. C)
      epsq = 1.e-12 ! threshold mixing ration (kg/kg)
      rhol=1000.    ! rho_liquid_water (kg/m**3)
      dbzmin=-20.
      zmin=10.**(0.1*dbzmin)
      rqr_drmin=n0r0*massr(mdrmin) ! rain content for mean drop dia. of .05 mm
      rqr_drmax=n0r0*massr(mdrmax) ! rain content for mean drop dia. of .45 mm
      c_n0r0=pi*rhol*n0r0
      cn0r0=1.e6/c_n0r0**.25
      cn0r_dmrmin=1./(pi*rhol*dmrmin**4)
      cn0r_dmrmax=1./(pi*rhol*dmrmax**4)
      cice=1.634e13
      nlimin=100.
      nlimax=10.e3
      rd=287.04
      fmw=18.015
      fmd=28.964
      epsfe=fmw/fmd
      onepsfe=1.-epsfe
c
      do k=1,mkzh
      do j=1,mjxec
      do i=1,miyec
c
         p1d=prs(i,j,k)*100. ! pressure (Pa)
         t1d=tmk(i,j,k)      ! temperature (K)
         q1d=qvp(i,j,k)*.001 ! specific humidity (kg/kg)
         c1d=qcw(i,j,k)      ! total condensate (CWM, kg/kg)
         fi1d=qci(i,j,k)     ! F_ICE (fraction of condensate in form of ice)
         fr1d=qra(i,j,k)     ! F_RAIN (fraction of liquid water in form of rain)
         fs1d=qsn(i,j,k)     ! F_RIMEF ("rime factor", ratio of total ice growth 
c                            !   to deposition growth)
c
         qw1=0.              ! cloud water mixing ratio (kg/kg)
         qi1=0.              ! cloud ice mixing ratio (kg/kg)
         qr1=0.              ! rain mixing ratio (kg/kg)
         qs1=0.              ! "snow" (precipitation ice) mixing ratio (kg/kg)
         dbz1=dbzmin         ! equivalent radar reflectivity factor in dbz;
c                            !   i.e., 10*log10(z)
         ztot=0.             ! equivalent radar reflectivity factor
c
         if (c1d .le. epsq) then
c
c         Skip rest of calculatiions if no condensate is present
c
            go to 1010
         else
            wc=c1d
         endif
c
c      Code below is from gsmdrive for determining:
c       qi1 - total ice (cloud ice & snow) mixing ratio
c       qw1 - cloud water mixing ratio
c       qr1 - rain mixing ratio
c
         tc=t1d-tfrz
         fice=fi1d
         frain=fr1d
         if (tc.le.t_ice .or. fice.ge.1.) then
            qi1=wc
         else if (fice .le. 0.) then
            qw1=wc
         else
            qi1=fice*wc
            qw1=wc-qi1
         endif   
         if (qw1.gt.0. .and. frain.gt.0.) then
            if (frain .ge. 1.) then 
               qr1=qw1
               qw1=0.
            else
               qr1=frain*qw1
               qw1=qw1-qr1
            endif 
         endif
         wv=q1d/(1.-q1d)
c
c      Saturation vapor pressure w/r/t water ( >=0C ) or ice ( <0C )
c
c         esat=1000.*fpvs(t1d)
c         (Above uses lookup table for esat, which I don't want to take
c          the trouble to import into RIP from WRF-POST, so I'll use RIP's
c          method that is based upon Bolton.  The WRF-POST version uses
c          esat_liquid for T>0C and esat_ice for T<0C.)
c
         if (t1d.lt.tfrz) then
            esat = ezero * exp( esicon1-esicon2/t1d )
         else
            esat = ezero * exp( eslcon1*(t1d-tfrz)/(t1d-eslcon2) )
         endif
         esat=esat*100.   ! hPa to Pa
         qsat=epsfe*esat/(p1d-esat)
         rho=p1d/(rgas*t1d*(1.+.608*q1d))  ! prs (in Pa)/( R * virtual T (in K)) 
         rrho=1./rho
c
c      Based on code from GSMCOLUMN in model to determine reflectivity from rain
c
         if (qr1 .gt. epsq) then
            rqr=rho*qr1
            if (rqr .le. rqr_drmin) then
               n0r=max(n0rmin, cn0r_dmrmin*rqr)
               indexr=mdrmin
            else if (rqr .ge. rqr_drmax) then
               n0r=cn0r_dmrmax*rqr
               indexr=mdrmax
            else
               n0r=n0r0
               indexr=max( xmrmin, min(cn0r0*rqr**.25, xmrmax) )
            endif
c
c         indexr is the mean drop size in microns; convert to mm
c
            drmm=1.e-3*real(indexr)
            ztot=ztot+0.72*n0r*drmm*drmm*drmm*drmm*drmm*drmm*drmm
         endif        !--- end if (qr1 .gt. epsq) block
c
c      Based on code from GSMCOLUMN in model to determine partition of 
c      total ice into cloud ice & snow (precipitation ice)
c
         if (qi1 .gt. epsq) then
            qice=qi1
            rho=p1d/(rd*t1d*(1.+onepsfe*q1d))
c           The above is actually incorrect (bug)--instead of onepsfe (which
c           is 1-epsfe=1-0.622=0.378), it should be 0.608 (which is 1/epsfe - 1).
            rrho=1./rho
c
c         ax is distance in km (converted from degrees) of the diagonal line
c         between two nearby H points on the E-grid
c
            ax=111.111111*sqrt(dxdeg*dxdeg+dydeg*dydeg)
            ax=min(100., max(5., ax) )
            rhgrd=0.90+.08*((100.-ax)/95.)**.5
            qsigrd=rhgrd*qsat
            wvqw=wv+qw1
c
c   Descriptioin of various things below:
c     flarge  - ratio of number of large ice to total (large & small) ice
c     fsmall  - ratio of number of small ice crystals to large ice particles
c      ->  small ice particles are assumed to have a mean diameter of 50 microns.
c     xsimass - used for calculating small ice mixing ratio
c     xlimass - used for calculating large ice mixing ratio
c     indexs  - mean size of snow to the nearest micron (units of microns)
c     rimef   - rime factor, which is the mass ratio of total (unrimed &
c               rimed) ice mass to the unrimed ice mass (>=1)
c     flimass - mass fraction of large ice
c     qtice   - time-averaged mixing ratio of total ice
c     qlice   - time-averaged mixing ratio of large ice
c     nlice   - time-averaged number concentration of large ice
c
            if (tc.ge.0. .or. wvqw.lt.qsigrd) then
               flarge=1.
            else
               flarge=.2
               if (tc.ge.-8. .and. tc.le.-3.) flarge=.5*flarge
            endif
            fsmall=(1.-flarge)/flarge
            xsimass=rrho*massi(mdimin)*fsmall
            dum=xmimax*exp(.0536*tc)
            indexs=min(mdimax, max(mdimin, int(dum) ) )
            rimef=amax1(1., fs1d )
            xlimass=rrho*rimef*massi(indexs)
            flimass=xlimass/(xlimass+xsimass)
            qlice=flimass*qice
            nlice=qlice/xlimass
            if (nlice.lt.nlimin .or. nlice.gt.nlimax) then
c
c            Force nlice to be between nlimin and nlimax
c
               dum=max(nlimin, min(nlimax, nlice) )
               xli=rho*(qice/dum-xsimass)/rimef
               if (xli .le. massi(mdimin) ) then
                  indexs=mdimin
               else if (xli .le. massi(450) ) then
                  dli=9.5885e5*xli**.42066         ! dli in microns
                  indexs=min(mdimax, max(mdimin, int(dli) ) )
               else if (xli .le. massi(mdimax) ) then
                  dli=3.9751e6*xli**.49870         ! dli in microns
                  indexs=min(mdimax, max(mdimin, int(dli) ) )
               else 
                  indexs=mdimax
c
c               8/22/01: increase density of large ice if maximum limits
c               are reached for number concentration (nlimax) and mean size
c               (mdimax).  done to increase fall out of ice.
c
                  if (dum .ge. nlimax)
     &              rimef=rho*(qice/nlimax-xsimass)/massi(indexs)
               endif             ! end if (xli .le. massi(mdimin) )
               xlimass=rrho*rimef*massi(indexs)
               flimass=xlimass/(xlimass+xsimass)
               qlice=flimass*qice
               nlice=qlice/xlimass
            endif               ! end if (nlice.lt.nlimin ...
            qs1=amin1(qi1, qlice)
            qi1=amax1(0., qi1-qs1)
c
c         Equation (c.8) in ferrier (1994, jas, p. 272), which when
c         converted from cgs units to mks units results in the same
c         value for cice, which is equal to the {} term below:
c
c         zi={.224*720*(10**18)/[(pi*rhol)**2]}*(rho*qlice)**2/nlice,
c         where rhol=1000 kg/m**3 is the density of liquid water
c
c         Valid only for exponential ice distributions
c
            ztot=ztot+cice*rho*rho*qlice*qlice/nlice 
         endif                 ! end if (qi1 .gt. 0.) then
c
c      Calculate total (convective + grid-scale) radar reflectivity
c
 1010    if (ztot .gt. zmin) dbz1=10.*alog10(ztot)
c
         qcw(i,j,k)=qw1*1000.   ! g/kg
         qci(i,j,k)=qi1*1000.
         qra(i,j,k)=qr1*1000.
         qsn(i,j,k)=qs1*1000.
         scr3(i,j,k)=dbz1
c
      enddo
      enddo
      enddo
c
      vardesc='Cloud water mixing ratio, g/kg'
      plchun='g kg~S~-1~N~'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(qcw,'qcw       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(qcw,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'qcw       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
      vardesc='Rain water mixing ratio, g/kg'
      plchun='g kg~S~-1~N~'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(qra,'qra       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(qra,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'qra       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
      vardesc='Cloud ice mixing ratio, g/kg'
      plchun='g kg~S~-1~N~'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(qci,'qci       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(qci,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'qci       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
      vardesc='Snow mixing ratio, g/kg'
      plchun='g kg~S~-1~N~'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(qsn,'qsn       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(qsn,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'qsn       ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
      vardesc='Ferrier equiv. refl. factor, dBZ'
      plchun='dBZ'
      icd=3
      ndim=3
      if (iinterp.eq.0) then
         call writefile_rdp(scr3,'dbz_fer   ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
      else
         call hinterp(scr3,miyec,mjxec,arref,miyef,mjxef,
     &      xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
         call writefile_rdp(arrib,'dbz_fer   ',ndim,icd,vardesc,
     &      plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
      endif
c
      elseif (iferrier.eq.0) then  !   ******** 2nd half of "iferrier" if-block
c
c   If available, get cloud water mixing ratio (QCLOUD), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QCLOUD', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QCLOUD.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            qcw(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
            if (qcw(i,j,k).lt.qdelta) qcw(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qcw,'qcw       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qcw,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qcw       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get rain water mixing ratio (QRAIN), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QRAIN', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QRAIN.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(087.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            qra(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
            if (qra(i,j,k).lt.qdelta) qra(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Rain water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qra,'qra       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qra,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qra       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get rain water number conc. (QNRAIN), do not convert
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
         do j=1,mjxec
         do i=1,miyec
            qnr(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
            if (qnr(i,j,k).lt.qdelta) qnr(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Rain water number concentration, /kg'
         plchun='kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qnr,'qnr       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qra,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qnr       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get cloud ice mixing ratio (QICE), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QICE', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QICE.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(088.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            qci(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
            if (qci(i,j,k).lt.qdelta) qci(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud ice mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qci,'qci       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qci,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qci       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get cloud ice number conc. (QNICE), do not convert,
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
         do j=1,mjxec
         do i=1,miyec
            qni(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
            if (qni(i,j,k).lt.qdelta) qni(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Cloud ice number concentration, /kg'
         plchun='kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qni,'qni       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qci,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qni       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get snow mixing ratio (QSNOW), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QSNOW', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QSNOW.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(089.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            qsn(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
            if (qsn(i,j,k).lt.qdelta) qsn(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Snow mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(qsn,'qsn       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(qsn,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qsn       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
c   If available, get graupel mixing ratio (QGRAUP), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QGRAUP', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QGRAUP.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(090.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
            if (scr3(i,j,k).lt.qdelta) scr3(i,j,k)=0.
         enddo
         enddo
         enddo
         vardesc='Graupel mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         icd=3
         ndim=3
         if (iinterp.eq.0) then
            call writefile_rdp(scr3,'qgr       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(scr3,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,'qgr       ',ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      endif
c
      endif  !   ******** end of "iferrier" if-block
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
      if (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(3).eq.dimid_bt.and.vardimids(2).eq.dimid_sn.and.
     &    vardimids(1).eq.dimid_we) then
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
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, ivar,
     &   'units', nf_att_text)
      call handle_err(094.,nf_status)
      plchun=nf_att_text
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, ivar,
     &   'description', nf_att_text)
      call handle_err(095.,nf_status)
      nf_status = nf_inq_attlen (ncid, ivar,
     &   'description', nf_att_len)
      call handle_err(096.,nf_status)
      if (nf_att_len .GT. LEN(nf_att_text)) then
         print*, 'NOT POSSIBLE TO CONTINUE'
         print*, ' MEMORY OVERWITE WILL RESULT BECAUSE'
         print*, ' nf_att_len is greater than declared'
         print*, ' size of nf_att_text.  Increase to a'
         print*, ' minimum value of ', nf_att_len
         STOP 'ABORT, UNABLE TO CONTINUE'
      endif

C..Added by G. Thompson to work around problem of potentially
C.. exceeding declared length of vardesc due to a very long
C.. description in the wrfout file, variable description.
      n1 = LEN(plchun)
      do n = n1, 1, -1
         ich = ichar(plchun(n:n))
         if (.not. ( (ich.ge.65 .and. ich.le.90) .or.     ! Letters A-Z
     &               (ich.ge.97 .and. ich.le.122) .or.    ! Letters a-z
     &               (ich.ge.45 .and. ich.le.58) .or.     ! digits 0-9, also [-./]
     &               (ich.eq.95) .or. ich.eq.32) ) then   ! underscore, blank
            iendn1 = n
            goto 44
         endif
      enddo
 44   continue
      if (iendn1.le.0) then
         iendn1 = n1
      endif
      n2 = LEN(nf_att_text)
      do n = n2, 1, -1
         ich = ichar(nf_att_text(n:n))
         if (.not. ( (ich.ge.65 .and. ich.le.90) .or.     ! Letters A-Z
     &               (ich.ge.97 .and. ich.le.122) .or.    ! Letters a-z
     &               (ich.ge.45 .and. ich.le.58) .or.     ! digits 0-9, also [-./]
     &               (ich.eq.95) .or. ich.eq.32) ) then   ! underscore, blank
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
      if (iendn1+iendn2 .lt. LEN(vardesc)) then
         iendn3 = nf_att_len
         iendt = iendn1+iendn2
      else
         iendn3 = MAX(0, nf_att_len - iendn1)
         iendt = LEN(vardesc)
      endif

      vardesc(1:iendt)=nf_att_text(1:iendn3)//', '//plchun(1:iendn1)

      if (itype.eq.1) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(097.,nf_status)
         do k=1,mkzh
         do j=1,mjxec
         do i=1,miyec
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
         ndim=3
         icd=3
         if (iinterp.eq.0) then
            call writefile_rdp(scr3,varname,ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(scr3,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,varname,ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      elseif (itype.eq.2) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(098.,nf_status)
         do j=1,mjxec
         do i=1,miyec
            scr2(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         ndim=2
         icd=3
         if (iinterp.eq.0) then
            call writefile_rdp(scr2,varname,ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
         else
            call hinterp(scr2,miyec,mjxec,arref,miyef,mjxef,
     &         xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
            call writefile_rdp(arrib,varname,ndim,icd,vardesc,
     &         plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
         endif
      elseif (itype.eq.3) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3s, nf_tarr3)
         call handle_err(099.,nf_status)
         isp=index(varname,' ')
         if (isp.eq.0.or.isp.eq.10) isp=9
         do k=1,nsoil
            do j=1,mjxec
            do i=1,miyec
               scr2(i,j)=nf_tarr3(j,i,k)
            enddo
            enddo
            write(varname(isp:isp+1),'(i2.2)') k
            vardesc=nf_att_text(1:nf_att_len)//
     &         ', layer '//varname(isp:isp+1)//', '//plchun
            ndim=2
            icd=3
            if (iinterp.eq.0) then
               call writefile_rdp(scr2,varname,ndim,icd,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miyec,mjxec,mkzh)
            else
               call hinterp(scr2,miyec,mjxec,arref,miyef,mjxef,
     &            xsrcc,ysrcc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
               call writefile_rdp(arrib,varname,ndim,icd,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miyib,mjxib,mkzh)
            endif
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
c     print*,' We''re outta here like Vladimir !! '
c     print*,'===================================='
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
         write(*,*)  'NetCDF error in ripdp_wrfnmm.  Marker = ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop
      endif
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine utodot(slab,maxiy,maxjx)
c
c   This routine converts data that is on the C-grid u-velocity
c   staggered grid to the B-grid velocity staggered grid (known in
c   MM5 lingo as "dot points") 
c
      dimension slab(maxiy,maxjx),bot(2000)
c
c   Extrapolate out to top and bottom edges.
c
      do j=1,maxjx
         bot(j)=(3.*slab(1,j)-slab(2,j))/2.
         slab(maxiy,j)=(3.*slab(maxiy-1,j)-slab(maxiy-2,j))/2.
      enddo
c
c   Interpolate in the interior.
c
      do j=maxjx,1,-1
      do i=maxiy-1,2,-1
         slab(i,j)=.5*(slab(i-1,j)+slab(i,j))
      enddo
      enddo
c
c   Put "bot" values into slab.
c
      do j=1,maxjx
         slab(1,j)=bot(j)
      enddo
c
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
c                                                                     c
c*********************************************************************c
c                                                                     c
      SUBROUTINE ETALL(IM,JM,TPH0D_in,TLM0D_in,DLMD_in,DPHD_in,
     & HLAT,HLON,VLAT,VLON)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM: ETALL         COMPUTE EARTH LATITUDE & LONIGTUDE OF
C                           ETA GRID POINTS
C   PRGMMR: ROGERS          ORG: W/NP22     DATE: 90-06-13
C
C ABSTRACT: COMPUTES THE EARTH LATITUDE AND LONGITUDE OF ETA GRID
C   POINTS (BOTH H AND V POINTS)
C
C PROGRAM HISTORY LOG:
C   90-06-13  E.ROGERS
C   98-06-09  M.BALDWIN - CONVERT TO 2-D CODE
C   01-01-03  T BLACK   - MODIFIED FOR MPI
C
C USAGE:    CALL ETALL(HLAT,HLON,VLAT,VLON)
C   INPUT ARGUMENT LIST:
c     IM       - NUMBER OF X POINTS (IN USUAL E-GRID SPECIFICATION)
c     JM       - NUMBER OF Y POINTS (IN USUAL E-GRID SPECIFICATION)
c     TPH0D_in - LATITUDE OF HIGH POINT OF ROTATED EQUATOR
c     TLM0D_in - LONGITUDE OF HIGH POINT OF ROTATED EQUATOR
c     DLMD_in  - LONGITUDINAL GRID SPACING IN DEGREES
c     DPHD_in  - LATITUDINAL GRID SPACING IN DEGREES
C
C   OUTPUT ARGUMENT LIST:
C     HLAT     - LATITUDE OF H GRID POINTS IN RADIANS (NEG=S)
C     HLON     - LONGITUDE OF H GRID POINTS IN RADIANS (E)
C     VLAT     - LATITUDE OF V GRID POINTS IN RADIANS (NEG=S)
C     VLON     - LONGITUDE OF V GRID POINTS IN RADIANS (E)
C
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM RS/6000 SP
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
       INTEGER, PARAMETER :: KNUM=SELECTED_REAL_KIND(13)
      REAL(KIND=KNUM) :: ONE=1.,R180=180.,TWO=2.
      REAL(KIND=KNUM) :: CLMH,CLMV,CTPH,CTPH0,CTPV,D2R,DLM,DLMD,DPH,DPHD     &
     &                  ,FACTR,PI,R2D,SB,SBD,SPHH,SPHV,STPH,STPH0,STPV   &
     &                  ,TDLM,TLMH,TLMV,TPH0,TPHH,TPHV,WB,WBD
      REAL(KIND=KNUM),DIMENSION(IM,JM) :: GLATH,GLATV,GLONH,GLONV
      REAL,DIMENSION(IM,JM) :: HLAT,HLON,VLAT,VLON

C-----------------------------------------------------------------------
C--------------DERIVED GEOMETRICAL CONSTANTS----------------------------
C----------------------------------------------------------------------

	TPH0D=TPH0D_in
	TLM0D=TLM0D_in
	DPHD=DPHD_in
	DLMD=DLMD_in

      WBD=-(IM-ONE)*DLMD
      SBD=-(JM-1)/2*DPHD
      PI=ACOS(-ONE)
      DTR = PI / R180
      TPH0 = TPH0D * DTR
      WB = WBD * DTR
      SB = SBD * DTR
      DLM = DLMD * DTR
      DPH = DPHD * DTR
      TDLM = DLM + DLM
      TDPH = DPH + DPH
C
      STPH0 = SIN(TPH0)
      CTPH0 = COS(TPH0)
C
C-----------------------------------------------------------------------
C---COMPUTE GEOGRAPHIC LAT AND LONG OF ETA GRID POINTS (H & V POINTS)---
C-----------------------------------------------------------------------
      DO 200 J = 1,JM
C
         TLMH = WB - TDLM + MOD(J+1,2) * DLM
         TPHH = SB+(J-1)*DPH
         TLMV = WB - TDLM + MOD(J,2) * DLM
         TPHV = TPHH
         STPH = SIN(TPHH)
         CTPH = COS(TPHH)
         STPV = SIN(TPHV)
         CTPV = COS(TPHV)
C----------------------------------------------------------------------
C---------- COMPUTE EARTH LATITUDE/LONGITUDE OF H POINTS --------------
C----------------------------------------------------------------------
         DO 201 I = 1,IM
           TLMH = TLMH + TDLM
           SPHH = CTPH0 * STPH + STPH0 * CTPH * COS(TLMH)
           GLATH(I,J) = ASIN(SPHH)
           CLMH = CTPH * COS(TLMH) / (COS(GLATH(I,J)) * CTPH0)
     1               - TAN(GLATH(I,J)) * TAN(TPH0)
           IF(CLMH .GT. ONE) CLMH = ONE
           IF(CLMH .LT. -ONE) CLMH = -ONE
           FACTH = ONE
           IF(TLMH .GT. 0.) FACTH = -ONE
           GLONH(I,J) = -TLM0D * DTR + FACTH * ACOS(CLMH)

           HLAT(I,J) = GLATH(I,J) / DTR
	   HLON(I,J)= -GLONH(I,J)/DTR
           IF(HLON(I,J) .GT. 180.) HLON(I,J) = HLON(I,J) - 360.
           IF(HLON(I,J) .LT. -180.) HLON(I,J) = HLON(I,J) + 360.
  201    CONTINUE
C----------------------------------------------------------------------
C---------- COMPUTE EARTH LATITUDE/LONGITUDE OF V POINTS --------------
C----------------------------------------------------------------------
         DO 202 I = 1,IM
           TLMV = TLMV + TDLM
           SPHV = CTPH0 * STPV + STPH0 * CTPV * COS(TLMV)
           GLATV(I,J) = ASIN(SPHV)
           CLMV = CTPV * COS(TLMV) / (COS(GLATV(I,J)) * CTPH0)
     1          - TAN(GLATV(I,J)) * TAN(TPH0)
           IF(CLMV .GT. 1.) CLMV = 1.
           IF(CLMV .LT. -1.) CLMV = -1.
           FACTV = 1.
           IF(TLMV .GT. 0.) FACTV = -1.
           GLONV(I,J) = -TLM0D * DTR + FACTV * ACOS(CLMV)
C
C    CONVERT INTO DEGREES AND EAST LONGITUDE
C
           VLAT(I,J) = GLATV(I,J) / DTR
           VLON(I,J) = -GLONV(I,J) / DTR
           IF(VLON(I,J) .GT. 180.) VLON(I,J) = VLON(I,J) - 360.
           IF(VLON(I,J) .LT. -180.) VLON(I,J) = VLON(I,J) + 360.

!	if (mod(I,10) .eq. 0 .and. mod(J,10) .eq. 0) then
!	write(6,*) 'I,J,HLAT,HLON,VLAT,VLON: ', I,J,HLAT(I,J),HLON(I,J)
!     &                                             ,VLAT(I,J),VLON(I,J)
!	endif
  202    CONTINUE
  200 CONTINUE
      RETURN
      END
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hinterp(arrec,miyec,mjxec,arref,miyef,mjxef,
     &   xsrc,ysrc,arrib,miyib,mjxib,mkzh,icd,ndim,rmsg)
c
      dimension arrec(miyec,mjxec,1+(mkzh-1)*(ndim-2))
      dimension arref(miyef,mjxef,1+(mkzh-1)*(ndim-2))
      dimension arrib(miyib,mjxib,1+(mkzh-1)*(ndim-2))
      dimension xsrc(miyib,mjxib),ysrc(miyib,mjxib)
c
      nlevs=1+(mkzh-1)*(ndim-2)
c
c     Copy the E-grid data array (arrec) to an array that contains
c     all E-grid points (staggered and unstaggered) (arreqexp).  Note,
c     data is still "scrunched" in the x direction.
c
      do k=1,nlevs
      do j=1,mjxec
      do i=1,miyec
         arref(i,j,k)=arrec(i,j,k)
      enddo
      enddo
      enddo
c
c     Now call egridfill, which will expand the compact data to the
c     proper points in the arref array, and will then fill in the
c     holes, resulting in an array that represents square grid boxes (at
c     all H and V points) instead of diamond-shaped
c
      call egridfill(arref,icd,miyef,mjxef,miyef,mjxef,
     &   nlevs)
c
c   Change icd.  Assume that if source array was on E-grid H (V)
c   points, data should be interpolated to cross (dot) points of
c   interpolation grid, and xsrc/ysrc are provided on the cross
c   (dot) point grid of the interpolation domain.  In either case,
c   xsrc/ysrc provide the grid location within the E-grid, defined
c   as lower-left H point being x=1,y=1 and upper-right H point
c   being x=mjxec,y=miy.
c
c   Note, this will change icd on output from this routine.
c
      icd=icd-2
c
c   Loop through interp domain points
c
      do jj=1,mjxib-icd
      do ii=1,miyib-icd
         jl=nint(xsrc(ii,jj)-.5)
         jr=jl+1
         ib=nint(ysrc(ii,jj)-.5)
         it=ib+1
         ratlr=xsrc(ii,jj)-jl
         ratbt=ysrc(ii,jj)-ib
         if (jl.lt.1.or.jl.gt.mjxef.or.jr.lt.1.or.jr.gt.mjxef.or.
     &       ib.lt.1.or.ib.gt.miyef.or.it.lt.1.or.it.gt.miyef)
     &       then
            do k=1,nlevs
               arrib(ii,jj,k)=rmsg
            enddo
         else
            fac1=(1.-ratlr)*(   ratbt)
            fac2=(   ratlr)*(   ratbt)
            fac3=(1.-ratlr)*(1.-ratbt)
            fac4=(   ratlr)*(1.-ratbt)
            do k=1,nlevs
               wk1=arref(it,jl,k)
               wk2=arref(it,jr,k)
               wk3=arref(ib,jl,k)
               wk4=arref(ib,jr,k)
               arrib(ii,jj,k)=
     &            fac1*wk1+fac2*wk2+fac3*wk3+fac4*wk4
            enddo
         endif
      enddo
      enddo
c
      return
      end
