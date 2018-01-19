      program rip
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIP is a post-processing program for the PSU/NCAR mesoscale       c
c   model, written by Mark Stoelinga at the University of             c
c   Washington and NCAR.  The name "RIP" stands for                   c
c   Read/Interpolate/Plot. The program is intended for use with       c
c   mesoscale model output, and has the flexibility to handle data    c
c   in any vertical coordinate system.  The data must first be        c
c   converted to the appropriate format expected by this              c
c   program, using the program RIPDP, which stands for RIP data       c
c   preparation.                                                      c
c                                                                     c
c   RIP Version 4.6 released March 2010.                              c
c   RIP Version 4.5 released March 2009.                              c
c   RIP Version 4.4 released May 2008.                                c
c   RIP Version 4.3 released June 2007.                               c
c   RIP Version 4.2 released December 2006.                           c
c   RIP Version 4.1 released Apr 2005.                                c
c   RIP Version 4.0 released Apr 2003.                                c
c   RIP Version 3 released Oct 2000.                                  c
c   RIP Version 2 released Dec 1997.                                  c
c   RIP was first coded informally by Mark Stoelinga in 1991.         c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The main program does very little "work".  It simply reads
c   a data header record in order to give dimensions to
c   subroutine driver, which does most of the "work" in RIP.
c   Some of the arguments to subroutine driver are used as
c   dimensions of large arrays, which are local to
c   subroutine driver.  The purpose of this set up is
c   for dynamic memory allocation.  This capability is
c   not standard to Fortran 77, but is allowable on some f77
c   compilers (such as Cray).  It is standard to Fortran 90.
c   Therefore, RIP may be compiled with f77 compilers that allow
c   adjustable dimensioning of local (non-argument) arrays in
c   subroutines, or with any f90 compiler.
c
c   The following parameters should not need to be changed for most
c   applications:
c
c      maxtraj: the maximum number of trajectories that can be
c         calculated in trajectory calculation mode
c      maxtavl: the maximum number of available times that can be
c         stored in the available times arrays (xtimeavl and cxtimeavl)
c
      parameter (maxtraj=7000,maxtavl=200)
c
      dimension xtimeavl(maxtavl)
c
      character argum(16)*256,fname*256,cxtimeavl(maxtavl)*10,
     &   casename*256
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes)
      character title*80,rootname*256,titlecolor*40,rip_root*256,
     &   ptimeunits*1, ncarg_type*10
      namelist/userin/ title,rip_root,rootname,flmin,frmax,fbmin,
     &   ftmax,ptimes,iptimes,ptimeunits,tacc,mdatebf,
     &   ntextq,ntextcd,ntextfn,idotser,idotitle,timezone,
     &   iusdaylightrule,inearesth,iinittime,ifcsttime,ivalidtime,
     &   noplots,itrajcalc,fcoffset,titlecolor,idescriptive,
     &   icgmsplit,maxfld,imakev5d,inewdom, istopmiss, ncarg_type
      dimension xjtraj(maxtraj),yitraj(maxtraj),zktraj(maxtraj),
     &   diag(maxtraj)
      character vctraj*1
      namelist/trajcalc/ rtim,ctim,dtfile,dttraj,
     &   vctraj,ihydrometeor,xjtraj,yitraj,zktraj
c
      dimension ptuse(maxptimes)
c
      include 'comconst'
c
c   Get command line arguments.
c
      nargum=iargc()
      do i=1,nargum
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:4).eq.'rip_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      iup=6         ! output unit# for standard rip printout
      if (nargum.eq.3.and.argum(1)(1:2).eq.'-f') then
         iup=10
         argum(1)=argum(2)
         argum(2)=argum(3)
      elseif (nargum.ne.2) then
         write(6,*) 'Usage:'
         write(6,*) '   rip [-f] model_case_name rip_case_name'
         write(6,*) 'Options:'
         write(6,*) '   -f: standard rip printout should go to a file'
         write(6,*) '       named rip_case_name.out rather than to'
         write(6,*) '       the screen.'
         write(6,*) 'Note: "rip_case_name" should be the root part of'
         write(6,*) 'the name of the ".in" file.  Do not inlcude ".in"'
         write(6,*) 'at the end of "rip_case_name" on the'
         write(6,*) 'rip command line.'
         stop
      endif
      iupt=iup
c
c   Get model case name from argument
c
      casename=argum(1)
      iendc=index(casename,' ')-1
      if (iendc.eq.-1) iendc=256
c
c   Get root name from argument
c
      rootname=argum(2)
      iendcr=index(rootname,' ')-1
      if (iendcr.eq.1) then
         write(iup,*) 'rip case name all blanks.'
         stop
      elseif (iendcr.gt.236) then
         write(iup,*) 'rip case name must be < or = 236 characters.'
         stop
      endif
      if (rootname(iendcr-2:iendcr).eq.'.in') then
         rootname(iendcr-2:iendcr)='   '
         iendcr=iendcr-3
      endif
c
c  make sure ncarg_root is set
c
      rip_root = 'ncarg_root_should_be_set'
      call getenv('NCARG_ROOT',rip_root)
      if (rip_root(1:10).eq.'          ') then
         write(iup,*)'The NCARG_ROOT environment variable does not seem 
     &to be set'
         stop 'ncarg_root'
      endif

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
c   Define unit numbers
c     (Note, unit numbers 25 through 30 are used for include files
c      in the plot specification table.  See readspec.f.)
c
      iuinput=7     ! input unit# for user input (namelists & plspecs)
      iustnlist=8   ! input unit# for the station list
      iutserstn=70  ! input unit# for the time series stations
      iutserdat=71  ! output unit# for the time series printout
      iunewdom=72   ! input unit# for new domain parameters
      iutrajin=74   ! input unit# for trajc positns (traj plot mode)
      iutrajout=75  ! output unit# for trajc positns (traj calc mode)
      iudiagout=76  ! output unit# for trajc diagnstcs (traj calc mode)
      iuv5dout=79   ! output unit# for vis5d files (make vis5d mode)
c
c   Open the print out file
c
      if (iup.eq.10) then
         fname=rootname(1:iendcr)//'.out'
         open (unit=iup,file=fname,form='formatted',status='unknown')
      endif
c
c   Read the namelist values.
c
      fname=rootname(1:iendcr)//'.in'
      open (unit=iuinput,file=fname,form='formatted',status='old')
      fname=rootname   ! save rootname temporarily
c
      title='auto'
      rip_root='/dev/null'
      ncarg_type='cgm'
      rootname=' '      ! no longer set in namelist
      mdatebf=99999999  ! no longer used
      flmin=.05
      frmax=.95
      fbmin=.10
      ftmax=.90
      do i=1,maxptimes
         ptimes(i)=9e9
         iptimes(i)=99999999
      enddo
      ptimeunits='h'
      tacc=1.
      ntextq=0
      ntextcd=0
      ntextfn=0
      idotser=0
      noplots=0
      idotitle=1
      timezone=-7.
      iusdaylightrule=1
      inearesth=0
      iinittime=1
      ifcsttime=1
      ivalidtime=1
      idescriptive=1
      itrajcalc=0
      imakev5d=0
      inewdom=0
      fcoffset=0.
      titlecolor=' '
      icgmsplit=0
      maxfld=10
      istopmiss=1
c
      read (iuinput,userin)
      rootname=fname

      if ( trim(ncarg_type) == 'x11') icgmsplit=0
c
c   RIP has a capability to define its own domain and gridded variables,
c   in order to do specialized tasks completely independently of any model
c   input/output files.  This capability is turned on with inewdom=1 in
c   the input namelist.  If inewdom=1, then here we open a file that
c   defines certain parameters of the new domain.
c
      if (inewdom.eq.1) then
         fname=casename(1:iendc)//'.newdom'
         open(unit=iunewdom,file=fname,form='formatted',
     &      status='old')
      endif
c
      if (itrajcalc.eq.1.and.imakev5d.eq.1) then
         write(iup,*) 'You cannot be in trajectory calculation mode',
     &          ' (itrajcalc=1)'
         write(iup,*) 'and "make vis5d data" mode (imakev5d=1) at the',
     &          ' same time.'
         write(iup,*) 'You must choose one or the other.'
      endif
c
      rtim=0.0
      ctim=0.0
      dtfile=3600.0   ! seconds
      dttraj=600.0    ! seconds
      vctraj='s'
      ihydrometeor=0
      do i=1,maxtraj
         xjtraj(i)=rmsg
         yitraj(i)=rmsg
         zktraj(i)=rmsg
      enddo
c
      if (itrajcalc.eq.1) then
         read (iuinput,trajcalc)
c
c      Check if a grid is defined
c
         if (xjtraj(1).lt.0.0) call mktrjpts(xjtraj,yitraj,zktraj,
     &                                        maxtraj,rmsg)
         do i=maxtraj,1,-1
            if (xjtraj(i).ne.rmsg) then
               ntraj=i
               goto 40
            endif
         enddo
 40      continue
      else
         ntraj=1
      endif
c
c   Get rip_root from environment variable RIP_ROOT
c   if rip_root='/dev/null '
c
      if (rip_root(1:10).eq.'/dev/null ') then
         call getenv('RIP_ROOT',rip_root)
      endif
      if (rip_root(1:10).eq.'          ') then
         write(iup,*)'Either RIP_ROOT environment variable or rip_root'
         write(iup,*)'is not properly set.'
         stop
      endif
c
c   Set up lookup table for getting temperature on a pseudoadiabat.
c   (Borrow the unit number for the stationlist, just for the moment.)
c
      iendrr=index(rip_root,' ')-1
      fname=rip_root(1:iendrr)//'/psadilookup.dat'
      open (unit=iustnlist,file=fname,form='formatted',status='old')
      do i=1,14
         read(iustnlist,*)
      enddo
      read(iustnlist,*) nthte,nprs
      if (nthte.ne.150.or.nprs.ne.150) then
         write(iup,*)
     &      'Number of pressure or theta_e levels in lookup table'
         write(iup,*) 'file not = 150.  Check lookup table file.'
         stop
      endif
      read(iustnlist,173) (psadithte(jt),jt=1,nthte)
      read(iustnlist,173) (psadiprs(ip),ip=1,nprs)
      read(iustnlist,173) ((psaditmk(ip,jt),ip=1,nprs),jt=1,nthte)
 173  format(5e15.7)
      close(iustnlist)
c
c   Obtain the available xtimes, in both character and
c   floating point arrays
c
      if (inewdom.eq.0) then
         call gettimes(casename,iendc,xtimeavl,cxtimeavl,nxtavl,
     &      maxtavl,ncxc,iup)
      else
         nxtavl=1
         ncxc=10
         read(iunewdom,*) cxtimeavl(1)
         read(cxtimeavl(1),'(f10.5)') xtimeavl(1)
         do i=nxtavl+1,maxtavl
            cxtimeavl(i)=' '
            xtimeavl(i)=0.
         enddo
      endif
c
c   Get grid dimensions
c
      if (inewdom.eq.0) then
         call getdims(casename,iendc,xtimeavl,cxtimeavl,
     &      maxtavl,ncxc,miy,mjx,mkzh,mabpl,morpl,iup)
      else
         read(iunewdom,*) miy,mjx,mkzh
         mabpl=2*max(miy,mjx)+1
         morpl=max(max(miy,mjx),mkzh)
      endif
c
      rewind (iuinput)
c
      call driver(miy,mjx,mkzh,mabpl,morpl,xtimeavl,
     &     cxtimeavl,ncxc,maxtavl,nxtavl,casename,iendc,
     &     title,rip_root,rootname,iendcr,ptimes,iptimes,
     &     ptuse,maxptimes,ptimeunits,tacc,ntextq,ntextcd,ntextfn,
     &     idotser,noplots,idotitle,timezone,iusdaylightrule,
     &     inearesth,iinittime,ifcsttime,ivalidtime,fcoffset,
     &     titlecolor,idescriptive,icgmsplit,maxfld,
     &     itrajcalc,rtim,ctim,dtfile,dttraj,
     &     vctraj,ihydrometeor,xjtraj,yitraj,zktraj,diag,
     &     ntraj,imakev5d,inewdom,istopmiss,ncarg_type)
c
      stop
      end
