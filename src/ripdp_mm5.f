      program ripdp_mm5
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIPDP_MM5 is a data preparation program that reads in data from   c
c   the PSU/NCAR mesoscale modeling system and makes data files       c
c   appropriate for input to the RIP data analysis and visulaization  c
c   program.  RIP stands for Read/Interpolate/Plot and RIPDP stands   c
c   for RIP data preparation.                                         c
c                                                                     c
c   FYI: Input unit numbers are 21,22,23,...                          c
c        Output unit number is 65.                                    c
c                                                                     c
c   Originally written by Mark Stoelinga, Univ. of Washington         c
c                                                                     c
c   Modified Feb 2003 for new, generalized vertical coordinate        c
c   version of RIP.  See code commented with "ccc".                   c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The main program does very little "work".  It simply reads
c   the model data header record in order to give information to
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
c   Model output header variables.
c
      integer   jyr(18),jmo(18),jdy(18),jhr(18)
      dimension sigf(200), plev(200), idumb(100), mif(30)
      real      mrf(10)
      logical   mlf(10)
      integer   mifv1(1000,20)
      real      mrfv1(1000,20)
      character*80 mifc(1000,20), mrfc(1000,20)
      integer bhi(50,20), flag
      real bhr(20,20)
      character*80 bhic(50,20),bhrc(20,20)
      integer ndim
      real time
      integer start_index(4), end_index(4)
      character staggering*4, ordering*4,
     &    current_date*24, name*9, units*25, description*46
c
      character dataform*8, argum(256)*256
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
      if (argum(1)(1:10).eq.'ripdp_mm5_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if ((argum(1)(1:3).eq.'-n '.and.nargum.lt.4).or.
     &    (argum(1)(1:3).ne.'-n '.and.nargum.lt.2)) then
         print*,'Usage:'
         print*,'  ripdp [-n namelist_file] casename',
     &          ' data_file_1 data_file_2 data_file_3 ...'
         print*
         print*,'Note: "namelist_file" can be any name, either with'
         print*,'an extension of your choosing or without an extension.'
         print*,'The ".in" extension is not assumed, as it is in rip.'
         stop
      endif
c
      if (argum(1)(1:3).eq.'-n ') then
         nnl=2
         ncn=3
         nsets=nargum-3
         nsetsbeg=4
      else
         nnl=0
         ncn=1
         nsets=nargum-1
         nsetsbeg=2
      endif
c
c   Read header record to get dimensions
c
      iudatin=21
      iexpanded=0
      open (unit=iudatin,file=argum(nsetsbeg),form='unformatted',
     &   status='old')
c
      read(iudatin,err=110,end=180)  mdateb,iprog,icoord,
     &              ksigf,(sigf(k),k=1,ksigf),
     &              mif,mrf,mlf,jyr,jmo,jdy,jhr,
     &              ioldprog,ioldcoord,ioldnlv,(plev(k),k=1,ioldnlv),
     &              ibltyp,  isfflx,  itgflg,  icdcon,
     &              isfpar,  ivmixm,  idry,    imoist,  icloud,
     &              iboudy,  imoiav,  ifsnow,  icustb,  itqpbl,
     &              ifrad
      if (mdateb.gt.99999999.or.mdateb.lt.999999.or.
     &    (icloud.ne.0.and.icloud.ne.1).or.
     &    (isfpar.ne.0.and.isfpar.ne.1)) goto 110
      dataform='mm4v8out'
      if (mlf(3)) then
         miy=mif(4)
         mjx=mif(5)
      else
         miy=mif(2)
         mjx=mif(3)
      endif
      mkzh=ksigf-1
      goto 200
 110  rewind (iudatin)
c
      read(iudatin,err=120,end=180)  mdateb,iprog,icoord,
     &   ksigf,(sigf(k),k=1,ksigf),mif,mrf,mlf,jyr,jmo,jdy,jhr,
     &   xtime,ndumb,(idumb(l), l = 1, ndumb),ioldprog,ioldcoord,
     &   ioldnlv,(plev(k),k=1,ioldnlv)
      if (mdateb.gt.99999999.or.mdateb.lt.999999) goto 120
      dataform='mm5v0out'
      if (mlf(3)) then
         miy=mif(4)
         mjx=mif(5)
      else
         miy=mif(2)
         mjx=mif(3)
      endif
      mkzh=ksigf-1
      goto 200
 120  rewind (iudatin)
c
      read(iudatin,err=130,end=180) mifv1,mrfv1,mifc,mrfc
      if (mifv1(1,1).eq.0) then
         print*,'In attempting to read the data as MM5V1, no error'
         print*,'occurred, but the first element of the integer'
         print*,'header array is 0.  RIPDP suspects the data is MM5V3,'
         print*,'and will proceed.'
         print*
         goto 130
      endif
      dataform='mm5v1   '
      iprog=mifv1(1,1)
      miy = mifv1(104,1)
      mjx = mifv1(105,1)
      if ((mifv1(103,1).eq.0) .and. (mifv1(5,1).eq.1) .and. 
     &     (iprog.le.3)) then
         miy = mifv1(6,1)
         mjx = mifv1(7,1)
      endif
      mkzh = mifv1(101,iprog)
      if ((iprog.eq.5).or.(iprog.eq.6))
     &   mkzh = nint(mrfv1(101,iprog))
      goto 200
 130  rewind (iudatin)
c
      read(iudatin,err=140,end=180) flag
      if (flag.ne.0) then
         print*,'MM5 V3 data does not begin with a big header.'
         print*,'Stopping.'
         stop
      endif
      read(iudatin) bhi, bhr, bhic, bhrc
c
c      do j=1,20
c         do i=1,50
c            if (bhi(i,j).ne.-999)
c     &         print 633, i,j,bhi(i,j),bhic(i,j)(1:50)
c         enddo
c         do i=1,20
c            if (abs(bhr(i,j)+999.).gt..01)
c     &         print 634, i,j,bhr(i,j),bhrc(i,j)(1:50)
c         enddo
c      enddo
c 633  format('bhi(',i2,',',i2,'): ',i6,' : ',a50)
c 634  format('bhr(',i2,',',i2,'): ',f9.2,' : ',a50)
c
      dataform='mm5v3   '
      iprog=bhi(1,1)
      miy = bhi(16,1)
      mjx = bhi(17,1)
      rewind (iudatin)
      if (iprog.eq.1) then
         mkzh=2  ! ccc "Fake" 3D data for TERRAIN output will have 2 vert. levs.
      elseif (iprog.eq.2.or.iprog.eq.3) then
c
c      Number of levels needs to be obtained from the header for a
c      3D array--can't trust bhi(12,iprog) to be correct.
c
         read(iudatin) flag
         if (flag.ne.0) then
            print*,'MM5 V3 data does not begin with a big header.'
            print*,'Stopping.'
            stop
         endif
         read(iudatin) bhi, bhr, bhic, bhrc
 703     read(iudatin) flag
         if (flag.eq.1) then
            read (iudatin) ndim,start_index,end_index,time,
     &         staggering,ordering,current_date,name,
     &         units,description
            if (name.ne.'U        ') then
               read(iudatin)
               goto 703
            else
               ntotlevels=end_index(3)-start_index(3)+1
               mkzh=ntotlevels-1   ! Don't include surface level in mkzh
            endif
         else
            print*,'Ran into a flag not =1 looking for U.'
            stop
         endif
         rewind (iudatin)
      elseif (iprog.eq.5.or.iprog.eq.11) then
         mkzh = bhi(12,iprog)
      else
         print*,'For MM5V3 data, can only read output from'
         print*,'Terrain, Regrid, Rawins/Little_r, Interp, or MM5.'
         stop
      endif
      if (iprog.ge.1.and.iprog.le.3) then
c
c      Need to determine if this is "expanded domain".
c      Look for 'TERRAIN' variable.
c
         read(iudatin) flag
         if (flag.ne.0) then
            print*,'MM5 V3 data does not begin with a big header.'
            print*,'Stopping.'
            stop
         endif
         read(iudatin) bhi, bhr, bhic, bhrc
 733     read(iudatin) flag
         if (flag.eq.1) then
            read (iudatin) ndim,start_index,end_index,time,
     &         staggering,ordering,current_date,name,
     &         units,description
            if (name.ne.'TERRAIN  ') then
               read(iudatin)
               goto 733
            else
               if (end_index(1).ne.miy.or.end_index(2).ne.mjx) then
c
c               Looks like this is an expanded domain.
c
                  miy=end_index(1)-start_index(1)+1
                  mjx=end_index(2)-start_index(2)+1
                  if (bhi(8,1).ne.1.or.bhi(13,1).ne.1.or.
     &                miy.ne.bhi(9,1).or.mjx.ne.bhi(10,1)) then
                     print*,'Looks like an expanded domain but ',
     &                  'some things in headers are not consistent.'
                     print*,'bhi(8,1),bhi(13,1),bhi(9,1),bhi(10,1)='
                     print*,bhi(8,1),bhi(13,1),bhi(9,1),bhi(10,1)
                     print*,'miy,mjx (dervied from "TERRAIN" ',
     &                      'variable header)='
                     print*,miy,mjx
                     stop
                  endif
                  iexpanded=1
               endif
            endif
         else
            print*,'Ran into a flag not =1 looking for TERRAIN.'
            stop
         endif
         rewind (iudatin)
      endif
      goto 200
 140  continue
c
      print*,'The model data header is not a format'//
     &   ' that RIPDP recognizes.  Stopping.'
      stop
c
 180  print*,'Unexpected EOF reached when trying to read'
      print*,'model data header.  Stopping.'
      stop
c
 200  continue
c
      close (iudatin)
      call process(miy,mjx,mkzh,dataform,iexpanded,
     &   argum,nnl,ncn,nsets,nsetsbeg,mif,mrf,mlf,mifv1,mrfv1,
     &   mifc,mrfc,bhi,bhr,bhic,bhrc)
      stop
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine process(miy,mjx,mkzh,dataform,iexpanded,
     &   argum,nnl,ncn,nsets,nsetsbeg,mif,mrf,mlf,mifv1,mrfv1,
     &   mifc,mrfc,bhi,bhr,bhic,bhrc)
c
c   This subroutine does most of the "work".
c
c   miy, and mjx are dot-point dimensions, in the y and x directions
c      respectively, of the domain to be analyzed.
c   mkzh is number of 1/2-sigma levels in the domain.
c   dataform tells what format the model data is.
c   iexpanded is a flag telling whether domain is expanded ot not
c   argum carries the names of the model data files.
c   nsets is the number of files in the model dataset.
c   nsetsbeg is the element of argum that holds the first file name
c      of the model dataset.
c   mif, mrf, mlf, mifv1, mrfv1, mifc, and mrfc are the model system
c      header variables
c
      character dataform*8, argum(256)*256
c
c   Other parameters (which shouldn't need to be changed
c      for most applications):
c
c   nvq, nvtq, and nvvq are the number of PV, temp., and u/v
c      partitions (set ipvdim to whatever ipv is in the data set).
c
      parameter (ipvdim=0, nvq=1, nvtq=5, nvvq=1)
c
ccc      dimension ter(miy,mjx),ter_tsf(miy,mjx),xmap(miy,mjx),
      real ter(miy,mjx),ter_pseudo(miy,mjx),xmap(miy,mjx),
     &   xlat(miy,mjx),xlon(miy,mjx),
     &   cor(miy,mjx),xlus(miy,mjx),sno(miy,mjx),pstx(miy,mjx),
     &   rtc(miy,mjx),rte(miy,mjx),pstd(miy,mjx),tgk(miy,mjx),
     &   dmap(miy,mjx),sst(miy,mjx),pblh(miy,mjx),
     &   regime(miy,mjx),sshflux(miy,mjx),slhflux(miy,mjx),
     &   ust(miy,mjx),swdown(miy,mjx),lwdown(miy,mjx),
     &   soil1(miy,mjx),soil2(miy,mjx),soil3(miy,mjx),soil4(miy,mjx),
     &   soil5(miy,mjx),soil6(miy,mjx),uuu_sfan(miy,mjx),
     &   vvv_sfan(miy,mjx),tmk_sfan(miy,mjx),qvp_sfan(miy,mjx),
     &   slp(miy,mjx),sigh(mkzh),sigf(mkzh+1),
     &   dotcor(miy,mjx),tmk(miy,mjx,mkzh),uuu(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),www(miy,mjx,mkzh+1),
     &   vvv(miy,mjx,mkzh),prs(miy,mjx,mkzh),sfp(miy,mjx),
     &   qvp(miy,mjx,mkzh),qcw(miy,mjx,mkzh),qra(miy,mjx,mkzh),
     &   qci(miy,mjx,mkzh),qsn(miy,mjx,mkzh),
     &   qgr(miy,mjx,mkzh),rnci(miy,mjx,mkzh),
     &   tke(miy,mjx,mkzh),radtnd(miy,mjx,mkzh),
     &   qqn(miy,1+ipvdim*(mjx-1),1+ipvdim*(mkzh-1),nvq),
     &   tqn(miy,1+ipvdim*(mjx-1),1+ipvdim*(mkzh-1),nvtq),
     &   uqn(miy,1+ipvdim*(mjx-1),1+ipvdim*(mkzh-1),nvvq),
     &   vqn(miy,1+ipvdim*(mjx-1),1+ipvdim*(mkzh-1),nvvq),
     &   scr3(miy,mjx,mkzh),scr3f(miy,mjx,mkzh+1),scr2(miy,mjx),
     &   rrbo(miy,mjx),clgo(miy,mjx),viso(miy,mjx),
     &   clgf(miy,mjx),visf(miy,mjx),
     &   xmav(miy,mjx)
c
ccccccccccccc temp tend
c      parameter (MIXttend =  142,   MJXttend =  151,   MKXttend = 37) 
c      dimension t3dten(mixttend,mjxttend,mkxttend)
c      character msg*3
ccccccccccccc temp tend
c
      character pvdfname*256,varname*10,fname*256,cxtime*10,
     &   cxtimeavl(200)*10
c      dimension tq4av(miy,mjx,mkzh)
c
c   Model output header variables.
c
      integer   jyr(18),jmo(18),jdy(18),jhr(18),plev(100),
     &   idumb(100), mif(30)
      real      mrf(10)
      logical   mlf(10)
      integer   mifv1(1000,20)
      real      mrfv1(1000,20)
      character*80 mifc(1000,20), mrfc(1000,20)
      integer bhi(50,20), flag
      real bhr(20,20)
      character*80 bhic(50,20),bhrc(20,20)
      integer ndim
      real time
      integer start_index(4), end_index(4)
      character staggering*4, ordering*4,
     &   current_date*24, name*9, units*25, description*46
      dimension prslvl(200)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c minfo variables  JFB 12/26/98
c
      character cupa(10)*10, bltyp(10)*10, mphys(10)*10, version*6
      data cupa /'No Cumulus','Anthes-Kuo',' Grell ','Ara-Schu',
     & 'Fritsch-Ch','Kain-Frsch','Betts-Mill','KF-2',' ',' '/
      data bltyp /'No frict','Bulk PBL','Blackadar','Burk-Thomp',
     & 'Eta PBL','MRF PBL',' ',' ',' ',' '/
      data mphys /' Dry','Stable','Warm rain','Simple ice','Reisner 1',
     & 'GSFC Graup','Reisner 2','Schultz',' ',' '/
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character discard(maxptimes)*9,ptimeunits*1
      namelist/userin/ ptimes,iptimes,ptimeunits,tacc,discard,
     &   iexpandedout,ipv,npvsets,iskpd1,iskppvd1,iobsprc,
     &   iobscnv,iftcnv
c
      print*,'Welcome to your friendly RIPDP (V4.6.5) output file !'    ! October 2013
c
c   Define some constants
c
      rgas=287.04  !J/K/kg
      grav=9.81           ! m/s**2
      sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
      eps=0.622
      ezero=6.112  ! hPa
      pvc=1.e6
      celkel=273.15
      eslcon1=17.67
      eslcon2=29.65
      ussalr=.0065      ! deg C per m
c
c   Define unit numbers
c
      iuinput=7     ! input unit# for namelist, color table,
c                        and plspec table
      iudatin=21    ! input unit# for the regular data.
      iupv=41       ! input unit# for the pv data.
      iuobsprc=91   ! input unit# for the observed precip data
      iuobsclg=92   ! input unit# for the observed ceiling data
      iuobsvis=93   ! input unit# for the observed visibility data
      iuftclg=94    ! input unit# for the FT ceiling data
      iuftvis=95    ! input unit# for the FT visibility data
c
      ifilecount=1
c
c   Read the namelist values.
c
      ipv=0
      npvsets=0
      iskpd1=0
      iskppvd1=0
      iobsprc=0
      iobscnv=0
      iftcnv=0
      do i=1,maxptimes
         ptimes(i)=9e9
         iptimes(i)=99999999
         discard(i)=' '
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
c
      iendc=index(argum(ncn),' ')-1
      if (iendc.eq.-1) iendc=256
      iendf1=iendc+12
      nxtavl=0
c
c   Open the observations files
c
      if (iobsprc.eq.1) open (unit=iuobsprc,
     &   file='obsprc.dat',form='unformatted',status='old')
      if (iobscnv.eq.1) then
         open (unit=iuobsclg,file='obsclg.dat',form='unformatted',
     &      status='old')
         open (unit=iuobsvis,file='obsvis.dat',form='unformatted',
     &      status='old')
      endif
      if (iftcnv.eq.1) then
         open (unit=iuftclg,file='ftclg.dat',form='unformatted',
     &      status='old')
         open (unit=iuftvis,file='ftvis.dat',form='unformatted',
     &      status='old')
      endif
c
c   Define iounit for reading conventional data.
c
      open (unit=iudatin,file=argum(nsetsbeg),
     &   form='unformatted',status='old')
c
      igotmm4v8h=0
c
      if (dataform(1:8).eq.'mm4v8out') then   ! mm4v8 output data
c
c   Read information from header record of mm4v8 ouput header file.
c     Note: plev in header is not used anywhere else in ripdp.
c
      print*
      print*,'*******   Reading mm4v8out header file   *******'
      read(iudatin,err=200,end=220)  mdateb,iprog,icoord,
     &              ksigf,(sigf(k),k=1,ksigf),
     &              mif,mrf,mlf,jyr,jmo,jdy,jhr,
     &              ioldprog,ioldcoord,ioldnlv,(plev(k),k=1,ioldnlv),
     &              ibltyp,  isfflx,  itgflg,  icdcon,
     &              isfpar,  ivmixm,  idry,    imoist,  icloud,
     &              iboudy,  imoiav,  ifsnow,  icustb,  itqpbl,
     &              ifrad
      goto 230
  200 print*,'   Error with reading header record of header file.'
      stop
  220 print*,'   End of file encountered.'
  230 continue
c
c   Read required 2-d arrays from header file.
c
      read(iudatin) ter                    ! terrain
      read(iudatin) xmap                   ! xmap
      read(iudatin) dmap                   ! dmap
      read(iudatin) dotcor                 ! dotcor
      if(itgflg.ne.3) read(iudatin)        ! reztmp
      read(iudatin) xlat                   ! xlat
      read(iudatin) xlon                   ! xlon
      read(iudatin) xlus                   ! xlus
      read(iudatin) sno                    ! sno
      igotmm4v8h=1
c
c   Read information from header file for pv output
c
      mjxqb=1
      miyqb=1
      if (ipv.eq.1) then
         iuqin=iupv
         write(pvdfname,'(a6,i1,1x)') 'pvdata',iuqin-iupv+1
         open (unit=iuqin,file=pvdfname,form='unformatted',
     &      status='old')
         print*
         print*,'*******   Reading mm4v8out pv header file   *******'
         read(iuqin) ixq,jxq,miyqb,mjxqb,nvarq,nvartq,nvarvq
         print *,' Header file says there should be ',
     &      nvarq,nvartq,nvarvq,' partitions'
         print *,'   of pv, theta, and velocity respectively.'
         if (nvarq.ne.nvq.or.nvartq.ne.nvtq.or.nvarvq.ne.nvvq) then
            print *,'   But you have them set to ',
     &               nvq,nvtq,nvvq
            stop
         else
            print*,'   You have set the pv dimensions correctly.'
         endif
         miyqset=miyqb-1
         mjxqset=mjxqb-1
         miyqe=miyqb+ixq-1
         mjxqe=mjxqb+jxq-1
c         ntq4av=0
c         do k=1,mkzh
c         do j=1,mjx
c         do i=1,miy
c            tq4av(i,j,k)=0.
c         enddo
c         enddo
c         enddo
      endif
c
c   Prepare to read timestep data.
c
c
      if (iskpd1.eq.1) then
         close (iudatin)
         ifilecount=ifilecount+1
         open (unit=iudatin,file=argum(ifilecount-1+nsetsbeg),
     &      form='unformatted',status='old')
      else
         call jumpendfile(iudatin,1)
      endif
      if (ipv.eq.1) then
         if (iskppvd1.eq.1) then
            close (iuqin)
            iuqin=iuqin+1
            write(pvdfname,'(a6,i1,1x)') 'pvdata',iuqin-iupv+1
            open (unit=iuqin,file=pvdfname,form='unformatted',
     &         status='old')
         else
            call jumpendfile(iuqin,1)
         endif
      endif
      ipvdo=0
c
      endif  ! end of special mm4v8 stuff
ccccccccccccc temp tend
c      open (unit=68,file='fort.68',form='unformatted',
c     &   status='old')
ccccccccccccc temp tend
c
c   Initialize "max" time level. This feature causes RIPDP to
c      ignore standard data which is out of chronological order.
c
      xtimemax=-1.
c
c   LOOP THROUGH TIME LEVELS.
c
c   Read header record of data file.
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
      ifirstread=1
c      ifile=0
  240 if (dataform(1:8).eq.'mm4v8out') then
         read(iudatin,end=254,err=250) xtime,idry,imoist,
     &      ibltyp,isfpar,itgflg
      elseif (dataform.eq.'mm5v0out') then
         read(iudatin,end=254,err=250)  mdateb,iprog,icoord,
     &      ksigf,(sigf(k),k=1,ksigf),mif,mrf,mlf,jyr,jmo,jdy,jhr,
     &      xtime,ndumb,(idumb(l), l = 1, ndumb),ioldprog,ioldcoord,
     &      ioldnlv,(plev(k),k=1,ioldnlv)
      elseif (dataform.eq.'mm5v1   ') then
         read(iudatin,end=254,err=250) mifv1,mrfv1,mifc,mrfc
         iprog=mifv1(1,1)
         if (iprog.eq.6) then
            xtime=mrfv1(1,6)/60.
            mdateb=mifv1(2,6)
            if (mdateb.lt.100000000) then
c              "8-digit" mdate, 20th century
               mdateb=mdateb   ! leave as is for RIP
            elseif (mdateb.gt.100000000) then
c              "8-digit" mdate with leading "1" (9th digit) -> 21st century
               mdateb=mdateb-100000000   ! drop leading 1
            endif
         elseif (iprog.eq.5) then
            mdateb=mifv1(4,5)
            if (mdateb.lt.100000000) then
c              "8-digit" mdate, 20th century
               mdateb=mdateb   ! leave as is for RIP
            elseif (mdateb.gt.100000000.and.mdateb.lt.200000000) then
c              "8-digit" mdate with leading "1" (9th digit) -> 21st century
               mdateb=mdateb-100000000   ! drop leading 1
c
c         The following are commented out due to problems with compilers
c         on some machines with IEEE-32-bit integers, but if your data
c         is V1 from a Cray run, you should uncomment these.
c
c            elseif (mdateb.gt.200000000.and.mdateb.lt.10000000000) then
cc              "10-digit" mdate, 20th century
c               mdateb=nint(mdateb/100.)   ! drop last two digits
c            elseif (mdateb.gt.10000000000) then
cc              "10-digit" mdate with leading "1" (11th digit) -> 21st century
c               mdateb=mdateb-10000000000   ! drop leading 1
c               mdateb=nint(mdateb/100.)   ! drop last two digits
            endif
            call mconvert(mdateb,mhourb,1,1940)
            mdate=mifv1(1,5)
            if (mdate.lt.100000000) then
c              "8-digit" mdate, 20th century
               mdate=mdate   ! leave as is for RIP
            elseif (mdate.gt.100000000.and.mdate.lt.200000000) then
c              "8-digit" mdate with leading "1" (9th digit) -> 21st century
               mdate=mdate-100000000   ! drop leading 1
c
c         The following are commented out due to problems with compilers
c         on some machines with IEEE-32-bit integers, but if your data
c         is V1 from a Cray run, you should uncomment these.
c
c            elseif (mdate.gt.200000000.and.mdate.lt.10000000000) then
cc              "10-digit" mdate, 20th century
c               mdate=nint(mdate/100.)   ! drop last two digits
c            elseif (mdate.gt.10000000000) then
cc              "10-digit" mdate with leading "1" (11th digit) -> 21st century
c               mdate=mdate-10000000000   ! drop leading 1
c               mdate=nint(mdate/100.)   ! drop last two digits
            endif
            call mconvert(mdate,mhour,1,1940)
            xtime=float(mhour-mhourb)
         else
            print*,'I only can digest model input or output.'
            print*,'This data is from program number ',iprog
            stop
         endif
      elseif (dataform.eq.'mm5v3   ') then
 796     read(iudatin,end=254,err=250) flag
         if (flag.eq.0) then ! big header
            read(iudatin) bhi, bhr, bhic, bhrc
            goto 796
         elseif (flag.eq.2) then
            goto 254
         endif
         read (iudatin) ndim,start_index,end_index,time,staggering,
     &      ordering,current_date,name,units,description
         backspace (iudatin)
         backspace (iudatin)
         iprog=bhi(1,1)
         if (iprog.ge.2) then
            mdateb=1000000*mod(bhi(5,iprog),100)+10000*bhi(6,iprog)+
     &         100*bhi(7,iprog)+bhi(8,iprog)
            rhourb=bhi(9,iprog)/60.+bhi(10,iprog)/3600.+
     &         bhi(11,iprog)/36000000.
            call mconvert(mdateb,mhourb,1,1940)
            read(current_date,'(1x,6(1x,i2),1x,i4)')
     &         iyr,imo,idy,ihr,imn,isc,itt
            mdate=1000000*iyr+10000*imo+100*idy+ihr
            rhour=imn/60.+isc/3600.+itt/36000000.
            call mconvert(mdate,mhour,1,1940)
            xtime=float(mhour-mhourb)+rhour-rhourb
            if (iprog.eq.11) then
               xtime2=time/60.
               if (abs(xtime-xtime2).gt..001) then
                  print*,'Seems to be an inconsistency between "time"'
                  print*,'and "current_date" in V3 header.'
                  print*,'current_date,time,xtime,xtime2='
                  print*,current_date,'  ',time,xtime,xtime2
                  stop
               endif
            endif
         elseif (iprog.eq.1) then ! Terrain data - date/time info irrelevant
            mdateb=651106   ! Mark Stoelinga's birthday
            call mconvert(mdateb,mhourb,1,1940)
            rhourb=0.
            mdate=mdateb
            mhour=mhourb
            rhour=0.
            xtime=0.
         endif
      else
         print*,'   Unrecognized dataform.'
         stop
      endif
c
c      ifile=ifile+1
c      if (ifile.lt.ifilebeg) then
c         if (dataform.eq.'mm5v1   ') then
c            call skpflds(iudatin,mifv1,1)
c         else
c            call jumpendfile(iudatin,1)
c         endif
c         goto 240
c      endif
c      if (ifile.gt.ifileend) goto 1000
c
      if (ifirstread.eq.1) then
         print*,'Data is recognized as ',dataform,','
         print*,'   from program number ',iprog,'.'
         print*
      endif
      if (dataform(1:5).eq.'mm4v8'.or.dataform(1:5).eq.'mm5v0'.or.
     &    dataform(1:5).eq.'mm5v1') then
         if (iprog.eq.6) then
            call mconvert(mdateb,mhourb,1,1940)
            mhour=mhourb+int(xtime)
            call mconvert(mdate,mhour,-1,1940)
         elseif (iprog.ne.5) then
            print*,'Can''t process these outputs prior to MM5V3'
         endif
         rhourb=0.00
         rhour=xtime-float(mhour-mhourb)
      endif
c
      secondspast=rhour*3600.
      if (iprog.eq.1) then
         print*,' ****  Reading terrain data.'
      elseif (iprog.eq.2.or.iprog.eq.3) then
         print*,' ****  Reading pressure-level analysis data at'
      elseif (iprog.eq.5) then
         print*,' ****  Reading sigma-level analysis data at'
      elseif (iprog.eq.6.or.iprog.eq.11) then
         print*,' ****  Reading sigma-level model output at'
         print*,'       forecast time=',xtime
      endif
      if (iprog.eq.1) goto 235
      if (secondspast.lt..2) then
         write(6,926) mdate
      else
         write(6,927) mdate,secondspast
      endif
 926  format('        (YYMMDDHH = ',i8.8,')')
 927  format('        (YYMMDDHH = ',i8.8,' plus ',f12.5,' seconds)')
 235  continue
c
      goto 252
  250 print*,'Error reading header record of output file.'
      stop
  252 continue
c
      goto 256
  254 close (iudatin)
      ifilecount=ifilecount+1
      if (ifilecount.eq.nsets+1) then
         print*
         print*,'*******   No more datasets   *******'
         goto 1000
      endif
      open (unit=iudatin,file=argum(ifilecount-1+nsetsbeg),
     &   form='unformatted',status='old')
      goto 240
  256 continue
c
      if (ifirstread.eq.1) then
         ifirstread=0
c
c      If using iptimes, convert the mdates in the iptimes array to
c      xtimes in the ptimes array.  Also, determine nptimes.
c
         if (ptimes(1).lt.0..or.iptimes(1).lt.0.or.(ptimes(1).eq.
     &       9e9.and.iptimes(1).eq.99999999)) then !user wants all times
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
 259        continue
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
 261        continue
         endif
c
c   Process time sequences in ptimes array.
c
         ii=0
         itime=0
         ptusemax=-9e9
  400    ii=ii+1
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
     &                 'Number of times requested exceeds maxptimes.'
                     print*,
     &                 'Increase maxptimes in ripdp code, recompile,'
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
  420    nptuse=itime
         if (nptuse.eq.0) ptusemax=9e9
c
c      Get some flags that are necessary for reading the right variables.
c
         if (dataform.eq.'mm5v1   ') then
            refslp=mrfv1(2,iprog)
            refslt=mrfv1(3,iprog)
            reflaps=mrfv1(4,iprog)
            refstratt=.1  ! just above absolute zero
            inhyd=mifv1(5,iprog)
            if (iprog.eq.6) then
               iice=mifv1(7,6)
            else
               iice=0
            endif
            if (iprog .eq. 6) then
              fname=argum(1)(1:iendc)//'.minfo'
              open (unit=58,file=fname,form='formatted',
     &             status='unknown')
              idx = index(mifc(1,1),'-')
              version = 'V2'
              if (idx .ne. 0) then
                version(2:2) = mifc(1,1)(idx+1:idx+1)
                version(3:6) = '.'//mifc(1,1)(idx+3:idx+5)
              endif
              itime = nint(mrfv1(309,6)*mrfv1(101,1)/mrfv1(1,1))
              write(58,1022) version, cupa(mifv1(318,6)),
     &          bltyp(mifv1(315,6)+1),
     &       mphys(mifv1(353,6)),nint(mrfv1(101,1)),nint(mrfv1(101,5))
     &         ,itime
 1022 format('Model info: ',a6,1x,3(a10,1x),i4,' km, ',i3,' levels, ',
     &   i4,' sec')
            endif
         elseif (dataform.eq.'mm5v3   ') then
            if (iprog.eq.5.or.iprog.eq.11) then
               refslp=bhr(2,5)
               refslt=bhr(3,5)
               reflaps=bhr(4,5)
               if (bhr(5,5).gt.100..and.bhr(5,5).lt.400.) then
                  refstratt=bhr(5,5)
               else
                  refstratt=.1  ! just above absolute zero
               endif
               inhyd=1
               if (iprog.eq.11) then
                  iice=bhi(18,11)
c
c               Check for possible lack of ice variables (in spite of iice
c               being =1 in V3 header), due to inconsistency in iice and
c               IMPHYS.
c
                  imphys=bhi(3,13)
                  iicemphys=0
                  if (imphys.ge.5) iicemphys=1
                  if (iice.eq.1.and.iicemphys.eq.0) then
                     print*,'Inconsistency in V3 header: iice=1 ',
     &                      'but imphys=',imphys
                     print*,'Imphys is more trustworthy, so we''ll',
     &                      ' set iice to 0.'
                     iice=0
                  endif
               else
                  iice=0
               endif
c
c            Find sigh (those @#$*&!$% NCAR people put it at the very END
c            of the output for a time period).
c
 603           read(iudatin) flag
               if (flag.eq.1) then
                  read (iudatin) ndim,start_index,end_index,time,
     &               staggering,ordering,current_date,name,
     &               units,description
                  if (name.ne.'SIGMAH   ') then
                     read(iudatin)
                     goto 603
                  else
                     if (end_index(1).ne.mkzh) then
                        print*,'SIGMAH has more elements than expected.'
                        print*,'end_index(1),mkzh=',end_index(1),mkzh
                        stop
                     endif
                     read(iudatin)sigh
                  endif
               else
                  print*,'Ran into a flag not =1 looking for SIGMAH.'
                  stop
               endif
               rewind (iudatin)
 199           read(iudatin) flag
               if (flag.eq.0) then ! big header
                  read(iudatin) bhi, bhr, bhic, bhrc
                  goto 199
               elseif (flag.eq.2) then
                  print*,'Something is wrong.  Ripdp has rewound'
                  print*,'the data file after finding SIGMAH, and'
                  print*,'is looking for the first big header, but'
                  print*,'flag=2, meaning "end of this time period".'
                  stop
               endif
               backspace (iudatin)
c
c            While we're here, write out model info to the Jim Bresch-inspired
c            ".minfo" file
c
               if (iprog .eq. 11) then
                  fname=argum(ncn)(1:iendc)//'.minfo'
                  open (unit=58,file=fname,form='formatted',
     &                  status='unknown')
                  write(version,'(a2,a1,i1,a1,i1)') 'V3','.',bhi(3,11),
     &               '.',bhi(4,11)
c                 print*,'version = ',version
c                 itime = nint(mrfv1(309,6)*mrfv1(101,1)/mrfv1(1,1))
                  itime = nint(bhr(2,12)/float(bhi(20,1)))
                  write(58,1021) version, cupa(bhi(2,13)),
     &               bltyp(bhi(4,13)+1),mphys(bhi(3,13)),
     &               nint(bhr(1,1)*.001/float(bhi(20,1))),
     &               bhi(12,5),itime
 1021             format('Model info: ',a6,1x,3(a10,1x),i4,' km, ',i3,
     &               ' levels, ',i4,' sec')
               endif
            elseif (iprog.eq.2.or.iprog.eq.3) then
               inhyd=0
               iice=0
c
c            Define pseudo-sigma levels.  First get pressure levels.
c
 713           read(iudatin) flag
               if (flag.eq.1) then
                  read (iudatin) ndim,start_index,end_index,time,
     &               staggering,ordering,current_date,name,
     &               units,description
                  if (name.ne.'PRESSURE ') then
                     read(iudatin)
                     goto 713
                  else
                     if (end_index(1).ne.mkzh+1) then
                        print*,'PRESSURE has more elements than',
     &                         ' expected.'
                        print*,'end_index(1),mkzh=',end_index(1),mkzh
                        stop
                     endif
c
c                  Read in pressure levels.  Reverse the order, so that
c                     prslvl(1) is the top (lowest prs).
c
                     read(iudatin)(prslvl(k),k=mkzh+1,1,-1)
c                     print*,'Surface given as p=',.01*prslvl(mkzh+1),
c     &                  ' hPa.'
                     do k=1,mkzh
                        prslvl(k)=.01*prslvl(k)  ! Pa to hPa
                     enddo
                  endif
               else
                  print*,'Ran into a flag not =1 looking for PRESSURE.'
                  stop
               endif
               rewind (iudatin)
 299           read(iudatin) flag
               if (flag.eq.0) then ! big header
                  read(iudatin) bhi, bhr, bhic, bhrc
                  goto 299
               elseif (flag.eq.2) then
                  print*,'Something is wrong.  Ripdp has rewound'
                  print*,'the data file after finding PRESSURE, and'
                  print*,'is looking for the first big header, but'
                  print*,'flag=2, meaning "end of this time period".'
                  stop
               endif
               backspace (iudatin)
c
c            Now define half sigma levels
c
c               dpavg=(prslvl(mkzh)-prslvl(1))/(mkzh-1)
c               ptop=max(prslvl(1)-.5*dpavg,10.)  ! min of 10 hPa for ptop
c               pbot=prslvl(mkzh)+.5*dpavg
               ptop=prslvl(1)-1.  ! 1 hPa above highest half sigma level
               pbot=prslvl(mkzh)+1.  ! 1 hPa below lowest half sigma level
               pstarconst=pbot-ptop
               do k=1,mkzh
                  sigh(k)=(prslvl(k)-ptop)/pstarconst
               enddo
            elseif (iprog.eq.1) then ! Terrain data
               inhyd=0
               iice=0
c
c            Define pseudo-sigma levels.
c
               ptop=50.
               pbot=1000.
               pstarconst=pbot-ptop
               sigh(1)=.5
            endif
         else
            refslp=1.e5  ! this is in Pascals
            refslt=mrf(9)
            reflaps=mrf(10)
            refstratt=.1  ! just above absolute zero
            if (dataform(1:8).eq.'mm4v8out') then
               iice=0
               inav=0
c
c   The following holerith conditional will have to be uncommented for
c   ripdp to process MM4V8 model output.  It is normally commented
c   out and replaced by a trivial conditional because some compilers
c   no longer recognize holleriths.
c
c               if (iprog.eq.5HMM4NH) then
               if (iprog.eq.0) then
                  inhyd=1
               else
                  inhyd=0
               endif
            elseif (dataform.eq.'mm5v0out') then
               idry =    idumb(1)
               imoist =  idumb(2)
               itgflg =  idumb(3)
               iice =    idumb(4)
               inav =    idumb(5)
               inhyd =   idumb(6)
            endif
         endif
c
      endif
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
c   (unless they are out of chronological order).  Also, for TERRAIN
c   output (iprog=1), the data are processed regardless of the values
c   in ptimes and iptimes.
c
      if (iprog.ne.1.and.iskipit.eq.0.and.nptuse.gt.0) then
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
c
         if (dataform.eq.'mm5v1   ') then
            call skpflds(iudatin,mifv1,1)
         elseif (dataform.eq.'mm5v3   ') then
 723        read(iudatin) flag
            if (flag.eq.1) then
               read(iudatin)
               read(iudatin)
               goto 723
            elseif (flag.eq.2) then
               continue
            endif
         else
            call jumpendfile(iudatin,1)
         endif
         goto 240
      endif
c
c   Get important info from headers.
c
      if (dataform(1:8).eq.'mm4v8out') then
         do k=1,mkzh
            sigh(k)=.5*(sigf(k)+sigf(k+1))
         enddo
         miycors=mif(2)
         mjxcors=mif(3)
         dskmc=mrf(1)
         xlatc=mrf(2)
         xlonc=mrf(3)
c
c   The following holerith conditional will have to be uncommented for
c   ripdp to process MM4V8 model output.  It is normally commented
c   out and replaced by a trivial conditional because some compilers
c   no longer recognize holleriths.
c
c         if (mif(17).eq.6HLAMCON) then
         if (mif(17).eq.0) then
            nproj=1
c         elseif (mif(17).eq.6HPOLSTR) then
         elseif (mif(17).eq.1) then
            nproj=2
c         elseif (mif(17).eq.6HMERCAT) then
         elseif (mif(17).eq.2) then
            nproj=3
         else
            print 298, mif(17)
            print*,'This is not a recognized map projection.'
            nproj=1
c            stop
         endif
 298     format('mif(17)=',a6)
         ilandset=1
         truelat1=91.
         truelat2=91.
         if (mif(25).gt.0.and.mif(25).lt.92.and.
     &       mif(26).gt.0.and.mif(26).lt.92) then
            truelat1=float(mif(25))
            truelat2=float(mif(26))
         endif
         ptop=mrf(6)
         if (mlf(3)) then
            yicorn=mif(6)
            xjcorn=mif(7)
            dskm=mrf(4)
         else
            yicorn=1.
            xjcorn=1.
            dskm=dskmc
         endif
      elseif (dataform.eq.'mm5v0out') then
         do k=1,mkzh
            sigh(k)=.5*(sigf(k)+sigf(k+1))
         enddo
         miycors=mif(2)
         mjxcors=mif(3)
         dskmc=mrf(1)
         xlatc=mrf(2)
         xlonc=mrf(3)
c
c   The following holerith conditional will have to be uncommented for
c   ripdp to process MM4V8 model output.  It is normally commented
c   out and replaced by a trivial conditional because some compilers
c   no longer recognize holleriths.
c
c         if (mif(17).eq.6HLAMCON) then
         if (mif(17).eq.0) then
            nproj=1
c         elseif (mif(17).eq.6HPOLSTR) then
         elseif (mif(17).eq.1) then
            nproj=2
c         elseif (mif(17).eq.6HMERCAT) then
         elseif (mif(17).eq.2) then
            nproj=3
         else
            write(6,'(''mif(17)='',a6)') mif(17)
            print*,'This is not a recognized map projection.'
            stop
         endif
         ilandset=1
         truelat1=91.
         truelat2=91.
         ptop=mrf(6)
         if (mlf(3)) then
            yicorn=mif(6)
            xjcorn=mif(7)
            dskm=mrf(4)
         else
            yicorn=1.
            xjcorn=1.
            dskm=dskmc
         endif
      elseif (dataform.eq.'mm5v1   ') then
         sigf(1)=0.
         do k=1,mkzh
            sigh(k)=mrfv1(101+k,iprog)
            sigf(k+1)=2.*sigh(k)-sigf(k)
         enddo
         miycors=mifv1(2,1)
         mjxcors=mifv1(3,1)
         dskmc=mrfv1(1,1)
         xlatc=mrfv1(2,1)
         xlonc=mrfv1(3,1)
         nproj=mifv1(4,1)
         if (nproj.gt.3.or.nproj.lt.1) then
            print*,'   Map proj. #',nproj,' is not recognized.'
            stop
         endif
         if (mifv1(11,1).lt.0.or.mifv1(11,1).gt.200) then
            ilandset=1
         elseif (mifc(11,1)(1:5).eq.'OLD  '.and.mifv1(11,1).eq.7) then
            ilandset=1
         elseif (mifc(11,1)(1:5).eq.'USGS '.and.mifv1(11,1).eq.16) then
            ilandset=2
         elseif (mifc(11,1)(1:5).eq.'SiB  '.and.mifv1(11,1).eq.15) then
            ilandset=3
         else
            print*,'RIPDP does not recognize land use data set'
            print*,'specified in header.  mifv1(11,1)=',mifv1(11,1)
            print*,'mifc(11,1)=',mifc(11,1)
            stop
         endif
         truelat1=mrfv1(5,1)
         truelat2=mrfv1(6,1)
         ptop=mrfv1(1,2)
         yicorn=mrfv1(102,1)
         xjcorn=mrfv1(103,1)
         dskm=mrfv1(101,1)
      elseif (dataform.eq.'mm5v3   ') then
         sigf(1)=0.
         if (iprog.eq.1) then
            sigf(2)=1.0
         elseif (iprog.eq.2.or.iprog.eq.3) then
            sigf(mkzh+1)=1.0
c
c         Note!! In the following, the full sigma levels are defined as
c         the midpoint between the two surrounding half sigma levels.
c         This is different from the standard definition, in which the
c         half sigma levels are defined as the midpoint between the two
c         surrounding full levels.  It must be done this was in order to
c         guarantee that pseudo-sigma levels can be defined so that half
c         sigma levels are exactly coincident with the pressure levels in
c         the pressure level output.  THEREFORE, NO CODE SHOULD BE WRITTEN
c         IN RIPDP OR RIP THAT ASSUMES THE STANDARD RELATIONSHIP BETWEEN
c         HALF AND FULL SIGMA LEVELS.
c
            do k=2,mkzh
               sigf(k)=.5*(sigh(k-1)+sigh(k))
            enddo
         elseif (iprog.eq.5.or.iprog.eq.11) then
            do k=1,mkzh
               sigf(k+1)=2.*sigh(k)-sigf(k)
            enddo
         endif
         miycors=bhi(5,1)
         mjxcors=bhi(6,1)
         if (iexpanded.eq.1) then
            miycors=miy
            mjxcors=mjx
            ioffexp=bhi(11,1)
            joffexp=bhi(12,1)
         endif
         dskmc=.001*bhr(1,1)
         xlatc=bhr(2,1)
         xlonc=bhr(3,1)
         nproj=bhi(7,1)
         if (nproj.gt.3.or.nproj.lt.1) then
            print*,'   Map proj. #',nproj,' is not recognized.'
            stop
         endif
         if (bhic(23,1)(1:5).eq.'OLD  '.and.bhi(23,1).eq.7) then
            ilandset=1
         elseif (bhic(23,1)(1:5).eq.'USGS '.and.bhi(23,1).eq.16) then
            ilandset=2
         elseif (bhic(23,1)(1:5).eq.'SiB  '.and.bhi(23,1).eq.15) then
            ilandset=3
         else
            print*,'RIPDP does not recognize land use data set'
            print*,'specified in header.  bhi(23,1)=',bhi(23,1)
            print*,'bhic(23,1)=',bhic(23,1)
            stop
         endif
         truelat1=bhr(5,1)
         truelat2=bhr(6,1)
         if (abs(bhr(7,1)).ne.90.) then
            print*,'Rip is only designed to deal with map backgrounds'
            print*,'that have pole=90 deg. or -90 deg.'
            stop
         endif
         if (iprog.eq.5.or.iprog.eq.11) ptop=.01*bhr(2,2) ! want it in hPa
         yicorn=bhr(10,1)
         xjcorn=bhr(11,1)
         dskm=.001*bhr(9,1)
      endif
c
      dsc=dskmc*1000.
      ds=dskm*1000.
      refrat=dsc/ds
c
c   Look for matching time in pv data file
c
      if (ipv.eq.1) then
  300    read(iuqin,end=305) xtimepv
  303    goto 306
  305    close (iuqin)
         iuqin=iuqin+1
         if (iuqin.eq.40+npvsets+1) then
            print*,'   Couldn''t find matching time in pv dataset.'
            goto 307
         endif
         write(pvdfname,'(a6,i1,1x)') 'pvdata',iuqin-iupv+1
         open (unit=iuqin,file=pvdfname,form='unformatted',
     &      status='old')
         goto 300
  306    if (xtimepv-xtime.lt.-.0001) then
            call jumpendfile(iuqin,1)
            goto 300
         elseif (xtimepv-xtime.gt..0001) then
            print*,'   pv dataset starts at later time.'
            backspace (iuqin)
            ipvdo=0
         else
            print*,'   Found matching time in pv dataset.'
            ipvdo=1
         endif
      endif
      goto 308
  307 continue
      ipv=0
      ipvdo=0
  308 continue
c
c   Look for matching time in observations data files
c
      call getobs(iobsprc,iobsprcfnd,mdate,iuobsprc)
      call getobs(iobscnv,iobsclgfnd,mdate,iuobsclg)
      call getobs(iobscnv,iobsvisfnd,mdate,iuobsvis)
      call getobs(iftcnv,iftclgfnd,mdate,iuftclg)
      call getobs(iftcnv,iftvisfnd,mdate,iuftvis)
c
c   Set up nonhydrostatic stuff.
c
      if (inhyd.eq.1) then
         if (refslp.lt.800e2.or.refslp.gt.1200e2
     &       .or.refslt.lt.240..or.refslt.gt.360.
     &       .or.reflaps.lt.-5..or.reflaps.gt.200.) then
            print*,'   Refslp, Refslt and/or Reflaps seem funny.'
            print*,'   Using "reasonable" values:'
            print*,'   Refslp=1000e2, Refslt=290., Reflaps=50.'
            refslp=1000e2
            refslt=290.
            reflaps=50.
            refstratt=.1  ! just above absolute zero
         endif
      else
         refslp=1000e2
         refslt=290.
         reflaps=50.
         refstratt=.1  ! just above absolute zero
      endif
c
c   Set all "gotit" flags to 0, except for arrays that have been
c   read in from the MM4V8 header
c
      igotituuu=0
      igotituuu_sfan=0   ! "sfan" means "surface analysis"
      igotitvvv=0
      igotitvvv_sfan=0
      igotittmk=0
      igotittmk_sfan=0
      igotitqvp=0
      igotitrh=0
      igotitrh_sfan=0
      igotitqcw=0
      igotitqra=0
      igotitqci=0
      igotitqsn=0
      igotitqgr=0
      igotitrnci=0
      igotittke=0
      igotitradtnd=0
      igotitwww=0
      igotitprs=0
      igotitght=0
      igotitslp=0
      igotitpstx=0
      igotittgk=0
      igotitsst=0
      igotitrtc=0
      igotitrte=0
      igotitter=0
ccc      igotitter_tsf=0   ! "tsf" means "at true surface"
      igotitxmap=0
      igotitdmap=0
      igotitxlat=0
      igotitxlon=0
      igotitdotcor=0
      igotitxlus=0
      igotitsno=0
      igotitpblh=0
      igotitregime=0
      igotitsshflux=0
      igotitslhflux=0
      igotitxmav=0     ! jfb
      igotitust=0
      igotitswdown=0
      igotitlwdown=0
      igotitsoil1=0
      igotitsoil2=0
      igotitsoil3=0
      igotitsoil4=0
      igotitsoil5=0
      igotitsoil6=0
c
      if (igotmm4v8h.eq.1) then
         igotitter=1
         igotitxmap=1
         igotitdmap=1
         igotitxlat=1
         igotitxlon=1
         igotitdotcor=1
         igotitxlus=1
         igotitsno=1
      endif
c
c   Initialize pressure array to zero
c
      do k=1,mkzh
      do j=1,mjx
      do i=1,miy
         prs(i,j,k)=0.
      enddo
      enddo
      enddo
c
c   Read required 2-d arrays from time-level file if mm5v0 output.
c
      if (dataform.eq.'mm5v0out') then     ! mm5v0 output data
         read(iudatin) ter                    ! terrain
         igotitter=1
         read(iudatin) xmap                   ! xmap
         igotitxmap=1
         read(iudatin) dmap                   ! dmap
         igotitdmap=1
         read(iudatin) dotcor                 ! dotcor
         igotitdotcor=1
         if(itgflg.ne.3) read(iudatin)        ! reztmp
         read(iudatin) xlat                   ! xlat
         igotitxlat=1
         read(iudatin) xlon                   ! xlon
         igotitxlon=1
         read(iudatin) xlus                   ! xlus
         igotitxlus=1
         read(iudatin) sno                    ! sno
         igotitsno=1
      endif
c
c   Set up default values of true latitudes if necessary.
c
      if (nproj.eq.3) then  ! Mercator
         if (truelat1.gt.90..or.truelat2.gt.90.) then
            truelat1=0.
            truelat2=truelat1   !not used
         endif
      elseif (nproj.eq.1) then  ! Lambert Conformal
         if (truelat1.gt.90..or.truelat2.gt.90.) then
            truelat1=sign(30.,xlatc)
            truelat2=sign(60.,xlatc)
         endif
      elseif (nproj.eq.2) then   ! Polar Stereographic
         if (truelat1.gt.90.) then
            truelat1=sign(60.,xlatc)
            truelat2=truelat1   !not used
         endif
      else
         write(6,*)'Unrecognized map projection.'
         stop
      endif
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
     &'map projection (1: Lam. Conf., 2: Pol. Ster., 3: Mercator)'
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
      chrip(10)=
     &'mdateb: YYMMDDHH (truncated hour) of hour-0 for this dataset'
      ihrip(10)=mdateb
      chrip(11)=
     &'mdate: YYMMDDHH (truncated hour) of this time'
      ihrip(11)=mdate
      chrip(12)=
     &'ice physics (1: sep. arrays for ice fields, 0: no sep. arrays)'
      ihrip(12)=iice
ccc      chrip(13)=
ccc     &'Program #: 1:TER. 2:DG/RG. 3:RAW. 5:INT. 6:MOD. 11:MOD.(MM5V3)'
      chrip(13)=
     &'ver. coord. type: <or=3: hgt. or prs.; >or=4: terrain-following'
      ihrip(13)=iprog
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
      rhrip(7)=yicorn
      chrip(ijmp+8)=
     &'coarse dom. x-position of lower left corner of this domain'
      rhrip(8)=xjcorn
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
      rhrip(14)=rhour
      chrip(ijmp+15)=
     &'xtime: exact time of this data relative to exact hour-0 (in h)'
      rhrip(15)=xtime
ccc      chrip(ijmp+16)=
ccc     &'reference stratospheric constant temperature (K)'
ccc      rhrip(16)=refstratt
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
c   Read in variables.
c
      if (dataform(1:8).eq.'mm4v8out'.or.
     &    dataform.eq.'mm5v0out') then
c
      read(iudatin) uuu
      igotituuu=1
      read(iudatin) vvv
      igotitvvv=1
      read(iudatin) tmk
      igotittmk=1
c
      if(idry.eq.0) then
         read(iudatin) qvp
         igotitqvp=1
      endif
c
      if (idry.eq.0.and.imoist.ge.2) then
         read(iudatin) qcw
         igotitqcw=1
         read(iudatin) qra
         igotitqra=1
      endif
c
      if (idry.eq.0.and.imoist.ge.2.and.iice.eq.1) then
         read(iudatin)
         read(iudatin) qci
         igotitqci=1
         read(iudatin) qsn
         igotitqsn=1
      endif
c
      read(iudatin) pstx
      igotitpstx=1
      if(itgflg.ne.3) then
         read(iudatin) tgk
         igotittgk=1
      endif
c
      if(idry.eq.0) then
         read(iudatin) rtc
         igotitrtc=1
         read(iudatin) rte
         igotitrte=1
      endif
c
      if (inhyd.eq.1) then
         read(iudatin) www
         igotitwww=1
         read(iudatin) prs
         igotitprs=1
      endif
c
      if (inav.eq.1) then
         read(iudatin)
      endif
      call jumpendfile(iudatin,1)
c
      elseif (dataform.eq.'mm5v1   ') then
c
      n3d = mifv1(201,iprog)
      n2d = mifv1(202,iprog)
      do ifld = 1, n3d
         icoupwarn=0
         idiscard=0
         do idis=1,ndiscard
            if (mifc(204+ifld,iprog)(1:9).eq.discard(idis)) then
               read(iudatin)
               idiscard=1
               print*,'   Discarding variable ',discard(idis)
               print*,'   because user does not want it.'
               goto 349
            endif
         enddo
         if (mifc(204+ifld,iprog)(1:2).eq.'U ') then
            read(iudatin) uuu
            igotituuu=1
         elseif (mifc(204+ifld,iprog)(1:2).eq.'V ') then
            read(iudatin) vvv
            igotitvvv=1
         elseif (mifc(204+ifld,iprog)(1:2).eq.'T ') then
            read(iudatin) tmk
            igotittmk=1
         elseif (mifc(204+ifld,iprog)(1:2).eq.'Q ') then
            read(iudatin) qvp
            igotitqvp=1
         elseif (mifc(204+ifld,iprog)(1:4).eq.'CLW ') then
            read(iudatin) qcw
            igotitqcw=1
         elseif (mifc(204+ifld,iprog)(1:4).eq.'RNW ') then
            read(iudatin) qra
            igotitqra=1
         elseif (mifc(204+ifld,iprog)(1:4).eq.'ICE ') then
            read(iudatin) qci
            igotitqci=1
         elseif (mifc(204+ifld,iprog)(1:5).eq.'SNOW ') then
            read(iudatin) qsn
            igotitqsn=1
         elseif (mifc(204+ifld,iprog)(1:8).eq.'GRAUPEL ') then
            read(iudatin) qgr
            igotitqgr=1
         elseif (mifc(204+ifld,iprog)(1:4).eq.'NCI ') then
            read(iudatin) rnci
            igotitrnci=1
         elseif (mifc(204+ifld,iprog)(1:4).eq.'TKE ') then
            read(iudatin) tke
            igotittke=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'RAD TEND ') then
            read(iudatin) radtnd
            igotitradtnd=1
         elseif (mifc(204+ifld,iprog)(1:2).eq.'W ') then
            read(iudatin) www
            igotitwww=1
         elseif (mifc(204+ifld,iprog)(1:3).eq.'PP ') then
            read(iudatin) prs
            igotitprs=1
         else
c
c         unknown 3d field
c
c         Note: V1 header: name is mifc(1:9) (9 chars.)
c                          units is mifc(10:26) (17 chars.)
c                          description is mifc(27:67) (41 chars.)
c                          coupled/uncoupled is mifc(68:75) (8 chars.)
c                          cross/dot is mifc(76:80) (5 chars.)
c              RIP header: varname*10,vardesc*64, plchun*24
c
            icoupwarn=1
            read(iudatin) scr3
            varname=mifc(204+ifld,iprog)(1:9)
            inname=0
            do ic=10,1,-1
               if (inname.eq.0) then
                  if (varname(ic:ic).ne.' ') inname=1
               else
                  if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
               endif
            enddo
            vardesc=mifc(204+ifld,iprog)(27:67)//', '//
     &         mifc(204+ifld,iprog)(10:26)
            plchun=mifc(204+ifld,iprog)(10:26)
            if (mifc(204+ifld,iprog)(76:80).eq.'CROSS') then
               icd=1
            else
               icd=0
            endif
            call writefile_rdp(scr3,varname,3,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
         endif
 349     continue
         if (idiscard.eq.0) then
            print*,'Processing MM5 variable ',mifc(204+ifld,iprog)(1:9)
            if (icoupwarn.eq.1) then
               print*,' Warning - this unexpected 3D array is probably'
               print*,' coupled to pstar (in kPa), but RIPDP cannot'
               print*,' decouple it.  You''ll have to do it yourself.'
            endif
         endif
      enddo
      do ifld = n3d+1, n3d+n2d
         idiscard=0
         do idis=1,ndiscard
            if (mifc(204+ifld,iprog)(1:9).eq.discard(idis)) then
               read(iudatin)
               idiscard=1
               print*,'   Discarding variable ',discard(idis)
               print*,'   because user does not want it.'
               goto 351
            endif
         enddo
         if (mifc(204+ifld,iprog)(1:9).eq.'RES TEMP '.or.
     &           mifc(204+ifld,iprog)(1:9).eq.'LATITDOT '.or.
     &           mifc(204+ifld,iprog)(1:9).eq.'LONGIDOT ') then
            read(iudatin)
            idiscard=1
            print*,'   Discarding variable ',
     &         mifc(204+ifld,iprog)(1:9)
            print*,'   because RIP does not need it.'
         elseif (mifc(204+ifld,iprog)(1:9).eq.'PSTARCRS ') then
            read(iudatin) pstx
            igotitpstx=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'GROUND T ') then
            read(iudatin) tgk
            igotittgk=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'RAIN CON ') then
            read(iudatin) rtc
            igotitrtc=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'RAIN NON ') then
            read(iudatin) rte
            igotitrte=1
         elseif (mifc(204+ifld,iprog)(1:8).eq.'TERRAIN ') then
            read(iudatin) ter
            igotitter=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'MAPFACCR ') then
            read(iudatin) xmap
            igotitxmap=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'MAPFACDT ') then
            read(iudatin) dmap
            igotitdmap=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'LATITCRS ') then
            read(iudatin) xlat
            igotitxlat=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'LONGICRS ') then
            read(iudatin) xlon
            igotitxlon=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'CORIOLIS ') then
            read(iudatin) dotcor
            igotitdotcor=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'LAND USE ') then
            read(iudatin) xlus
            igotitxlus=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'SNOWCOVR ') then
            read(iudatin) sno
            igotitsno=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'PBL HGT  ') then
            read(iudatin) pblh
            igotitpblh=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'REGIME   ') then
            read(iudatin) regime
            igotitregime=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'SHFLUX   ') then
            read(iudatin) sshflux
            igotitsshflux=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'LHFLUX   ') then
            read(iudatin) slhflux
            igotitslhflux=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'MAVAIL   ') then
            read(iudatin) xmav
            igotitxmav=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'UST      ') then
            read(iudatin) ust
            igotitust=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'SWDOWN   ') then
            read(iudatin) swdown
            igotitswdown=1
         elseif (mifc(204+ifld,iprog)(1:9).eq.'LWDOWN   ') then
            read(iudatin) lwdown
            igotitlwdown=1
         else
c
c         unknown 2d field
c
c         Note: V1 header: name is mifc(1:9) (9 chars.)
c                          units is mifc(10:26) (17 chars.)
c                          description is mifc(27:67) (41 chars.)
c                          coupled/uncoupled is mifc(68:75) (8 chars.)
c                          cross/dot is mifc(76:80) (5 chars.)
c              RIP header: varname*10,vardesc*64, plchun*24
c
            read(iudatin) scr2
            varname=mifc(204+ifld,iprog)(1:9)
            inname=0
            do ic=10,1,-1
               if (inname.eq.0) then
                  if (varname(ic:ic).ne.' ') inname=1
               else
                  if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
               endif
            enddo
            vardesc=mifc(204+ifld,iprog)(27:67)//', '//
     &         mifc(204+ifld,iprog)(10:26)
            plchun=mifc(204+ifld,iprog)(10:26)
            if (mifc(204+ifld,iprog)(76:80).eq.'CROSS') then
               icd=1
            else
               icd=0
            endif
            call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
         endif
 351     continue
         if (idiscard.eq.0) then
            print*,'Processing MM5 variable ',mifc(204+ifld,iprog)(1:9)
         endif
      enddo
c
      elseif (dataform.eq.'mm5v3   ') then
c
 421  read(iudatin) flag
      if (flag.eq.1) then
         read (iudatin) ndim,start_index,end_index,time,
     &      staggering,ordering,current_date,name,
     &      units,description
c         print*
c         print*,'ndim=',ndim
c         print*,'start_index=',start_index
c         print*,'end_index=',end_index
c         print*,'time=',time
c         print*,'staggering=',staggering
c         print*,'ordering=',ordering
c         print*,'current_date=',current_date
c         print*,'name=',name
c         print*,'units=',units
c         print*,'description=',description
         idiscard=0
         do idis=1,ndiscard
            if (name.eq.discard(idis)) then
               read(iudatin)
               idiscard=1
               print*,'   Discarding variable ',discard(idis)
               print*,'   because user does not want it.'
               goto 353
            endif
         enddo
         if (name.eq.'PSEALVLD '.or.
     &           name.eq.'LATITDOT '.or.
     &           name.eq.'LONGIDOT '.or.
     &           name.eq.'RES TEMP '.or.
     &           name.eq.'RES TMP  '.or.
     &           name.eq.'HSFC     ') then
            read(iudatin)
            print*,'   Discarding variable ',name
            print*,'   because RIP does not need it.'
            idiscard=1
         elseif (name.eq.'U        ') then
            if (iprog.eq.2.or.iprog.eq.3) then
               read(iudatin) uuu_sfan,
     &            (((uuu(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
               igotituuu_sfan=1
            else
               read(iudatin) uuu
            endif
            igotituuu=1
         elseif (name.eq.'V        ') then
            if (iprog.eq.2.or.iprog.eq.3) then
               read(iudatin) vvv_sfan,
     &            (((vvv(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
               igotitvvv_sfan=1
            else
               read(iudatin) vvv
            endif
            igotitvvv=1
         elseif (name.eq.'T        ') then
            if (iprog.eq.2.or.iprog.eq.3) then
               read(iudatin) tmk_sfan,
     &            (((tmk(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
               igotittmk_sfan=1
            else
               read(iudatin) tmk
            endif
            igotittmk=1
         elseif (name.eq.'Q        ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qvp in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qvp
            igotitqvp=1
         elseif (name.eq.'RH       ') then
            if (iprog.ne.2.and.iprog.ne.3) then
               print*,'RIP only expects RH in prs-level data.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qvp_sfan,
     &            (((qvp(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
            igotitrh_sfan=1
            igotitrh=1
         elseif (name.eq.'CLW      ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qcw in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qcw
            igotitqcw=1
         elseif (name.eq.'RNW      ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qra in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qra
            igotitqra=1
         elseif (name.eq.'ICE      ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qci in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qci
            igotitqci=1
         elseif (name.eq.'SNOW     ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qsn in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qsn
            igotitqsn=1
         elseif (name.eq.'GRAUPEL  ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects qgr in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) qgr
            igotitqgr=1
         elseif (name.eq.'NCI      ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects rnci in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) rnci
            igotitrnci=1
         elseif (name.eq.'W        ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,'RIP only expects www in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) www
            igotitwww=1
         elseif (name.eq.'PP       ') then
            if (iprog.ne.5.and.iprog.ne.11) then
               print*,
     &          'RIP only expects prs pert. in model input or output.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) prs
            igotitprs=1
         elseif (name.eq.'H        ') then
            if (iprog.ne.2.and.iprog.ne.3) then
               print*,'RIP only expects geop. hgt. in prs-level data.'
               print*,'Can''t handle it here.  Stopping.'
               stop
            endif
            read(iudatin) scr2,    ! surface ght is same as ter
     &            (((ght(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
            igotitght=1
         elseif (name.eq.'PSTARCRS ') then
            read(iudatin) pstx
            igotitpstx=1
         elseif (name.eq.'PSEALVLC ') then
            read(iudatin) slp
            igotitslp=1
         elseif (name.eq.'GROUND T ') then
            read(iudatin) tgk
            igotittgk=1
         elseif (name.eq.'TSEASFC  ') then
            read(iudatin) sst
            igotitsst=1
         elseif (name.eq.'RAIN CON ') then
            read(iudatin) rtc
            igotitrtc=1
         elseif (name.eq.'RAIN NON ') then
            read(iudatin) rte
            igotitrte=1
         elseif (name.eq.'TERRAIN  ') then
ccc            if (iprog.eq.2.or.iprog.eq.3) then
ccc               read(iudatin) ter_tsf
ccc               igotitter_tsf=1
ccc            else
ccc               read(iudatin) ter
ccc               igotitter=1
ccc            endif
            read(iudatin) ter
            igotitter=1
         elseif (name.eq.'MAPFACCR ') then
            read(iudatin) xmap
            igotitxmap=1
         elseif (name.eq.'MAPFACDT '.or.name.eq.'MAPFADOT ') then
            read(iudatin) dmap
            igotitdmap=1
         elseif (name.eq.'LATITCRS ') then
            read(iudatin) xlat
            igotitxlat=1
         elseif (name.eq.'LONGICRS ') then
            read(iudatin) xlon
            igotitxlon=1
         elseif (name.eq.'CORIOLIS ') then
            read(iudatin) dotcor
            igotitdotcor=1
         elseif (name.eq.'LAND USE ') then
            read(iudatin) xlus
            igotitxlus=1
         elseif (name.eq.'SNOWCOVR ') then
            read(iudatin) sno
            igotitsno=1
         elseif (name.eq.'PBL HGT  ') then
            read(iudatin) pblh
            igotitpblh=1
         elseif (name.eq.'REGIME   ') then
            read(iudatin) regime
            igotitregime=1
         elseif (name.eq.'SHFLUX   ') then
            read(iudatin) sshflux
            igotitsshflux=1
         elseif (name.eq.'LHFLUX   ') then
            read(iudatin) slhflux
            igotitslhflux=1
         elseif (name.eq.'MAVAIL   ') then
            read(iudatin) xmav
            igotitxmav=1
         elseif (name.eq.'UST      ') then
            read(iudatin) ust
            igotitust=1
         elseif (name.eq.'SWDOWN   ') then
            read(iudatin) swdown
            igotitswdown=1
         elseif (name.eq.'LWDOWN   ') then
            read(iudatin) lwdown
            igotitlwdown=1
         elseif (name.eq.'SOIL T 1 ') then
            read(iudatin) soil1
            igotitsoil1=1
         elseif (name.eq.'SOIL T 2 ') then
            read(iudatin) soil2
            igotitsoil2=1
         elseif (name.eq.'SOIL T 3 ') then
            read(iudatin) soil3
            igotitsoil3=1
         elseif (name.eq.'SOIL T 4 ') then
            read(iudatin) soil4
            igotitsoil4=1
         elseif (name.eq.'SOIL T 5 ') then
            read(iudatin) soil5
            igotitsoil5=1
         elseif (name.eq.'SOIL T 6 ') then
            read(iudatin) soil6
            igotitsoil6=1
         elseif (ordering.eq.'YXS '.or.ordering.eq.'YXW ') then
c
c         unknown 3d field
c
c         Note: V3 header: name*9, units*25, description*46
c              RIP header: varname*10,vardesc*64, plchun*24
c
            if (iprog.eq.2.or.iprog.eq.3) then
               read(iudatin) scr2,
     &            (((scr3(i,j,k),i=1,miy),j=1,mjx),k=mkzh,1,-1)
            elseif (end_index(3).eq.mkzh+1) then
               read(iudatin) scr3f
            else
               read(iudatin) scr3
            endif
            varname=name
            inname=0
            do ic=10,1,-1
               if (inname.eq.0) then
                  if (varname(ic:ic).ne.' ') inname=1
               else
                  if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
               endif
            enddo
            vardesc=description//', '//units ! Last 9 chars. of "units"
c                                              get cut off (oh well)
            plchun=units(1:24) ! last character of "units" gets cut off,
c                                and it's not in plotchar format (oh well)
            icd=1
            if (staggering.eq.'D   ') icd=0
            if (end_index(3).eq.mkzh+1) then ! interpolate to half sigma levels
               do k=1,mkzh
               do j=1,mjx-icd
               do i=1,miy-icd
                  scr3(i,j,k)=.5*(scr3f(i,j,k)+scr3f(i,j,k+1))
               enddo
               enddo
               enddo
               end_index(3)=mkzh
            endif
            if (end_index(3).eq.mkzh) then
               call writefile_rdp(scr3,varname,3,icd,vardesc,plchun,
     &            fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            endif
c
c         Write out surface part from pressure level data.
c
            if (iprog.eq.2.or.iprog.eq.3) then
               iendvarname=index(varname,' ')-1
               if (iendvarname.eq.-1) iendvarname=10
               iendvarname=min(iendvarname,8)
               varname=varname(1:iendvarname)//'_s'
c
c            Last 13 chars. of "units" get cut off (oh well)
c
               vardesc='SFC '//description//', '//units
               plchun=units(1:24) ! last character of "units" gets cut off,
c                                and it's not in plotchar format (oh well)
               icd=1
               if (staggering.eq.'D   ') icd=0
               call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &            fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            endif
         elseif (ordering.eq.'YX  ') then ! unknown 2d field
            read(iudatin) scr2
            varname=name
            inname=0
            do ic=10,1,-1
               if (inname.eq.0) then
                  if (varname(ic:ic).ne.' ') inname=1
               else
                  if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
               endif
            enddo
            vardesc=description//', '//units ! Last 9 chars. of "units"
c                                              get cut off (oh well)
            plchun=units(1:24) ! last character of "units" gets cut off,
c                                and it's not in plotchar format (oh well)
            icd=1
            if (staggering.eq.'D   ') icd=0
            call writefile_rdp(scr2,varname,2,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
         elseif (name.eq.'SIGMAH   ') then
            read(iudatin)
            idiscard=1
            print*,'   Discarding 1-D SIGMAH array --',
     &         ' RIPDP read it earlier.'
         elseif (name.eq.'PRESSURE ') then
            read(iudatin)
            idiscard=1
            if (iprog.le.3) then
               print*,'   Discarding 1-D PRESSURE array --',
     &            ' RIPDP read it earlier.'
            else
               print*,'   Discarding 1-D PRESSURE array --',
     &            ' not needed for sigma-level data.'
            endif
         else
            read(iudatin)
            idiscard=1
            print*,'   Discarding ',name,
     &             '   because it is not 2- or 3-d data.'
         endif
 353     continue
         if (idiscard.eq.0) then
            print*,'Processing MM5 variable ',name
         endif
         goto 421
      elseif (flag.eq.2) then
         continue
      else
         print*,'Ran into a flag not =1 or 2 while reading data.'
         stop
      endif
c
ccccccccccccc temp tend
c      read(68,end=229) msg,inest,xtime_ttend
c      read(68) t3dten
c      print*,'temperature tendency from CPS read in, time=',
c     &   xtime_ttend,' minutes.'
c      do k=1,mkzh
c      do j=1,mjx-icd
c      do i=1,miy-icd
c         scr3(i,j,k)=t3dten(i,j,k)*3600.
c      enddo
c      enddo
c      enddo
c      vardesc='T tendency from CPS, K/h'
c      plchun='K h~S~-1~N~'
c      call writefile_rdp(scr3,'ttendcps  ',3,1,vardesc,plchun,
c     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
cc
c      read(68) msg,inest,xtime_ttend
c      read(68) t3dten
c      print*,'temperature tendency from RAD read in, time=',
c     &   xtime_ttend,' minutes.'
c      do k=1,mkzh
c      do j=1,mjx-icd
c      do i=1,miy-icd
c         scr3(i,j,k)=t3dten(i,j,k)*3600.
c      enddo
c      enddo
c      enddo
c      vardesc='T tendency from RAD, K/h'
c      plchun='K h~S~-1~N~'
c      call writefile_rdp(scr3,'ttendrad  ',3,1,vardesc,plchun,
c     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
cc
c      read(68) msg,inest,xtime_ttend
c      read(68) t3dten
c      print*,'temperature tendency from BMP read in, time=',
c     &   xtime_ttend,' minutes.'
c      do k=1,mkzh
c      do j=1,mjx-icd
c      do i=1,miy-icd
c         scr3(i,j,k)=t3dten(i,j,k)*3600.
c      enddo
c      enddo
c      enddo
c      vardesc='T tendency from BMP, K/h'
c      plchun='K h~S~-1~N~'
c      call writefile_rdp(scr3,'ttendbmp  ',3,1,vardesc,plchun,
c     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
c     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c 229  continue
ccccccccccccc temp tend
c
      endif
c
c   Read in der. q. data variables.
c
      if (ipvdo.eq.1) then
         call fillarray(qqn,miy*mjx*mkzh*nvarq,0.)
         call fillarray(tqn,miy*mjx*mkzh*nvartq,0.)
         call fillarray(uqn,miy*mjx*mkzh*nvarvq,0.)
         call fillarray(vqn,miy*mjx*mkzh*nvarvq,0.)
         do ivq=1,nvarq
            read(iuqin) (((qqn(i,j,k,ivq),i=miyqb,miyqe),
     &                     j=mjxqb,mjxqe),k=1,mkzh)
            do k=1,mkzh
               do j=mjxqb,mjxqe
                  qqn(miyqe,j,k,ivq)=0.
               enddo
               do i=miyqb,miyqe-1
                  qqn(i,mjxqe,k,ivq)=0.
               enddo
            enddo
         enddo
         do ivq=1,nvartq
            read(iuqin) (((tqn(i,j,k,ivq),i=miyqb,miyqe),
     &                     j=mjxqb,mjxqe),k=1,mkzh)
            do k=1,mkzh
               do j=mjxqb,mjxqe
                  tqn(miyqe,j,k,ivq)=0.
               enddo
               do i=miyqb,miyqe-1
                  tqn(i,mjxqe,k,ivq)=0.
               enddo
            enddo
         enddo
c         if (xtime.gt.6.5.and.xtime.lt.30.5) then
c            ntq4av=ntq4av+1
c            do k=1,mkzh
c            do j=mjxqb,mjxqe-1
c            do i=miyqb,miyqe-1
c               tq4av(i,j,k)=tq4av(i,j,k)+tqn(i,j,k,5)
c            enddo
c            enddo
c            enddo
c         endif
         do ivq=1,nvarvq
            read(iuqin) (((uqn(i,j,k,ivq),i=miyqb,miyqe),
     &                     j=mjxqb,mjxqe),k=1,mkzh)
         enddo
         do ivq=1,nvarvq
            read(iuqin) (((vqn(i,j,k,ivq),i=miyqb,miyqe),
     &                     j=mjxqb,mjxqe),k=1,mkzh)
         enddo
         call jumpendfile(iuqin,1)
      endif
c
c   Read in observations data
c
      if (iobsprc.eq.1.and.iobsprcfnd.eq.1) then
         read(iuobsprc) rrbo
         vardesc='Obs. precip. since last output time, mm'
         plchun='mm'
         call writefile_rdp(rrbo,'rrbo      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (iobscnv.eq.1.and.iobsclgfnd.eq.1) then
         read(iuobsclg) clgo
         vardesc='Obs. cloud ceiling, feet'
         plchun='ft'
         call writefile_rdp(clgo,'clgo      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (iobscnv.eq.1.and.iobsvisfnd.eq.1) then
         read(iuobsvis) viso
         vardesc='Obs. visibility, miles'
         plchun='mi'
         call writefile_rdp(viso,'viso      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (iftcnv.eq.1.and.iftclgfnd.eq.1) then
         read(iuftclg) clgf
         vardesc='FT-fcst. cloud ceiling, feet'
         plchun='ft'
         call writefile_rdp(clgf,'clgf      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (iftcnv.eq.1.and.iftvisfnd.eq.1) then
         read(iuftvis) visf
         vardesc='FT-fcst. visibility, miles'
         plchun='mi'
         call writefile_rdp(visf,'visf      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
      print*,'   Finished reading data for this time.'
c
c   Write out some fields that are OK as is.
c
      if (igotitter.eq.1) then
         plchun='m'
         vardesc='Terrain height AMSL, m'
         call writefile_rdp(ter,'ter       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
ccc      if (igotitter_tsf.eq.1) then
ccc         plchun='m'
ccc         vardesc='Terrain height (of true sfc) AMSL, m'
ccc         call writefile_rdp(ter_tsf,'ter       ',2,1,vardesc,plchun,
ccc     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
ccc     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
ccc      endif
      if (igotitdmap.eq.1) then
         vardesc='Map factor on dot points'
         plchun='none'
         call writefile_rdp(dmap,'dmap      ',2,0,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitxmap.eq.1) then
         vardesc='Map factor on cross points'
         plchun='none'
         call writefile_rdp(xmap,'xmap      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitxlat.eq.1) then
         vardesc='Latitude on cross points'
         plchun='degrees'
         call writefile_rdp(xlat,'xlat      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitxlon.eq.1) then
         vardesc='Longitude on cross points'
         plchun='degrees'
         call writefile_rdp(xlon,'xlon      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   Fill outer row and column of xlus with values from next
c      inner row or column.
c
      if (igotitxlus.eq.1) then
         do j=2,mjx-2
            xlus(1,j)=xlus(2,j)
            xlus(miy-1,j)=xlus(miy-2,j)
         enddo
         do i=1,miy-1
            xlus(i,1)=xlus(i,2)
            xlus(i,mjx-1)=xlus(i,mjx-2)
         enddo
         vardesc='Land use category'
         plchun='none'
         call writefile_rdp(xlus,'xlus      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   Create cor.
c
      if (igotitdotcor.eq.1) then
         do j=1,mjx
         do i=1,miy
            cor(i,j)=0.
         enddo
         enddo
         do j=1,mjx-1
         do i=1,miy-1
            cor(i,j)=.25*(dotcor(i,j)+dotcor(i+1,j)+
     &                    dotcor(i,j+1)+dotcor(i+1,j+1))
         enddo
         enddo
         vardesc='Coriolis parameter, per s'
         plchun='s~S~-1~N~'
         call writefile_rdp(cor,'cor       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   May need to combine ground temperature and SST
c
      if (igotittgk.eq.1.or.igotitsst.eq.1) then
         if (igotittgk.eq.0) then
            do j=1,mjx-1
            do i=1,miy-1
               tgk(i,j)=sst(i,j)
            enddo
            enddo
         elseif (igotittgk.eq.1.and.igotitsst.eq.1) then
            if (igotitxlus.eq.1) then
               if (ilandset.eq.1) then
                  iwater=7
               elseif (ilandset.eq.2) then
                  iwater=16
               elseif (ilandset.eq.3) then
                  iwater=15
               endif
               do j=1,mjx-1
               do i=1,miy-1
                  if (nint(xlus(i,j)).eq.iwater) tgk(i,j)=sst(i,j)
               enddo
               enddo
            else
               print*,'Cannot discriminate between ground temp. and'
               print*,'SST without land use information.'
               stop
            endif
         endif
         vardesc='Ground/sea-surface temperature, K'
         plchun='K'
         call writefile_rdp(tgk,'tgk       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   SNO might either be either snow cover (1=snow,0=no snow) or
c   snow depth in mm.
c
      if (igotitsno.eq.1) then
         inotzeroorone=0
         do j=2,mjx-2
         do i=2,miy-2
            icheck=nint(sno(i,j)*10.)
            if (icheck.gt.0.and.icheck.ne.10.and.icheck.ne.5)
     &          inotzeroorone=inotzeroorone+1
         enddo
         enddo
         if (inotzeroorone.gt.0) then ! it's snow depth
            varname='snod      '
            vardesc='Snow depth, mm'
            plchun='mm'
         else                         ! it's snow cover
            do j=2,mjx-2
            do i=2,miy-2
               if (sno(i,j).ne.1.) sno(i,j)=0.
            enddo
            enddo
            varname='sno       '
            vardesc='Snow cover (0.0 or 1.0)'
            plchun='none'
         endif
         call writefile_rdp(sno,varname,2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   If we got SLP, convert it to hPa and write out      
c
      if (igotitslp.eq.1) then
         vardesc='Sea-level pressure, hPa'
         plchun='hPa'
         do j=1,mjx-1
         do i=1,miy-1
            slp(i,j)=.01*slp(i,j)
         enddo
         enddo
         call writefile_rdp(slp,'slp_sfan  ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitpblh.eq.1) then
         vardesc='PBL height, m'
         plchun='m'
         call writefile_rdp(pblh,'pblh      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitregime.eq.1) then
         vardesc='PBL regime (1,2,3,or 4)'
         plchun='none'
         call writefile_rdp(regime,'regime    ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsshflux.eq.1) then
         vardesc='Surface sensible heat flux, W/m**2'
         plchun='W m~S~-2~N~'
         call writefile_rdp(sshflux,'sshflux   ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitslhflux.eq.1) then
         vardesc='Surface latent heat flux, W/m**2'
         plchun='W m~S~-2~N~'
         call writefile_rdp(slhflux,'slhflux   ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitxmav.eq.1) then
         vardesc='Soil moisture availability, %'
         plchun='%'
         call writefile_rdp(xmav,'mava      ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitust.eq.1) then
         vardesc='Friction velocity, m/s'
         plchun='m s~S~-1~N~'
         call writefile_rdp(ust,'ust       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitswdown.eq.1) then
         vardesc='Surface downward shortwave radiation, W/m**2'
         plchun='W m~S~-2~N~'
         call writefile_rdp(swdown,'swdown    ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitlwdown.eq.1) then
         vardesc='Surface downward longwave radiation, W/m**2'
         plchun='W m~S~-2~N~'
         call writefile_rdp(lwdown,'lwdown    ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil1.eq.1) then
         vardesc='Soil temperature in layer 1, K'
         plchun='K'
         call writefile_rdp(soil1,'soil1     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil2.eq.1) then
         vardesc='Soil temperature in layer 2, K'
         plchun='K'
         call writefile_rdp(soil2,'soil2     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil3.eq.1) then
         vardesc='Soil temperature in layer 3, K'
         plchun='K'
         call writefile_rdp(soil3,'soil3     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil4.eq.1) then
         vardesc='Soil temperature in layer 4, K'
         plchun='K'
         call writefile_rdp(soil4,'soil4     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil5.eq.1) then
         vardesc='Soil temperature in layer 5, K'
         plchun='K'
         call writefile_rdp(soil5,'soil5     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitsoil6.eq.1) then
         vardesc='Soil temperature in layer 6, K'
         plchun='K'
         call writefile_rdp(soil6,'soil6     ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   Create pstd
c
      if (iprog.eq.1.or.iprog.eq.2.or.iprog.eq.3) then
c
c      Make pstx in Pa, to be consistent with pstx that is read in
c      from V3 sigma-level output
c
         do j=1,mjx
         do i=1,miy
            pstx(i,j)=100.*pstarconst
         enddo
         enddo
         igotitpstx=1
      endif
      if (igotitpstx.eq.1) then
         do j=1,mjx
         do i=1,miy
            iph=min(i,miy-1)
            jph=min(j,mjx-1)
            imh=max(i-1,1)
            jmh=max(j-1,1)
            pstd(i,j)=.25*(pstx(iph,jph)+pstx(imh,jph)+
     &                     pstx(iph,jmh)+pstx(imh,jmh))
         enddo
         enddo
      else
         print*,'Didn''t find pstx.  Stopping'
         stop
      endif
c
c   Process and write out other variables
c
      if (igotituuu.eq.1) then
         if (dataform.ne.'mm5v3   ')
     &      call decouple(uuu,pstd,0,1.,miy,mjx,mkzh)
         vardesc='Horizontal wind (x-comp.), m/s'
         plchun='m s~S~-1~N~'
c diff test
c         do k=1,mkzh
c         do j=1,mjx
c         do i=1,miy
c            uuu(i,j,k)=uuu(i,j,k)+(i+j)/10.
c         enddo
c         enddo
c         enddo
         call writefile_rdp(uuu,'uuu       ',3,0,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotituuu_sfan.eq.1) then
         vardesc='Hor. wind (sfc. anal.) (x-comp.), m/s'
         plchun='m s~S~-1~N~'
         call writefile_rdp(uuu_sfan,'uuu_sfan  ',2,0,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitvvv.eq.1) then
         if (dataform.ne.'mm5v3   ')
     &      call decouple(vvv,pstd,0,1.,miy,mjx,mkzh)
         vardesc='Horizontal wind (y-comp.), m/s'
         plchun='m s~S~-1~N~'
         call writefile_rdp(vvv,'vvv       ',3,0,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitvvv_sfan.eq.1) then
         vardesc='Hor. wind (sfc. anal.) (y-comp.), m/s'
         plchun='m s~S~-1~N~'
         call writefile_rdp(vvv_sfan,'vvv_sfan  ',2,0,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotittmk.eq.1) then
         if (dataform.ne.'mm5v3   ')
     &      call decouple(tmk,pstx,1,1.,miy,mjx,mkzh)
         vardesc='Temperature, K'
         plchun='K'
         call writefile_rdp(tmk,'tmk       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotittmk_sfan.eq.1) then
         vardesc='Temperature (sfc. anal.), K'
         plchun='K'
         call writefile_rdp(tmk_sfan,'tmk_sfan  ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitqvp.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qvp,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qvp(i,j,k)=1000.*qvp(i,j,k)! convert to g/kg
            enddo
            enddo
            enddo
         endif
         vardesc='Water vapor mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qvp,'qvp       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitrh.eq.1) then
         vardesc='Relative humidity (w.r.t. water), %'
         plchun='%'
         call writefile_rdp(qvp,'rhu       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c      Convert to mixing ratio (g/kg) and write out again
c
         if ((iprog.eq.2.or.iprog.eq.3).and.igotittmk.eq.1) then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               es=ezero*exp(eslcon1*(tmk(i,j,k)-celkel)/
     &            (tmk(i,j,k)-eslcon2))
               ws=eps*es/(prslvl(k)-es)  ! in kg/kg
c
c            Convert qvp from RH in % to mix rat in g/kg
c
               qvp(i,j,k)=10.*qvp(i,j,k)*ws
c diff test
c               qvp(i,j,k)=qvp(i,j,k)+2.
            enddo
            enddo
            enddo
         endif
         vardesc='Water vapor mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qvp,'qvp       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
         igotitqvp=1
      endif
      if (igotitqcw.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qcw,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qcw(i,j,k)=1000.*qcw(i,j,k)! convert to g/kg
            enddo
            enddo
            enddo
         endif
         vardesc='Cloud water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qcw,'qcw       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitqra.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qra,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qra(i,j,k)=1000.*qra(i,j,k)! convert to g/kg
            enddo
            enddo
            enddo
         endif
         vardesc='Rain water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qra,'qra       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitwww.eq.1) then
c
c      Decouple www, interpolate it to half sigma levels
c      and convert it to cm/s.
c
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (dataform.ne.'mm5v3   ') then
               www(i,j,k)=50.*(www(i,j,k)+www(i,j,k+1))/pstx(i,j)
            else
               www(i,j,k)=50.*(www(i,j,k)+www(i,j,k+1))
            endif
         enddo
         enddo
         enddo
         vardesc='Vertical velocity, cm/s'
         plchun='cm s~S~-1~N~'
         call writefile_rdp(www,'www       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
      if (iprog.eq.5.or.iprog.eq.6.or.iprog.eq.11) then
c
c   Note, in the above "if" statement,
c   that pressure will not be processed if iprog=1,2 or 3.  If
c   iprog=2 or 3, pressure can be easily created in rip, since
c   data levels are on const. pressure surfaces.
ccc
ccc  New note: 3D pressure array is now created and written for
ccc  pressure-level data, below.
ccc
c
      if ((inhyd.eq.1.and.igotitprs.eq.1.and.igotitpstx.eq.1).or.
     &    (inhyd.eq.0.and.igotitprs.eq.0.and.igotitpstx.eq.1)) then
c
c      Process pressure, and write it out.
c
         igotitprs=1
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (dataform.ne.'mm5v3   ') then
c
c            pstx is still in kPa here.
c            prs is either zero (if inhyd=0), or is pressure pert.
c            (in Pa), coupled to pstx (in kPa).
c
               prs(i,j,k)=sigh(k)*pstx(i,j)*10.+ptop+
     &            .01*prs(i,j,k)/pstx(i,j)  ! want in hPa
            else
c
c            pstx is still in Pa here.
c            prs is pressure pert. (in Pa), uncoupled
c
               prs(i,j,k)=sigh(k)*pstx(i,j)*.01+ptop+
     &            .01*prs(i,j,k)  ! want in hPa
            endif
         enddo
         enddo
         enddo
         vardesc='Pressure, hPa'
         plchun='hPa'
         call writefile_rdp(prs,'prs       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
ccc
ccc      Process surface pressure, and write it out.
ccc
         igotitsfp=1
         do j=1,mjx-1
         do i=1,miy-1
c           pstx is still in Pa here.
            sfp(i,j)=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)*.01
         enddo
         enddo
         vardesc='Surface pressure, hPa'
         plchun='hPa'
         call writefile_rdp(sfp,'sfp       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      elseif (inhyd.eq.1) then
         print*,'That is odd:  data is supposedly nonhydrostatic,'
         print*,'but no pressure perturbation array was found.'
         stop
      elseif (inhyd.eq.0) then
         print*,'That is odd:  data is supposedly hydrostatic,'
         print*,'but a pressure perturbation array was found.'
         stop
      else
         print*,'Something odd about pressure situation.'
         print*,'inhyd,igotitprs,igotitpstx=',inhyd,igotitprs,igotitpstx
         print*,'dataform=',dataform
         stop
      endif
ccc
ccc   Make vertical velocity if hyrdostatic data.
ccc
      if (inhyd.eq.0) then
         call wcalchyd(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,miy,mjx,mkzh,ds,rgas,grav,ptop)
         vardesc='Vertical velocity, cm/s'
         plchun='cm s~S~-1~N~'
         call writefile_rdp(www,'www       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
      elseif (iprog.eq.2.or.iprog.eq.3) then
c
      igotitprs=1
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         prs(i,j,k)=prslvl(k)
      enddo
      enddo
      enddo
      vardesc='Pressure, hPa'
      plchun='hPa'
      call writefile_rdp(prs,'prs       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      endif
c
      if (igotitqci.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qci,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qci(i,j,k)=1000.*qci(i,j,k)
            enddo
            enddo
            enddo
         endif
         vardesc='Cloud ice mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qci,'qci       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitqsn.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qsn,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qsn(i,j,k)=1000.*qsn(i,j,k)
            enddo
            enddo
            enddo
         endif
         vardesc='Snow mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qsn,'qsn       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitqgr.eq.1) then
         if (dataform.ne.'mm5v3   ') then
            call decouple(qgr,pstx,1,1000.,miy,mjx,mkzh)! convert to g/kg
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qgr(i,j,k)=1000.*qgr(i,j,k)
            enddo
            enddo
            enddo
         endif
         vardesc='Graupel mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile_rdp(qgr,'qgr       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitrnci.eq.1) then
         if (dataform.ne.'mm5v3   ')
     &      call decouple(rnci,pstx,1,1.,miy,mjx,mkzh)
         vardesc='Number concentration of cloud ice, #/m**3'
         plchun='m~S~-3~N~'
         call writefile_rdp(rnci,'rnci      ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotittke.eq.1) then
c        Note: no need to decouple - this var. not coupled with pstar
         vardesc='Turbulent kinetic energy, J/kg'
         plchun='J kg~S~-1~N~'
         call writefile_rdp(tke,'tke       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitradtnd.eq.1) then
c        Note: no need to decouple - this var. not coupled with pstar
         vardesc='Radiation tendency, K/day'
         plchun='K day~S~-1~N~'
         call writefile_rdp(radtnd,'radtnd    ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
      if (ipvdo.eq.1) then
c
c      Decouple pv-related variables and convert pv to PVU,
c      and convert temp. partitions to theta.
c
         rgas=287.04
         pvc=1.e6
         if (nvq.gt.1) then
            varname='qq'
            do ivq=1,nvq
            do k=1,mkzh
            do j=mjxqb,mjxqe-1
            do i=miyqb,miyqe-1
               qqn(i,j,k,ivq)=qqn(i,j,k,ivq)/pstx(i,j)*pvc
            enddo
            enddo
            enddo
               if (ivq.eq.1) vardesc=
     &            'Conserved PV partition, PVU'
               if (ivq.eq.2) vardesc=
     &            'PV partition due to thermal diffusion, PVU'
               if (ivq.eq.3) vardesc=
     &            'PV partition due to explicit latent heating, PVU'
               if (ivq.eq.4) vardesc=
     &            'PV partition due to cumulus parameterization, PVU'
               if (ivq.eq.5) vardesc=
     &            'PV partition due to surface heat fluxes, PVU'
               if (ivq.eq.6) vardesc=
     &            'PV partition due to momentum diffusion, PVU'
               if (ivq.eq.7) vardesc=
     &            'PV partition due to surface momentum fluxes'//
     &            ' (friction), PVU'
               write(varname(3:3),'(i1)') ivq-1
               plchun='PVU'
               call writefile_rdp(qqn(1,1,1,ivq),varname,3,1,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            enddo
         endif
         if (nvtq.gt.1) then
            varname='tq'
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3(i,j,k)=(1000./prs(i,j,k))**.286
            enddo
            enddo
            enddo
            do ivq=1,nvtq
            do k=1,mkzh
            do j=mjxqb,mjxqe-1
            do i=miyqb,miyqe-1
               tqn(i,j,k,ivq)=tqn(i,j,k,ivq)/pstx(i,j)*
     &            scr3(i,j,k)
            enddo
            enddo
            enddo
               if (ivq.eq.1) vardesc=
     &            'Conserved theta partition, K'
               if (ivq.eq.2) vardesc=
     &            'Theta partition due to thermal diffusion, K'
               if (ivq.eq.3) vardesc=
     &            'Theta partition due to explicit latent heating, K'
               if (ivq.eq.4) vardesc=
     &            'Theta partition due to cumulus parameterization, K'
               if (ivq.eq.5) vardesc=
     &            'Theta partition due to surface heat fluxes, K'
               write(varname(3:3),'(i1)') ivq-1
               plchun='K'
               call writefile_rdp(tqn(1,1,1,ivq),varname,3,1,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            enddo
         endif
         if (nvvq.gt.1) then
            varname='uq'
            do ivq=1,nvvq
            do k=1,mkzh
            do j=mjxqb,mjxqe
            do i=miyqb,miyqe
               uqn(i,j,k,ivq)=uqn(i,j,k,ivq)/pstd(i,j)
            enddo
            enddo
            enddo
               if (ivq.eq.1) vardesc=
     &            'Conserved u-momentum partition, m/s'
               if (ivq.eq.2) vardesc=
     &            'U-momentum partition due to momentum diffusion, m/s'
               if (ivq.eq.3) vardesc=
     &            'U-momentum partition due to surf. mom. fluxes'//
     &            ' (friction), m/s'
               write(varname(3:3),'(i1)') ivq-1
               plchun='m s~S~-1~N~'
               call writefile_rdp(uqn(1,1,1,ivq),varname,3,0,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            enddo
            varname='vq'
            do ivq=1,nvvq
            do k=1,mkzh
            do j=mjxqb,mjxqe
            do i=miyqb,miyqe
               vqn(i,j,k,ivq)=vqn(i,j,k,ivq)/pstd(i,j)
            enddo
            enddo
            enddo
               if (ivq.eq.1) vardesc=
     &            'Conserved v-momentum partition, m/s'
               if (ivq.eq.2) vardesc=
     &            'V-momentum partition due to momentum diffusion, m/s'
               if (ivq.eq.3) vardesc=
     &            'V-momentum partition due to surf. mom. fluxes'//
     &            ' (friction), m/s'
               write(varname(3:3),'(i1)') ivq-1
               plchun='m s~S~-1~N~'
               call writefile_rdp(vqn(1,1,1,ivq),varname,3,0,vardesc,
     &            plchun,fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &            iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            enddo
         endif
      endif
      if (igotitpstx.eq.1) then
         do j=1,mjx-1
         do i=1,miy-1
            if (dataform.ne.'mm5v3   ') then
               pstx(i,j)=pstx(i,j)*10.  ! convert to hPa before writing
            else
               pstx(i,j)=pstx(i,j)*.01  ! convert to hPa before writing
            endif
         enddo
         enddo
         vardesc='P-star on cross points, hPa'
         plchun='hPa'
ccc         call writefile_rdp(pstx,'pstx      ',2,1,vardesc,plchun,
ccc     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
ccc     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitpstx.eq.1) then
         do j=1,mjx
         do i=1,miy
            if (dataform.ne.'mm5v3   ') then
               pstd(i,j)=pstd(i,j)*10.  ! convert to hPa before writing
            else
               pstd(i,j)=pstd(i,j)*.01  ! convert to hPa before writing
            endif
         enddo
         enddo
         vardesc='P-star on dot points, hPa'
         plchun='hPa'
ccc         call writefile_rdp(pstd,'pstd      ',2,0,vardesc,plchun,
ccc     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
ccc     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (iprog.gt.3.and.igotitpstx.eq.1.and.igotitprs.eq.1.and.
     &    igotittmk.eq.1.and.igotitter.eq.1) then
c
c      Calculate ght.
c
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            ght(i,j,k)=0.0
         enddo
         enddo
         enddo
         if (igotitqvp.ne.1) then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               qvp(i,j,k)=0.0
            enddo
            enddo
            enddo
         endif
         call ghtcalc(sigh,pstx,prs,qvp,tmk,ter,ght,rgas,ussalr,
     &      grav,refslp,refslt,reflaps,refstratt,ptop,inhyd,
     &      miy,mjx,mkzh)
         igotitght=1
      endif
      if (igotitght.eq.1) then
         vardesc='Geopotential height, m'
         plchun='m'
         call writefile_rdp(ght,'ght       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
ccc
      expon = rgas*ussalr/grav
      exponi = 1./expon
c
      if (iprog.eq.1) then
ccc
ccc   Make "fake" fields for TERRAIN output (based on standard atmosphere).
ccc   Note, there are 2 model levels in the 3D "fake" fields.  These fields
ccc   are created just so RIP will function properly using this data.
ccc
      do j=1,mjx-1
      do i=1,miy-1
         prs(i,j,1)=1000.
         prs(i,j,2)=900.
         ght(i,j,1)=288.15/ussalr*(1.-(prs(i,j,1)/1013.25)**expon)
         ght(i,j,2)=288.15/ussalr*(1.-(prs(i,j,2)/1013.25)**expon)
         sfp(i,j)=1013.25*(1.-ussalr/288.15*ter(i,j))**exponi
         tmk(i,j,1)=288.15-ussalr*ght(i,j,1)
         tmk(i,j,2)=288.15-ussalr*ght(i,j,2)
      enddo
      enddo
      do j=1,mjx
      do i=1,miy
         uuu(i,j,1)=1.
         uuu(i,j,2)=1.
         vvv(i,j,1)=1.
         vvv(i,j,2)=1.
      enddo
      enddo
c
      vardesc='Pressure, hPa'
      plchun='hPa'
      call writefile_rdp(prs,'prs       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Surface pressure, hPa'
      plchun='hPa'
      call writefile_rdp(sfp,'sfp       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Geopotential height, m'
      plchun='m'
      call writefile_rdp(ght,'ght       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Temperature, K'
      plchun='K'
      call writefile_rdp(tmk,'tmk       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Horizontal wind (x-comp.), m/s'
      plchun='m s~S~-1~N~'
      call writefile_rdp(uuu,'uuu       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Horizontal wind (y-comp.), m/s'
      plchun='m s~S~-1~N~'
      call writefile_rdp(vvv,'vvv       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      elseif (iprog.eq.2.or.iprog.eq.3) then
c
c      Special stuff to do for pressure level data.
c      First, get desired sigma level to use temperature from.
c
         sigc=(850.-ptop)/(1000.-ptop)
         do k=1,mkzh
            if (sigc.ge.sigf(k).and.sigc.le.sigf(k+1)) kupper=k
         enddo
c
         do j=1,mjx-1
         do i=1,miy-1
c
c      Create ter_pseudo (i.e., height at the pseudo-surface defined by p=pbot).
c
         ter_pseudo(i,j)=ght(i,j,mkzh)-rgas/grav*
     *      virtual(tmk(i,j,mkzh),.001*qvp(i,j,mkzh))*
     &      log(pbot/prslvl(mkzh))
         igotitter=1
c
ccc      Create sfp (i.e., true surface pressure, at hgt=ter).
c
         zlhsl=ght(i,j,mkzh)   ! lhsl => "lowest half sigma level"
         ezlhsl=exp(-zlhsl/sclht)
         plhsl=prslvl(mkzh)
         zpsf=ter_pseudo(i,j)   ! psf => "at the pseudo-surface"
         ezpsf=exp(-zpsf/sclht)
         ppsf=pbot
         ztsf=ter(i,j)   ! ccc tsf => "at the true surface"
         eztsf=exp(-ztsf/sclht)
c        ptsf = sfp(i,j) = ????  This is what we're after
c
c      First check if true surface is above top-most half sigma level.
c
         if (ztsf.gt.ght(i,j,1)) then
            print*,'True surface is above highest level of data.'
            print*,'i,j=',i,j
            print*,'This seems very strange.'
            print*,'ztsf,ght(i,j,1)=',ztsf,ght(i,j,1)
            stop
         endif
c
c      Check if true surface is somewhere within sigma levels.
c
         do k=mkzh,2,-1
            if (ztsf.ge.ght(i,j,k).and.
     &          ztsf.le.ght(i,j,k-1)) then
               ezk=exp(-ght(i,j,k)/sclht)
               ezkm1=exp(-ght(i,j,k-1)/sclht)
               sfp(i,j)=((ezk-eztsf)*prslvl(k-1)+
     &                       (eztsf-ezkm1)*prslvl(k))/(ezk-ezkm1)
               goto 716
            endif
         enddo
c
c      Otherwise, true surface is below the lowest half sigma level
c
         if (ztsf.ge.zpsf) then
c
c         True surface is below the lowest half sigma level
c         but above the pseudo-surface.
c
            sfp(i,j)=((ezpsf-eztsf)*plhsl+(eztsf-ezlhsl)*ppsf)/
     &         (ezpsf-ezlhsl)
c
         else
c
c         True surface is below the pseudo-surface.
c         Extrapolate using USSALR
c
            tpsfbg=tmk(i,j,kupper)*(ppsf/prslvl(kupper))**expon
            tpsfbgv=virtual(tpsfbg,.001*qvp(i,j,mkzh))
            sfp(i,j)=ppsf*(1.+ussalr/tpsfbgv*(zpsf-ztsf))**exponi
c
         endif
c
 716     continue
c
         enddo
         enddo
c
ccc         vardesc='Terrain height AMSL, m'
ccc         plchun='m'
ccc         call writefile_rdp(ter,'ter       ',2,1,vardesc,plchun,
ccc     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
ccc     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
ccc         vardesc='Pressure (at true sfc), hPa'
         vardesc='Surface pressure, hPa'
         plchun='hPa'
         call writefile_rdp(sfp,'sfp       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c      Also, now that we have sfp, we can process rhu_sfan and
c      qvp_sfan
c
         if (igotitrh_sfan.eq.1) then
            vardesc='Relative humidity (sfc. anal.), %'
            plchun='%'
            call writefile_rdp(qvp_sfan,'rhu_sfan  ',2,1,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c         Convert to mixing ratio (g/kg) and write out again
c
            if ((iprog.eq.2.or.iprog.eq.3).and.igotittmk_sfan.eq.1) then
               do j=1,mjx-1
               do i=1,miy-1
                  es=ezero*exp(eslcon1*(tmk_sfan(i,j)-celkel)/
     &               (tmk_sfan(i,j)-eslcon2))
                  ws=eps*es/(sfp(i,j)-es)  ! in kg/kg
c
c               Convert qvp from RH in % to mix rat in g/kg
c
                  qvp_sfan(i,j)=10.*qvp_sfan(i,j)*ws
               enddo
               enddo
            endif
            vardesc='Water vapor mixing ratio (sfc. anal.), g/kg'
            plchun='g kg~S~-1~N~'
            call writefile_rdp(qvp_sfan,'qvp_sfan  ',2,1,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
            igotitqvp_sfan=1
         endif
      endif
c
      if (igotitrtc.eq.1) then
         do j=1,mjx-1
         do i=1,miy-1
            rtc(i,j)=10.*rtc(i,j) ! cm to mm
         enddo
         enddo
         vardesc='Cumulus precip. since h 0, mm'
         plchun='mm'
         call writefile_rdp(rtc,'rtc       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
      if (igotitrte.eq.1) then
         do j=1,mjx-1
         do i=1,miy-1
            rte(i,j)=10.*rte(i,j) ! cm to mm
         enddo
         enddo
         vardesc='Explicit precip. since h 0, mm'
         plchun='mm'
         call writefile_rdp(rte,'rte       ',2,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
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
      goto 240       ! End of time loop.
c
 1000 continue
c
c      if (ipv.eq.1) then
c         print*,'ntq4av=',ntq4av
c         do k=1,mkzh
c         do j=mjxqb,mjxqe-1
c         do i=miyqb,miyqe-1
c            tq4av(i,j,k)=tq4av(i,j,k)/ntq4av
c         enddo
c         enddo
c         enddo
c         open(unit=67,file='tq4av.dat',form='unformatted',status='new')
c         write(67)(((tq4av(i,j,k),
c     &      i=miyqb,miyqe),j=mjxqb,mjxqe),k=1,mkzh)
c      endif
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
      subroutine decouple(array,pstar,icd,factor,miy,mjx,mkzh)
      dimension array(miy,mjx,mkzh),pstar(miy,mjx)
      if (factor.eq.1) then
         do k=1,mkzh
         do j=1,mjx-icd
         do i=1,miy-icd
            array(i,j,k)=array(i,j,k)/pstar(i,j)
         enddo
         enddo
         enddo
      else
         do k=1,mkzh
         do j=1,mjx-icd
         do i=1,miy-icd
            array(i,j,k)=array(i,j,k)/pstar(i,j)*factor
         enddo
         enddo
         enddo
      endif
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getobs(iobs,iobsfnd,mdate,iunit)
c
      if (iobs.eq.1) then
  280    read(iunit,end=285) mdategp
  283    goto 286
  285    print*,'   Couldn''t find matching obs dataset.'
         iobs=0
         iobsfnd=0
         goto 287
  286    if (mdategp.gt.mdate) then
            backspace (iunit)
            iobsfnd=0
            print*,'   Model output time is earlier than'
            print*,'      first obs data time.'
         elseif (mdategp.lt.mdate) then
            read(iunit)
            goto 280
         else
            iobsfnd=1
            print*,'   Found matching time in obs dataset.'
         endif
      endif
  287 continue
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ghtcalc(sigh,pstx,prs,qvp,tmk,ter,ght,rgas,ussalr,
     &      grav,refslp,refslt,reflaps,refstratt,ptop,inhyd,
     &      miy,mjx,mkzh)
c
      dimension sigh(mkzh),pstx(miy,mjx),
     &   qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),ter(miy,mjx),
     &   ght(miy,mjx,mkzh),prs(miy,mjx,mkzh)
c
      dimension tv(100)
c
      expon=rgas*ussalr/grav
      exponi=1./expon
c
      if (inhyd.eq.0) then  ! hydrostatic sigma-coord. form
c
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
c
c   Calculate tv
c
      do k=1,mkzh
         tv(k)=virtual(tmk(i,j,k),.001*qvp(i,j,k))
      enddo
c
c   Calculate geopotential height at lowest half sigma level,
c      assuming tv folows standard lapse rate below sigh(mkzh),
c      and using altimeter equation.
c
      psurf=pstx(i,j)+ptop
      ght(i,j,mkzh)=ter(i,j)+tv(mkzh)/ussalr*
     &   ((psurf/prs(i,j,mkzh))**expon - 1.)
c
      do k = mkzh-1,1,-1
         tvavg=.5*(tv(k)+tv(k+1))
         ght(i,j,k) = ght(i,j,k+1) + (rgas*tvavg/grav)*
     &      log(prs(i,j,k+1)/prs(i,j,k))
      enddo
c
 1000 continue
c
      else  ! nonhydrostatic sigma-coord. form
c
      cc1=rgas/grav*(-.5)*reflaps
      cc2=rgas/grav*(reflaps*log(.01*refslp)-refslt)
      cc3=rgas/grav*(refslt-.5*reflaps*log(.01*refslp))*log(.01*refslp)
c
      alnpreftpause=(refstratt-refslt)/reflaps+log(.01*refslp)
      ztpause=cc1*alnpreftpause*alnpreftpause+
     &   cc2*alnpreftpause+cc3
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         alnpref=log(sigh(k)*pstx(i,j)+ptop)
         if (alnpref.gt.alnpreftpause) then
            ght(i,j,k)=cc1*alnpref*alnpref+cc2*alnpref+cc3
         else
            ght(i,j,k)=ztpause+rgas*refstratt/grav*
     &         (alnpreftpause-alnpref)
         endif
      enddo
      enddo
      enddo
c
      endif
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine jumpendfile(iunit,ijump)
      if (ijump.le.0) return
      isofar=0
   10 continue
   20 continue
      read(iunit,end=100)
      goto 20
  100 continue
      isofar=isofar+1
      if (isofar.lt.ijump) goto 10
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine skpflds(iunit,mifv1,nskip)
      integer mifv1(1000,20)
      n3d = mifv1(201,mifv1(1,1))
      n2d = mifv1(202,mifv1(1,1))
      do iskp=1,nskip
         do iskpflds = 1, n3d
            read(iunit,end=30)
         enddo
         do iskpflds=n3d+1, n3d+n2d
            read(iunit,end=30)
         enddo
      enddo
   30 continue
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine wcalchyd(pstd,dmap,sigf,xmap,pstx,sigh,
     &   uuu,vvv,qvp,tmk,verv,miy,mjx,mkzh,ds,rgas,grav,ptop)
c
      dimension pstd(miy,mjx),dmap(miy,mjx),sigf(mkzh+1),
     &   xmap(miy,mjx),pstx(miy,mjx),sigh(mkzh),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   verv(miy,mjx,mkzh)
c
      dimension hmfd(100),sgd(100)
c
      ds8=8.*ds
      ds2=2.*ds
c
      do 1000 j=1,mjx-1
c
      jp1=min(j+1,mjx-1)
      jm1=max(j-1,1)
      xjfact=2./(jp1-jm1)
c
      do 1000 i=1,miy-1
c
      ip1=min(i+1,miy-1)
      im1=max(i-1,1)
      yifact=2./(ip1-im1)
c
c   Calculate pstar tendency and horizontal mass flux div.
c
      pstten=0.
      do 200 k=1,mkzh
       hmfd(k)=
     &  (pstd(i+1,j+1)*( uuu(i+1,j+1,k)+vvv(i+1,j+1,k))/dmap(i+1,j+1)+
     &   pstd(i  ,j+1)*( uuu(i  ,j+1,k)-vvv(i  ,j+1,k))/dmap(i  ,j+1)+
     &   pstd(i+1,j  )*(-uuu(i+1,j  ,k)+vvv(i+1,j  ,k))/dmap(i+1,j  )+
     &   pstd(i  ,j  )*(-uuu(i  ,j  ,k)-vvv(i  ,j  ,k))/dmap(i  ,j  ))/
     &   ds2
         pstten=pstten-hmfd(k)*(sigf(k+1)-sigf(k))*xmap(i,j)*xmap(i,j)
  200 continue
c
c   Calculate sigmadot
c
      sgd(1)=0.
      do 400 k=2,mkzh+1
         sgd(k)=sgd(k-1)-(pstten+xmap(i,j)*xmap(i,j)*hmfd(k-1))*
     &      (sigf(k)-sigf(k-1))/pstx(i,j)
  400 continue
c
c   Check sigmadot at ground is less than ~1 cm/s
c
      sgdcms=sgd(mkzh+1)*1.e6
      if (abs(sgdcms).gt.1.) then
         write(6,*)'Surface sigmadot too big!  i,j,sgdcms= ',
     &      i,j,sgdcms
         stop
      endif
c
c   Interpolate sigmadot to half levels.
c
      do 420 k=1,mkzh
         verv(i,j,k)=.5*(sgd(k)+sgd(k+1))
  420 continue
c
c   Calculate omega.
c
      do 700 k=1,mkzh
         verv(i,j,k)=(pstx(i,j)*verv(i,j,k)+sigh(k)*(pstten+(
     &    (uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+uuu(i+1,j+1,k))*
     &    (pstx(i,jp1)-pstx(i,jm1))*xjfact+
     &    (vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+vvv(i+1,j+1,k))*
     &    (pstx(ip1,j)-pstx(im1,j))*yifact)/ds8*xmap(i,j)))*
     &    1000.  !dPa/sec
  700 continue
c
c   Calculate approximate w (ignoring slope and vertical velocity
c      of pressure surface). No pert. needed here (hydrostatic).
c
      do 710 k=1,mkzh
         verv(i,j,k)=-rgas*virtual(tmk(i,j,k),.001*qvp(i,j,k))/
     &      (grav*(sigh(k)*pstx(i,j)+ptop))*verv(i,j,k)*
     &      .1  ! verv was in dPa/sec, want it in cm/sec
  710 continue
c
 1000 continue
      return
      end
