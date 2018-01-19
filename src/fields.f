c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine fields(cfeld,wk,indwk,icdwk,rlevl,rlavl,unwk,
     &   uuu,vvv,tmk,qvp,prs,ght,www,sfp,sfpsm,dmap,xmap,ter,
     &   cor,unorth,vnorth,rstmv,rrfst,pslab1,pslab2,
     &   incwk,ipl,iplstrt,idimn,rcrag,ismcp,
     &   rcrbg,cptyp,mdate,rhour,ydist,
     &   xdist,xseclen,nscrs,raddf,csave,lredo,ihrip,rhrip,rsepa,
     &   chrip,vardesc,lgrad,llapl,lhadv,
     &   igdir,iqgsm,plchun,casename,iendc,engplttl,ctjfl,rtjst,
     &   rtjen,cdiff,rdiff,iovly,ldfrl,xtime,tacch,ccalb,
     &   xtimeavl,cxtimeavl,ncxc,maxtavl,nxtavl,nxt,ldiffsuccess,
     &   rip_root,maxslab,maxlev,maxpl,miy,mjx,mkzh,mabpl,morpl,
     &   istopmiss,ixwin,iywin)
c
c   Compute fields to plot. When adding a new field, the first call to
c   getpt must be for the output field (pl2 or pl3), then call getpt
c   for the scratch arrays.
c
      dimension uuu(miy,mjx,mkzh), vvv(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh), qvp(miy,mjx,mkzh), prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh), www(miy,mjx,mkzh),
     &   dmap(miy,mjx), xmap(miy,mjx), ter(miy,mjx),
     &   cor(miy,mjx), sfp(miy,mjx), sfpsm(miy,mjx)
c
      dimension wk(miy,mjx,maxslab),indwk(3,maxpl),icdwk(maxpl),
     &   pslab1(mabpl,morpl),pslab2(mabpl,morpl),
     &   unorth(miy,mjx),vnorth(miy,mjx),rstmv(2,maxpl),
     &   rcrag(2,maxpl),rcrbg(2,maxpl),
     &   idimn(maxpl),igdir(maxpl),rsepa(32,maxpl),
     &   rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),rrfst(4,maxpl),
     &   raddf(maxpl),ismcp(maxpl),iqgsm(maxpl),
     &   rtjst(maxpl),rtjen(maxpl),rdiff(maxpl),iovly(maxpl),
     &   xtimeavl(maxtavl),xtimeavl_sv(maxtavl)
c
      character cfeld(3,maxpl)*10,cptyp(maxpl)*2,
     &   csave(maxpl)*10,ctjfl(maxpl)*256,engplttl(maxpl)*36,
     &   cdiff(maxpl)*256,unwk(maxpl)*24,varname*10,
     &   casename*256,cxtimeavl(maxtavl)*10,ccalb(maxpl)*256,
     &   casename_sv*256,cxtimeavl_sv(maxtavl)*10,rip_root*256
      logical lredo(maxpl),lgrad(maxpl),
     &   llapl(maxpl),lhadv(maxpl),ldfrl(maxpl),ldiffsuccess
      dimension uscratch(1000),vscratch(1000)
      dimension calb(1000)
      dimension ixwin(2,maxpl),iywin(2,maxpl)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      character getfromcore(20)*10, hydrometeors(20)*10
c
      logical numeric
      external numeric
c
      include 'comconst'
c
      include 'pointers'
c
c   Following are fields that are already in core memory (i.e.,
c   in arrays by the same name) and therefore do not need to be
c   calculated or read from a file.
c
      data (getfromcore(i),i=1,12) /'uuu       ','vvv       ',
     &    'tmk       ','prs       ','qvp       ','ght       ',
     &    'ter       ','dmap      ','xmap      ',
     &    'cor       ','sfp       ','www       '/
c
      data (hydrometeors(i),i=1,5) /'qcw       ','qra       ',
     &    'qci       ','qsn       ','qgr       '/
c
      do 2000 ifld=1,3
c
      if (ifld.eq.2.and.
     &    (cptyp(ipl)(2:2).eq.'c'.or.cptyp(ipl)(2:2).eq.'h')) goto 2000
      if (ifld.eq.3.and.cptyp(ipl).ne.'vv') goto 2000
c
c   Don't read or calculate field if already in work array
c
      do iplchk=iplstrt,ipl-1
         if ((cfeld(ifld,ipl).eq.cfeld(ifld,iplchk)).and.
     &       (ismcp(ipl).eq.ismcp(iplchk)).and.
     &       (lgrad(ipl).eqv.lgrad(iplchk)).and.
     &       (llapl(ipl).eqv.llapl(iplchk)).and.
     &       (lhadv(ipl).eqv.lhadv(iplchk)).and.
     &       igdir(ipl).eq.igdir(iplchk).and.
     &       iqgsm(ipl).eq.iqgsm(iplchk).and.
     &       rdiff(ipl).eq.rdiff(iplchk).and.
     &       cdiff(ipl).eq.cdiff(iplchk)) then
            if (iplchk.gt.iplstrt.and.raddf(iplchk-1).ne.0..and.
     &          raddf(iplchk).eq.0.) goto 291
            if (rsepa(1,ipl).eq.rsepa(1,iplchk).and.
     &          rsepa(2,ipl).eq.rsepa(2,iplchk).and.
     &          rsepa(3,ipl).eq.rsepa(3,iplchk).and.
     &          rsepa(4,ipl).eq.rsepa(4,iplchk).and.
     &          rsepa(5,ipl).eq.rsepa(5,iplchk).and.
     &          rsepa(6,ipl).eq.rsepa(6,iplchk).and.
     &          rsepa(7,ipl).eq.rsepa(7,iplchk).and.
     &          rsepa(8,ipl).eq.rsepa(8,iplchk).and.
     &          rsepa(9,ipl).eq.rsepa(9,iplchk)) then
               indwk(ifld,ipl)=indwk(ifld,iplchk)
               icdwk(ipl)=icdwk(iplchk)
               unwk(ipl)=unwk(iplchk)
               engplttl(ipl)=engplttl(iplchk)
               idimn(ipl)=3
               if (idimn(iplchk).eq.2) idimn(ipl)=2
               goto 1980
            endif
         endif
      enddo
 291  continue
c
c   Some preliminary stuff
c
      idiffpass=0
      idimn(ipl)=3
      engplttl(ipl)=' '
 35   ifree=incwk
      idiffpass=idiffpass+1
c
c   If the field is one that is not in core memory, it may be in
c   a file, either because it was processed by ripdp, or because
c   it was calculated by RIP and saved.  Therefore, for fields
c   not in core memory (i.e., not in the "getfromcore" list),
c   check to see if a file exists and if so, read it from the
c   file (unless the user has set "redo" because he/she specifically
c   wants the field to be re-calculated; or unless it is one of the
c   "hydrometeor" fields, in which case it should be processed by the
c   approapriate "if" block, because special things are done with the
c   hydrometeor fields).
c
      if (lredo(ipl)) goto 9
      do i=1,20
         if (getfromcore(i).eq.'          ') goto 7
         if (cfeld(ifld,ipl).eq.getfromcore(i)) goto 9
      enddo
 7    continue
      do i=1,20
         if (hydrometeors(i).eq.'          ') goto 8
         if (cfeld(ifld,ipl).eq.hydrometeors(i)) goto 9
      enddo
 8    continue
      call getvarinfo(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,cfeld(ifld,ipl),maxtavl,ndim,icd,vardesc,plchun,
     &   0,istat,iup)
      if (istat.ne.1) goto 9
      write(iup,*)'    Reading field ',cfeld(ifld,ipl),
     &   ' from a file.'
      if (ndim.eq.2) then
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &      ncxc,cfeld(ifld,ipl),miy,mjx,mkzh,maxtavl,2,
     &      0,pl2,istat)
      elseif (ndim.eq.3) then
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &      ncxc,cfeld(ifld,ipl),miy,mjx,mkzh,maxtavl,3,
     &      0,pl3,istat)
      else
         write(iup,*)'Dimension of data is not 2 or 3.  Stopping.'
         stop
      endif
      idimn(ipl)=ndim
      indwk(ifld,ipl)=incwk
      icdwk(ipl)=icd
      unwk(ipl)=plchun
      do i=35,1,-1
         if (vardesc(i:i+1).eq.', ') then
            iendv=i-1
            goto 71
         endif
      enddo
      if (vardesc(36:36).eq.',') then
         iendv=35
      else
         iendv=36
      endif
 71   continue
      engplttl(ipl)=vardesc(1:iendv)
      goto 1970
c
 9    continue
c
c   Load or calculate new field.
c
      if (cfeld(ifld,ipl)(1:4).eq.'uuu ') then! u-velocity, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(uuu,pl3,miy,mjx,mkzh,3,0,1.,0.)
         if (rstmv(2,ipl).ne.0.0) then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               pl3(i,j,k)=uuu(i,j,k)-rstmv(2,ipl)
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvv ') then! v-velocity, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(vvv,pl3,miy,mjx,mkzh,3,0,1.,0.)
         if (rstmv(1,ipl).ne.0.0) then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               pl3(i,j,k)=vvv(i,j,k)-rstmv(1,ipl)
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmc ') then! temperature, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=tmk(i,j,k)-celkel
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:8).eq.'k-index ') then! model k (v. lev.) index
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do j=1,mjx
         do i=1,miy
            do k=2,mkzh-1
               pl3(i,j,k)=float(k)
            enddo
            pl3(i,j,1)=.9999
            pl3(i,j,mkzh)=float(mkzh)+.0001
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Model k-index'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tptslr ') then
c        Temperature pert. from US standard lapse rate, K (or deg. C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            trefer=celkel+max(-55.,15.-.0065*ght(i,j,k))
            pl3(i,j,k)=tmk(i,j,k)-trefer
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. pert. (from US std. atm.)'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmk ') then! temperature, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(tmk,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmf ') then! temperature, deg F
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=(tmk(i,j,k)-celkel)*1.8 + 32.0
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='~S~o~N~F'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tpt ') then! tmp. prtrb., K (or deg C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         alnslp=log(.01*rrfst(1,ipl))
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            alnp=log(prs(i,j,k))
            reftmk=rrfst(2,ipl)+rrfst(3,ipl)*(alnp-alnslp)
            reftmk=max(reftmk,rrfst(4,ipl))
            pl3(i,j,k)=tmk(i,j,k)-reftmk
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. pert. (from MM5 std. atm.)'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvp ') then! vir. tmp. prb., K (or deg C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         alnslp=log(.01*rrfst(1,ipl))
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            alnp=log(prs(i,j,k))
            reftmk=rrfst(2,ipl)+rrfst(3,ipl)*(alnp-alnslp)
            reftmk=max(reftmk,rrfst(4,ipl))
            pl3(i,j,k)=virtual(tmk(i,j,k),qvp(i,j,k))-reftmk
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vir. temp. prt. (from MM5 std. atm.)'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:6).eq.'tvphm ') then! vir. tmp. prb., incl.
c                                                   ! hydroms., K (or deg C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &      3,1,1.,1.)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3a,istat)
            if (istat.ge.0) call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
         endif
         alnslp=log(.01*rrfst(1,ipl))
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            hydrometmr=.001*pl3(i,j,k) ! g/kg to kg/kg
            alnp=log(prs(i,j,k))
            reftmk=rrfst(2,ipl)+rrfst(3,ipl)*(alnp-alnslp)
            reftmk=max(reftmk,rrfst(4,ipl))
            pl3(i,j,k)=virtualhyd(tmk(i,j,k),qvp(i,j,k),hydrometmr)
     &         -reftmk
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vir. temp. prt. incl. hydroms.'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdp ' .or. 
     &        cfeld(ifld,ipl)(1:4).eq.'tdc ') then! dewpoint, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sdpc'  .or. ! surface dewpoint, deg C
     &        cfeld(ifld,ipl)(1:4).eq.'sdpf'  .or. ! surface dewpoint, deg F
     &        cfeld(ifld,ipl)(1:4).eq.'sdpk') then ! surface dewpoint, deg K
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         if (iplevdata.le.3) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qvp_sfan  ',miy,mjx,mkzh,maxtavl,2,1,
     &           scr2a,istat)
            do j = 1, mjx-1
            do i = 1, miy-1
               scr2a(i,j)= .001 * scr2a(i,j)
            enddo
            enddo
         else
            istat = -1
            if ( xtime .gt. 0.) then    ! Check if we have Q2
              call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'Q2        ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat)
            endif
            if ( istat .lt. 0 ) then
              do j = 1, mjx-1
              do i = 1, miy-1
                 scr2a(i,j)=qvp(i,j,mkzh)
              enddo
              enddo
            endif
         endif
         call tdpcalc(scr2a,sfp,pl2,miy,mjx,1)
         unwk(ipl)='~S~o~N~C'
         if (cfeld(ifld,ipl)(1:4).eq.'sdpf') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j) = pl2(i,j) * 1.8 +32.
            enddo
            enddo
            unwk(ipl)='~S~o~N~F'
         else if (cfeld(ifld,ipl)(1:4).eq.'sdpk') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j) = pl2(i,j) + celkel
            enddo
            enddo
            unwk(ipl)='~S~o~N~K'
         endif
         idimn(ipl)=2
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Surface dewpoint temperature'
      elseif (cfeld(ifld,ipl)(1:4).eq.'twb ') then! wet-blb tmp., deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wetbulbcalc(prs,tmk,qvp,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Wetbulb temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvk ') then! virtual temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
            pl3(i,j,k)=virtual(tmk(i,j,k),qvp(i,j,k))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvc ') then! virtual temp., deg. C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
            pl3(i,j,k)=virtual(tmk(i,j,k),qvp(i,j,k))-celkel
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvp ') then! vapor mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(qvp,pl3,miy,mjx,mkzh,3,1,1000.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Water vapor mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qra') then! rain mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         if (iqgsm(ipl).ne.0)
     &      call smoothrain(pl3,miy,mjx,mkzh,pslab1,pslab2,
     &         mabpl,morpl,iqgsm(ipl))
         if (iice.eq.0.and.cfeld(ifld,ipl)(4:4).ne.'b') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).eq.'b') then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            do j=1,mjx-1
            do i=1,miy-1
               probsnow=scr3b(i,j,3)
            do k=1,mkzh
               if (probsnow.ge.50..or.(probsnow.lt.50..and.
     &             k.ne.mkzh.and.tmk(i,j,k).lt.celkel)) 
     &            pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Rain water mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qsn') then! snow mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).ne.'b') then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).ge.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).eq.'b') then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            do j=1,mjx-1
            do i=1,miy-1
               probsnow=scr3b(i,j,3)
            do k=1,mkzh
               if (.not.(probsnow.ge.50..or.(probsnow.lt.50..and.
     &             k.ne.mkzh.and.tmk(i,j,k).lt.celkel))) 
     &            pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Snow mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qgr ') then! graupel mixing ratio, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &        scr3a,istat)
         if (istat.ge.0) then
            call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Graupel mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qpr ') then! qsn+qra+qgra, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         if (iqgsm(ipl).ne.0)
     &      call smoothrain(pl3,miy,mjx,mkzh,pslab1,pslab2,
     &         mabpl,morpl,iqgsm(ipl))
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3a,istat)
            if (istat.ge.0)
     &         call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Total precipitation mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'snf ') then! % snow=(qsn+qgra)/qpr
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         if (iice.eq.0) then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (pl3(i,j,k).le..0001e-3) then
                  pl3(i,j,k)=rmsg
               elseif (tmk(i,j,k).ge.celkel) then
                  pl3(i,j,k)=0.
               else
                  pl3(i,j,k)=100.
               endif
            enddo
            enddo
            enddo
         elseif (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3b,istat)
            if (istat.ge.0)
     &         call addorfill(scr3b,scr3a,miy,mjx,mkzh,3,1,1.,1.)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qtot=pl3(i,j,k)+scr3a(i,j,k)
               if (qtot.le..0001e-3) then
                  pl3(i,j,k)=rmsg
               else
                  pl3(i,j,k)=100.*scr3a(i,j,k)/qtot
               endif
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Frozen precipitation fraction'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:8).eq.'pr_stand') then! prcp rate (stand), mm/h
         wfac=1.
         if (cfeld(ifld,ipl)(9:9).eq.'p') wfac=0.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat) 
            if (iqgsm(ipl).ne.0)
     &         call smoothrain(scr3b,miy,mjx,mkzh,pslab1,pslab2,
     &            mabpl,morpl,iqgsm(ipl))
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)      
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=(3.6*((prs(i,j,k)*100./
c                      mm/h
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3a(i,j,k))*
     &         (((100.*19.3*(((1.18452108874)/((prs(i,j,k)*100./
c                      a_fall_g    rho_0                              
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.5)
c                                                        exp for dens corr fac
     &         *9.7309*(1./6.)*((((prs(i,j,k)*100./
c          gamma(4+b_fall_g) 
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3a(i,j,k))/(PI*400.*amax1(1.0E4,amin1(4.0E6,
c                                rho_g          
     &         (2.38*(PI*400/((prs(i,j,k)*100./(rgas*
     &         virtual(tmk(i,j,k),qvp(i,j,k))))*
     &         scr3a(i,j,k)))**0.92)))))**0.0925))
c               3 lines above are N_0_g     b_fall_g*.25              
     &         -(wfac*www(i,j,k)))/100.))+
c
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3b(i,j,k))*
     &         (((100.*842.*(((1.18452108874)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.5)
     &         *17.8379*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3b(i,j,k))/(PI*1000.*(4.9E8*tanh((.0002-
     &         (scr3b(i,j,k)*.001))/1.E-4)+5.1E8)))**0.20))
     &         -(wfac*www(i,j,k)))/100.))+
c
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3c(i,j,k))*
     &         (((11.72*100.*(((1.18452108874)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.5)
     &         *10.2754*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*.001*
     &         scr3c(i,j,k))/(PI*100.*(amin1(2.0E8,2.0E6*exp(-0.12*
     &         (tmk(i,j,k)-273.15))))))**0.1025))
     &         -(wfac*www(i,j,k)))/100.))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Precip rate (standard)'
         unwk(ipl)='mm h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'pr_coldt') then! prcp rate (coldt), mm/h
         wfac=1.
         if (cfeld(ifld,ipl)(9:9).eq.'p') wfac=0.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat) 
            if (iqgsm(ipl).ne.0)
     &         call smoothrain(scr3b,miy,mjx,mkzh,pslab1,pslab2,
     &            mabpl,morpl,iqgsm(ipl))
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)      
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=(3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3a(i,j,k))*
     &         (((100.*19.3*(((1.13571616794)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.33333)
     &         *9.7309*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3a(i,j,k))/(PI*400.*amax1(1.0E4,amin1(4.0E6,
     &         (2.38*(PI*400/((prs(i,j,k)*100./(rgas*
     &         virtual(tmk(i,j,k),qvp(i,j,k))))*
     &         scr3a(i,j,k)))**0.92)))))**0.0925))
     &         -(wfac*www(i,j,k)))/100.))+
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3b(i,j,k))*
     &         (((100.*842.*(((1.20443751401)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.4)
     &         *17.8379*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3b(i,j,k))/(PI*1000.*(4.9E8*tanh((.0002-
     &         (scr3b(i,j,k)*.001))/1.E-4)+5.1E8)))**0.20))
     &         -(wfac*www(i,j,k)))/100.))+
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3c(i,j,k))*
     &         (((11.72*100.*(((1.13571616794)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.33333)
     &         *2.7114*(1./1.8274)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*.001*
     &         scr3c(i,j,k))/(0.01854*1.8274
     &         *(amin1(2.0E8,2.0E6*exp(-0.12*
     &         (tmk(i,j,k)-273.15))))))**(.41/2.9)))
     &         -(wfac*www(i,j,k)))/100.))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Precip rate (cold-type)'
         unwk(ipl)='mm h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'pr_dendr') then! prcp rate (dendr), mm/h
         wfac=1.
         if (cfeld(ifld,ipl)(9:9).eq.'p') wfac=0.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
                   call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat) 
            if (iqgsm(ipl).ne.0)
     &         call smoothrain(scr3b,miy,mjx,mkzh,pslab1,pslab2,
     &            mabpl,morpl,iqgsm(ipl))
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)      
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=(3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3a(i,j,k))*
     &         (((100.*19.3*(((1.13571616794)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.33333)
     &         *9.7309*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3a(i,j,k))/(PI*400.*amax1(1.0E4,amin1(4.0E6,
     &         (2.38*(PI*400/((prs(i,j,k)*100./(rgas*
     &         virtual(tmk(i,j,k),qvp(i,j,k))))*
     &         scr3a(i,j,k)))**0.92)))))**0.0925))
     &         -(wfac*www(i,j,k)))/100.))+
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3b(i,j,k))*
     &         (((100.*842.*(((1.20443751401)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.4)
     &         *17.8379*(1./6.)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*0.001*
     &         scr3b(i,j,k))/(PI*1000.*(4.9E8*tanh((.0002-
     &         (scr3b(i,j,k)*.001))/1.E-4)+5.1E8)))**0.20))
     &         -(wfac*www(i,j,k)))/100.))+
     &         (3.6*((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*scr3c(i,j,k))*
     &         (((1.68*100.*(((1.13571616794)/((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))))**0.33333)
     &         *3.0036*(1./2.3999)*((((prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k))))*.001*
     &         scr3c(i,j,k))/(0.05239*2.3999
     &         *(amin1(2.0E8,2.0E6*exp(-0.12*
     &         (tmk(i,j,k)-273.15))))))**(.217/3.19)))
     &         -(wfac*www(i,j,k)))/100.))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Precip rate (dendrite)'
         unwk(ipl)='mm h~S~-1~N~'
c
c KWM -- Max reflectivity from REFL_10CM field
c
      elseif (cfeld(ifld,ipl)(1:9).eq.'maxrefl10') then
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'REFL_10CM ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
         refmax = -30.
         refmin = 100.
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j) = -30.
            do k=1,mkzh
               pl2(i,j) = max(pl2(i,j),scr3b(i,j,k))
            enddo
            if (pl2(i,j).gt.refmax) refmax = pl2(i,j)
            if (pl2(i,j).lt.refmin) refmin = pl2(i,j)
         enddo
         enddo
         idimn(ipl)=2
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='maxdbz_refl_10cm'
         unwk(ipl)='dBz'
c
c    compute max reflectivity in a column: WW, April 29, 2005
c
      elseif (cfeld(ifld,ipl)(1:6).eq.'maxdbz') then!max reflectivity, dBZ
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,scr3a,
     &        istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3c,istat)
            if (istat.eq.-1) call fillarray(scr3c,miy*mjx*mkzh,0.)
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) then
                  scr3b(i,j,k)=scr3a(i,j,k)
                  scr3a(i,j,k)=0.
               endif
               scr3c(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=.001*scr3a(i,j,k)   ! Convert to kg/kg
            scr3b(i,j,k)=.001*scr3b(i,j,k)   ! Convert to kg/kg
            scr3c(i,j,k)=.001*scr3c(i,j,k)   ! Convert to kg/kg
         enddo
         enddo
         enddo
         if (cfeld(ifld,ipl)(7:10).eq.'    ') then
            in0r=0
            in0s=0
            in0g=0
            iliqskin=0
            engplttl(ipl)='Max Reflectivity'
         else
            read(cfeld(ifld,ipl)(7:7),'(i1)') in0r
            read(cfeld(ifld,ipl)(8:8),'(i1)') in0s
            read(cfeld(ifld,ipl)(9:9),'(i1)') in0g
            read(cfeld(ifld,ipl)(10:10),'(i1)') iliqskin
          engplttl(ipl)='Max Reflectivity ('//cfeld(ifld,ipl)(7:10)//')'
         endif
         call dbzcalc(qvp,scr3a,scr3b,scr3c,tmk,prs,scr3d,miy,mjx,mkzh,
     &      in0r,in0s,in0g,iliqskin)
         refmax = -30.
         refmin = 100.
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j) = -30.
            do k=1,mkzh
               pl2(i,j) = max(pl2(i,j),scr3d(i,j,k))
            enddo
            if (pl2(i,j).gt.refmax) refmax = pl2(i,j)
            if (pl2(i,j).lt.refmin) refmin = pl2(i,j)
         enddo
         enddo
         idimn(ipl)=2
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='dBZ'
c
c    compute max ferrier reflectivity in a column
c
      elseif (cfeld(ifld,ipl)(1:6).eq.'maxdbf') then!max ferrier refl, dBZ
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'dbz_fer   ',miy,mjx,mkzh,maxtavl,3,1,scr3a,
     &        istat)
         refmax = -30.
         refmin = 100.
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j) = -30.
            do k=1,mkzh
               pl2(i,j) = max(pl2(i,j),scr3a(i,j,k))
            enddo
            if (pl2(i,j).gt.refmax) refmax = pl2(i,j)
            if (pl2(i,j).lt.refmin) refmin = pl2(i,j)
         enddo
         enddo
         idimn(ipl)=2
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Max Refl (Ferrier)'
         unwk(ipl)='dBZ'
      elseif (cfeld(ifld,ipl)(1:3).eq.'dbz') then   ! reflectivity, dBZ
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,scr3a,
     &        istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3c,istat)
            if (istat.eq.-1) call fillarray(scr3c,miy*mjx*mkzh,0.)
         else
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) then
                  scr3b(i,j,k)=scr3a(i,j,k)
                  scr3a(i,j,k)=0.
               endif
               scr3c(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=.001*scr3a(i,j,k)   ! Convert to kg/kg
            scr3b(i,j,k)=.001*scr3b(i,j,k)   ! Convert to kg/kg
            scr3c(i,j,k)=.001*scr3c(i,j,k)   ! Convert to kg/kg
         enddo
         enddo
         enddo
         if (cfeld(ifld,ipl)(4:7).eq.'    ') then
            in0r=0
            in0s=0
            in0g=0
            iliqskin=0
            engplttl(ipl)='Reflectivity'
         else
            read(cfeld(ifld,ipl)(4:4),'(i1)') in0r
            read(cfeld(ifld,ipl)(5:5),'(i1)') in0s
            read(cfeld(ifld,ipl)(6:6),'(i1)') in0g
            read(cfeld(ifld,ipl)(7:7),'(i1)') iliqskin
          engplttl(ipl)='Reflectivity ('//cfeld(ifld,ipl)(4:7)//')'
         endif
         call dbzcalc(qvp,scr3a,scr3b,scr3c,tmk,prs,pl3,miy,mjx,mkzh,
     &      in0r,in0s,in0g,iliqskin)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='dBZ'
      elseif (cfeld(ifld,ipl)(1:5).eq.'pcpwv'.or.  ! precip'bl water (vapor), mm
     &        cfeld(ifld,ipl)(1:5).eq.'pcptw'.or. ! precip'bl water (vap+hydrometeors), mm
     &        cfeld(ifld,ipl)(1:5).eq.'pcpw '.or. ! for backward compatability
     &        cfeld(ifld,ipl)(1:7).eq.'intcld '.or. ! integ. cld. hyds., mm
     &        cfeld(ifld,ipl)(1:7).eq.'intclq '.or. ! integ. cld. lq. wat., mm
     &        cfeld(ifld,ipl)(1:7).eq.'intpcp ') then! integ. pcp. hyds., mm
         if (cfeld(ifld,ipl)(1:5).eq.'pcptw' .or. 
     &       cfeld(ifld,ipl)(1:5).eq.'pcpw ') then
            iv=1
            icliq=1
            icice=1
            ip=1
            engplttl(ipl)='Total precipitable water'
         elseif (cfeld(ifld,ipl)(1:5).eq.'pcpwv') then
            iv=1
            icliq=0
            icice=0
            ip=0
            engplttl(ipl)='Precipitable water vapor'
         elseif (cfeld(ifld,ipl)(1:7).eq.'intcld ') then
            iv=0
            icliq=1
            icice=1
            ip=0
            engplttl(ipl)='Column-integ. cloud hydrometeors'
         elseif (cfeld(ifld,ipl)(1:7).eq.'intclq ') then
            iv=0
            icliq=1
            icice=0
            ip=0
            engplttl(ipl)='Column-integ. cloud liq. water'
         elseif (cfeld(ifld,ipl)(1:7).eq.'intpcp ') then
            iv=0
            icliq=0
            icice=0
            ip=1
            engplttl(ipl)='Column-integ. precip. hydrometeors'
         endif
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call fillarray(scr3a,miy*mjx*mkzh,0.)
         call fillarray(scr3b,miy*mjx*mkzh,0.)
         call fillarray(scr3c,miy*mjx*mkzh,0.)
         call addorfill(qvp,scr3a,miy,mjx,mkzh,3,1,1.,1.)
         if (iv.eq.1) call addorfill(qvp,scr3b,miy,mjx,mkzh,
     &                               3,1,1.,1.)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3c,istat)
         call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
         if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                               3,1,.001,1.)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3c,istat)
         call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
         if (iice.eq.1) then
            if (icliq.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                  3,1,.001,1.)
         elseif (icliq.eq.1.or.icice.eq.1) then
            do j=1,mjx-1
            do i=1,miy-1
            do k=1,mkzh
               if (icliq.eq.1.and.tmk(i,j,k).ge.celkel) then
                  scr3b(i,j,k)=scr3b(i,j,k)+.001*scr3c(i,j,k)
               elseif (icice.eq.1.and.tmk(i,j,k).lt.celkel) then
                  scr3b(i,j,k)=scr3b(i,j,k)+.001*scr3c(i,j,k)
               endif
            enddo
            enddo
            enddo
         endif
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
            if (icice.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                  3,1,.001,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
            if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                  3,1,.001,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3c,istat)
            if (istat.ge.0) then
               call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
               if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                     3,1,.001,1.)
            endif
         endif
         call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=0.
            do k=1,mkzh
               if (k.eq.1) then
                  pabove=2.*prs(i,j,k)-scr3c(i,j,k)
               else
                  pabove=scr3c(i,j,k-1)
               endif
               pl2(i,j)=pl2(i,j)+scr3b(i,j,k)/(1.+scr3a(i,j,k))*
     &            100.*(scr3c(i,j,k)-pabove) !mix rats are kg/kg, prs's are hPa
            enddo
            pl2(i,j)=pl2(i,j)/grav  ! dens of water cancels with conv to mm
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='mm'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qcw') then! cld wat mix r., g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         if (iice.eq.0) then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud water mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'lwc ') then! liquid water content, g/m^3
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))   !air density
            pl3(i,j,k)=pl3(i,j,k)*scr3a(i,j,k)
            if (iice.eq.0.and.tmk(i,j,k).lt.celkel) pl3(i,j,k)=0.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud liquid water content'
         unwk(ipl)='g m~S~-3~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qci') then! cld ice mix r., g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
         elseif (iice.eq.0) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).ge.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud ice mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qcl ') then! qcw+qci, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Total cloud mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ght ') then! geop. height, m
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(ght,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Geopotential height'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:7).eq.'ghtagl ') then! height above ground
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(ght,pl3,miy,mjx,mkzh,3,1,1.,0.)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)-ter(i,j)
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Geopotential height above ground'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:4).eq.'prs ') then! pressure, hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(prs,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Pressure'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rho ') then! density, kg/m**3
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Density'
         unwk(ipl)='kg m~S~-3~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'phydp') then! hyd. prs. pert., hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c
         read(cfeld(ifld,ipl)(6:7),'(i2)') klev
         cc1=rgas/grav*(-.5)*rrfst(3,ipl)
         cc2=rgas/grav*(rrfst(3,ipl)*log(.01*rrfst(1,ipl))-
     &      rrfst(2,ipl))
         cc3=rgas/grav*(rrfst(2,ipl)-.5*rrfst(3,ipl)*
     &      log(.01*rrfst(1,ipl)))*log(.01*rrfst(1,ipl))
         alnpreftpause=(rrfst(4,ipl)-rrfst(2,ipl))/rrfst(3,ipl)+
     &      log(.01*rrfst(1,ipl))
         ztpause=cc1*alnpreftpause*alnpreftpause+
     &      cc2*alnpreftpause+cc3
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (ght(i,j,k).le.ztpause) then
               alnpref=
     &          (-cc2-sqrt(cc2*cc2-4.*cc1*(cc3-ght(i,j,k))))/(2.*cc1)
            else
               alnpref=alnpreftpause-grav/(rgas*rrfst(4,ipl))*
     &          (ght(i,j,k)-ztpause)
            endif
            scr3a(i,j,k)=exp(alnpref)   ! ref. pressure in hPa
         enddo
         enddo
         enddo
c
         do j=1,mjx-1
         do i=1,miy-1
c
         pl3(i,j,klev)=prs(i,j,klev)-scr3a(i,j,klev)
         alogpold=alog(prs(i,j,klev))
         do k=klev+1,mkzh
            tvbar=.5*(virtual(tmk(i,j,k-1),qvp(i,j,k-1))
     &            +virtual(tmk(i,j,k),qvp(i,j,k)))
            alogp=alogpold-grav*(ght(i,j,k)-ght(i,j,k-1))/(rgas*tvbar)
            pl3(i,j,k)=exp(alogp)-scr3a(i,j,k)
            alogpold=alogp
         enddo
         alogpold=alog(prs(i,j,klev))
         do k=klev-1,1,-1
            tvbar=.5*(virtual(tmk(i,j,k),qvp(i,j,k))
     &            +virtual(tmk(i,j,k+1),qvp(i,j,k+1)))
            alogp=alogpold-grav*(ght(i,j,k)-ght(i,j,k+1))/(rgas*tvbar)
            pl3(i,j,k)=exp(alogp)-scr3a(i,j,k)
            alogpold=alogp
         enddo
c
         enddo
         enddo
c
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Hydrostatic pressure pert.'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ppt ') then! press pertbn, hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cc1=rgas/grav*(-.5)*rrfst(3,ipl)
         cc2=rgas/grav*(rrfst(3,ipl)*log(.01*rrfst(1,ipl))-
     &      rrfst(2,ipl))
         cc3=rgas/grav*(rrfst(2,ipl)-.5*rrfst(3,ipl)*
     &      log(.01*rrfst(1,ipl)))*log(.01*rrfst(1,ipl))
         alnpreftpause=(rrfst(4,ipl)-rrfst(2,ipl))/rrfst(3,ipl)+
     &      log(.01*rrfst(1,ipl))
         ztpause=cc1*alnpreftpause*alnpreftpause+
     &      cc2*alnpreftpause+cc3
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            if (ght(i,j,k).le.ztpause) then
               alnpref=
     &          (-cc2-sqrt(cc2*cc2-4.*cc1*(cc3-ght(i,j,k))))/(2.*cc1)
            else
               alnpref=alnpreftpause-grav/(rgas*rrfst(4,ipl))*
     &          (ght(i,j,k)-ztpause)
            endif
            pl3(i,j,k)=prs(i,j,k)-exp(alnpref)
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Pressure pert. (from MM5 std. atm.)'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'the ') then! pot. temperature, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call thecalc(prs,tmk,qvp,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:8).eq.'refmos3 ') then ! obs 3D refl mosaic, dbz
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
c
c        ctjfl holds the file name of the 3D mosaic data
c
         nxmos=1201
c         nymos=1001    ! tile 3
         nymos=501     ! tile5
         nzmos=21
         call refmos3calc(ctjfl(ipl),pl3,ght,miy,mjx,mkzh,
     &      nxmos,nymos,nzmos,1)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Observed 3D reflectivity mosaic'
         unwk(ipl)='dBZ'
      elseif (cfeld(ifld,ipl)(1:10).eq.'refmos3mx ') then ! obs 3Dmx refl mosaic, dbz
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c
c        ctjfl holds the file name of the 3D mosaic data
c
         nxmos=1201
         nymos=1001
c         nymos=501
         nzmos=21
         call refmos3calc(ctjfl(ipl),scr3a,ght,miy,mjx,mkzh,
     &      nxmos,nymos,nzmos,2)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=scr3a(i,j,1)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Observed 3Dmx reflectivity mosaic'
         unwk(ipl)='dBZ'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:7).eq.'refmos2') then ! obs 2D refl mosaic, dbz
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
c
c        ctjfl holds the file name of the 3D mosaic data
c
         nxmos=1201
         nymos=1001
         read(cfeld(ifld,ipl)(8:8),'(i1)') ichoice
         call refmos2calc(ctjfl(ipl),pl2,miy,mjx,nxmos,nymos,ichoice)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Observed 2D refl mosaic (ch  )'
         write(engplttl(ipl)(29:29),'(i1)') ichoice
         unwk(ipl)='dBZ'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:8).eq.'profuuu '.or.   ! profiler vel
     &        cfeld(ifld,ipl)(1:8).eq.'profvvv ') then ! comps, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
c
c        ctjfl holds the file name of the profiler data
c        pl3 will hold the u or v field
c        rsepa(1,ipl) holds the assumed storm motion direction (deg. compass)
c        rsepa(2,ipl) holds the assumed storm speed (m/s)
c
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call profvelcalc(cfeld(ifld,ipl)(5:5),ctjfl(ipl),mdate,rhour,
     &      xtime,cosa,sina,rsepa(1,ipl),rip_root,
     &      unorth,vnorth,ght,pl3,miy,mjx,mkzh)
c
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         if (cfeld(ifld,ipl)(5:5).eq.'u') then
            engplttl(ipl)='Profiler wind (x-comp.)'
         elseif (cfeld(ifld,ipl)(5:5).eq.'v') then
            engplttl(ipl)='Profiler wind (y-comp.)'
         endif
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tgradx ') then ! grad(T) (x-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            pl3(i,j,k)=pl3(i,j,k)*1e5  ! convert from K/m to K/100km
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. grad. (x-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tgrady ') then ! grad(T) (y-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'y',miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            pl3(i,j,k)=pl3(i,j,k)*1e5  ! convert from K/m to K/100km
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. grad. (y-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tgradm ') then ! grad(T) (mag) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            pl3(i,j,k)=1e5*sqrt(pl3(i,j,k)*pl3(i,j,k)+
     &         scr3a(i,j,k)*scr3a(i,j,k)) ! compute mag, conv from K/m to K/100km
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. grad. (mag.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'tadv ') then ! temp adv, K/day
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivc(tmk,1,prs,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            ucross=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                       uuu(i+1,j+1,k))
            vcross=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                       vvv(i+1,j+1,k))
            pl3(i,j,k)=-(ucross*pl3(i,j,k)+vcross*scr3a(i,j,k))
     &         *86400.   ! convert from K/s to K/day
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temp. adv.'
         unwk(ipl)='K day~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'tgradrx ') then ! retrieved
c                                                      ! grad(T) (x-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               vscratch(k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                          vvv(i+1,j+1,k))
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               pl3(i,j,k)=-1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &            (vscratch(kp1)-vscratch(km1))/
     &            (prs(i,j,kp1)-prs(i,j,km1))
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Retr. temp. grad. (x-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'tgradry ') then ! retrieved
c                                                      ! grad(T) (y-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               uscratch(k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                          uuu(i+1,j+1,k))
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               pl3(i,j,k)= 1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &            (uscratch(kp1)-uscratch(km1))/
     &            (prs(i,j,kp1)-prs(i,j,km1))
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Retr. temp. grad. (y-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'tgradrm ') then ! retrieved
c                                                      ! grad(T) (mag) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               uscratch(k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                          uuu(i+1,j+1,k))
               vscratch(k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                          vvv(i+1,j+1,k))
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               dprs=prs(i,j,kp1)-prs(i,j,km1)
               ushear=(uscratch(kp1)-uscratch(km1))/dprs
               vshear=(vscratch(kp1)-vscratch(km1))/dprs
               pl3(i,j,k)= 1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &            sqrt(ushear*ushear+vshear*vshear)
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Retr. temp. grad. (mag.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'tadvr ') then ! retrieved temp adv, K/day
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               uscratch(k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                          uuu(i+1,j+1,k))
               vscratch(k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                          vvv(i+1,j+1,k))
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               dprs=prs(i,j,kp1)-prs(i,j,km1)
               ushear=(uscratch(kp1)-uscratch(km1))/dprs
               vshear=(vscratch(kp1)-vscratch(km1))/dprs
               pl3(i,j,k)= -86400.*cor(i,j)*prs(i,j,k)/rgas*
     &            (uscratch(k)*(-vshear)+vscratch(k)*(ushear))
c                (converted from K/s to K/day)
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Retr. temp. adv.'
         unwk(ipl)='K day~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'tgradrpx ') then ! retrieved (from prof)
c                                                      ! grad(T) (x-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profvvv   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               n=0
               vscratch(k)=0.
               if (scr3a(i,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3a(i,j,k)
                  n=n+1
               endif
               if (scr3a(i+1,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3a(i+1,j,k)
                  n=n+1
               endif
               if (scr3a(i,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3a(i,j+1,k)
                  n=n+1
               endif
               if (scr3a(i+1,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3a(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  vscratch(k)=vscratch(k)/n
               else
                  vscratch(k)=rmsg
               endif
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               if (vscratch(kp1).ne.rmsg.and.vscratch(km1).ne.rmsg) then
                  pl3(i,j,k)=-1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &               (vscratch(kp1)-vscratch(km1))/
     &               (prs(i,j,kp1)-prs(i,j,km1))
               else
                  pl3(i,j,k)=rmsg
               endif
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Prf-retr. temp. grad. (x-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'tgradrpy ') then ! retrieved (from prof)
c                                                      ! grad(T) (y-comp) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profuuu   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               n=0
               uscratch(k)=0.
               if (scr3a(i,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j,k)
                  n=n+1
               endif
               if (scr3a(i+1,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j,k)
                  n=n+1
               endif
               if (scr3a(i,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j+1,k)
                  n=n+1
               endif
               if (scr3a(i+1,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  uscratch(k)=uscratch(k)/n
               else
                  uscratch(k)=rmsg
               endif
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               if (uscratch(kp1).ne.rmsg.and.uscratch(km1).ne.rmsg) then
                  pl3(i,j,k)= 1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &               (uscratch(kp1)-uscratch(km1))/
     &               (prs(i,j,kp1)-prs(i,j,km1))
               else
                  pl3(i,j,k)=rmsg
               endif
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Prf-retr. temp. grad. (y-comp.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'tgradrpm ') then ! retrieved (from prof)
c                                                      ! grad(T) (mag) K/100km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profuuu   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profvvv   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               n=0
               uscratch(k)=0.
               if (scr3a(i,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j,k)
                  n=n+1
               endif
               if (scr3a(i+1,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j,k)
                  n=n+1
               endif
               if (scr3a(i,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j+1,k)
                  n=n+1
               endif
               if (scr3a(i+1,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  uscratch(k)=uscratch(k)/n
               else
                  uscratch(k)=rmsg
               endif
               n=0
               vscratch(k)=0.
               if (scr3b(i,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i,j,k)
                  n=n+1
               endif
               if (scr3b(i+1,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i+1,j,k)
                  n=n+1
               endif
               if (scr3b(i,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i,j+1,k)
                  n=n+1
               endif
               if (scr3b(i+1,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  vscratch(k)=vscratch(k)/n
               else
                  vscratch(k)=rmsg
               endif
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               if (uscratch(kp1).ne.rmsg.and.uscratch(km1).ne.rmsg.and.
     &             vscratch(kp1).ne.rmsg.and.vscratch(km1).ne.rmsg) then
                  dprs=prs(i,j,kp1)-prs(i,j,km1)
                  ushear=(uscratch(kp1)-uscratch(km1))/dprs
                  vshear=(vscratch(kp1)-vscratch(km1))/dprs
                  pl3(i,j,k)= 1.e5*cor(i,j)*prs(i,j,k)/rgas*
     &               sqrt(ushear*ushear+vshear*vshear)
               else
                  pl3(i,j,k)=rmsg
               endif
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Prf-retr. temp. grad. (mag.)'
         unwk(ipl)='K (100 km)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tadvrp ') then ! retrieved (from prof)
c                                                       temp adv, K/day
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profuuu   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'profvvv   ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
         do j=1,mjx-1
         do i=1,miy-1
            do k=1,mkzh
               n=0
               uscratch(k)=0.
               if (scr3a(i,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j,k)
                  n=n+1
               endif
               if (scr3a(i+1,j,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j,k)
                  n=n+1
               endif
               if (scr3a(i,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i,j+1,k)
                  n=n+1
               endif
               if (scr3a(i+1,j+1,k).ne.rmsg) then
                  uscratch(k)=uscratch(k)+scr3a(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  uscratch(k)=uscratch(k)/n
               else
                  uscratch(k)=rmsg
               endif
               n=0
               vscratch(k)=0.
               if (scr3b(i,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i,j,k)
                  n=n+1
               endif
               if (scr3b(i+1,j,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i+1,j,k)
                  n=n+1
               endif
               if (scr3b(i,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i,j+1,k)
                  n=n+1
               endif
               if (scr3b(i+1,j+1,k).ne.rmsg) then
                  vscratch(k)=vscratch(k)+scr3b(i+1,j+1,k)
                  n=n+1
               endif
               if (n.gt.0) then
                  vscratch(k)=vscratch(k)/n
               else
                  vscratch(k)=rmsg
               endif
            enddo
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
               if (uscratch(kp1).ne.rmsg.and.uscratch(km1).ne.rmsg.and.
     &             vscratch(kp1).ne.rmsg.and.vscratch(km1).ne.rmsg.and.
     &             uscratch(k).ne.rmsg.and.vscratch(k).ne.rmsg) then
                  dprs=prs(i,j,kp1)-prs(i,j,km1)
                  ushear=(uscratch(kp1)-uscratch(km1))/dprs
                  vshear=(vscratch(kp1)-vscratch(km1))/dprs
                  pl3(i,j,k)= -86400.*cor(i,j)*prs(i,j,k)/rgas*
     &               (uscratch(k)*(-vshear)+vscratch(k)*(ushear))
c                   (converted from K/s to K/day)
               else
                  pl3(i,j,k)=rmsg
               endif
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Prf-retr. temp. adv.'
         unwk(ipl)='K day~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'eth ') then! eqv. pot. temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Equivalent potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:7).eq.'sateth ') then! sat. eqv. pot. temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call sateqthecalc(tmk,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Saturated equiv. pot. temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rcum'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'rexp'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'rfra'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'stot'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'mtot'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'rtot') then  ! rainfall, mm
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2c,wk,maxslab)
         if (cfeld(ifld,ipl)(1:4).eq.'rcum') then
            icum=1
            iexp=0
            engplttl(ipl)='Cumulus precip.'
            iendeng=15
         elseif (cfeld(ifld,ipl)(1:4).eq.'rexp') then
            icum=0
            iexp=1
            engplttl(ipl)='Explicit precip.'
            iendeng=16
         elseif (cfeld(ifld,ipl)(1:4).eq.'rfra') then
            icum=1
            iexp=1
            engplttl(ipl)='Cumulus fraction'
            iendeng=16
         elseif (cfeld(ifld,ipl)(1:4).eq.'rtot') then
            icum=1
            iexp=1
            engplttl(ipl)='Total precip.'
            iendeng=13
! total snow and total mixed precip requires the WRF SR field (Snow Ratio)
         elseif (cfeld(ifld,ipl)(1:4).eq.'stot') then
            icum=1
            iexp=1
            engplttl(ipl)='Total snow precip.'
            iendeng=18
         elseif (cfeld(ifld,ipl)(1:4).eq.'mtot') then
            icum=1
            iexp=1
            engplttl(ipl)='Total mixed precip.'
            iendeng=19
         endif
         iendcf=4
         if (cfeld(ifld,ipl)(iendcf+1:iendcf+2).eq.'sh') then
c
c         Given time is interpreted as "since hour x"
c
            irel=0
            ispos=iendcf+3
            iepos=index(cfeld(ifld,ipl),' ')-1
            if (iepos.eq.-1) iepos=10
         else
c
c         Given time is interpreted as "in past x hours"
c
            irel=1
            ispos=iendcf+1
            iepos=iendcf+index(cfeld(ifld,ipl)(iendcf+1:),'h')-1
         endif
         if (.not.numeric(cfeld(ifld,ipl)(ispos:iepos))) then
            write(iup,*)'Processing rainfall specifier in fields.f.'
            write(iup,*)'Non-numeric characters found where time (in h)'
            write(iup,*)'is supposed to be specified.'
            write(iup,*)'feld was given as ',cfeld(ifld,ipl)
            stop
         endif
         read(cfeld(ifld,ipl)(ispos:iepos),fmt=*)raintime
         if (irel.eq.0) then
            xtime_get=raintime
            engplttl(ipl)(iendeng+1:)=' since h '//
     &         cfeld(ifld,ipl)(ispos:iepos)
         else
            xtime_get=max(0.0,xtime-raintime)
            engplttl(ipl)(iendeng+1:)=' in past '//
     &         cfeld(ifld,ipl)(ispos:iepos)//' h'
         endif
c
c      Get current rainfall
c
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=0.
         enddo
         enddo
         if (icum.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'rtc       ',miy,mjx,mkzh,maxtavl,2,1,
     &           scr2a,istat)
            call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,1.,1.)
            if (cfeld(ifld,ipl)(1:4).eq.'rfra') 
     &           call addorfill(scr2a,scr2b,miy,mjx,mkzh,2,1,1.,0.)
         endif
         if (iexp.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'rte       ',miy,mjx,mkzh,maxtavl,2,1,
     &           scr2a,istat)
            call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,1.,1.)
         endif
         if (cfeld(ifld,ipl)(1:4).eq.'stot' .or. 
     &       cfeld(ifld,ipl)(1:4).eq.'mtot') then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'SR        ',miy,mjx,mkzh,maxtavl,2,1,
     &           scr2c,istat)
         endif
c
c      Subtract previous rainfall if xtime_get is greater than ~.1 seconds
c
         if (xtime_get.ge.3e-5) then
c
c         Check if requested time is available
c
            iavail=0
            do i=1,nxtavl
               if (abs(xtimeavl(i)-xtime_get).le.tacch) then
                  iavail=1
                  nxt_get=i
                  goto 57
               endif
            enddo
 57         continue
            if (iavail.eq.0) then
               write(iup,*)'   Requested previous time ',
     &            xtime_get,' for'
               write(iup,*)'   acccumulated rainfall is not available.'
               write(iup,*)'   RIP will plot requested rainfall type'
               write(iup,*)'   since h 0.'
               goto 59
            endif
c
            if (icum.eq.1) then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_get,
     &              ncxc,'rtc       ',miy,mjx,mkzh,
     &              maxtavl,2,1,scr2a,istat)
               call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,-1.,1.)
            if (cfeld(ifld,ipl)(1:4).eq.'rfra') 
     &           call addorfill(scr2a,scr2b,miy,mjx,mkzh,2,1,-1.,1.)
            endif
            if (iexp.eq.1) then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt_get,
     &              ncxc,'rte       ',miy,mjx,mkzh,
     &              maxtavl,2,1,scr2a,istat)
               call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,-1.,1.)
            endif

	   if (cfeld(ifld,ipl)(1:4).eq.'stot' .or.
     &         cfeld(ifld,ipl)(1:4).eq.'mtot') then
	      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &             'SR        ',miy,mjx,mkzh,maxtavl,2,1,
     &             scr2a,istat)
	 ! simply average the SR over the period of interest - almost every approach
	 ! is going to have problems. scr2c contains the final SR.
               call addorfill(scr2a,scr2c,miy,mjx,mkzh,2,1,0.5,0.5)
	   endif
         endif
 59      continue
	 if (cfeld(ifld,ipl)(1:4).eq.'stot' .or.
     &       cfeld(ifld,ipl)(1:4).eq.'mtot') then
	   do j=1,mjx-1
	   do i=1,miy-1
	    if (cfeld(ifld,ipl)(1:4).eq.'stot' ) then
	      if (scr2c(i,j) .le. 0.7 ) pl2(i,j) = 0.
	    else if (cfeld(ifld,ipl)(1:4).eq.'mtot' ) then
	      if (scr2c(i,j) .gt. 0.7 .or. 
     &            scr2c(i,j) .lt. 0.3 ) pl2(i,j) = 0.
	     endif
	   enddo
	   enddo
	 endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='mm'
	 domtot = 0.
         if (cfeld(ifld,ipl)(1:4).eq.'rfra') then
           do j=1,mjx-1
           do i=1,miy-1
             if (pl2(i,j).gt..0001e-3) then
              pl2(i,j) = 100.*scr2b(i,j)/pl2(i,j)
             else
               pl2(i,j) = rmsg
             endif
           enddo
           enddo
           unwk(ipl)='%'
	 else
           do j=ixwin(1,ipl),ixwin(2,ipl)-1
           do i=iywin(1,ipl),iywin(2,ipl)-1
	     domtot = domtot + pl2(i,j)
	   enddo
	   enddo
           if (ixwin(1,ipl) .ne. 1 .or. ixwin(2,ipl) .ne. mjx .or.
     &         iywin(1,ipl) .ne. 1 .or. iywin(2,ipl) .ne. miy) then
            write(6,*) 'area total ',cfeld(ifld,ipl)(1:4),' = ',domtot
           else
            write(6,*) 'domain total ',cfeld(ifld,ipl)(1:4),' = ',domtot
           endif
         endif
      elseif (cfeld(ifld,ipl)(1:4).eq.'ter ') then! terrain, m
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(ter,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Terrain height AMSL'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xmap ') then! map fac (x-points)
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(xmap,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Map factor on cross points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:5).eq.'dmap ') then! map fac (dot points)
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(dmap,pl2,miy,mjx,mkzh,2,0,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Map factor on dot points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'cor ') then! Coriolis, per s
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(cor,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Coriolis parameter'
         unwk(ipl)='s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sfp ') then ! surface pressure, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=sfp(i,j)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Surface pressure'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'frgm'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'frem') then
c
c      Frontogenesis (Miller form), K per 100 km per day
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
c
c      Put theta or theta-e in scr3a
c
         if (cfeld(ifld,ipl)(3:3).eq.'g') then
            call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(3:3).eq.'e') then
            call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         endif
c
         call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3d(i,j,k)=sqrt(scr3b(i,j,k)*scr3b(i,j,k)+
     &         scr3c(i,j,k)*scr3c(i,j,k))
         enddo
         enddo
         enddo
         if (cfeld(ifld,ipl)(5:7).eq.'had') then
            call derivc(scr3d,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(scr3d,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               ucross=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))
               vcross=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))
               pl3(i,j,k)=-(ucross*scr3e(i,j,k)+vcross*scr3f(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (hor. adv.)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'vad') then
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=
     &            -.01*www(i,j,k)*(scr3d(i,j,kp1)-scr3d(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (ver. adv.)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'div') then
            call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-.5/scr3d(i,j,k)*(scr3b(i,j,k)*
     &               scr3b(i,j,k)+scr3c(i,j,k)*scr3c(i,j,k))*
     &               (scr3e(i,j,k)+scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (div. term)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'def') then
            call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-.5/scr3d(i,j,k)*(scr3b(i,j,k)*
     &               scr3b(i,j,k)-scr3c(i,j,k)*scr3c(i,j,k))*
     &               (scr3e(i,j,k)-scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'y',miy,mjx,mkzh)
            call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'x',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=pl3(i,j,k)-1./scr3d(i,j,k)*scr3b(i,j,k)*
     &               scr3c(i,j,k)*(scr3e(i,j,k)+scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (def. term)'
         elseif (cfeld(ifld,ipl)(5:8).eq.'dia ') then
            if (rsepa(1,ipl).eq.200..and.rsepa(2,ipl).eq.1..and.
     &          rsepa(3,ipl).eq.1.6) then
               irhthtype=0  ! use RH wrt liquid (0) or ice (1) for thresholding
               rhthresh=99. ! RH threshold
               ilrtype=0    ! pseudoadi. lapse rate is wrt liq (0) or ice (1)
            else
               irhthtype=rsepa(1,ipl)
               rhthresh=rsepa(2,ipl)
               ilrtype=rsepa(3,ipl)
            endif
            call condheat(tmk,qvp,www,irhthtype,ilrtype,0,prs,
     &         miy,mjx,mkzh,pl3,rhthresh)
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
c                 Factor of 2.7778e-4 needed because heating is in K/hour
c                 and we want per second.
                  pl3(i,j,k)=1./scr3d(i,j,k)*2.778e-4*(scr3b(i,j,k)*
     &               scr3e(i,j,k)+scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (dia. term)'
         elseif (cfeld(ifld,ipl)(5:8).eq.'dia2') then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'ttendbmp  ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &           istat)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               gammam=gamma*(1.+gammamd*qvp(i,j,k))
               pl3(i,j,k)=pl3(i,j,k)*(1000./prs(i,j,k))**gammam
            enddo
            enddo
            enddo
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
c              Factor of 2.7778e-4 needed because heating is in K/hour
c              and we want per second.
               pl3(i,j,k)=1./scr3d(i,j,k)*2.778e-4*(scr3b(i,j,k)*
     &            scr3e(i,j,k)+scr3c(i,j,k)*scr3f(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (dia2. term)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'dlg') then
            if (rsepa(1,ipl).eq.200..and.rsepa(2,ipl).eq.1..and.
     &          rsepa(3,ipl).eq.1.6) then
               irhthtype=0  ! use RH wrt liquid (0) or ice (1) for thresholding
               rhthresh=99. ! RH threshold
               ilrtype=0    ! pseudoadi. lapse rate is wrt liq (0) or ice (1)
            else
               irhthtype=rsepa(1,ipl)
               rhthresh=rsepa(2,ipl)
               ilrtype=rsepa(3,ipl)
            endif
            call condheat(tmk,qvp,www,irhthtype,ilrtype,1,prs,
     &         miy,mjx,mkzh,pl3,rhthresh)
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(pl3,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=www(i,j,k)*1.e-5/scr3d(i,j,k)*(
     &              scr3b(i,j,k)*scr3e(i,j,k)+scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller fgen. (dia. term, MALR grad)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'dwg') then
            if (rsepa(1,ipl).eq.200..and.rsepa(2,ipl).eq.1..and.
     &          rsepa(3,ipl).eq.1.6) then
               irhthtype=0  ! use RH wrt liquid (0) or ice (1) for thresholding
               rhthresh=99. ! RH threshold
               ilrtype=0    ! pseudoadi. lapse rate is wrt liq (0) or ice (1)
            else
               irhthtype=rsepa(1,ipl)
               rhthresh=rsepa(2,ipl)
               ilrtype=rsepa(3,ipl)
            endif
            call condheat(tmk,qvp,www,irhthtype,ilrtype,1,prs,
     &         miy,mjx,mkzh,pl3,rhthresh)
            call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=pl3(i,j,k)*1.e-5/scr3d(i,j,k)*(
     &              scr3b(i,j,k)*scr3e(i,j,k)+scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller fgen. (dia. term, w grad)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'til') then
            call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-1./scr3d(i,j,k)*.01*(scr3a(i,j,kp1)-
     &               scr3a(i,j,km1))/(ght(i,j,kp1)-ght(i,j,km1))*
     &               (scr3b(i,j,k)*scr3e(i,j,k)+
     &                scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (til. term)'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)*3.6e8 ! K/m/s to K/100km/hour
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='K (100 km h)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'frg2d'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'fre2d') then
c
c      Frontogenesis (uni-directional form), K per 100 km per day
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
c
c      Put theta or theta-e in scr3a
c
         if (cfeld(ifld,ipl)(3:3).eq.'g') then
            call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(3:3).eq.'e') then
            call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         endif
         cosa=xdist/xseclen
         sina=ydist/xseclen
c
c      For some terms, put thermal gradient into scr3a instead
c
         if (cfeld(ifld,ipl)(6:8).ne.'dia'.and.
     &       cfeld(ifld,ipl)(6:8).ne.'dii'.and.
     &       cfeld(ifld,ipl)(6:8).ne.'til') then
c
         call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (cfeld(ifld,ipl)(6:8).ne.'shr') then
c              thermal grad. along cross-sec:
               scr3a(i,j,k)=cosa*scr3b(i,j,k)+sina*scr3c(i,j,k)
            else
c              thermal grad. into cross-sec:
               scr3a(i,j,k)=-sina*scr3b(i,j,k)+cosa*scr3c(i,j,k) !cor!
            endif
         enddo
         enddo
         enddo
c
         endif
c
c      Smooth scr3a (on constant pressure), since it gets used in
c         almost every term.
c
         call smoothcp(scr3a,1,scr3b,prs,pslab1,pslab2,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
c
         if (cfeld(ifld,ipl)(6:8).eq.'aad') then
            call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
            call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3d(i,j,k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))-rstmv(2,ipl)
               scr3e(i,j,k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))-rstmv(1,ipl)
            enddo
            enddo
            enddo
            call smoothcp(scr3d,1,scr3a,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call smoothcp(scr3e,1,scr3a,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-(cosa*scr3d(i,j,k)+sina*scr3e(i,j,k))*
     &                     (cosa*scr3b(i,j,k)+sina*scr3c(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (adv. along c.s.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'iad') then
            call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
            call derivc(scr3a,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3d(i,j,k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))-rstmv(2,ipl)
               scr3e(i,j,k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))-rstmv(1,ipl)
            enddo
            enddo
            enddo
            call smoothcp(scr3d,1,scr3a,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call smoothcp(scr3e,1,scr3a,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-(-sina*scr3d(i,j,k)+cosa*scr3e(i,j,k))*
     &                     (-sina*scr3b(i,j,k)+cosa*scr3c(i,j,k))  !cor!
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (adv. into c.s.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'vad') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3b(i,j,k)=www(i,j,k)
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=
     &            -.01*scr3b(i,j,k)*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (ver. adv.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'cnf') then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               scr3b(i,j,k)=cosa*uuu(i,j,k)+sina*vvv(i,j,k) !wind along x-sec.
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivc(scr3b,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivc(scr3b,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               duadxa=cosa*scr3c(i,j,k)+sina*scr3d(i,j,k)
               pl3(i,j,k)=-scr3a(i,j,k)*duadxa
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (cnf. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'shr') then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               scr3b(i,j,k)=-sina*uuu(i,j,k)+cosa*vvv(i,j,k) !wind into x-sec.!cor!
            enddo
            enddo
            enddo
            call smoothcp(scr3b,0,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivc(scr3b,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivc(scr3b,0,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               duidxa=cosa*scr3c(i,j,k)+sina*scr3d(i,j,k)  !cor!
               pl3(i,j,k)=-scr3a(i,j,k)*duidxa
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (shr. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'til') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3b(i,j,k)=www(i,j,k)
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivc(scr3b,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivc(scr3b,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-.01*(scr3a(i,j,kp1)-
     &            scr3a(i,j,km1))/(ght(i,j,kp1)-ght(i,j,km1))*
     &            (cosa*scr3c(i,j,k)+sina*scr3d(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (til. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'dia') then
            if (rsepa(1,ipl).eq.200..and.rsepa(2,ipl).eq.1..and.
     &          rsepa(3,ipl).eq.1.6) then
               irhthtype=0  ! use RH wrt liquid (0) or ice (1) for thresholding
               rhthresh=99. ! RH threshold
               ilrtype=0    ! pseudoadi. lapse rate is wrt liq (0) or ice (1)
            else
               irhthtype=rsepa(1,ipl)
               rhthresh=rsepa(2,ipl)
               ilrtype=rsepa(3,ipl)
            endif
            call condheat(tmk,qvp,www,irhthtype,ilrtype,0,prs,
     &         miy,mjx,mkzh,scr3b,rhthresh)
            call smoothcp(scr3b,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivc(scr3b,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivc(scr3b,1,ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=2.778e-4*(cosa*scr3c(i,j,k)+
     &            sina*scr3d(i,j,k)) ! K/m/h to K/m/s
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (dia. term)'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)*3.6e8 ! K/m/s to K/100km/hour
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='K (100 km h)~S~-1~N~'
         if (cfeld(ifld,ipl)(3:3).eq.'e') then
            engplttl(ipl)(4:16)='theta-e frgn.'
         endif
      elseif (cfeld(ifld,ipl)(1:5).eq.'ugeo '.or.
     &        cfeld(ifld,ipl)(1:6).eq.'uageo ') then ! u_geos. or
c                                                      u_ageos., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'y',miy,mjx,mkzh)
         fbar=cor(miy/2,mjx/2)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               pl3(i,j,k)=-rgas*tv*pl3(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
            if (cfeld(ifld,ipl)(1:6).eq.'uageo ') then
               do j=1,mjx
               do i=1,miy
                  pl3(i,j,k)=uuu(i,j,k)-pl3(i,j,k)
               enddo
               enddo
            else
               if (rstmv(2,ipl).ne.0.0) then
                  do j=1,mjx
                  do i=1,miy
                     pl3(i,j,k)=pl3(i,j,k)-rstmv(2,ipl)
                  enddo
                  enddo
               endif
            endif
         enddo
         if (cfeld(ifld,ipl)(1:6).eq.'uageo ') then
            engplttl(ipl)='Ageostrophic wind (x-comp.)'
         else
            engplttl(ipl)='Geostrophic wind (x-comp.)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'vgeo '.or.
     &        cfeld(ifld,ipl)(1:6).eq.'vageo ') then ! v_geos. or
c                                                      v_ageos., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         fbar=cor(miy/2,mjx/2)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               pl3(i,j,k)=rgas*tv*pl3(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
            if (cfeld(ifld,ipl)(1:6).eq.'vageo ') then
               do j=1,mjx
               do i=1,miy
                  pl3(i,j,k)=vvv(i,j,k)-pl3(i,j,k)
               enddo
               enddo
            else
               if (rstmv(1,ipl).ne.0.0) then
                  do j=1,mjx
                  do i=1,miy
                     pl3(i,j,k)=pl3(i,j,k)-rstmv(1,ipl)
                  enddo
                  enddo
               endif
            endif
         enddo
         if (cfeld(ifld,ipl)(1:6).eq.'vageo ') then
            engplttl(ipl)='Ageostrophic wind (y-comp.)'
         else
            engplttl(ipl)='Geostrophic wind (y-comp.)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'xptgeo') then ! geo.wind prlel. to x-sec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)-rstmv(2,ipl)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)-rstmv(1,ipl)
               pl3(i,j,k)=cosa*ugeo+sina*vgeo
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Geostrophic wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'xptageo') then !ageo. wind prll. to xsec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
               pl3(i,j,k)=cosa*(uuu(i,j,k)-ugeo)+sina*(vvv(i,j,k)-vgeo)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Ageostrophic wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'xntgeo') then ! geo.wind norm. to x-sec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)-rstmv(2,ipl)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)-rstmv(1,ipl)
               pl3(i,j,k)=-sina*ugeo+cosa*vgeo
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Geostrophic wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'xntageo') then !ageo. wind norm. to xsec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
               pl3(i,j,k)=-sina*(uuu(i,j,k)-ugeo)+cosa*(vvv(i,j,k)-vgeo)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Ageostrophic wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'slp') then! sea-lev press, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (cfeld(ifld,ipl)(1:6).eq.'slpbm '.or.
     &       cfeld(ifld,ipl)(1:4).eq.'slp ') then   ! B&M SLP, hPa
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
	    do k=1,mkzh
	    do j=1,mjx-1
	    do i=1,miy-1
	       scr3a(i,j,k)=exp(-ght(i,j,k)/sclht)
	       scr3b(i,j,k)=prs(i,j,k)
	    enddo
	    enddo
	    enddo
            call vinterp('z',0.,1,1,1,scr3a,tmk,qvp,
     &         prs,ght,ter,sfp,sfpsm,.false.,0,'prs       ',
     &         scr3b,pslab1,mabpl,morpl,mjx-1,miy-1,miy,mjx,mkzh)
	    do j=1,mjx-1
	    do i=1,miy-1
	       pl2(i,j)=pslab1(j,i)
	    enddo
	    enddo
         elseif (cfeld(ifld,ipl)(1:6).eq.'slpgr ') then! GRAPH SLP, hPa
c
c         This is the version in Graph (added back to rip4 on 7/22/06 jfb)
c
            call seaprs(tmk,prs,ter,sfp,miy,mjx,mkzh,pl2,iup)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Sea-level pressure'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'dpbh') then! pressure diff.
c                                                   betw. htl levels, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=exp(-ght(i,j,k)/sclht)
            scr3b(i,j,k)=prs(i,j,k)
         enddo
         enddo
         enddo
         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') h1,h2
         h1=h1*.1     ! hm to km
         h2=h2*.1     ! hm to km
         call vinterp('z',h1,1,1,1,scr3a,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,.false.,0,'prs       ',
     &      scr3b,pslab1,mabpl,morpl,mjx-1,miy-1,miy,mjx,mkzh)
         call vinterp('z',h2,1,1,1,scr3a,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,.false.,0,'prs       ',
     &      scr3b,pslab2,mabpl,morpl,mjx-1,miy-1,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=abs(pslab2(j,i)-pslab1(j,i))
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Press. diff. bet. '//cfeld(ifld,ipl)(5:6)//
     &      '.'//cfeld(ifld,ipl)(7:7)//' and '//cfeld(ifld,ipl)(8:9)//
     &      '.'//cfeld(ifld,ipl)(10:10)//' km'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'thck') then! thickness
c                                                   betw. prs levels, dam
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=ght(i,j,k)
            scr3b(i,j,k)=prs(i,j,k)
         enddo
         enddo
         enddo
         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') p1,p2
         p1=p1*10.    ! kPa to hPa
         p2=p2*10.    ! kPa to hPa
         call vinterp('p',p1,1,1,1,scr3b,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,.false.,0,'ght       ',
     &      scr3a,pslab1,mabpl,morpl,mjx-1,miy-1,miy,mjx,mkzh)
         call vinterp('p',p2,1,1,1,scr3b,tmk,qvp,
     &      prs,ght,ter,sfp,sfpsm,.false.,0,'ght       ',
     &      scr3a,pslab2,mabpl,morpl,mjx-1,miy-1,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=0.1*abs(pslab2(j,i)-pslab1(j,i))
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)=cfeld(ifld,ipl)(5:7)//'0'//' to '//
     &      cfeld(ifld,ipl)(8:10)//'0 hPa thickness'
         unwk(ipl)='dam'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ctt ') then! cld-top temp, deg C
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
         endif
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
         call cttcalc(prs,scr3c,tmk,scr3a,scr3b,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud-top temperature'
         unwk(ipl)='~S~o~N~C'
c      elseif (cfeld(ifld,ipl)(1:9).eq.'rad_elev ') then! radar elev. angle, deg.
c         idimn(ipl)=2
c         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
c         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c         open (unit=iutrajin,file=ctjfl(ipl),form='unformatted',
c     &      status='old')
c         read (iutrajin,*)
c         read (iutrajin,*)
c         read (iutrajin,*) radar_lat
c         read (iutrajin,*) radar_lon
c         read (iutrajin,*) radar_elev
c         close (iutrajin)
cc
cc      Constants
cc
c         rke = 4./3.  ! Four-thirds earth approximation
c         r_earth = 6.37e6   ! radius of earth in meters
c         r_eff = rke * r_earth
c
c
c
c         indwk(ifld,ipl)=incwk
c         icdwk(ipl)=1
c         engplttl(ipl)='Radar elevation angle'
c         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:5).eq.'cap3 '.or.   ! 3D CAPE, J/kg
     &        cfeld(ifld,ipl)(1:5).eq.'cin3 ') then ! 3D conv. inhib., J/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,1)
         if (cfeld(ifld,ipl)(1:4).eq.'cap3') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=scr3a(i,j,k)
            enddo
            enddo
            enddo
            engplttl(ipl)='CAPE'
         elseif (cfeld(ifld,ipl)(1:4).eq.'cin3') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=scr3b(i,j,k)
            enddo
            enddo
            enddo
            engplttl(ipl)='Convective inhibition'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='J kg~S~-1~N~'
c
c      In the following: if "save" was requested, it will be done at the end
c      of this routine, but we might as well also save the other field
c      that was calculated (either cap3 or cin3), because they take so
c      long to calculate. 
c
         if (csave(ipl).ne.'dontsave  '.and..not.lgrad(ipl).and..not.
     &       llapl(ipl).and..not.lhadv(ipl).and.ismcp(ipl).eq.0) then
            if (cfeld(ifld,ipl)(1:5).eq.'cap3 ') then
               vardesc='Conv. Inh., J/kg'
               plchun=unwk(ipl)
               ndimen=3
               call writefile (scr3b,'cin3      ',0,
     &            ndimen,icdwk(ipl),vardesc,plchun,ihrip,rhrip,
     &            chrip,casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &            maxtavl,miy,mjx,mkzh)
            elseif (cfeld(ifld,ipl)(1:5).eq.'cin3 ') then
               vardesc='CAPE, J/kg'
               plchun=unwk(ipl)
               ndimen=3
               call writefile (scr3a,'cap3      ',0,
     &            ndimen,icdwk(ipl),vardesc,plchun,ihrip,rhrip,
     &            chrip,casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &            maxtavl,miy,mjx,mkzh)
            endif
         endif
      elseif (cfeld(ifld,ipl)(1:5).eq.'mcap '.or.   ! Max CAPE, J/kg
     &        cfeld(ifld,ipl)(1:5).eq.'mcin ') then ! assoc. CIN, J/kg
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,0)
         if (cfeld(ifld,ipl)(1:4).eq.'mcap') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j)=scr3a(i,j,mkzh)
            enddo
            enddo
            engplttl(ipl)='CAPE (for parcel with max theta-e)'
         elseif (cfeld(ifld,ipl)(1:4).eq.'mcin') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j)=scr3b(i,j,mkzh)
            enddo
            enddo
            engplttl(ipl)='CIN (for parcel with max theta-e)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='J kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'lcl '.or.
     &        cfeld(ifld,ipl)(1:4).eq.'lfc ') then! LCL or LFC, meters AGL
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
c
c      Put LCL (LFC) in the mkzh-1 (mkzh-2) slab of scr3b
c
         call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,0)
c
c      Now put it in pl2
c
         if (cfeld(ifld,ipl)(1:4).eq.'lfc ') then
            kget=mkzh-2
            engplttl(ipl)='LFC (for parcel with max theta-e)'
         else
            kget=mkzh-1
            engplttl(ipl)='LCL (for parcel with max theta-e)'
         endif
         do j = 1, mjx-1
         do i = 1, miy-1
            pl2(i,j)=scr3b(i,j,kget)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='m (AGL)'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmcl') then!tmp. of lifted parcel, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         read(cfeld(ifld,ipl)(5:6),'(i2)') kpar !nth level from bottom
         kpar=mkzh-kpar+1   ! k-level
         call liftparcel(prs,tmk,qvp,ght,pl3,kpar,0,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Lifted T from k=    '
         write(engplttl(ipl)(18:20),'(i3)') kpar
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvcl') then! T_v of lifted parc., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         read(cfeld(ifld,ipl)(5:6),'(i2)') kpar !nth level from bottom
         kpar=mkzh-kpar+1   ! k-level
         call liftparcel(prs,tmk,qvp,ght,pl3,kpar,1,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Lifted Tv from k=    '
         write(engplttl(ipl)(19:21),'(i3)') kpar
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:5).eq.'ethmx') then
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2c,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            eqthetamax=0.
            do k=1,mkzh
               if (ght(i,j,k)-ter(i,j).lt.3000..and.
     &             scr3a(i,j,k).gt.eqthetamax) then
                  kmax=k
                  eqthetamax=scr3a(i,j,k)
               endif
            enddo
            pmax=prs(i,j,kmax)
            if (cfeld(ifld,ipl)(6:6).eq.'v') then
               pl2(i,j)=eqthetamax
            elseif (cfeld(ifld,ipl)(6:6).eq.'p') then
               pl2(i,j)=pmax
            endif
         enddo
         enddo
         if (cfeld(ifld,ipl)(6:6).eq.'v') then! max eth, K
            engplttl(ipl)='Max theta-e below 3000 m AGL'
            unwk(ipl)='K'
         elseif (cfeld(ifld,ipl)(6:6).eq.'p') then! p-lev. of max eth, hPa
            engplttl(ipl)='Prs. at max theta-e below 3000 m AGL'
            unwk(ipl)='hPa'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)=' '
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'freg') then! flight regul.
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         if (cfeld(ifld,ipl)(6:6).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2c,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3b,istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3d,istat)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'c'.or.
     &       cfeld(ifld,ipl)(5:5).eq.'b') then! get ceiling
            clgfac=1.
            call pfcalc(prs,sfp,scr3e,miy,mjx,mkzh)
            call ceilingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr3e,scr2c,cfeld(ifld,ipl)(6:6),clgfac,scr2a,
     &         miy,mjx,mkzh)
         else
            call fillarray(scr2a,miy*mjx,rmsg)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'v'.or.
     &       cfeld(ifld,ipl)(5:5).eq.'b') then! get visibility
            call viscalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,prs,
     &         scr2c,cfeld(ifld,ipl)(6:6),scr2b,
     &         miy,mjx,mkzh)
         else
            call fillarray(scr2b,miy*mjx,rmsg)
         endif
         call fregcalc(scr2a,scr2b,cfeld(ifld,ipl)(5:5),pl2,miy,mjx)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Flight regulation category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:3).eq.'vis'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'o'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'f') then! hor. vis., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3b,istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3d,istat)
         endif
         call viscalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,prs,
     &      scr2a,cfeld(ifld,ipl)(4:4),pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Visibility'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:3).eq.'clg'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'o'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'f') then! cloud ceiling, m
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3b,istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3d,istat)
         endif
         call pfcalc(prs,sfp,scr3e,miy,mjx,mkzh)
         clgfac=1.
         call ceilingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &      prs,scr3e,scr2a,cfeld(ifld,ipl)(4:4),clgfac,pl2,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud ceiling'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:3).eq.'ext') then! extinc. coef, per km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3b,istat)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3c,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3d,istat)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'c') then! due to cloud water
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),1,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'r') then! due to rain
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),2,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'i') then! due to cloud ice
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),3,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'s') then! due to snow
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),4,pl2,
     &         miy,mjx,mkzh)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Extinction coefficient due to '//
     &      cfeld(ifld,ipl)(5:5)
         unwk(ipl)='km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'boc') then! Bocchieri prec prob
c
c      Fourth character determines desired prob:
c        'l':liquid; 'f':freezing; 'i':ice
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
         if (cfeld(ifld,ipl)(4:4).eq.'l') then
            n=1
            engplttl(ipl)='Prob. of liquid precip.'
         elseif (cfeld(ifld,ipl)(4:4).eq.'f') then
            n=2
            engplttl(ipl)='Prob. of freezing/mixed precip.'
         elseif (cfeld(ifld,ipl)(4:4).eq.'i') then
            n=3
            engplttl(ipl)='Prob. of frozen precip.'
         endif
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=scr3b(i,j,n)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tgc '.or.   ! ground temp., deg.C
     &        cfeld(ifld,ipl)(1:4).eq.'tgk ') then ! ground temp., deg.C
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (iplevdata.le.3) then
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'tgk       ',miy,mjx,mkzh,maxtavl,2,0,pl2,istat)
           if ( istat .lt. 0 ) then
             call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'SKINTEMP  ',miy,mjx,mkzh,maxtavl,2,1,pl2,istat)
           endif
         else
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'tgk       ',miy,mjx,mkzh,maxtavl,2,1,pl2,
     &        istat)
         endif
         subtract=0.
         if (cfeld(ifld,ipl)(1:4).eq.'tgc ') then
            subtract=celkel
            unwk(ipl)='~S~o~N~C'
         else
            subtract=0.
            unwk(ipl)='K'
         endif
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=pl2(i,j)-subtract
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Ground/sea-surface temperature'
      elseif (cfeld(ifld,ipl)(2:8).eq.'llmptf ') then
c
c      x or y value (on fine grid) based on conversion of lat/lon
c      array to x/y using maptform
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'xlat      ',miy,mjx,mkzh,maxtavl,2,1,pl2,
     &        istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'xlon      ',miy,mjx,mkzh,maxtavl,2,1,scr2a,
     &        istat)
         do j=1,mjx-1
         do i=1,miy-1
            call maptform(riy,rjx,pl2(i,j),scr2a(i,j),-1)
            rjxn=1.+(rjx-xjcorn)*refrat
            riyn=1.+(riy-yicorn)*refrat
            pl2(i,j)=rjxn
            if (cfeld(ifld,ipl)(1:1).eq.'y') pl2(i,j)=riyn
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='xlat/xlon-to-x-value'
         if (cfeld(ifld,ipl)(1:1).eq.'y') 
     &      engplttl(ipl)='xlat/xlon-to-y-value'
         unwk(ipl)='grid points'
      elseif (cfeld(ifld,ipl)(1:10).eq.'ivalcross ') then ! "i" (dot grid) value of cross points
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=float(i)+.5
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='"i" value of cross points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:10).eq.'jvalcross ') then ! "j" (dot grid) value of cross points
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=float(j)+.5
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='"j" value of cross points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xlat ') then ! latitude, degrees
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'xlat      ',miy,mjx,mkzh,maxtavl,2,1,pl2,
     &        istat)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Latitude'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xlon ') then ! longitude, degrees
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'xlon      ',miy,mjx,mkzh,maxtavl,2,1,pl2,
     &        istat)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Longitude'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:9).eq.'xlatmptf '.or.
     &        cfeld(ifld,ipl)(1:9).eq.'xlonmptf ') then
c
c      Maptform-calculated lat or lon, degrees
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx-1
            rj1=xjcorn+(j-.5)/refrat
         do i=1,miy-1
            ri1=yicorn+(i-.5)/refrat
            call maptform(ri1,rj1,rlat,rlon,1)
            rl=rlat
            if (cfeld(ifld,ipl)(1:9).eq.'xlonmptf ') rl=rlon
            pl2(i,j)=rl
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Latitude (maptform)'
         if (cfeld(ifld,ipl)(1:9).eq.'xlonmptf ')
     &      engplttl(ipl)='Longitude (maptform)'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:8).eq.'xlaterr '.or.
     &        cfeld(ifld,ipl)(1:8).eq.'xlonerr ') then
c
c      Error in maptform-calculated lat or lon,
c      compared to model output arrays
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        cfeld(ifld,ipl)(1:4)//'      ',
     &        miy,mjx,mkzh,maxtavl,2,1,pl2,istat)
         do j=1,mjx-1
            rj1=xjcorn+(j-.5)/refrat
         do i=1,miy-1
            ri1=yicorn+(i-.5)/refrat
            call maptform(ri1,rj1,rlat,rlon,1)
            rl=rlat
            if (cfeld(ifld,ipl)(1:8).eq.'xlonerr ') rl=rlon
            pl2(i,j)=pl2(i,j)-rl
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Latitude error'
         if (cfeld(ifld,ipl)(1:8).eq.'xlonerr ')
     &      engplttl(ipl)='Longitude error'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xlus ') then ! land use category
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'xlus      ',miy,mjx,mkzh,maxtavl,2,1,pl2,
     &        istat)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Land use category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'pvo ') then! pot. vorticity, PVU
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call pvocalc(xmap,uuu,vvv,cor,scr3a,prs,
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Potential vorticity'
         unwk(ipl)='PVU'
      elseif (cfeld(ifld,ipl)(1:4).eq.'pvm ') then! moist PV, PVU
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         call pvocalc(xmap,uuu,vvv,cor,scr3a,prs,
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Moist potential vorticity'
         unwk(ipl)='PVU'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vor '.or.! relative vort., per s
     &    cfeld(ifld,ipl)(1:4).eq.'avo ') then! absolute vort., per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         if (cfeld(ifld,ipl)(1:1).eq.'v') then
            icor=0
            engplttl(ipl)='Relative vorticity'
         elseif (cfeld(ifld,ipl)(1:1).eq.'a') then
            icor=1
            engplttl(ipl)='Absolute vorticity'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)-scr3a(i,j,k)+icor*cor(i,j)
c
c         Scale the vorticity by 1.e5
c
            pl3(i,j,k)=  pl3(i,j,k) * 1.e5
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'dil') then ! comp of dilat axis, s^-1
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
c
         call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'x',miy,mjx,mkzh)
         call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'y',miy,mjx,mkzh)
         call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
         call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            f1=scr3a(i,j,k)-scr3b(i,j,k)
            f2=scr3c(i,j,k)+scr3d(i,j,k)
            fff=sqrt(f1*f1+f2*f2)
            psi=0.5*atan2(f2,f1)
            if (cfeld(ifld,ipl)(4:4).eq.'x') then
               pl3(i,j,k)=1.e5*fff*cos(psi)
            elseif (cfeld(ifld,ipl)(4:4).eq.'y') then
               pl3(i,j,k)=1.e5*fff*sin(psi)
            endif
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         if (cfeld(ifld,ipl)(4:4).eq.'x') then
            engplttl(ipl)='Dilatation axis (x-comp.)'
         elseif (cfeld(ifld,ipl)(4:4).eq.'y') then
            engplttl(ipl)='Dilatation axis (y-comp.)'
         endif
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vox ') then ! x-comp. of vort., per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'y',miy,mjx,mkzh)
c        pl3 now holds (dw/dy)_z=const, in cm/s/m
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            dvdz=.25*(vvv(i,j,kp1)+vvv(i+1,j,kp1)+
     &                vvv(i,j+1,kp1)+vvv(i+1,j+1,kp1)-
     &                vvv(i,j,km1)-vvv(i+1,j,km1)-
     &                vvv(i,j+1,km1)-vvv(i+1,j+1,km1))/
     &             (ght(i,j,kp1)-ght(i,j,km1))
            pl3(i,j,k)=(.01*pl3(i,j,k)-dvdz)*1.e5
c           In the above, the .01 is to convert dw/dy from cm/s/m to /m, and
c           the 1.e5 scales the result by the commonly used 10^5 factor.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vorticity (x-comp.)'
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'voy ') then ! y-comp. of vort., per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
c        pl3 now holds (dw/dx)_z=const, in cm/s/m
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            dudz=.25*(uuu(i,j,kp1)+uuu(i+1,j,kp1)+
     &                uuu(i,j+1,kp1)+uuu(i+1,j+1,kp1)-
     &                uuu(i,j,km1)-uuu(i+1,j,km1)-
     &                uuu(i,j+1,km1)-uuu(i+1,j+1,km1))/
     &             (ght(i,j,kp1)-ght(i,j,km1))
            pl3(i,j,k)=(dudz-.01*pl3(i,j,k))*1.e5
c           In the above, the .01 is to convert dw/dx from cm/s/m to /m, and
c           the 1.e5 scales the result by the commonly used 10^5 factor.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vorticity, (y-comp.)'
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'vor3d ') then ! mag of 3-D rel vort, per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         call derivc(www,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            vvor=pl3(i,j,k)-scr3a(i,j,k)
            dudz=.25*(uuu(i,j,kp1)+uuu(i+1,j,kp1)+
     &                uuu(i,j+1,kp1)+uuu(i+1,j+1,kp1)-
     &                uuu(i,j,km1)-uuu(i+1,j,km1)-
     &                uuu(i,j+1,km1)-uuu(i+1,j+1,km1))/
     &             (ght(i,j,kp1)-ght(i,j,km1))
            dvdz=.25*(vvv(i,j,kp1)+vvv(i+1,j,kp1)+
     &                vvv(i,j+1,kp1)+vvv(i+1,j+1,kp1)-
     &                vvv(i,j,km1)-vvv(i+1,j,km1)-
     &                vvv(i,j+1,km1)-vvv(i+1,j+1,km1))/
     &             (ght(i,j,kp1)-ght(i,j,km1))
            xvor=(scr3c(i,j,k)-dvdz)
            yvor=(dudz-scr3b(i,j,k))
            pl3(i,j,k)=sqrt(vvor*vvor+xvor*xvor+yvor*yvor)*1.e5  ! scale by 1.e5
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Mag. of 3-D Rel. Vorticity'
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'div ') then! divergence, per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call derivc(uuu,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivc(vvv,0,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         call addorfill(pl3,pl3,miy,mjx,mkzh,3,1,1.e5,0.)  ! scale by 1.e5
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Divergence'
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rhu ') then! RH wrt liq, percent
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call rhucalc(qvp,tmk,prs,0,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Relative humidity (w.r.t. water)'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rhi ') then! RH wrt ice, percent
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call rhucalc(qvp,tmk,prs,1,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Relative humidity (w.r.t. ice)'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'www ') then
c
c      vert velocity (w in cm/s)
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(www,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vertical velocity'
         unwk(ipl)='cm s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'omg ') then
c
c      vert velocity (omega in dPa/s)
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Omega'
         unwk(ipl)='dPa s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'netasc ') then
c
c      Net ascent of trajectories (hPa).  Must have a trajectory
c      position file from a properly organized trajectory run.
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call netasc(casename,iendc,cxtimeavl,xtimeavl,
     &      ncxc,nxtavl,maxtavl,scr3a,scr3b,scr3c,scr3d,
     &      ctjfl(ipl),rtjst(ipl),rtjen(ipl),tacch,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Net ascent of trajectories'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvx ') then! Q_vector_x,
c                                                    10**-6 m/s**3/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         call bvfricalc(prs,ght,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,scr2a,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,1,1,0,0,0,0,100.,
     &      ihrip,rhrip,chrip,
     &      vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &      nxt,ncxc,maxtavl,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q (x-comp.)'
         unwk(ipl)='m s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-6
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvy ') then! Q_vector_y,
c                                                    10**-6 m/s**3/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         call bvfricalc(prs,ght,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,scr2a,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,2,1,0,0,0,0,100.,
     &      ihrip,rhrip,chrip,
     &      vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &      nxt,ncxc,maxtavl,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q (y-comp.)'
         unwk(ipl)='m s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-6
      elseif (cfeld(ifld,ipl)(1:6).eq.'qvdiv ') then !div-Q,
c                                                  10**-12(s**3 hPa)**-1
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         call bvfricalc(prs,ght,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,scr2a,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         ivar2=igdir(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,3,ivar2,0,0,0,0,100.,
     &      ihrip,rhrip,
     &      chrip,vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &      nxt,ncxc,maxtavl,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q-vector divergence'
         unwk(ipl)='s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-12
      elseif (cfeld(ifld,ipl)(1:5).eq.'qgomf') then !full omega, ubar/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         call bvfricalc(prs,ght,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,scr2a,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,0,1,0,0,0,0,100.,
     &      ihrip,rhrip,
     &      chrip,vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &      nxt,ncxc,maxtavl,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Full omega'
         unwk(ipl)='dPa s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'qgomg'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'qmomg') then! QG dp/dt, ubar/s
         read(cfeld(ifld,ipl)(6:6),'(i1)') iqvecforc
         read(cfeld(ifld,ipl)(7:7),'(i1)') itopobc
         read(cfeld(ifld,ipl)(8:8),'(i1)') iekmnbc
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call bvfricalc(prs,ght,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,scr2a,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         imo=0
         rhithresh=100.
         engplttl(ipl)='QG omega (dry)'
         if (cfeld(ifld,ipl)(2:2).eq.'m') then
            engplttl(ipl)='QG omega (moist)'
            imo=1
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
            call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
            if (iice.eq.1) then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qci       ',miy,mjx,mkzh,
     &              maxtavl,3,1,scr3b,istat)
               call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qsn       ',miy,mjx,mkzh,
     &              maxtavl,3,1,scr3b,istat)
               call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qgr       ',miy,mjx,mkzh,
     &              maxtavl,3,0,scr3b,istat)
               if (istat.ne.-1) call addorfill(scr3b,pl3,
     &            miy,mjx,mkzh,3,1,.001,1.)
            endif
            call bvfricalc(prs,ght,tmk,scr3c,qvp,pl3,
     &         uuu,vvv,scr2a,0.,0.,0.,'bvfsqi',1,scr3b,miy,mjx,mkzh)
            read(cfeld(ifld,ipl)(9:10),'(i2)') irhithresh
            rhithresh=float(irhithresh)
            call rhucalc(qvp,tmk,prs,1,scr3c,miy,mjx,mkzh)
         endif
         call omgcalc(qvp,tmk,www,prs,pl3,miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,4,1,iqvecforc,itopobc,iekmnbc,imo,rhithresh,
     &      ihrip,rhrip,
     &      chrip,vardesc,plchun,casename,iendc,cxtimeavl,xtimeavl,
     &      nxt,ncxc,maxtavl,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='dPa s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:2).eq.'se'.or.
     &        cfeld(ifld,ipl)(1:2).eq.'sm') then ! Sawyer-Eliassen
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3g,wk,maxslab)
c
c   Put v-geos. in 3a, v-ageos. in 3b, u on cross points in 3c,
c   omega in 3d, and d(ter)/dx in 2a.
c   Note, y axis is assumed to be toward left of cross sec.,
c   x axis is assumed to be into page.
c
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            tv=virtual(tmk(i,j,k),qvp(i,j,k))
            ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
            vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
            scr3a(i,j,k)=-cosa*ugeo-sina*vgeo
            utot=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &         uuu(i+1,j+1,k))
            vtot=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &         vvv(i+1,j+1,k))
            scr3b(i,j,k)=-cosa*utot-sina*vtot-scr3a(i,j,k)
            scr3c(i,j,k)=-sina*utot+cosa*vtot
         enddo
         enddo
         enddo
         do j=1,mjx-1
            jp1=min(j+1,mjx-1)
            jm1=max(j-1,1)
         do i=1,miy-1
            ip1=min(i+1,miy-1)
            im1=max(i-1,1)
            dxtr=ds/xmap(i,j)
            if (nproj.ne.4) then
               dytr=ds/xmap(i,j)
            else
               dytr=ds
            endif
            dx=dxtr*(jp1-jm1)
            dy=dytr*(ip1-im1)
            dterdx=(ter(i,jp1)-ter(i,jm1))/dx
            dterdy=(ter(ip1,j)-ter(im1,j))/dy
            scr2a(i,j)=-sina*dterdx+cosa*dterdy
         enddo
         enddo
         call omgcalc(qvp,tmk,www,prs,scr3d,miy,mjx,mkzh)
c
c      Get stability fields
c
         call bvfricalc(prs,ght,tmk,scr3f,qvp,scr3g,
     &      uuu,vvv,scr2b,0.,0.,0.,'bvfsqd    ',0,scr3e,
     &      miy,mjx,mkzh)
c
c   Set Saw.-El. parameters
c
         imaxlit=nint(rsepa(1,ipl))
         errmin=rsepa(2,ipl)
         alphor=rsepa(3,ipl)
         ixaverage=nint(rsepa(4,ipl))
         smfac=rsepa(5,ipl)
         imaxbig=nint(rsepa(6,ipl))
         rhithresh=rsepa(7,ipl)
         smfstb=rsepa(8,ipl)
         ilhs=nint(rsepa(9,ipl))
         if (cfeld(ifld,ipl)(2:2).eq.'m') then
            imo=1
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &           pl3,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qra       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3f,istat)
            call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
            if (iice.eq.1) then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qci       ',miy,mjx,mkzh,
     &              maxtavl,3,1,scr3f,istat)
               call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qsn       ',miy,mjx,mkzh,
     &              maxtavl,3,1,scr3f,istat)
               call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              'qgr       ',miy,mjx,mkzh,
     &              maxtavl,3,0,scr3f,istat)
               if (istat.ne.-1) call addorfill(scr3f,pl3,
     &            miy,mjx,mkzh,3,1,.001,1.)
            endif
            call bvfricalc(prs,ght,tmk,scr3g,qvp,pl3,
     &         uuu,vvv,scr2b,0.,0.,0.,'bvfsqi    ',1,scr3f,
     &         miy,mjx,mkzh)
            call rhucalc(qvp,tmk,prs,1,scr3g,miy,mjx,mkzh)
         else
            imo=0
            imaxbig=1
         endif
c
c   Calculate number of pressure levels.  We want grid
c   aspect ratio to be similar to typical frontal slope,
c   i.e. dz/dx ~ .01
c
         pbotse=-9e9
         ptopse=-9e9
         do j=1,mjx-1
         do i=1,miy-1
            pbotse=max(pbotse,prs(i,j,mkzh))
            ptopse=max(ptopse,prs(i,j,1))
         enddo
         enddo
         ptopse=ptopse+1.
         dpapprox=ds*xseclen/(nscrs-1.)*.01/13.
         mkp=nint((pbotse-ptopse)/dpapprox)+1
c
         call fillarray(pl3,miy*mjx*mkzh,0.)
         call saweli(ght,tmk,qvp,prs,scr3a,scr3b,scr3c,scr3d,
     &      scr3e,scr3f,scr3g,imo,ilhs,rhithresh,imaxlit,imaxbig,
     &      errmin,alphor,ter,scr2a,cor,xmap,rrfst(1,ipl),
     &      pl3,rcrag(1,ipl),rcrbg(1,ipl),ixaverage,smfac,smfstb,
     &      nscrs,nscrs,xdist,ydist,xseclen,ptopse,pbotse,
     &      cfeld(ifld,ipl),mkp,miy,mjx,mkzh)
c
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='unknown'
         if (cfeld(ifld,ipl)(3:5).eq.'vvv'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vge'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vag'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vab'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vtb') then
            unwk(ipl)='m s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'omf'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'omb') then
            unwk(ipl)='dPa s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'pvo') then
            unwk(ipl)='PVU'
         elseif (cfeld(ifld,ipl)(3:5).eq.'con'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'she'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'tot') then
            unwk(ipl)='m s~S~-2~N~ hPa~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'std'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'stm'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'ste'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'saj') then
            unwk(ipl)='m~S~2~N~/(s hPa)~S~2~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'rhi') then
            unwk(ipl)='%'
         elseif (cfeld(ifld,ipl)(3:5).eq.'ghp') then
            unwk(ipl)='m'
         elseif (cfeld(ifld,ipl)(3:5).eq.'alp') then
            unwk(ipl)='m~S~2~N~/(hPa s~S~2~N~)'
         elseif (cfeld(ifld,ipl)(3:5).eq.'rsw'.or.
     &           (cfeld(ifld,ipl)(3:4).eq.'fr'.and.
     &            cfeld(ifld,ipl)(5:5).ne.'6')) then
            unwk(ipl)='none'
         elseif (cfeld(ifld,ipl)(3:5).eq.'vor'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'vaj') then
            unwk(ipl)='s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'bcl'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'baj') then
            unwk(ipl)='m s~S~-1~N~ hPa~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'fst'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'faj') then
            unwk(ipl)='m s~S~-1~N~ hPa~S~-2~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'fr6'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'psi') then
            unwk(ipl)='m hPa s~S~-1~N~'
         endif
         engplttl(ipl)=cfeld(ifld,ipl)
      elseif (cfeld(ifld,ipl)(1:5).eq.'bvfsq'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'richn'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'spsq') then
c
c      Brunt-Vaisala frequency squared (per sec squared), Richardson
c      number (dim'less), or Scorer parameter squared (per km squared)
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call fillarray(scr3b,miy*mjx*mkzh,0.)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,0,pl3,
     &        istat)
         if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &      3,1,.001,1.)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,0,pl3,
     &        istat)
         if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &      3,1,.001,1.)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,0,
     &           pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,0,
     &           pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
         endif
         cosa=xdist/xseclen
         sina=ydist/xseclen
         cphase=cosa*rstmv(2,ipl)+sina*rstmv(1,ipl)
c         cphase=sqrt(rstmv(2,ipl)*rstmv(2,ipl)+
c     &      rstmv(1,ipl)*rstmv(1,ipl))
c         if (cphase.gt.0.) then
c            cosa=rstmv(2,ipl)/cphase
c            sina=rstmv(1,ipl)/cphase
c         else
c            cosa=1.
c            sina=0.
c         endif
         call bvfricalc(prs,ght,tmk,scr3a,qvp,scr3b,
     &      uuu,vvv,scr2a,cosa,sina,cphase,cfeld(ifld,ipl),0,pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         if (cfeld(ifld,ipl)(1:5).eq.'bvfsq') then
            engplttl(ipl)='Brunt-Vaisala freq. (squared)'
            unwk(ipl)='s~S~-2~N~'
         elseif (cfeld(ifld,ipl)(1:5).eq.'richn') then
            engplttl(ipl)='Richardson number'
            unwk(ipl)='none'
         elseif (cfeld(ifld,ipl)(1:4).eq.'spsq') then
            engplttl(ipl)='Scorer parameter (squared)'
            unwk(ipl)='km~S~-2~N~'
         endif
      elseif (cfeld(ifld,ipl)(1:3).eq.'rim') then
c
c      Various forms of Richardson number from MRF/HIRPBL schemes
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call fillarray(scr3a,miy*mjx*mkzh,0.)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,0,pl3,
     &        istat)
         if (istat.ne.-1) call addorfill(pl3,scr3a,miy,mjx,mkzh,
     &      3,1,.001,1.)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,0,
     &           pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3a,miy,mjx,mkzh,
     &         3,1,.001,1.)
         endif
         if (cfeld(ifld,ipl)(4:4).eq.'d') then
            im=1
         elseif (cfeld(ifld,ipl)(4:4).eq.'m') then
            im=2
         elseif (cfeld(ifld,ipl)(4:4).eq.'b') then
            im=3
         endif
         read(cfeld(ifld,ipl)(5:5),'(i1)') itype
         call ricalc(itype,im,prs,ght,tmk,qvp,scr3a,
     &      uuu,vvv,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Richardson number, '
     &      //cfeld(ifld,ipl)(4:4)//', '//cfeld(ifld,ipl)(5:5)
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rib ') then
c
c      Near-surface Richardson number (dim'less), used in HIRPBL scheme
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'tgk       ',miy,mjx,mkzh,maxtavl,2,1,
     &        scr2a,istat)
         do j=1,mjx-1
         do i=1,miy-1
            za=ght(i,j,mkzh)-ter(i,j)
            gammam=gamma*(1.+gammamd*qvp(i,j,mkzh))
            pfaca=(1000./prs(i,j,mkzh))**gammam
            psfc=sfp(i,j)
            pfacg=(1000./psfc)**gammam
            tha=tmk(i,j,mkzh)*pfaca
            thva=virtual(tmk(i,j,mkzh),qvp(i,j,mkzh))*pfaca
            thvg=virtual(scr2a(i,j),qvp(i,j,mkzh))*pfacg
            ucross=.25*(uuu(i,j,mkzh)+uuu(i+1,j,mkzh)+
     &                  uuu(i,j+1,mkzh)+uuu(i+1,j+1,mkzh))
            vcross=.25*(vvv(i,j,mkzh)+vvv(i+1,j,mkzh)+
     &                  vvv(i,j+1,mkzh)+vvv(i+1,j+1,mkzh))
            vsq=ucross*ucross+vcross*vcross
            pl2(i,j)=grav*za*(thva-thvg)/(tha*vsq)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='s~S~-2~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'condheat ') then !cond. htg., K/h
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (rsepa(1,ipl).eq.200..and.rsepa(2,ipl).eq.1..and.
     &       rsepa(3,ipl).eq.1.6) then
            irhthtype=0  ! use RH wrt liquid (0) or ice (1) for thresholding
            rhthresh=99. ! RH threshold
            ilrtype=0    ! pseudoadi. lapse rate is wrt liq (0) or ice (1)
         else
            irhthtype=rsepa(1,ipl)
            rhthresh=rsepa(2,ipl)
            ilrtype=rsepa(3,ipl)
         endif
         call condheat(tmk,qvp,www,irhthtype,ilrtype,0,prs,
     &      miy,mjx,mkzh,pl3,rhthresh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Condensational heating'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'condheat2') then !cond. htg., K/h
c                                                       !(direct from model)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'ttendbmp  ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &           istat)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            gammam=gamma*(1.+gammamd*qvp(i,j,k))
            pl3(i,j,k)=pl3(i,j,k)*(1000./prs(i,j,k))**gammam
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Expl. heating (from model)'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'vadvtheta') then !v.adv. of theta, K/hr
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call smoothcp(scr3a,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3b(i,j,k)=www(i,j,k)
         enddo
         enddo
         enddo
         call smoothcp(scr3b,1,scr3c,prs,pslab1,pslab2,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
         do k=1,mkzh
            kp1=min(mkzh,k+1)
            km1=max(1,k-1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=-.01*scr3b(i,j,k)*
     &         (scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &         (ght(i,j,kp1)-ght(i,j,km1))*3600.  ! K/s to K/hr
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vert. adv. of potential temperature'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'wsp ' .or. 
     &        cfeld(ifld,ipl)(1:5).eq.'wspk ') then! Hor wind speed, m/s or kt
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wspcalc(uuu,vvv,0.,0.,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind speed'
         if (cfeld(ifld,ipl)(1:4).eq.'wsp ') then
           unwk(ipl)='m s~S~-1~N~'
         else
           unwk(ipl)='kt'
           call addorfill(pl3,pl3,miy,mjx,mkzh,3,0,1.94,0.)
         endif
      elseif (cfeld(ifld,ipl)(1:6).eq.'wspsf ' .or.
     &        cfeld(ifld,ipl)(1:7).eq.'wspksf ') then! Surface Hor wind speed, m/s or kt
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'U10       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat2)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'V10       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2b,istat2)
c        print *, 'Reading U10 and V10, istat2 = ', istat2
         call wspcalc(scr2a,scr2b,0.,0.,rstmv(1,ipl),pl2,
     &      miy,mjx,1)
c        print *, 'U10,V10,wsp ',scr2a(60,80),scr2b(60,80),pl2(60,80)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Surface wind speed'
         if (cfeld(ifld,ipl)(1:5).eq.'wspsf') then
           unwk(ipl)='m s~S~-1~N~'
         else
           unwk(ipl)='kt'
           call addorfill(pl2,pl2,miy,mjx,1,3,0,1.94,0.)
         endif
      elseif (cfeld(ifld,ipl)(1:6).eq.'wsptr ' .or.
     &        cfeld(ifld,ipl)(1:7).eq.'wspktr '.or.     ! trop wind speed, m/s or kt
     &        cfeld(ifld,ipl)(1:6).eq.'wspmw ' .or.     
     &        cfeld(ifld,ipl)(1:7).eq.'wspkmw ' ) then ! max wind speed, m/s or kt
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         if ( index(cfeld(ifld,ipl)(1:7),'tr') .gt. 0 ) then
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'UUT       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat2)
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'VVT       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2b,istat2)
           engplttl(ipl)='Tropopause wind speed'
         else
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'UUW       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat2)
           call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'VVW       ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2b,istat2)
           engplttl(ipl)='Maximum wind speed'
         endif
         call wspcalc(scr2a,scr2b,0.,0.,rstmv(1,ipl),pl2,
     &      miy,mjx,1)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         if (cfeld(ifld,ipl)(1:5).eq.'wsptr' .or.
     &       cfeld(ifld,ipl)(1:5).eq.'wspmw' ) then
           unwk(ipl)='m s~S~-1~N~'
         else
           unwk(ipl)='kt'
           call addorfill(pl2,pl2,miy,mjx,1,3,0,1.94,0.)
         endif
      elseif (cfeld(ifld,ipl)(1:4).eq.'wdr ') then! Hor wind dir, deg.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wdircalc(uuu,vvv,unorth,vnorth,rstmv(1,ipl),
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind direction'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:4).eq.'unor') then! westerly wind comp., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=vnorth(i,j)*(uuu(i,j,k)-rstmv(2,ipl))-
     &                 unorth(i,j)*(vvv(i,j,k)-rstmv(1,ipl))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Westerly wind component'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vnor') then! southerly wind comp., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=unorth(i,j)*(uuu(i,j,k)-rstmv(2,ipl))+
     &                 vnorth(i,j)*(vvv(i,j,k)-rstmv(1,ipl))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Southerly wind component'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'xnt ') then! X-nrm tot wind, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,-sina,cosa,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'amt ') then! Abs. mom., m/s
c      Here amt is defined as total wind velocity into the page,
c      plus f (at middle of x-sec.) times left-to-right distance
c      along cross-section.  If the x-axis is taken to be parallel
c      to the cross section and pointing toward the right, then this
c      definition of abs. mom. corresponds to v+fx.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,-sina,cosa,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
         caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
         caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
         cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
         cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
         do j=1,mjx
         do i=1,miy
            dltr=ds*((j-caxgn)*cosa+(i-caygn)*sina)
            do k=1,mkzh
               pl3(i,j,k)=pl3(i,j,k)+fbar*dltr
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Absolute momentum'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'amtg ') then! Geost. abs. mom., m/s
c      Here amt is defined as geostrophic wind velocity into the page,
c      plus f (at middle of x-sec.) times left-to-right distance
c      along cross-section.  If the x-axis is taken to be parallel
c      to the cross section and pointing toward the right, then this
c      definition of abs. mom. corresponds to v_g+fx.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
c
         icen=nint(1.+(.5*(rcrag(1,ipl)+rcrbg(1,ipl))-xjcorn)*refrat)
         jcen=nint(1.+(.5*(rcrag(2,ipl)+rcrbg(2,ipl))-xjcorn)*refrat)
         fbar=cor(icen,jcen)
c
c      Put dp/dy in scr3a, then convert to u_geo
c
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               scr3a(i,j,k)=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(scr3a(1,1,k),miy,mjx)
         enddo
c
c      Put dp/dx in scr3b, then convert to v_geo
c
         call derivc(prs,1,ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               scr3b(i,j,k)=rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(scr3b(1,1,k),miy,mjx)
         enddo
c
c      Make geost. abs. mom.
c
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(scr3a,scr3b,-sina,cosa,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
         caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
         cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
         cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
         do j=1,mjx
         do i=1,miy
            dltr=ds*((j-caxgn)*cosa+(i-caygn)*sina)
            do k=1,mkzh
               pl3(i,j,k)=pl3(i,j,k)+fbar*dltr
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Absolute geostrophic momentum'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'xpt ') then! X-paral. wind, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,cosa,sina,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'stb ') then! -d(theta)/dp, K/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            scr3a(i,j,k)=-scr3a(i,j,k)
         enddo
         enddo
         enddo
         call ddpcalc(prs,scr3a,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='-d(theta)/dp '
         unwk(ipl)='K hPa~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'stbe ') then! -d(theta_e)/dp, K/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            scr3a(i,j,k)=-scr3a(i,j,k)
         enddo
         enddo
         enddo
         call ddpcalc(prs,scr3a,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='-d(theta_e)/dp '
         unwk(ipl)='K hPa~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'stbz ') then! d(theta)/dz, K/km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=1000.*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='d(theta)/dz'
         unwk(ipl)='K km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'dthtedz') then! d(theta_e)/dz, K/km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=1000.*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='d(theta-e)/dz'
         unwk(ipl)='K km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'brnshr ') then
c        Bulk Richardson Number Shear (a la Stensrud et al. 1997), m**2/s**2
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call pfcalc(prs,sfp,scr3a,miy,mjx,mkzh)
         call brnshr(uuu,vvv,ght,scr3a,ter,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m~S~2~N~ s~S~-2~N~'
c
c   Following are new fields added by J.F. Bresch between Dec/97 and Mar/00
c
      elseif (cfeld(ifld,ipl)(1:4).eq.'tsfc' .or.    ! surface temp, deg C
     &        cfeld(ifld,ipl)(1:4).eq.'tsfk' .or.    ! surface temp, K
     &        cfeld(ifld,ipl)(1:4).eq.'tsff' .or.    ! surface temp, deg F
     &        cfeld(ifld,ipl)(1:4).eq.'thsf') then   ! surface theta, K
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         if (iplevdata.le.3) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'tmk_sfan  ',miy,mjx,mkzh,maxtavl,2,0,
     &           pl2,istat)
            if (istat .lt. 0) then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &              't2        ',miy,mjx,mkzh,maxtavl,2,1,
     &              pl2,istat)    ! use t2 from metgrid
            endif
            if (cfeld(ifld,ipl)(1:4).eq.'thsf') then
               call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc, ! g/kg,
     &              'qvp_sfan  ',miy,mjx,mkzh,
     &              maxtavl,2,0,scr2a,istat) 
               if (istat .lt. 0) then
                 call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &               'q2        ',miy,mjx,mkzh,maxtavl,2,1,
     &               scr2a,istat)    ! use q2 from metgrid
               endif
               do j=1,mjx-1
               do i=1,miy-1
                  scr2a(i,j)=scr2a(i,j)*.001 ! g/kg to kg/kg
               enddo
               enddo
            endif
         else
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'tgk       ',miy,mjx,mkzh,maxtavl,2,1,
     &           pl2,istat)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'T2        ',miy,mjx,mkzh,maxtavl,2,0,
     &           scr2b,istat2)
            print *, 'Reading T2, istat2 = ', istat2
            if ( istat2 .gt. 0 ) print *, 'Found T2 and use it as sfc T'
            do j=1,mjx-1
            do i=1,miy-1
c
c            Added by WW: if T2 exists in model data, use it.
c            If air is warmer than ground, use air temp.
c            If air is colder than ground, use avg. of air temp.
c               and ground. temp.
c
               if ( istat2 .gt. 0 .and. xtime .gt. 0. ) then
                  pl2(i,j)=scr2b(i,j)
               else
                  if (pl2(i,j) .le. tmk(i,j,mkzh)) then
                     pl2(i,j)= tmk(i,j,mkzh)
                  else
                     pl2(i,j)=0.5*(tmk(i,j,mkzh)+pl2(i,j))
                  endif
               endif
               if (cfeld(ifld,ipl)(1:4).eq.'thsf')
     &            scr2a(i,j)=qvp(i,j,mkzh)
            enddo
            enddo
         endif
         if (cfeld(ifld,ipl)(1:4).eq.'tsfk') then
            engplttl(ipl)='Surface air temperature'
            unwk(ipl)='K'
         elseif (cfeld(ifld,ipl)(1:4).eq.'tsfc') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j) = pl2(i,j) - celkel
            enddo
            enddo
            engplttl(ipl)='Surface air temperature'
            unwk(ipl)='~S~o~N~C'
         elseif (cfeld(ifld,ipl)(1:4).eq.'tsff') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j) = (pl2(i,j)-celkel) * 1.8 +32.
            enddo
            enddo
            engplttl(ipl)='Surface air temperature'
            unwk(ipl)='~S~o~N~F'
         elseif (cfeld(ifld,ipl)(1:4).eq.'thsf') then
            do j=1,mjx-1
            do i=1,miy-1
               gammam=gamma*(1.+gammamd*scr2a(i,j))
               pl2(i,j)=pl2(i,j)*(1000./sfp(i,j))**gammam
            enddo
            enddo
            engplttl(ipl)='Surface air potential temperature'
            unwk(ipl)='K'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdf ') then! dewpoint, deg F
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           pl3(i,j,k) = pl3(i,j,k) * 1.8 + 32.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint temperature'
         unwk(ipl)='~S~o~N~F'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdk ') then! dewpoint, deg K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           pl3(i,j,k) = pl3(i,j,k) + celkel
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint temperature'
         unwk(ipl)='~S~o~N~F'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdd ') then! dewpoint depression, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           pl3(i,j,k) = tmk(i,j,k) - celkel - pl3(i,j,k) 
           if (pl3(i,j,k) .lt. 0. ) pl3(i,j,k) = 0.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint depression'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tfp ') then! frostpoint, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tfpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Frostpoint temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'thv ') then! virtual pot. temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
           tv=virtual(tmk(i,j,k),qvp(i,j,k))
           gammam=gamma*(1.+gammamd*qvp(i,j,k))
           pl3(i,j,k)=tv*(1000./prs(i,j,k))**gammam
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:6).eq.'thvhm ') then! vir. pot. tmp., incl.
c                                                   ! hydroms., K (or deg C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qra       ',miy,mjx,mkzh,maxtavl,3,1,pl3,
     &        istat)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'qcw       ',miy,mjx,mkzh,maxtavl,3,1,
     &        scr3a,istat)
         call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &      3,1,1.,1.)
         if (iice.eq.1) then
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qsn       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qci       ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &           'qgr       ',miy,mjx,mkzh,maxtavl,3,0,
     &           scr3a,istat)
            if (istat.ge.0) call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
         endif
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
            hydrometmr=.001*pl3(i,j,k) ! g/kg to kg/kg
            tv=virtualhyd(tmk(i,j,k),qvp(i,j,k),hydrometmr)
            gammam=gamma*(1.+gammamd*qvp(i,j,k))
            pl3(i,j,k)=tv*(1000./prs(i,j,k))**gammam
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vir. pot. temp. incl. hydrom.'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sreh') then! storm-rel helicity
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (cfeld(ifld,ipl)(5:5).eq.'1') then
           call relhl(uuu,vvv,ght,ter,1000.,pl2,miy,mjx,mkzh)
           engplttl(ipl)='Sfc-1 km Storm-Rel Helicity 75%:30R'
         else if (cfeld(ifld,ipl)(5:5).eq.'3') then
           call relhl(uuu,vvv,ght,ter,3000.,pl2,miy,mjx,mkzh)
           engplttl(ipl)='Sfc-3 km Storm-Rel Helicity 75%:30R'
         else
           call relhl(uuu,vvv,ght,ter,3000.,pl2,miy,mjx,mkzh)
           engplttl(ipl)='Sfc-3 km Storm-Rel Helicity 75%:30R'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='m~S~2~N~ s~S~-2~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:3).eq.'ehi') then! energy-helicity index
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
c
c      Attempt to read cape from a file, into scr2a.
c      If it's not there, then calculate it.
c
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'mcap      ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat)
         if (istat.ge.0) then
            write(iup,*) 'read in mcap'
         else
            call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
            call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &         miy,mjx,mkzh,0)
            do j = 1, mjx-1
            do i = 1, miy-1
               scr2a(i,j)=scr3a(i,j,mkzh)
            enddo
            enddo
         endif
c
         if (cfeld(ifld,ipl)(4:4).eq.'1') then
            write(iup,*) 'setting top to 1km'
            call relhl(uuu,vvv,ght,ter,1000.,pl2,miy,mjx,mkzh)
         else if (cfeld(ifld,ipl)(4:4).eq.'3') then
            call relhl(uuu,vvv,ght,ter,3000.,pl2,miy,mjx,mkzh)
         else
            call relhl(uuu,vvv,ght,ter,3000.,pl2,miy,mjx,mkzh)
         endif
         do j = 1, mjx-1
         do i = 1, miy-1
            pl2(i,j)= (scr2a(i,j) * pl2(i,j)) / 160000.
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Energy-helicity index'
         unwk(ipl)='none'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vgp') then! vorticity-gen. potential
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
c
c      Attempt to read cape from a file, into scr2a.
c      If it's not there, then calculate it.
c
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'mcap      ',miy,mjx,mkzh,maxtavl,2,0,
     &        scr2a,istat)
         if (istat.ge.0) then
            write(iup,*) 'read in mcap'
         else
            call pfcalc(prs,sfp,scr3c,miy,mjx,mkzh)
            call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &         miy,mjx,mkzh,0)
            do j = 1, mjx-1
            do i = 1, miy-1
               scr2a(i,j)=scr3a(i,j,mkzh)
            enddo
            enddo
         endif
c
         call vgp (uuu,vvv,ght,ter,scr2a,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vorticity generation potential'
         unwk(ipl)='m s~S~-2~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:6).eq.'refalt') then! 10dBZ altitude
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &        'REFL_10CM ',miy,mjx,mkzh,maxtavl,3,1,
     &           scr3b,istat)
         call reflalt (scr3b, prs, pl2, miy, mjx, mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Ten-dBZ altitude'
         unwk(ipl)='hPa'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:5).eq.'lifti') then! lifted index
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call lifted_index (tmk,qvp,prs,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Lifted index'
         unwk(ipl)='~S~o~N~C'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:5).eq.'tropp') then! tropopause pressue
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call trop (tmk,ght,prs,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Tropopause pressure'
         unwk(ipl)='hPa'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'srfl') then! storm-rel low-level inflow
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo(uuu,vvv,ght,ter,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Low-level storm-relative flow'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'srfh') then! storm-rel mid-level flow
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo(uuu,vvv,ght,ter,1,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Mid-level storm-relative flow'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'uubs') then! bulk shear, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (cfeld(ifld,ipl)(5:5).ne.' ') then
           read(cfeld(ifld,ipl),'(4x,i1)') llimit
         else
           llimit = 6
         endif
         call bshear (llimit,0,uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         write(engplttl(ipl),'("0-",i1," km shear (x-comp.)")') llimit
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
c        engplttl(ipl)='0-6 km shear (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvbs') then! bulk shear, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (cfeld(ifld,ipl)(5:5).ne.' ') then
           read(cfeld(ifld,ipl),'(4x,i1)') llimit
         else
           llimit = 6
         endif
         llimit = -1 * llimit
         call bshear (llimit,0,uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         write(engplttl(ipl),'("0-",i1," km shear (y-comp.)")')
     &             abs(llimit)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
c        engplttl(ipl)='0-6 km shear (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:3).eq.'bsh') then! Hor bulk wind shear, m/s or kt
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).ne.' ' .and.
     &       cfeld(ifld,ipl)(4:4).ne.'k') then
           read(cfeld(ifld,ipl),'(3x,i1)') llimit
         else
           llimit = 6
         endif
         call bshear (llimit,0,uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         llimit = -1 * llimit
         call bshear (llimit,0,uuu,vvv,ght,ter,scr2a,miy,mjx,mkzh)
         call wspcalc(pl2,scr2a,0.,0.,rstmv(1,ipl),pl2,
     &      miy,mjx,1)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         write(engplttl(ipl),'("0-",i1," km shear")')
     &               abs(llimit)
         if (cfeld(ifld,ipl)(4:4).eq.'k' .or.
     &       cfeld(ifld,ipl)(5:5).eq.'k') then
           unwk(ipl)='kt'
           call addorfill(pl2,pl2,miy,mjx,1,3,0,1.94,0.)
         else
           unwk(ipl)='m s~S~-1~N~'
         endif
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'uusr') then! suprcell mvmnt x-comp., m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,0,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Supercell motion (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvsr') then! suprcell mvmnt y-comp., m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,0,1,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Supercell motion (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:3).eq.'sr9') then! storm-rel. flow at 9km, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,1,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Supercell type (9-10 km rel. flow)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'cat ') then! Turbulence index
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call turb(prs,uuu,vvv,ght,xmap,dmap,tmk,qvp,
     &     scr2a,scr2b,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Clear air turbulence index'
         unwk(ipl)='s~S~-2~N~'
c
c   End of new fields added by J.F. Bresch between Dec '97 and Mar '00
c
c   20 April 2012: unity field for use with the addf option.
      elseif (cfeld(ifld,ipl)(1:5).eq.'unity') then   ! constant value of 1 
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx
         do i=1,miy
            pl2(i,j)=1.
         enddo
         enddo
         engplttl(ipl)='    '
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='  '
      else
        if ( istopmiss .eq. 1 ) then
          write(iup,*)'   Stopping in FIELDS: I don''t recognize',
     &      ' this field:'
          write(iup,*)'   <',cfeld(ifld,ipl),'>'
          stop
        else
          write(iup,*)'   In FIELDS: I don''t recognize',
     &      ' this field:'
          write(iup,*)'   <',cfeld(ifld,ipl),'>','  continuing anyway'
        endif
      endif
c
c   End of big "if" block
c
 1970 continue
c
      if (idimn(ipl).eq.3) then
         incwk=incwk+mkzh
      elseif (idimn(ipl).eq.2) then
         incwk=incwk+1
      else
         write(iup,*)'For ipl = ',ipl,'  idimn not 2 or 3.  idimn= ',
     &      idimn(ipl)
         stop
      endif
c
c   Do calibration, if asked for
c
      if (ccalb(ipl)(1:5).ne.'junk ') then
         open(unit=77,file=ccalb(ipl),form='formatted',status='old')
         read(77,*)vvals
         nvals=nint(vvals)
         read(77,*)xx1,xx2
         dxx=xx2-xx1
         read(77,*)(calb(i),i=1,nvals)
         close(77)
         do j=1,mjx-icdwk(ipl)
         do i=1,miy-icdwk(ipl)
            if (idimn(ipl).eq.3) then
               do k=1,mkzh
                  findx = 1.+(pl3(i,j,k)-xx1)/dxx
                  indx = max(1,min(nvals-1,nint(findx-0.5)))
                  rem=findx-float(indx)
                  pl3(i,j,k)=calb(indx)+rem*(calb(indx+1)-calb(indx))
               enddo
            else
               findx = 1.+(pl2(i,j)-xx1)/dxx
               indx = max(1,min(nvals-1,nint(findx-0.5)))
               rem=findx-float(indx)
               pl2(i,j)=calb(indx)+rem*(calb(indx+1)-calb(indx))
            endif
         enddo
         enddo
      endif
c
c   Do horizontal gradient, if asked for, and if field is 3D
c
      if (lgrad(ipl).or.llapl(ipl)) then
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         ipass=0
 543     ipass=ipass+1
         if (idimn(ipl).eq.3) then
            call derivc(pl3,icdwk(ipl),ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3a,icdwk(ipl),'x',miy,mjx,mkzh)
            if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &         (scr3b,pl3,miy,mjx,mkzh,3,icdwk(ipl),1.,0.)
            call derivc(pl3,icdwk(ipl),ght,xmap,dmap,qvp,tmk,
     &         scr2a,scr2b,scr3b,icdwk(ipl),'y',miy,mjx,mkzh)
         else
            if (icdwk(ipl).eq.0) then
               call ddx(pl2,icdwk(ipl),scr2a,xmap,dmap,
     &            icdwk(ipl),miy,mjx)
               if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &            (scr2b,pl2,miy,mjx,mkzh,2,icdwk(ipl),1.,0.)
               call ddy(pl2,icdwk(ipl),scr2b,xmap,dmap,
     &            icdwk(ipl),miy,mjx)
            elseif (icdwk(ipl).eq.1) then
               call ddx(pl2,icdwk(ipl),scr2a,xmap,dmap,
     &            icdwk(ipl),miy,mjx)
               if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &            (scr2b,pl2,miy,mjx,mkzh,2,icdwk(ipl),1.,0.)
               call ddy(pl2,icdwk(ipl),scr2b,xmap,dmap,
     &            icdwk(ipl),miy,mjx)
            endif
         endif
         if (igdir(ipl).ge.0.and.igdir(ipl).le.360) then ! gradient or
c                                            2nd deriv. in specified dir
c            cosa=cos(rpd*(90-igdir(ipl)))
c            sina=sin(rpd*(90-igdir(ipl)))
            print*,'Setting gdir to a compass direction has been'
            print*,'disabled because it doesn''t work properly.'
            print*,'The way the code is written, gdir between'
            print*,'0 and 360 specifies a direction like compass'
            print*,'direction, but relative to the y-axis of the grid'
            print*,'rather than to north.'
            stop
         elseif (igdir(ipl).eq.361) then  ! mag. of grad., or total lapl.
            cosa=0.
            sina=0.
         elseif (igdir(ipl).eq.362) then  ! grad. or 2nd drv along cross sec.
            cosa=xdist/xseclen
            sina=ydist/xseclen
         elseif (igdir(ipl).eq.363) then  ! grad. or 2nd drv perp to cross sec.
            cosa=-ydist/xseclen
            sina=xdist/xseclen
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (idimn(ipl).eq.3) then
               if (cosa.eq.0..and.sina.eq.0.) then
                  if (llapl(ipl).and.ipass.eq.1) then
                     pl3(i,j,k)=scr3a(i,j,k)
                  else
                     pl3(i,j,k)=sqrt(scr3a(i,j,k)*scr3a(i,j,k)+
     &                  scr3b(i,j,k)*scr3b(i,j,k))
                  endif
               else
                  pl3(i,j,k)=cosa*scr3a(i,j,k)+sina*scr3b(i,j,k)
               endif
            elseif (idimn(ipl).eq.2) then
               if (cosa.eq.0..and.sina.eq.0.) then
                  if (llapl(ipl).and.ipass.eq.1) then
                     pl2(i,j)=scr2a(i,j)
                  else
                     pl2(i,j)=sqrt(scr2a(i,j)*scr2a(i,j)+
     &                  scr2b(i,j)*scr2b(i,j))
                  endif
               else
                  pl2(i,j)=cosa*scr2a(i,j)+sina*scr2b(i,j)
               endif
            endif
         enddo
         enddo
         enddo
         if (llapl(ipl).and.ipass.eq.1) goto 543
      endif
c
c   Do advection, if asked for, and if field is 3D
c
      if (lhadv(ipl)) then
         if (idimn(ipl).ne.3) then
            write(iup,*)'Can only do advection for 3-D fields.'
         endif
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivc(pl3,icdwk(ipl),ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3a,1,'x',miy,mjx,mkzh)
         call derivc(pl3,icdwk(ipl),ght,xmap,dmap,qvp,tmk,
     &      scr2a,scr2b,scr3b,1,'y',miy,mjx,mkzh)
         if (igdir(ipl).ge.0.and.igdir(ipl).le.360) then ! adv. in 
c                                                          specified dir.
c            cosa=cos(rpd*(90-igdir(ipl)))
c            sina=sin(rpd*(90-igdir(ipl)))
            print*,'Setting gdir to a compass direction has been'
            print*,'disabled because it doesn''t work properly.'
            print*,'The way the code is written, gdir between'
            print*,'0 and 360 specifies a direction like compass'
            print*,'direction, but relative to the y-axis of the grid'
            print*,'rather than to north.'
            stop
         elseif (igdir(ipl).eq.361) then  ! total adv.
            cosa=0.
            sina=0.
         elseif (igdir(ipl).eq.362) then  ! adv. along cross sec.
            cosa=xdist/xseclen
            sina=ydist/xseclen
         elseif (igdir(ipl).eq.363) then  ! adv. perp. to cross sec.
            cosa=-ydist/xseclen
            sina=xdist/xseclen
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            uuucross=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &         uuu(i+1,j+1,k))
            vvvcross=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &         vvv(i+1,j+1,k))
            if (cosa.eq.0..and.sina.eq.0.) then
                  pl3(i,j,k)=
     &               -uuucross*scr3a(i,j,k)-vvvcross*scr3b(i,j,k)
            else
               pl3(i,j,k)=-(cosa*uuucross+sina*vvvcross)*
     &            (cosa*scr3a(i,j,k)+sina*scr3b(i,j,k))
            endif
         enddo
         enddo
         enddo
      endif
c
c   Do constant-pressure smoothing if asked for.
c
      if (ismcp(ipl).gt.0) then
         if (idimn(ipl).ne.3) then
            write(iup,*)'You can only do constant-pres.',
     &         ' smoothing of a 3-d field.'
            write(iup,*)'  ipl,ifld,cfeld=',ipl,ifld,cfeld(ifld,ipl)
            stop
         endif
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call smoothcp(pl3,icdwk(ipl),
     &      scr3a,prs,pslab1,pslab2,ismcp(ipl),miy,mjx,mkzh,
     &      mabpl,morpl)
      endif
c
c   Do time/model run differencing, if requested.
c
      if (cdiff(ipl)(1:5).ne.'none '.or.rdiff(ipl).ne.rmsg) then
c
      if (idiffpass.eq.2) goto 779
c
      ldiffsuccess=.true.
c
c   Save various parameters from present data set and time.
c
      xtime_sv=xtime
      casename_sv=casename
      iendc_sv=iendc
      nxtavl_sv=nxtavl
      ncxc_sv=ncxc
      do i=1,maxtavl
         cxtimeavl_sv(i)=cxtimeavl(i)
         xtimeavl_sv(i)=xtimeavl(i)
      enddo
      nxt_sv=nxt
c
c  Set new xtime
c
      if (rdiff(ipl).ne.rmsg) then
         if (ldfrl(ipl)) then
            xtime=xtime+rdiff(ipl)
         else
            xtime=rdiff(ipl)
         endif
      endif
c
c   Set new case name, load new available times
c
      if (cdiff(ipl)(1:5).ne.'none ') then
         casename=cdiff(ipl)
         iendc=index(casename,' ')-1
         call gettimes(casename,iendc,xtimeavl,cxtimeavl,nxtavl,
     &      maxtavl,ncxc,iup)
      endif
c
c   Determine if requested time is available in requested data set.
c
      iavail=0
      do i=1,nxtavl
         if (abs(xtimeavl(i)-xtime).le.tacch) then
            iavail=1
            nxt=i
            goto 640
         endif
      enddo
 640  continue
      if (iavail.eq.0) then
         write(iup,*)'   Requested time ',xtime,' is not available'
         write(iup,*)'   in the requested difference data set.'
         write(iup,*)'   RIP will plot undifferenced field.'
         ldiffsuccess=.false.  ! This is so difference info won't be
c                              ! printed in routine pltitle
c
c      If time was not available, restore some things and get out
c
         xtime=xtime_sv
         nxt=nxt_sv
         if (cdiff(ipl)(1:5).ne.'none ') then
            casename=casename_sv
            iendc=iendc_sv
            nxtavl=nxtavl_sv
            ncxc=ncxc_sv
            do i=1,maxtavl
               cxtimeavl(i)=cxtimeavl_sv(i)
               xtimeavl(i)=xtimeavl_sv(i)
            enddo
         endif
         goto 781
      endif
c
c   Get new header information and check domain compatibility
c      (only if doing case differencing)
c
      if (cdiff(ipl)(1:5).ne.'none ') then
c
      call getdims(casename,iendc,xtimeavl,cxtimeavl,
     &   maxtavl,ncxc,miy_ch,mjx_ch,mkzh_ch,mabpl_ch,morpl_ch,iup)
c
c      As far as compatibility between the original domain and the domain to
c      be subtracted, the only thing that we will require is that the
c      dimensions be the same in all three directions.  All other aspects
c      of the subtraction data set can (in principal) be different from the
c      original data set.
c
         if(miy.ne.miy_ch.or.mjx.ne.mjx_ch.or.mkzh.ne.mkzh_ch) then
            write(iup,*)'   Subtraction domain dimensions are ',
     &            miy_ch,mjx_ch,mkzh_ch,'.'
            write(iup,*)'   This is incompatible with original domain'
            write(iup,*)'   whose dimensions are ',miy,mjx,mkzh
            write(iup,*)'   RIP will plot undifferenced field.'
            ldiffsuccess=.false.  ! This is so difference info won't be
c                                 ! printed in routine pltitle
c
c         If domain was incompatible, restore some things and get out
c
            xtime=xtime_sv
            casename=casename_sv
            iendc=iendc_sv
            nxtavl=nxtavl_sv
            ncxc=ncxc_sv
            do i=1,maxtavl
               cxtimeavl(i)=cxtimeavl_sv(i)
               xtimeavl(i)=xtimeavl_sv(i)
            enddo
            nxt=nxt_sv
            goto 781
         endif
c
c      Get information from new header record
c
         call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &      nxt,maxtavl,ncxc,nproj,miycors,mjxcors,mdateb,mhourb,iice,
     &      iplevdata,true1,true2,xlatc,
     &      xlonc,dskmc,dskm,yicorn,xjcorn,rhourb,dsc,ds,refrat,iup)
c
      endif
c
c   Read the basic fields.
c
      call getbasicvars(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,maxtavl,uuu,vvv,tmk,qvp,www,prs,ght,sfp,sfpsm,ter,
     &   dmap,xmap,cor,dskm,miy,mjx,mkzh)
c
c   Load/calculate field from new time/model run
c
      goto 35
c
 779  continue
c
c   Subtract the fields
c
      if (idimn(ipl).eq.3) then
c
c      New field is at slab marker incwk-mkzh, and is already in
c      pointed array pl3.  Old field is at slab marker incwk-2*mkzh.
c      Set up pointer so that old field is in pointed array scr3a.
c
         ifree=incwk-2*mkzh
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c
c      Put subtracted field in scr3a
c
         if (iovly(ipl) .eq. 0) then   ! difference field (default)
         do k=1,mkzh
         do j=1,mjx-icdwk(ipl)
         do i=1,miy-icdwk(ipl)
            scr3a(i,j,k)=scr3a(i,j,k)-pl3(i,j,k)
         enddo
         enddo
         enddo
         else                          ! overlay field
         do k=1,mkzh
         do j=1,mjx-icdwk(ipl)
         do i=1,miy-icdwk(ipl)
            scr3a(i,j,k)=pl3(i,j,k)
         enddo
         enddo
         enddo
         endif
c
c      Adjust slab markers
c
         indwk(ifld,ipl)=incwk-2*mkzh  ! where subtracted data is now located
         incwk=incwk-mkzh ! Reset to end of subtracted data.
c
c      Ready for next field.  Data from new time/model run will be overwritten
c
      elseif (idimn(ipl).eq.2) then
c
c      New field is at slab marker incwk-1, and is already in
c      pointed array pl2.  Old field is at slab marker incwk-2.
c      Set up pointer so that old field is in pointed array scr2a.
c
         ifree=incwk-2
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
c
c      Put subtracted field in scr2a
c
         if (iovly(ipl) .eq. 0) then   ! difference field (default)
           do j=1,mjx-icdwk(ipl)
           do i=1,miy-icdwk(ipl)
              scr2a(i,j)=scr2a(i,j)-pl2(i,j)
           enddo
           enddo
         else                          ! overlay field
           do j=1,mjx-icdwk(ipl)
           do i=1,miy-icdwk(ipl)
              scr2a(i,j)=pl2(i,j)
           enddo
           enddo
         endif
c
c      Adjust slab markers
c
         indwk(ifld,ipl)=incwk-2  ! where subtracted data is now located
         incwk=incwk-1 ! Reset to end of subtracted data.
c
c      Data from new time/model run will be overwritten by next field.
c
      endif
c
c   Restore original parameters and data.
c
      xtime=xtime_sv
      nxt=nxt_sv
      if (cdiff(ipl)(1:5).ne.'none ') then
         casename=casename_sv
         iendc=iendc_sv
         nxtavl=nxtavl_sv
         ncxc=ncxc_sv
         do i=1,maxtavl
            cxtimeavl(i)=cxtimeavl_sv(i)
            xtimeavl(i)=xtimeavl_sv(i)
         enddo
         call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &      nxt,maxtavl,ncxc,nproj,miycors,mjxcors,mdateb,mhourb,iice,
     &      iplevdata,true1,true2,xlatc,
     &      xlonc,dskmc,dskm,yicorn,xjcorn,rhourb,dsc,ds,refrat,iup)
      endif
      call getbasicvars(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,maxtavl,uuu,vvv,tmk,qvp,www,prs,ght,sfp,sfpsm,ter,
     &   dmap,xmap,cor,dskm,miy,mjx,mkzh)
c
 781  continue
c
      endif  ! end of time/model run differencing
c
 1980 continue
c
c   Add the previous field to this field if asked for.
c
      if (ipl.gt.iplstrt.and.raddf(ipl).eq.0.0) then
      if (raddf(ipl-1).ne.0.0) then
         if (cfeld(ifld,ipl)(1:2).ne.'se'.and.
     &       cfeld(ifld,ipl)(1:2).ne.'sm') then
            icda=icdwk(ipl)
         else
            icda=0
         endif
         do ipla=ipl-1,iplstrt,-1
            if (raddf(ipla).eq.0.0) then
               goto 1900
            elseif (raddf(ipla).eq.1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)+
     &               wk(i,j,kpla)
               enddo
               enddo
               enddo
            elseif (raddf(ipla).eq.-1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)-
     &               wk(i,j,kpla)
               enddo
               enddo
               enddo
            else
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)+
     &               raddf(ipla)*wk(i,j,kpla)
               enddo
               enddo
               enddo
            endif
         enddo
 1900    continue
      endif
      endif
c
c   Save the field to a file, if asked for.
c
      if (csave(ipl).ne.'dontsave  ') then
         vardesc=engplttl(ipl)
         if (csave(ipl).eq.'sameasfeld') then
            varname=cfeld(ifld,ipl)
         else
            varname=csave(ipl)
         endif
         plchun=unwk(ipl)
         call writefile (wk(1,1,indwk(ifld,ipl)),varname,0,
     &      idimn(ipl),icdwk(ipl),vardesc,plchun,ihrip,rhrip,
     &      chrip,casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &      maxtavl,miy,mjx,mkzh)
      endif
c
 2000 continue
c
c   Set all levels to 1 for 2-d variables.
c
      if (idimn(ipl).eq.2) then
         do ilev=1,maxlev
            rlevl(ilev,ipl)=1
            rlavl(ilev,ipl)=1
         enddo
      endif
c
      return
      end
