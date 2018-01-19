c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pltitle(ctitl,cfeld,cptyp,rlevl,rlavl,cvcor,
     &   rstmv,lgrad,llapl,lhadv,rcrag,rcrbg,rslcg,ismth,
     &   ismcp,toptextclg,ilev,ixavg,icomg,csout,idimn,raddf,
     &   rtjst,rtjen,rtjti,titlstr,iovly,
     &   cdiff,rdiff,ldfrl,xtime,itjns,rtim,ctim,lnttl,lnsmm,lnvlb,
     &   ldiffsuccess,idescriptive,engplttl,maxlev,maxpl,mkzh,ipl,
     &   noplots)
c
      dimension rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),
     &   rcrag(2,maxpl),rcrbg(2,maxpl),ismth(maxpl),iovly(maxpl),
     &   ixavg(maxpl),icomg(maxpl),rstmv(2,maxpl),idimn(maxpl),
     &   raddf(maxpl),ismcp(maxpl),rslcg(2,maxpl),rdiff(maxpl),
     &   rtjst(maxpl),rtjen(maxpl),rtjti(maxpl),itjns(maxpl)
      character cfeld(3,maxpl)*10,cptyp(maxpl)*2,cvcor(maxpl)*1,
     &   csout(maxpl)*58,str*82,field1*10,field2*10,field3*10,
     &   ctitl(maxpl)*82,engplttl(maxpl)*36,cdiff(maxpl)*256,
     &   titlstr*82
      logical lnttl(maxpl),lnsmm(maxpl),lnvlb(maxpl),
     &   ldfrl(maxpl),lgrad(maxpl),llapl(maxpl),lhadv(maxpl),
     &   ldiffsuccess
c
      include 'comconst'
c
      chwd=.011
      chht=1.155*chwd
      chgap=.6*chht
c
c   Make set call, set color
c
      if (noplots.eq.0) then
        call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
        call gsplci(icomg(ipl))
        call gstxci(icomg(ipl))
      endif
c
      if (ctitl(ipl)(1:20).ne.'auto                ') then
         str=ctitl(ipl)
         do i=1,82
            if (str(i:i).eq.'_') str(i:i)=' '
         enddo
         goto 60
      endif
c
      idescr=idescriptive
      if (engplttl(ipl)(1:10).eq.'          ') idescr=0
c
      field1=cfeld(1,ipl)
      field2=cfeld(2,ipl)
      field3=cfeld(3,ipl)
      str=' '
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
c
c   Write title for individual plot
c
      if (cptyp(ipl)(2:2).eq.'t') then
         iend=0
         if ((itjns(ipl).gt.0.or.rstmv(1,ipl).ne.0..or.
     &        rstmv(2,ipl).ne.0.).and.field1.ne.'circle    ') then
            str(iend+1:iend+15)='Storm-relative '
            iend=iend+15
         endif
         if (rtjst(ipl).eq.rmsg) then
            starttim=min(rtim,ctim)
         else
            starttim=max(rtjst(ipl),min(rtim,ctim))
         endif
         if (rtjen(ipl).eq.rmsg) then
            endtim=max(rtim,ctim)
         else
            endtim=min(rtjen(ipl),max(rtim,ctim))
         endif
         if (field1.eq.'ribbon    '.or.
     &       field1.eq.'arrow     ') then
            str(iend+1:iend+41)='Trajectories from hour '//
     &         'xxx.xxx to xxx.xxx'
            write(str(iend+24:iend+30),'(f7.3)') starttim
            write(str(iend+35:iend+41),'(f7.3)') endtim
            iend=iend+41
         elseif (field1.eq.'circle    ') then
            str(iend+1:iend+50)='Net trajectory ascent from hour '//
     &         'xxx.xxx to xxx.xxx'
            write(str(iend+33:iend+39),'(f7.3)') starttim
            write(str(iend+44:iend+50),'(f7.3)') endtim
            iend=iend+50
         elseif (field1.eq.'swarm     '.or.
     &           field1.eq.'gridswarm ') then
            if (endtim.ge.starttim) then
               isign=1
            else
               isign=-1
            endif
            str(iend+1:iend+16)='Swarms at hours '
            iend=iend+16
            xtimetrj=starttim
 43         continue
            if (xtimetrj.lt.10.) then
               write(str(iend+1:iend+6),'(f5.3,a1)') xtimetrj,','
               iend=iend+6
            elseif (xtimetrj.ge.10..and.xtimetrj.lt.100.) then
               write(str(iend+1:iend+7),'(f6.3,a1)') xtimetrj,','
               iend=iend+7
            elseif (xtimetrj.ge.100..and.xtimetrj.lt.1000.) then
               write(str(iend+1:iend+8),'(f7.3,a1)') xtimetrj,','
               iend=iend+8
            endif
            xtimetrj=xtimetrj+isign*rtjti(ipl)
            if ((isign.eq. 1.and.xtimetrj.le.endtim).or.
     &          (isign.eq.-1.and.xtimetrj.ge.endtim)) goto 43
            str(iend:iend)=' '
            iend=iend-1
         endif
         if (cptyp(ipl)(1:1).eq.'v') then
            write(str(iend+2:iend+31),50)caxgn,caygn,cbxgn,cbygn
         endif
         goto 60
      endif
c
      ie1=lennonblank(field1)
      ie2=lennonblank(field2)
      ie3=lennonblank(field3)
c
      str=' '
      if (cptyp(ipl)(1:1).eq.'s') then   ! sounding
         if (cptyp(ipl)(2:2).eq.'c') then   ! contours
            str(1:ie1)=field1(1:ie1)
         elseif (cptyp(ipl)(2:2).eq.'v') then    ! vectors
            str(1:ie1+ie2+3)='<'//field1(1:ie1)//','//
     &         field2(1:ie2)//'>'
         endif
         if (idescr.eq.1) then
            if (cptyp(ipl)(2:2).eq.'c') then
               str(1:23)=engplttl(ipl)(1:23)
            elseif (cptyp(ipl)(2:2).eq.'v') then
               isv=index(engplttl(ipl),'-comp.')-2
               if (isv.gt.0) then
                  str(1:23)=engplttl(ipl)(1:23)
                  str(isv:23)='vectors'
               endif
            endif
         endif
         sxgn=1.+(rslcg(2,ipl)-xjcorn)*refrat
         sygn=1.+(rslcg(1,ipl)-yicorn)*refrat
         write(csout(ipl)(5:10),'(f6.2)') sxgn
         write(csout(ipl)(12:17),'(f6.2)') sygn
         write(str(25:82),'(a58)') csout(ipl)
         goto 60
      endif
c
      if (cptyp(ipl)(2:2).eq.'c'.or.
     &    cptyp(ipl)(2:2).eq.'h'   ) then   ! contours or characters
         str(1:ie1)=field1(1:ie1)
         if (idescr.eq.1) then
            str(1:36)=engplttl(ipl)
         endif
      elseif (cptyp(ipl).eq.'vv') then    ! vertical vector in x-sec.
         str(1:ie1+ie2+ie3+12)='<'//field1(1:ie1)//','//
     &         field2(1:ie2)//','//field3(1:ie3)//'> Vectors'
         if (idescr.eq.1) then
            if (field1(1:4).eq.'uuu '.and.field2(1:4).eq.'vvv '.and.
     &          (field3(1:4).eq.'sgd '.or.field3(1:4).eq.'omg '.or.
     &           field3(1:4).eq.'www ')) then
               str(1:36)='Circulation vectors'
            elseif (field1(1:6).eq.'uageo '.and.
     &              field2(1:6).eq.'vageo '.and.
     &          (field3(1:4).eq.'sgd '.or.field3(1:4).eq.'omg '.or.
     &           field3(1:4).eq.'www ')) then
               str(1:36)='Ageost. circulation vectors'
            endif
         endif
      elseif (cptyp(ipl).eq.'hv') then    ! horiz. winds on hor. plot
         str(1:ie1+ie2+11)='<'//field1(1:ie1)//','//
     &         field2(1:ie2)//'> Vectors'
         if (idescr.eq.1) then
            isv=index(engplttl(ipl),'-comp.')-2
            if (isv.gt.0) then
               str(1:36)=engplttl(ipl)
               str(isv:36)='vectors'
            endif
         endif
      elseif (cptyp(ipl).eq.'vw') then    ! horiz. winds in x-sec
         str(1:ie1+ie2+18)='<'//field1(1:ie1)//','//
     &         field2(1:ie2)//'> Horiz. Vectors'
         if (idescr.eq.1) then
            isv=index(engplttl(ipl),'-comp.')-2
            if (isv.gt.0) then
               str(1:36)=engplttl(ipl)
               str(isv:36)='vectors'
            endif
         endif
      elseif (cptyp(ipl)(2:2).eq.'s') then   ! streamlines
         str(1:ie1+ie2+15)='<'//field1(1:ie1)//','//
     &         field2(1:ie2)//'> Streamlines'
         if (idescr.eq.1) then
            isv=index(engplttl(ipl),'-comp.')-2
            if (isv.gt.0) then
               str(1:36)=engplttl(ipl)
               str(isv:36)='streamlines'
            endif
         endif
      endif
c
      if (raddf(ipl).ne.0.0) then
         write(str(39:70),'(a25,f6.2,a1)')
     &      ' (Added field, factor of ',raddf(ipl),')'
         goto 300
      endif
c
      if (cptyp(ipl)(1:1).eq.'h'.and.
     &    idimn(ipl).eq.2) goto 300 ! hor. plot, 2D field
c
      if (cptyp(ipl)(1:1).eq.'h'.and.
     &    rlevl(ilev,ipl).eq.rlavl(ilev,ipl)) then ! hor. plts, no avg.

         if(field1(1:6).eq.'track ') then
           write(str(1:20),'(a20)') 'Typhoon Track       '

         else if (cvcor(ipl).eq.'s') then   ! k index
            write(str(39:57),'(a16,i3)')    '   at k-index = ',
     &                    nint(rlevl(ilev,ipl))
         elseif (cvcor(ipl).eq.'z') then    ! height (geop.)
            heit=rlevl(ilev,ipl)
            write(str(39:60),'(a14,f5.2,a3)') '  at height = ',
     &                                   heit,' km'
         elseif (cvcor(ipl).eq.'f') then    ! height (geop.) above frz. lev.
            heit=rlevl(ilev,ipl)
            write(str(39:60),'(a14,f5.2,a3)') '  at h(AFL) = ',
     &                                   heit,' km'
         elseif (cvcor(ipl).eq.'p'.or.
     &           cvcor(ipl).eq.'l'.or.
     &           cvcor(ipl).eq.'x') then   ! pressure
            if (nint(rlevl(ilev,ipl)).eq.1001) then
               write(str(39:55),'(a10)')   'at surface'
            else
               write(str(39:60),'(a14,i4,a4)')   'at pressure = ',
     &                          nint(rlevl(ilev,ipl)),' hPa'
            endif
         elseif (cvcor(ipl).eq.'t') then   ! theta
            write(str(39:58),'(a14,i3,a3)')   '   at theta = ',
     &                          nint(rlevl(ilev,ipl)),' K '
         elseif (cvcor(ipl).eq.'m') then   ! temperature
            write(str(39:56),'(a10,i3,a5)')   '   at T = ',
     &                          nint(rlevl(ilev,ipl)),' dg C'
         elseif (cvcor(ipl).eq.'e') then   ! theta_e
            write(str(39:58),'(a14,i3,a3)')   ' at theta_e = ',
     &                          nint(rlevl(ilev,ipl)),' K '
         elseif (cvcor(ipl).eq.'q') then   ! PV
            write(str(39:60),'(a14,f4.1,a4)') '      at PV = ',
     &                                rlevl(ilev,ipl),' PVU'
         endif
      elseif (cptyp(ipl)(1:1).eq.'h'.and.
     &    rlevl(ilev,ipl).ne.rlavl(ilev,ipl).and.
     &    rlavl(ilev,ipl).ge.0.) then  ! hor. plts w/ avg.
         if (cvcor(ipl).eq.'s') then   ! k-index
            write(str(39:63),'(a15,i3,a4,i3)')    'Avg, k-index = ',
     &           nint(rlevl(ilev,ipl)),' to ',
     &           nint(rlavl(ilev,ipl))
c         elseif (cvcor(ipl).eq.'z') then   ! height (geop.)
c            heit=rlevl(ilev,ipl)
c            heita=rlavl(ilev,ipl)
c            write(str(39:66),'(a11,f5.2,a4,f5.2,a3)') 'Avg, hgt = ',
c     &                                   heit,' to ',heita,' km'
c         elseif (cvcor(ipl).eq.'p'.or.
c     &           cvcor(ipl).eq.'l'.or.
c     &           cvcor(ipl).eq.'x') then   !pressure
c            write(str(39:65),'(a11,i4,a4,i4,a4)')     'Avg, prs = ',
c     &         nint(rlevl(ilev,ipl)),' to ',nint(rlavl(ilev,ipl)),' hPa'
c         elseif (cvcor(ipl).eq.'t') then   ! theta
c            write(str(39:62),'(a11,i3,a4,i3,a3)')     'Avg, tht = ',
c     &        nint(rlevl(ilev,ipl)),' to ',nint(rlavl(ilev,ipl)),' K '
c         elseif (cvcor(ipl).eq.'e') then   ! theta_e
c            write(str(39:62),'(a11,i3,a4,i3,a3)')     'Avg, th_e= ',
c     &        nint(rlevl(ilev,ipl)),' to ',nint(rlavl(ilev,ipl)),' K '
         endif
      elseif (cptyp(ipl)(1:1).eq.'h'.and.
     &    rlevl(ilev,ipl).ne.rlavl(ilev,ipl).and.
     &    rlavl(ilev,ipl).lt.0.) then  ! hor. plts w/ vert.-level differencing
         if (cvcor(ipl).eq.'s') then   !k-index
            write(str(39:69),'(a20,i3,a5,i3)')  'Diff betw k-index = ',
     &           nint(rlevl(ilev,ipl)),' and ',
     &           nint(-rlavl(ilev,ipl))
         else
            write(iup,*) 'RIP is only set up to do averaging when'
            write(iup,*) 'the vertical coordinate is k-index.'
            stop
         endif
      elseif (cptyp(ipl)(1:1).eq.'v'.and.
     &      ixavg(ipl).eq.0) then      ! vert. plot w/out averaging
         write(str(39:68),50)caxgn,caygn,cbxgn,cbygn
      elseif (cptyp(ipl)(1:1).eq.'v'.and.
     &      ixavg(ipl).ne.0) then      ! vert. plot w/ averaging
         write(str(39:75),51)caxgn,caygn,cbxgn,cbygn,ixavg(ipl)
      endif
 300  continue
   50 format( 'XY= ',f5.1,',',f5.1,' to ',f5.1,',',f5.1)
   51 format( 'XY= ',f5.1,',',f5.1,' to ',f5.1,',',f5.1,',av=',i3 )
      if (ismth(ipl).gt.0.and.ismcp(ipl).gt.0) then
         write(iup,*)'In pltitle, RIP noticed you used values of both'
         write(iup,*)'smth and smcp. RIP doesn''t like it.  ipl=',ipl
         write(iup,*)'Choose one or the other.'
         stop
      elseif (ismth(ipl).gt.0) then
         ism=mod(ismth(ipl),100)
      elseif (ismcp(ipl).gt.0) then
         ism=mod(ismcp(ipl),100)
      else
         ism=0
      endif
      if ((cptyp(ipl)(2:2).eq.'c'.or.
     &     cptyp(ipl)(2:2).eq.'v'.or.
     &     cptyp(ipl)(2:2).eq.'s'    ).and.
     &     ism.gt.0 .and. .not.lnsmm(ipl)) then
         write(str(78:82),'(a3,i2)') 'sm=',ism
      endif
 60   continue
      write(iup,*)str(1:82)
      if (lnvlb(ipl)) write(str(39:60),'(a21)') '                     '
      if (.not.lnttl(ipl)) then
         ploc=toptextclg-.5*chht
         if (noplots.eq.0) call plchhq(.01,ploc,str,.011,0.,-1.)
         toptextclg=toptextclg-(chht+chgap)
      endif
c
c   Add some additional information in a second line
c
      titlstr = str
      str=' '
      iendstr=1
      if (rstmv(1,ipl).ne.0..or.rstmv(2,ipl).ne.0.) then
         speed=sqrt(rstmv(1,ipl)*rstmv(1,ipl)+
     &              rstmv(2,ipl)*rstmv(2,ipl))
         write(str(iendstr+1:),'(a15,f6.1,a6)')
     &      '(relative to V=',speed,' m/s) '
         iendstr=iendstr+27
      endif
      if (lgrad(ipl)) then
         write(str(iendstr+1:),'(a)') '(gradient) '
         iendstr=iendstr+11
      endif
      if (llapl(ipl)) then
         write(str(iendstr+1:),'(a)') '(laplacian) '
         iendstr=iendstr+12
      endif
      if (lhadv(ipl)) then
         write(str(iendstr+1:),'(a)') '(hor. advec.) '
         iendstr=iendstr+14
      endif
      if ((cdiff(ipl)(1:5).ne.'none '.or.rdiff(ipl).ne.rmsg).and.
     &    ldiffsuccess) then
         if (rdiff(ipl).eq.rmsg) then
            xtime_df=xtime
         else
            if (ldfrl(ipl)) then
               xtime_df=xtime+rdiff(ipl)
            else
               xtime_df=rdiff(ipl)
            endif
         endif
         if (cdiff(ipl)(1:5).ne.'none ') then
            iendcdiff=index(cdiff(ipl),' ')-1
            ibegcdiff=1
            do i=256,2,-1
               if (cdiff(ipl)(i-1:i-1).eq.'/') then
                  ibegcdiff=i
                  goto 75
               endif
            enddo
 75         continue
            ilencdiff=min(iendcdiff-ibegcdiff+1,15)
            if (iovly(ipl) .eq. 0) then
              write(str(iendstr+1:),'(a17,a,a7,f6.2,a2)')
     &         '(diff. from case=',
     &         cdiff(ipl)(ibegcdiff:ibegcdiff+ilencdiff-1),
     &         ', time=',xtime_df,') '
            else
              write(str(iendstr+1:),'(a17,a,a7,f6.2,a2)')
     &         '(      from case=',
     &         cdiff(ipl)(ibegcdiff:ibegcdiff+ilencdiff-1),
     &         ', time=',xtime_df,') '
            endif
            iendstr=iendstr+32+iendcdiff
         else
            if (iovly(ipl) .eq. 0) then
              write(str(iendstr+1:),'(a17,f6.2,a2)')
     &         '(diff. from time=',xtime_df,') '
            else
              write(str(iendstr+1:),'(a17,f6.2,a2)')
     &         '(      from time=',xtime_df,') '
            endif
            iendstr=iendstr+25
         endif
      endif
      if (iendstr.gt.1) then
         write(iup,*)str(1:79)
         if (.not.lnttl(ipl)) then
            ploc=toptextclg-.5*chht
            if (noplots.eq.0) call plchhq(.01,ploc,str,.011,0.,-1.)
            toptextclg=toptextclg-(chht+chgap)
         endif
      endif
c
c     call flush(iup)
      if (noplots.eq.0) then
        call gsplci(1)
        call gstxci(1)
      endif
      return
      end
