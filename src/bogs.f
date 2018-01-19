      subroutine bogs(uuu,vvv,tmk,qvp,prs,ght,unorth,vnorth,ter,
     &   casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,maxtavl,
     &   vardesc,plchun,ihrip,rhrip,chrip,
     &   mdate,rhour,rlat,rlon,rslcg,ipl,maxpl,miy,mjx,mkzh)

      real uuu(miy,mjx,mkzh), vvv(miy,mjx,mkzh),
     & tmk(miy,mjx,mkzh), qvp(miy,mjx,mkzh), prs(miy,mjx,mkzh)
     & ,rslcg(2,maxpl), unorth(miy,mjx), vnorth(miy,mjx), ter(miy,mjx)
     & ,ght(miy,mjx,mkzh)
c
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10, casename*256
c
c arrays local to bogs
      real tdp(miy,mjx,mkzh), p(mkzh), u(mkzh), v(mkzh), td(mkzh),
     &   t(mkzh), z(mkzh)
      character *40 string1,string2,string3,string4
      data iseq_num /0/
c
      include 'comconst'
c
c  make mdate rounded to the nearest hour  JFB 7/27/05
      call jdate (mdate, rhour, ndate)
c
c   Interpolate data to sounding location
c
c     write(6,*) 'in bogs, ipl = ',ipl
c     write(6,*) 'rslcg(2,ipl) = ',rslcg(2,ipl),' rslcg(1,ipl) = ',
c    &  rslcg(1,ipl)
      sxgn=1.+(rslcg(2,ipl)-xjcorn)*refrat
      sygn=1.+(rslcg(1,ipl)-yicorn)*refrat
      if (sxgn.le..5.or.sxgn.ge.mjx-.5.or.
     &    sygn.le..5.or.sygn.ge.miy-.5) then
         write(iup,*)'I don''t do soundings outside the'
         write(iup,*)'cross-point domain.'
         write(iup,*) 'sxgn = ',sxgn,' sygn = ',sygn
         write(iup,*) 'mjx = ',mjx,' miy = ',miy
         stop
      endif
c
c   Make p
c
      posx=sxgn-.5
      posy=sygn-.5
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      do 2 k=1,mkzh
        p(k)= (    (1.-ratlr)*(   ratbt)*prs(it,jl,k)+
     &             (   ratlr)*(   ratbt)*prs(it,jr,k)+
     &             (1.-ratlr)*(1.-ratbt)*prs(ib,jl,k)+
     &             (   ratlr)*(1.-ratbt)*prs(ib,jr,k) )
     &   * 100.
    2 continue
      tern = (   (1.-ratlr)*(   ratbt)*ter(it,jl)+
     &           (   ratlr)*(   ratbt)*ter(it,jr)+
     &           (1.-ratlr)*(1.-ratbt)*ter(ib,jl)+
     &           (   ratlr)*(1.-ratbt)*ter(ib,jr) )
      if(abs(tern) .lt. .01 ) tern = 0.   ! change sea level point to elev of 0.
c
c   Make velocity components
c
      posx=sxgn
      posy=sygn
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
      do 3 k=1,mkzh
         uuus = (
     &                (1.-ratlr)*(   ratbt)*uuu(it,jl,k)+
     &                (   ratlr)*(   ratbt)*uuu(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*uuu(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*uuu(ib,jr,k) )
         vvvs = (
     &                (1.-ratlr)*(   ratbt)*vvv(it,jl,k)+
     &                (   ratlr)*(   ratbt)*vvv(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*vvv(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*vvv(ib,jr,k) )
         v(k)= (unorths*uuus+vnorths*vvvs)
         u(k)= (vnorths*uuus-unorths*vvvs)
    3 continue
c
c   Make temperature array
c
      posx=sxgn-.5
      posy=sygn-.5
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      do 4 k=1,mkzh
         if (tmk(it,jl,k).eq.rmsg.or.tmk(it,jr,k).eq.rmsg.or.
     &       tmk(ib,jl,k).eq.rmsg.or.tmk(ib,jr,k).eq.rmsg) then
            t(k)=rmsg
         else
            t(k)= (
     &                (1.-ratlr)*(   ratbt)*tmk(it,jl,k)+
     &                (   ratlr)*(   ratbt)*tmk(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*tmk(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*tmk(ib,jr,k) )
         endif
    4 continue

c  get dewpoint
c  calling tdpcalc is costly. So, save it to a file.  JFB 8/2/05
      call getvarinfo(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,'tdp       ',maxtavl,ndim,icd,vardesc,plchun,
     &   0,istat,iup)
      if (istat.eq.1) then
         call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &      ncxc,'tdp       ',miy,mjx,mkzh,maxtavl,3,
     &      0,tdp,istat)
      else
        call tdpcalc(qvp,prs,tdp,miy,mjx,mkzh)
        call writefile (tdp,'tdp       ',0,
     &      3,1,'Dewpoint              ','C       ',ihrip,rhrip,
     &      chrip,casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &      maxtavl,miy,mjx,mkzh)
      endif
c
c   Make dewpoint array
c
      posx=sxgn-.5
      posy=sygn-.5
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      do 5 k=1,mkzh
         if (tdp(it,jl,k).eq.rmsg.or.tdp(it,jr,k).eq.rmsg.or.
     &       tdp(ib,jl,k).eq.rmsg.or.tdp(ib,jr,k).eq.rmsg) then
            td(k)=rmsg
         else
            td(k)= (
     &                (1.-ratlr)*(   ratbt)*tdp(it,jl,k)+
     &                (   ratlr)*(   ratbt)*tdp(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*tdp(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*tdp(ib,jr,k) )
     &          + 273.16   ! convert to kelvin
         endif
    5 continue

c   Make geopotential height array
c
      posx=sxgn-.5
      posy=sygn-.5
      jl=int(posx)
      jr=jl+1
      ib=int(posy)
      it=ib+1
      ratlr=posx-jl
      ratbt=posy-ib
      do 6 k=1,mkzh
         if (ght(it,jl,k).eq.rmsg.or.ght(it,jr,k).eq.rmsg.or.
     &       ght(ib,jl,k).eq.rmsg.or.ght(ib,jr,k).eq.rmsg) then
            z(k)=rmsg
         else
            z(k)= (
     &                (1.-ratlr)*(   ratbt)*ght(it,jl,k)+
     &                (   ratlr)*(   ratbt)*ght(it,jr,k)+
     &                (1.-ratlr)*(1.-ratbt)*ght(ib,jl,k)+
     &                (   ratlr)*(1.-ratbt)*ght(ib,jr,k) )
         endif
    6 continue

      do k = mkzh, 1, -1
	vv = sqrt( u(k)*u(k) + v(k)*v(k) ) 
	v(k) = atan2(-u(k),-v(k))*180./3.14159
	if(v(k).lt.0.) v(k) = v(k) + 360.
	u(k) = vv
c       z(k) = -888888.
c u contains speed, v contains dir
c       write(6,100) p(k), u(k), v(k), t(k), td(k)
      enddo
  100 format (3x,5f12.5)
      slp = -888888.
      string1 = '                                        '
      iseq_num = iseq_num + 1
      write(string1,'(i6.6)') iseq_num
      string2 = 'RIP BOGUS SOUNDING                      '
      string3 = 'FM-35 TEMP (BOGUS)                      '
      string4 = '                                        '
      iunit = 66
      call write_bog ( p, z, t, td, u, v, 
     * slp, tern, rlat, rlon, ndate, mkzh, 
     * string1, string2, string3, string4, .true., iseq_num,
     * iunit )
      return
      end
C--------------------------------------------------------------------------
      SUBROUTINE write_bog ( p, z, t, td, spd, dir, 
     * slp, ter, xlat, xlon, mdate, mkzh, 
     * string1, string2, string3, string4, bogus, iseq_num,
     * iunit )
c  Write out a sounding in little-r format. p in Pa, t and td in K,
c  speed in m/s.
c
      real p(mkzh), z(mkzh),t(mkzh),td(mkzh),spd(mkzh),dir(mkzh)

      CHARACTER*20 date_char
      CHARACTER*40 string1, string2 , string3 , string4
      CHARACTER*84 rpt_format 
      CHARACTER*22 meas_format 
      CHARACTER*14 end_format

      LOGICAL bogus

      rpt_format =  ' ( 2f20.5 , 2a40 , ' 
     *             // ' 2a40 , 1f20.5 , 5i10 , 3L10 , ' 
     *             // ' 2i10 , a20 ,  13( f13.5 , i7 ) ) '
      meas_format =  ' ( 10( f13.5 , i7 ) ) '
      end_format = ' ( 3 ( i7 ) ) ' 

      WRITE (date_char(9:16),fmt='(i8.8)') mdate
      IF (mdate/1000000 .GT. 70 ) THEN
         date_char(7:8)='19'
      ELSE
         date_char(7:8)='20'
      ENDIF
      date_char(17:20)='0000'
      date_char(1:6)='      '

      WRITE ( UNIT = iunit , ERR = 19 , FMT = rpt_format ) 
     *        xlat,xlon, string1 , string2 , 
     *        string3 , string4 , ter, mkzh*6, 0,0,iseq_num,0, 
     *        .true.,bogus,.false., 
     *         -888888, -888888, date_char , 
     *         slp,0,-888888.,0, -888888.,0, -888888.,0, -888888.,0, 
     *               -888888.,0, 
     *               -888888.,0, -888888.,0, -888888.,0, -888888.,0, 
     *               -888888.,0, 
     *               -888888.,0, -888888.,0
   
c      z(mkzh) = ter    ! set surface height to lowest model level height
      DO 100 k = mkzh, 1, -1
         WRITE ( UNIT = iunit , ERR = 19 , FMT = meas_format ) 
     *          p(k), 0, z(k),0, t(k),0, td(k),0, 
     *          spd(k),0, dir(k),0, 
     *          -888888.,0, -888888.,0,-888888.,0, -888888.,0
100   CONTINUE
      WRITE ( UNIT = iunit , ERR = 19 , FMT = meas_format ) 
     *         -777777.,0, -777777.,0,float(mkzh),0,
     *         -888888.,0, -888888.,0, -888888.,0, 
     *         -888888.,0, -888888.,0, -888888.,0, 
     *         -888888.,0
      WRITE ( UNIT = iunit , ERR = 19 , FMT = end_format )  mkzh, 0, 0

      RETURN
19    CONTINUE
      PRINT *,'troubles writing a sounding'
      STOP 'write_obs'
      END
c-----------------------------------------------------------
      subroutine jdate (mdate, rhour, ndate)
      character t8*8
      integer nd(12)
      data nd /31,28,31,30,31,30,31,31,30,31,30,31/
      write(t8,'(i8)') mdate
      read(t8,'(4i2)') iy, im, id, ih
      if (mod(iy,4) .eq. 0 ) nd(2) = 29
      if (rhour .gt. 0.5) then
        ih = ih + 1
        if ( ih .ge. 24 ) then
          ih = 0
          id = id + 1
          if ( id .gt. nd(im)) then
            id = 1
            im = im + 1
            if ( im .gt. 12 ) then
              im = 1
              iy = iy + 1
            endif
          endif
        endif
      endif
      write(t8,'(4i2.2)') iy, im, id, ih
      read(t8,'(i8)') ndate
      return
      end
