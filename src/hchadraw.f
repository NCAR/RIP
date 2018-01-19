c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hchadraw(ixwin,iywin,rlevl,work,icdwk,ilev,
     &         iintv,icosq,rcosq,lchfl,icolr,incsq,
     &         pslab1,pslabt,ipslab,
     &         maxcosq,mabpl,morpl,maxlev,maxpl,miy,mjx,mkzh,ipl,irota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),
     &   icosq(maxcosq,maxpl),rlevl(maxlev,maxpl),iintv(maxpl),
     &   work(miy,mjx,mkzh),icdwk(maxpl),
     &   rcosq(maxcosq,maxpl),incsq(maxpl),icolr(maxpl),
     &   pslab1(mabpl,morpl),pslabt(mabpl,morpl),ipslab(mabpl,morpl),
     &   irota(maxpl)
      logical lchfl(maxpl)
c
      character pchar*1, pchar2*2
c
      include 'comconst'
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
      xintervs=ixwin(2,ipl)-ixwin(1,ipl)
      yintervs=iywin(2,ipl)-iywin(1,ipl)
      if (irota(ipl).ne.90.and.irota(ipl).ne.-90) then
         aspect=yintervs/xintervs
      else
         aspect=xintervs/yintervs
      endif
      faspect=(ftmax-fbmin)/(frmax-flmin)
      if (aspect.lt.faspect) then
         fl=flmin
         fr=frmax
         fextra=.5*((ftmax-fbmin)-aspect*(frmax-flmin))
         fb=fbmin+fextra
         ft=ftmax-fextra
      else
         fb=fbmin
         ft=ftmax
         fextra=.5*((frmax-flmin)-1./aspect*(ftmax-fbmin))
         fl=flmin+fextra
         fr=frmax-fextra
      endif
      if (icdwk(ipl).eq.0) then
         ul=1.
         ur=1.+xintervs
         ub=1.
         ut=1.+yintervs
         niy=nint(ut)
         njx=nint(ur)
      else
         ul=.5
         ur=.5+xintervs
         ub=.5
         ut=.5+yintervs
         niy=int(ut)
         njx=int(ur)
      endif
c
c   Transpose ul,ur,ub,ut if view is rotated +/- 90 degrees.
c   niy, njx, and pslab arrays will not be adjusted for rotation
c   until a later point.
c
      if (irota(ipl).eq.90.or.irota(ipl).eq.-90) then
         usv=ul
         ul=ub
         ub=usv
         usv=ur
         ur=ut
         ut=usv
      endif
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
c
c   Put appropriate data into horizontal slab.
c
      do 90 j=1,njx
         jj=j+ixwin(1,ipl)-1
         do 90 i=1,niy
            ii=i+iywin(1,ipl)-1
            pslab1(j,i)=work(ii,jj,nint(rlevl(ilev,ipl)))
   90 continue
c
c   Rotate pslab1 array, and switch njx,niy, if needed.
c
      call rotpslab(pslab1,pslabt,mabpl,morpl,njx,niy,irota(ipl))
c
c   Plot the characters
c
      if (.not.lchfl(ipl)) then
c
      setfrac=(fr-fl)/(ur-ul)*50.
      if (iintv(ipl).gt.0) then
         chsiz=max(.004,.012*setfrac*iintv(ipl))
      else
         chsiz=max(.004,.012*setfrac*iintv(ipl)*.5)
      endif
      chsiz2=max(.004,.7*chsiz)
c
      intp=iintv(ipl)
      if (intp.gt.0) then  ! regular staggering
         intph=0
      else   ! "diamond pattern" staggering
         intp=-intp
         intph=intp/2
      endif
      do j=1,njx
      do i=1,niy
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       pslab1(j,i).ne.rmsg) then
            ips=nint(min(pslab1(j,i),100.4))
            if (ips.gt.99.or.ips.lt.0) goto 100
            rj=j
            ri=i
            do 95 ir = 1,incsq(ipl)
               if (nint(rcosq(ir,ipl)).eq.ips) then
                  call gstxci(icosq(ir,ipl))
                  call gsplci(icosq(ir,ipl))
                  goto 96
               endif
   95       continue
            call gstxci(icolr(ipl))
            call gsplci(icolr(ipl))
   96       continue        
            if (ips.lt.10) then
               write(pchar,'(i1)') ips
               call plchhq(rj,ri,pchar,chsiz,0.,0.)
            else
               write(pchar2,'(i2)') ips
               call plchhq(rj,ri,pchar2,chsiz2,0.,0.)
            endif
         endif
  100    continue
      enddo
      enddo
c
      else
c
      do j=1,njx
      do i=1,niy
         ips=nint(max(min(pslab1(j,i),101.),-1.))
         if (ips.gt.99.or.ips.lt.0) then
            ipslab(j,i)=1
         else
            do ir = 1,incsq(ipl)
               if (nint(rcosq(ir,ipl)).eq.ips) then
                  ipslab(j,i)=icosq(ir,ipl)
                  goto 106
               endif
            enddo
            ipslab(j,i)=1
  106       continue        
         endif
      enddo
      enddo
c
      call gca (.5,.5,float(njx)+.5,float(niy)+.5,mabpl,morpl,1,1,
     1     njx,niy,ipslab)
      call sflush
c
      endif
c
      call gstxci(1)
      call gsplci(1)
      return
      end
