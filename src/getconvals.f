c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getconvals(cbeg,cend,cint,ncon,mult,cval,ncvl,mjsk,
     &   pslab,mabpl,njx,niy,rmsg,maxcon,valcon,majcon,
     &   cintuse,numcon,iup)
c
c
c   cbeg is the beginning contour value.  rmsg means pick a
c   nice value close to the max value in the field. -rmsg means
c   pick a nice value close to the min value in the field.
c
c   cend is the ending contour value.  rmsg means pick a nice
c   value close to the max value in the field. -rmsg means pick
c   a nice value close to the min value in the field.
c
c   cint is the contour interval. rmsg means pick a nice value
c   that generates close to ncon contours.
c
c   mult means use a multiplicative contour interval instead of
c   an additive one.  If mult is true, cbeg and cend must be of
c   the same sign (otherwise cend will be changed to a value
c   that is of the same sign as cbeg), and all contour values
c   will also be the same sign.
c
c   cval is an array holding user-specified contour values.  ncvl is
c   the number of intended values in cval.  If ncvl>0, it is assumed
c   that the user wants to use their own specified values for the
c   contours, and this overrides any information supplied in cint,
c   cbeg, cend, and mult.
c
c   mjsk is the number of minor (unlabeled) contours desired
c   between major (labeled) contours.
c
      dimension pslab(mabpl,niy),valcon(maxcon),majcon(maxcon)
      dimension cval(maxcon)
      logical mult, choosecb, chooseci
c
c   First deal with the case where the user has specified the
c   desired contour values with cval.
c
      if (ncvl.gt.0) then
         imaj=-999
         do i=1,ncvl
            valcon(i)=cval(i)
            if (valcon(i).eq.0.0) then
               imaj=i
            endif
         enddo
         numcon=ncvl
         if (imaj.eq.-999) then
            imaj=min((mjsk+3)/2,numcon)
         endif
         imajrel=mod(imaj,mjsk+1)-(mjsk+1)
         do i=1,numcon
            if (valcon(i).eq.0.0) then    ! zero contour
               majcon(i)=0
            elseif (mod(i-imajrel,mjsk+1).eq.0) then
               if (valcon(i).gt.0) then ! positive major (labeled) contour
                  majcon(i)=2
               else                     ! negative major (labeled) contour
                  majcon(i)=-2
               endif
            else
               if (valcon(i).gt.0) then   ! positive unlabeled contour
                  majcon(i)=1
               else                       ! negative unlabeled contour
                  majcon(i)=-1
               endif
            endif
         enddo
         cintuse=0.0
         return
      endif
c
c   First get max and min in field.
c
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).ne.rmsg) then
            fv=pslab(j,i)
            goto 20
         endif
      enddo
      enddo
      write(iup,*)'   In getconvals: not generating any contours'
      write(iup,*)'    because all values are special value.'
      numcon=0
      return
 20   continue
c
      rmax=fv
      rmin=fv
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).ne.rmsg) then
            rmax=max(pslab(j,i),rmax)
            rmin=min(pslab(j,i),rmin)
         endif
      enddo
      enddo
      if (rmax.eq.rmin) then
         write(iup,*)'   In getconvals: not generating any contours'
         write(iup,*)'    because field is constant.  Value = ',rmax
         numcon=0
         return
      endif
c
      choosecb=.false.
      if (cbeg.eq.rmsg) then
         cb=rmax
         choosecb=.true.
      elseif (cbeg.eq.-rmsg) then
         cb=rmin
         choosecb=.true.
      else
         cb=cbeg
      endif
c
      if (cend.eq.rmsg) then
         ce=rmax
      elseif (cend.eq.-rmsg) then
         ce=rmin
      else
         ce=cend
      endif
c
      chooseci=.false.
      if (cint.eq.rmsg.or.cint.eq.0.0) then
         chooseci=.true.
      elseif (cint.lt.0.0) then
         ci=-cint
      else
         ci=cint
      endif
c
      if (mult) then
         if (cb.eq.0.0.or.ce.eq.0.0) then
            write(iup,*)'In getconvals: not generating any contours'
            write(iup,*)'   because mult is true and cbeg or cend'
            write(iup,*)'   are zero. cbeg,cend=',cbeg,cend
            numcon=0
            return
         endif
         if (cb*ce.lt.0.0) then
            write(iup,*)'In getconvals: not generating any contours'
            write(iup,*)'   because mult is true and cbeg and cend are'
            write(iup,*)'   of opposite sign.  cbeg,cend=',cbeg,cend
            numcon=0
            return
         endif
         if (ci.eq.1.0) chooseci=.true.
         if (cb.lt.0.0) then
            iswitch=1
            cb=-cb
            ce=-ce
            temp=rmax
            rmax=-rmin
            rmin=-temp
         else
            iswitch=0
         endif
         valm=1e-8*cb
         rmax=max(rmax,valm)
         rmin=max(rmin,valm)
         ce=max(ce,valm)
      endif
c
c   Set up start, end, and interval.
c
      if (.not.mult) then
         if (chooseci) then
            crange=max(abs(ce-cb),1.e-10)
            ci = crange/ncon
            p = 10.**(int(alog10(ci)+50000.)-50000)
            ci=max(ci,p)
            ii=int(ci/p)
            if (ii.ge.7) then
               ci = 10.*p
            else
               ci = ii*p
            endif
         endif
         if (choosecb) then
            cb=nint(cb/ci)*ci
            if (cb.ge.rmax) cb=cb-ci
            if (cb.le.rmin) cb=cb+ci
         endif
      else
         if (choosecb) then
            cmax=10.**(int(alog10(rmax)+50000.)-50000)
            cmin=10.**(int(alog10(rmin)+50000.)-50000+1)
            cb=10.**(nint(alog10(cb)+50000.)-50000)
            cb=min(max(cb,cmin),cmax)
         endif
         if (chooseci) then
            ci=(max(ce,cb)/min(ce,cb))**(1./(max(ncon,2)-1))
            if (ci.le.1.6) then
               ci=sqrt(2.)
            elseif (ci.le.3.5) then
               ci=2.
            elseif (ci.le.7.5) then
               ci=5.
            else
               ci=10.**nint(alog10(ci))
            endif
         endif
      endif
      cintuse=ci
c
c   Generate contour levels.
c
      if (.not.mult) then
         if (ce.lt.cb) then
            ci=-abs(ci)
         else
            ci=abs(ci)
         endif
      else
         alci=alog(ci)
         if (abs(ce).lt.abs(cb)) then
            alci=-abs(alci)
         else
            alci=abs(alci)
         endif
         ci=exp(alci)
      endif
      numcon=1
      valcon(numcon)=cb
 50   numcon=numcon+1
      if (numcon.gt.maxcon) goto 100
      if (.not.mult) then
         valcon(numcon)=valcon(numcon-1)+ci
      else
         valcon(numcon)=valcon(numcon-1)*ci
      endif
      if ((ce.ge.cb.and.valcon(numcon).gt.ce).or.
     &    (ce.lt.cb.and.valcon(numcon).lt.ce)) goto 100
      goto 50
 100  numcon=numcon-1
c
c   Return mult contours to opposite sign, if iswitch=1
c
      if (mult) then
         if (iswitch.eq.1) then
            do i=1,numcon
               valcon(i)=-valcon(i)
            enddo
         endif
      endif
c
c   Determine major contours
c
      do i=1,numcon
         if (.not.mult) then
            fac=valcon(i)/((mjsk+1)*abs(ci))
         else
            fac=alog10(abs(valcon(i)))
         endif
         diff=abs(fac-float(nint(fac)))
         if (diff.lt..001) then
            imaj=i
            goto 200
         endif
      enddo
      imaj=1
 200  continue
      imajrel=mod(imaj,mjsk+1)-(mjsk+1)
      do i=1,numcon
         if (.not.mult) then
            fac=valcon(i)/((mjsk+1)*abs(ci))
         else
            fac=1.
         endif
         if (abs(fac).lt.1e-4) then    ! zero contour
            valcon(i)=0.0
            majcon(i)=0
         elseif (mod(i-imajrel,mjsk+1).eq.0) then
            if (valcon(i).gt.0) then ! positive major (labeled) contour
               majcon(i)=2
            else                     ! negative major (labeled) contour
               majcon(i)=-2
            endif
         else
            if (valcon(i).gt.0) then   ! positive unlabeled contour
               majcon(i)=1
            else                       ! negative unlabeled contour
               majcon(i)=-1
            endif
         endif
      enddo
c
c   Reorder if valcon(1) > valcon(numcon)
c
      if (valcon(1).gt.valcon(numcon)) then
         do i=1,numcon/2
            ii=numcon+1-i
            tmp=valcon(i)
            valcon(i)=valcon(ii)
            valcon(ii)=tmp
            itmp=majcon(i)
            majcon(i)=majcon(ii)
            majcon(ii)=itmp
         enddo
      endif
c
      return
      end
