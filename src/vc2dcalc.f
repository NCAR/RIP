c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vc2dcalc(vc3d,prs,cvcor,rcrag,
     &   rcrbg,nscrs,sfp,vcground,mabpl,miy,mjx,mkzh)
c
      dimension vc3d(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   rcrag(2),rcrbg(2),vcground(mabpl),sfp(miy,mjx)
      character cvcor*1
c
      include 'comconst'
      include 'comvctran'
c
c   Interpolate vc3d to x-section.
c
      caxgn=1.+(rcrag(2)-xjcorn)*refrat
      caygn=1.+(rcrag(1)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1)-yicorn)*refrat
      xj1t=caxgn
      xj2t=cbxgn
      yi1t=caygn
      yi2t=cbygn
      if (xj1t.le.1.5.or.xj1t.ge.mjx-.5.or.
     &    xj2t.le.1.5.or.xj2t.ge.mjx-.5.or.
     &    yi1t.le.1.5.or.yi1t.ge.miy-.5.or.
     &    yi2t.le.1.5.or.yi2t.ge.miy-.5) then
         write(iup,*)'In vc2dclac: Cross sec. endpoints must be greater'
         write(iup,'(2a,f6.1,a,f6.1,a)')' than 1.5 and less than',
     &  ' (miy-.5) (',miy-.5,') or (mjx-.5) (',mjx-.5,').'
         stop
      endif
c
      do ls=1,nscrs
         posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t)-.5
         posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t)-.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
c
c      Interpolate vc3d and pressure to x-sec
c
         do k=1,mkzh
            if (vc3d(it,jl,k).eq.rmsg.or.
     &          vc3d(it,jr,k).eq.rmsg.or.
     &          vc3d(ib,jl,k).eq.rmsg.or.
     &          vc3d(ib,jr,k).eq.rmsg) then
               write(iup,*)'vc3d has an rmsg value.'
               stop
            else
               wk1=vc3d(it,jl,k)
               wk2=vc3d(it,jr,k)
               wk3=vc3d(ib,jl,k)
               wk4=vc3d(ib,jr,k)
            endif
            vc2d(ls,k)=
     &         (1.-ratlr)*(   ratbt)*wk1+
     +         (   ratlr)*(   ratbt)*wk2+
     +         (1.-ratlr)*(1.-ratbt)*wk3+
     +         (   ratlr)*(1.-ratbt)*wk4 
            if (prs(it,jl,k).eq.rmsg.or.
     &          prs(it,jr,k).eq.rmsg.or.
     &          prs(ib,jl,k).eq.rmsg.or.
     &          prs(ib,jr,k).eq.rmsg) then
               write(iup,*)'prs has an rmsg value.'
               stop
            else
               wk1=prs(it,jl,k)
               wk2=prs(it,jr,k)
               wk3=prs(ib,jl,k)
               wk4=prs(ib,jr,k)
            endif
            prs2d(ls,k)=
     &         (1.-ratlr)*(   ratbt)*wk1+
     +         (   ratlr)*(   ratbt)*wk2+
     +         (1.-ratlr)*(1.-ratbt)*wk3+
     +         (   ratlr)*(1.-ratbt)*wk4 
         enddo
c
c      Interpolate surface pressure to x-sec, put in vcground
c
         if (sfp(it,jl).eq.rmsg.or.
     &       sfp(it,jr).eq.rmsg.or.
     &       sfp(ib,jl).eq.rmsg.or.
     &       sfp(ib,jr).eq.rmsg) then
            write(iup,*)'sfp has an rmsg value.'
            stop
         else
            sfp1=sfp(it,jl)
            sfp2=sfp(it,jr)
            sfp3=sfp(ib,jl)
            sfp4=sfp(ib,jr)
         endif
         vcground(ls)=
     &      (1.-ratlr)*(   ratbt)*sfp1+
     +      (   ratlr)*(   ratbt)*sfp2+
     +      (1.-ratlr)*(1.-ratbt)*sfp3+
     +      (   ratlr)*(1.-ratbt)*sfp4 
      enddo
c
c   Calculate the vertical coordinate at the ground surface,
c   by interpolation/extrapolation, linear in pressure
c
      do ls=1,nscrs
         if (vcground(ls).le.prs2d(ls,1)) then ! above model top - extrap.
            vcground(ls)=vc2d(ls,1)-
     &         (vc2d(ls,2)-vc2d(ls,1))*
     &         (prs2d(ls,1)-vcground(ls))/(prs2d(ls,2)-prs2d(ls,1))
         elseif (vcground(ls).ge.prs2d(ls,mkzh)) then ! below bottom - extrap.
            vcground(ls)=vc2d(ls,mkzh)+
     &         (vc2d(ls,mkzh)-vc2d(ls,mkzh-1))*
     &         (vcground(ls)-prs2d(ls,mkzh))/(prs2d(ls,mkzh)-
     &            prs2d(ls,mkzh-1))
         else  ! somewhere within model levels
            do k=mkzh,2,-1
               if (vcground(ls).le.prs2d(ls,k).and.
     &             vcground(ls).ge.prs2d(ls,k-1)) then
                  vcground(ls)=
     &               ((vcground(ls)-prs2d(ls,k-1))*vc2d(ls,k)+
     &                (prs2d(ls,k)-vcground(ls))*vc2d(ls,k-1))/
     &               (prs2d(ls,k)-prs2d(ls,k-1))
                  goto 107
               endif
            enddo
 107        continue
         endif
      enddo
c
c   Special things to be done for particular
c      vertical coordinates
c
      do ls=1,nscrs
         if (cvcor.eq.'z'.or.cvcor.eq.'f') then ! convert to height (km)
            vcground(ls)=-.001*sclht*alog(vcground(ls))
         elseif (cvcor.eq.'l') then ! convert to log pressure
            vcground(ls)=alog(vcground(ls))
         elseif (cvcor.eq.'x') then ! convert to exner function
            vcground(ls)=(vcground(ls))**gamma
         endif
         do k=1,mkzh
            if (cvcor.eq.'z'.or.cvcor.eq.'f') then ! convert to height (km)
               vc2d(ls,k)=-.001*sclht*alog(vc2d(ls,k))

            elseif (cvcor.eq.'l') then ! convert to log pressure
               vc2d(ls,k)=alog(vc2d(ls,k))
            elseif (cvcor.eq.'x') then ! convert to exner function
               vc2d(ls,k)=(vc2d(ls,k))**gamma
            endif
         enddo
      enddo
c
      return
      end
