c                                                                     c
c*********************************************************************c
c                                                                     c
      real function finterp(varr,arr,icd,miy,mjx,mkzh,yi,xj,vval,
     &   refrat,yicorn,xjcorn,rmsg,iup)
c
c   This function returns an interpolated value of the field contained
c   in the 3D array "arr", given an x value, a y value, a vertical
c   coordinate value, and a 3D array of that vertical coordinate
c   ("varr").  Note, this function assumes that yi and xj are with
c   respect to the coarse domain dot-point grid.  The input source
c   array ("arr") can be either dot or cross (icd must be set
c   accordingly).  The input vertical coordinate array must be on
c   cross points, and must be monotonically either increasing or
c   decreasing with k index.
c
      dimension varr(miy,mjx,mkzh), arr(miy,mjx,mkzh)
c
      dimension varrcol(1000), arrcol(1000)
c
c   Get x/y values for this domain (instead of coarse) and make sure
c   point is within domain.
c
      xjthisdom=1.+(xj-xjcorn)*refrat
      yithisdom=1.+(yi-yicorn)*refrat
      if ((icd.eq.0.and.
     &      (xjthisdom.lt.1..or.xjthisdom.gt.float(mjx).or.
     &       yithisdom.lt.1..or.yithisdom.gt.float(miy)))   .or.
     &    (icd.eq.1.and.
     &      (xjthisdom.lt.1.5.or.xjthisdom.gt.float(mjx)-.5.or.
     &       yithisdom.lt.1.5.or.yithisdom.gt.float(miy)-.5))) then
         finterp=rmsg
         return
      endif
c
c   Use bilinear interpolation in x and y to make a column array for
c   vertical coordinate.
c
      jarrg=int(xjthisdom-.5)
      jarrgp=jarrg+1
      iarrg=int(yithisdom-.5)
      iarrgp=iarrg+1
      ratx=xjthisdom-.5-jarrg
      raty=yithisdom-.5-iarrg
      do k=1,mkzh
         if (varr(iarrgp,jarrg ,k).ne.rmsg.and.
     &       varr(iarrg ,jarrg ,k).ne.rmsg.and.
     &       varr(iarrgp,jarrgp,k).ne.rmsg.and.
     &       varr(iarrg ,jarrgp,k).ne.rmsg) then
            varrcol(k)=((1.-ratx)*(   raty)*varr(iarrgp,jarrg ,k)+
     &                  (1.-ratx)*(1.-raty)*varr(iarrg ,jarrg ,k)+
     &                  (   ratx)*(   raty)*varr(iarrgp,jarrgp,k)+
     &                  (   ratx)*(1.-raty)*varr(iarrg ,jarrgp,k))
         else
            write(iup,*)'Found missing value in vertical cordinate'
            write(iup,*)'array in routine finterp.  This is highly'
            write(iup,*)'unexpected, so RIP must terminate.'
            stop
         endif
      enddo
c
c     Use bilinear interpolation in x and y to make a column array for
c     source array ("arr").
c
      jarrg=int(xjthisdom-icd*.5)
      jarrgp=jarrg+1
      iarrg=int(yithisdom-icd*.5)
      iarrgp=iarrg+1
      ratx=xjthisdom-icd*.5-jarrg
      raty=yithisdom-icd*.5-iarrg
      do k=1,mkzh
         if (arr(iarrgp,jarrg ,k).ne.rmsg.and.
     &       arr(iarrg ,jarrg ,k).ne.rmsg.and.
     &       arr(iarrgp,jarrgp,k).ne.rmsg.and.
     &       arr(iarrg ,jarrgp,k).ne.rmsg) then
            arrcol(k)=((1.-ratx)*(   raty)*arr(iarrgp,jarrg ,k)+
     &                 (1.-ratx)*(1.-raty)*arr(iarrg ,jarrg ,k)+
     &                 (   ratx)*(   raty)*arr(iarrgp,jarrgp,k)+
     &                 (   ratx)*(1.-raty)*arr(iarrg ,jarrgp,k))
         else
            arrcol(k)=rmsg
         endif
      enddo
c
c     Rearrange column arrays to make sure that varrcol is monotonically
c     increasing with k index.
c
      diff=varrcol(mkzh)-varrcol(1)
      if (diff.lt.0.) then
         do k=1,mkzh/2
            save=varrcol(k)
            varrcol(k)=varrcol(mkzh+1-k)
            varrcol(mkzh+1-k)=save
            save=arrcol(k)
            arrcol(k)=arrcol(mkzh+1-k)
            arrcol(mkzh+1-k)=save
         enddo
      endif
c
c   If point is below bottom (or above top) of domain, assign value at
c   bottom (or top).  Otherwise, perform vertical linear interpolation.
c
      if (vval.le.varrcol(1)) then
         finterp=arrcol(1)
      elseif (vval.ge.varrcol(mkzh)) then
         finterp=arrcol(mkzh)
      else
         do k=2,mkzh
            if (vval.ge.varrcol(k-1).and.vval.le.varrcol(k)) then
               if (arrcol(k-1).ne.rmsg.and.arrcol(k).ne.rmsg) then
                  finterp=(arrcol(k-1)*(varrcol(k)-vval)+
     &                     arrcol(k)*(vval-varrcol(k-1)))/
     &                    (varrcol(k)-varrcol(k-1))
               else
                  finterp=rmsg
               endif
               goto 55
            endif
         enddo
 55      continue
      endif
c
      return
      end
