c                                                                     c
c*********************************************************************c
c                                                                     c
      real function finterp2d(arr,icd,miy,mjx,yi,xj,
     &   refrat,yicorn,xjcorn,rmsg,iup)
c
c     This function returns an interpolated value of the field contained
c     in the 2D array "arr", given an x and y value.
c     Note, this function assumes that yi and xj are with
c     respect to the coarse domain dot-point grid.  The input source
c     array ("arr") can be either dot or cross (icd must be set
c     accordingly).
c
      dimension arr(miy,mjx)
c
c     Get x/y values for this domain (instead of coarse) and make sure
c     point is within domain.
c
      xjthisdom=1.+(xj-xjcorn)*refrat
      yithisdom=1.+(yi-yicorn)*refrat
      if ((icd.eq.0.and.
     &      (xjthisdom.lt.1..or.xjthisdom.gt.float(mjx).or.
     &       yithisdom.lt.1..or.yithisdom.gt.float(miy)))   .or.
     &    (icd.eq.1.and.
     &      (xjthisdom.lt.1.5.or.xjthisdom.gt.float(mjx)-.5.or.
     &       yithisdom.lt.1.5.or.yithisdom.gt.float(miy)-.5))) then
         finterp2d=rmsg
         return
      endif
c
c     Use bilinear interpolation in x and y to get the value.
c
      jarrg=int(xjthisdom-icd*.5)
      jarrgp=jarrg+1
      iarrg=int(yithisdom-icd*.5)
      iarrgp=iarrg+1
      ratx=xjthisdom-icd*.5-jarrg
      raty=yithisdom-icd*.5-iarrg
      if (arr(iarrgp,jarrg ).ne.rmsg.and.
     &    arr(iarrg ,jarrg ).ne.rmsg.and.
     &    arr(iarrgp,jarrgp).ne.rmsg.and.
     &    arr(iarrg ,jarrgp).ne.rmsg) then
         finterp2d=((1.-ratx)*(   raty)*arr(iarrgp,jarrg )+
     &              (1.-ratx)*(1.-raty)*arr(iarrg ,jarrg )+
     &              (   ratx)*(   raty)*arr(iarrgp,jarrgp)+
     &              (   ratx)*(1.-raty)*arr(iarrg ,jarrgp))
      else
         finterp2d=rmsg
      endif
c
      return
      end
