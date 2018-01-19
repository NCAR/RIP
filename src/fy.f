c                                                                     c
c*********************************************************************c
c                                                                     c
      function fy(xinp,yinp)
c
      include 'comconst'
      include 'comvctran'
c
      if (ivcs.eq.1) then
         l=min(nscross-1,max(1,int(xinp)))
         k=min(mkzhcross-1,max(1,int(yinp)))
         lpl=l+1
         kpl=k+1
         reml=xinp-float(l)
         remk=yinp-float(k)
         fy=(   reml)*(   remk)*vc2d(lpl,kpl)+
     &      (1.-reml)*(   remk)*vc2d(l,kpl)+
     &      (   reml)*(1.-remk)*vc2d(lpl,k)+
     &      (1.-reml)*(1.-remk)*vc2d(l,k)
      else
         fy=yinp
      endif
      return
      end
