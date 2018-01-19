c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cpmpxy(imap,xinp,yinp,xotp,yotp)
c
c   This is a user-supplied subroutine that transforms each grid
c   point's vertical level index to the chosen vertical coordinate.
c
      include 'comconst'
      include 'comvctran'
c
      if (imap.eq.0) then
         if (int(xinp).eq.1) then
           yinp=3.
         else
           yinp=0.
         endif
      elseif (imap.eq.1) then
         l=min(nscross-1,max(1,int(xinp)))
         k=min(mkzhcross-1,max(1,int(yinp)))
         lpl=l+1
         kpl=k+1
         reml=xinp-float(l)
         remk=yinp-float(k)
         yotp=(   reml)*(   remk)*vc2d(lpl,kpl)+
     &      (1.-reml)*(   remk)*vc2d(l,kpl)+
     &      (   reml)*(1.-remk)*vc2d(lpl,k)+
     &      (1.-reml)*(1.-remk)*vc2d(l,k)
         xotp=xinp
      elseif (imap.eq.-1) then
         l=min(nscross-1,max(1,int(xinp)))
         lpl=l+1
         reml=xinp-float(l)
         do k=1,mkzhcross
            tmpcol(k)=(   reml)*vc2d(lpl,k)+
     &         (1.-reml)*vc2d(l,k)
         enddo
         if ((tmpcol(mkzhcross).gt.tmpcol(1).and.
     &       yinp.ge.tmpcol(mkzhcross)).or.
     &       (tmpcol(mkzhcross).lt.tmpcol(1).and.
     &       yinp.le.tmpcol(mkzhcross))) then
            k=mkzhcross-1
            kpl=mkzhcross
         elseif ((tmpcol(mkzhcross).gt.tmpcol(1).and.
     &       yinp.le.tmpcol(1)).or.
     &       (tmpcol(mkzhcross).lt.tmpcol(1).and.
     &       yinp.ge.tmpcol(1))) then
            k=1
            kpl=2
         else
            do kl=1,mkzhcross-1
               if ((yinp.ge.tmpcol(kl).and.yinp.le.tmpcol(kl+1)).or.
     &             (yinp.le.tmpcol(kl).and.yinp.ge.tmpcol(kl+1))) then
                  k=kl
                  kpl=kl+1
                  goto 50
               endif
            enddo
 50         continue
         endif
         yotp=((tmpcol(kpl)-yinp)*float(k)+
     &         (yinp-tmpcol(k))*float(kpl))/
     &        (tmpcol(kpl)-tmpcol(k))
         xotp=xinp
      else
         yotp=yinp
         xotp=xinp
      endif
      return
      end
