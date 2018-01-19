c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cpchll (iflag)
c
c   This is a user-supplied conpack subroutine which allows the user
c   to alter the characteristics of the contour line labels.  This
c   particular version is used to alter the fill color of label boxes.
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
c
      if (iflag .eq. 2) then   ! a cont. lab. box is about to be filled
         call cpgeti('PAI',i)
         if (mconcp(i).eq.2) then
            call gsfaci (icpfclb) ! Set up fill color for pos. contour
         elseif (mconcp(i).eq.-2) then
            call gsfaci (icpfcnl) ! Set up fill color for neg. contour
         elseif (mconcp(i).eq.0) then
            call gsfaci (icpfczr) ! Set up fill color for zero contour
         endif
      endif
c
      return
      end
