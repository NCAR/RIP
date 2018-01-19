c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cpchhl (iflag)
c
c   This is a user-supplied conpack subroutine which allows the user
c   to alter the characteristics of the high and low labels.  This
c   particular version is used to alter the fill color of label boxes.
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
c
      if (iflag .eq. 2) then   ! a box for a high is about to be filled
         call gsfaci (icpfchl)   ! Set up high label box fill color
      elseif (iflag .eq. 6) then   ! a low box is about to be filled
        call gsfaci (icpfclo)   ! Set up low label box fill color index
      endif
c
      return
      end
