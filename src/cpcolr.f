c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cpcolr (xcra,ycra,ncra,iaia,igia,naia)
c
c   This is a user-supplied areas subroutine which allows the user
c   to color areas in a particular way.
c
      dimension xcra(*),ycra(*),iaia(*),igia(*)
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
c
c   Assume polygon will be filled until we find otherwise.
c
      ifll=1
c
c   If any area identifier is negative, don't fill the polygon
c
      do 101 i=1,naia
         if (iaia(i).lt.0.) ifll=0
 101  continue
c
c   Otherwise, fill the polygon in the color implied by its area
c   identifier relative to edge group 3 (the contour-line group).
c
      if (ifll.ne.0) then
         ifll=0
         do 102 i=1,naia
            if (igia(i).eq.3) ifll=iaia(i)
 102     continue
c
c      Note that if icoindcp(ifll) is negative, that means this
c      polygon should remain transparent (i.e. not filled).
c
         if ( ifll .ge. 1 .and. ifll .le. nconarea ) then
           if ( icoindcp(ifll) .ge. 0 ) then
             call gsfaci(icoindcp(ifll))
             call gfa (ncra-1,xcra,ycra)
           endif
         endif
      endif
c
      return
      end
