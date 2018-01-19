c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine drawcl (xcs,ycs,ncs,iai,iag,nai)
c
c   This is a user-supplied conpack subroutine which allows the user
c   to draw contour lines in a particular way.  This particular version
c   is used to draw a polyline if and only if none of the area
c   identifiers for the area containing the polyline are negative.
c
      dimension xcs(*),ycs(*),iai(*),iag(*)
c
c   Use iag, just to avoid the "unused dummy" warning from
c      some compilers
c
      iuse=iag(1)
c
c   Turn on drawing.
c
      idr=1
c
c   If any area identifier is negative, turn off drawing.
c
      do 101 i=1,nai
         if (iai(i).lt.0.) idr=0
 101  continue
c
c   If drawing is turned on, draw the polyline.
c
      if (idr.eq.1) call curved (xcs,ycs,ncs)
c
      return
      end
