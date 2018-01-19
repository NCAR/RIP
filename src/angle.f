c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real function angle(xor,yor,xpt,ypt,ftouhor,ftouver)
c
c   This function returns the angle between an x-axis with origin at
c   (xor,yor) and the line joining (xor,yor) and (xpt,ypt).
c
      xd=(xpt-xor)/ftouhor
      yd=(ypt-yor)/ftouver
      divinsure=.000001/ftouhor
      hypot = ( xd**2 + yd**2 ) ** .5 + divinsure
      if (yd.ge.0) then
         angle = acos(xd/hypot)
      else
         angle = 2.*3.14159 - acos(xd/hypot)
      endif
      return
      end
