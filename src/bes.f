c                                                                     c
c*********************************************************************c
c                                                                     c
      function bes(x)
      rint=0.
      do i=1,1000
         u=i*.001-.0005
         rint=rint+sqrt(1.-u*u)*cos(x*u)*.001
      enddo
      bes=2.*x*rint/(4.*atan(1.))
      return
      end
