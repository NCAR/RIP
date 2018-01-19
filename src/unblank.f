c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine unblank(string)
      character string*(*)
c
      include 'comconst'
c
      ilastorig=lennonblank(string)
      if (ilastorig.eq.0) then
         write(iup,*)'String is all blanks.'
         stop
      endif
      ilast=ilastorig
      do 100 j=ilastorig,1,-1
         if (string(j:j).eq.' ') then
            do 30 k=j,ilast-1
               string(k:k)=string(k+1:k+1)
   30       continue
            string(ilast:ilast)=' '
            ilast=ilast-1
         endif
  100 continue
      return
      end
