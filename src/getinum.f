c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getinum(string,ipos,ival)
      character string*(*)
c
      include 'comconst'
c
      lenstr=len(string)
      do 10 i=ipos,lenstr
         if (string(i:i).eq.';'.or.string(i:i).eq.','.or.
     &         string(i:i).eq.' ') then
            goto 20
         endif
   10 continue
   20 ilast=i-1
      read(string(ipos:ilast),fmt=*,err=50,end=51) rval
      ival=nint(rval)
      ipos=ilast+2
      return
   50 write(iup,*)'   Can''t read list-directed format from "',
     &   string(ipos:ilast),'"'
      stop
   51 write(iup,*)'   End of internal file found when reading '
      write(iup,*) 'ipos = ',ipos,' ilast = ',ilast
      write(iup,*) 'string = ',string
      write(iup,*) 'offending character string(',ipos,') is ',
     &string(ipos:ipos),'\n which appears after string(',ilast,
     &') which is ', string(ilast:ilast)
      end
