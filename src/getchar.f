c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getchar(string,ipos,chval,ignorecomma)
      character string*(*),chval*(*)
c
      include 'comconst'
c
      lenstr=len(string)
      lenchv=len(chval)
      do 10 i=ipos,lenstr
         if (string(i:i).eq.';'.or.
     &       (ignorecomma.eq.0.and.string(i:i).eq.',').or.
     &       string(i:i).eq.' ') then
            goto 20
         endif
   10 continue
   20 ilast=i-1
      ilf=ilast-ipos+1
      if (ilf.gt.lenchv) then
         write(iup,*)'In getchar, chval is not long enough.'
         write(iup,*)'lenchv,ilf,string(ipos:ilast)='
         write(iup,*)lenchv,ilf,' #',string(ipos:ilast),'#'
         stop
      endif
      chval(1:ilf)=string(ipos:ilast)
      if (ilf.lt.lenchv) then
         do i=ilf+1,lenchv
            chval(i:i)=' '
         enddo
      endif
      ipos=ilast+2
      return
      end
