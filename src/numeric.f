c                                                             c
c*************************************************************c
c                                                             c
      logical function numeric(string)
      character string*(*), ch*1
      lenstr=len(string)
      il=0
 10   il=il+1
      if (il.gt.lenstr) goto 200
      do 100 i=0,9
         write(ch,'(i1)')i
         if (string(il:il).eq.ch) goto 10
 100  continue
      if (string(il:il).eq.'.') goto 10
      if (string(il:il).eq.'-') goto 10
      numeric=.false.
      return
 200  numeric=.true.
      return
      end
