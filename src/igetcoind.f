c                                                                     c
c*********************************************************************c
c                                                                     c
      function igetcoind(cospec,conam,nco,string)
      character cospec*(*), conam(0:nco)*(*), string*(*)
c
      include 'comconst'
c
      if (cospec(1:11).eq.'transparent') then
         igetcoind=-1
         return
      endif
      do 50 i=0,nco
         if (cospec.eq.conam(i)) goto 100
   50 continue
      write(iup,*)'   Couldn''t find color named [',cospec,']'
      write(iup,*)'   There are ',nco,' colors defined. They are:'
      do 80 i=0,nco
         write(iup,*)'     ',i,' is [',conam(i),']'
   80 continue
         write(iup,*)' Offending string is ',string
      stop
  100 igetcoind=i
      return
      end
