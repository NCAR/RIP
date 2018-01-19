c                                                                     c
c*********************************************************************c
c                                                                     c
      integer function mlen(fname)
c
c   This routine finds the location of the date within the
c   rip outout file name, based on the location of the last
c   underscore character ("_").
c
      character*(*) fname
      n = len(fname)
      do i = n,1,-1
        if (fname(i:i) .eq. '_') then
          mlen = i - 9
          return
        endif
      enddo
      stop 'error in mlen'
      end
