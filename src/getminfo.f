c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getminfo(casename,iendc,minfo,iup)
c
      character casename*(*),fname*256
      character minfo(5)*256
c
c   Read in the minfo string from the .minfo file.
c
      igotminfo=0
      do iline=1,5
         minfo(iline)=' '
      enddo
      fname=casename(1:iendc)//'.minfo'
      open(unit=25,file=fname,err=1021,form='formatted',status='old')
      do iline=1,5
         read(25,'(a)',end=901,err=901) minfo(iline)
         iendm=lennonblank(minfo(iline))
         if (iendm.gt.0) then
            igotminfo=1
            write(iup,*) minfo(iline)(1:iendm)
         endif
 901     continue
      enddo
      close (25)
 1021 continue
c
      if (igotminfo.eq.0) then
         write(iup,*) 'Model info file corrupted, empty, or'
         write(iup,*) 'unavailable, so model info won''t be plotted'
         write(iup,*) 'at bottom of frame.  Proceeding.'
      endif
c
      return
      end
