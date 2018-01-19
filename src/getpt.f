c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getpt(miy,mjx,mkzh,ifree,ndim,ipoint,wk,maxslab)
c
      dimension wk(miy,mjx,maxslab)
      pointer(ipoint,dumdum)
c
c   Note, the variable dumdum is never used.  Its pointer is defined
c   just to insure that ipoint, being a pointer rather than an integer
c   variable, will always be compatible in size (32 or 64 bit) with
c   the corresponding pointer variable in the call to getpt.  This way
c   we don't have to worry about declaring ipoint as integer*4 or
c   integer*8 depending on machine type.
c
      include 'comconst'
c
      nslab=1+(ndim-2)*(mkzh-1)
c     write(6,*) 'nslab = ',nslab
      ifreenew=ifree+nslab
c     write(6,*) 'ifreenew = ',ifreenew
      if (ifreenew.gt.maxslab) then
         write(iup,*)'The program is trying to allocate more work'
         write(iup,*)'space than is available.  Try increasing maxfld'
         write(iup,*)'and run rip again.'
         stop
      endif
c     write(6,*) 'ifree = ',ifree
      ipoint=loc(wk(1,1,ifree))
c     write(6,*) 'ipoint = ',ipoint
      ifree=ifreenew
      return
      end
