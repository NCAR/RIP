c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hiresmap(directory,filename,iunit)
c
      character directory*256,filename*32
      character fullpath*256
c
      dimension xlatmap(10000),xlonmap(10000)
c
c   This routine was adapted from GRAPH, for making custom map plots.
c   It assumes all the line width, color, map position, etc. has already
c   been set.  Jim Bresch 5/6/99, Mark Stoelinga 12/9/2002
c
c   The ncarg routines mapit and mapiq plot out the lines and flush
c   the buffers after each line segment is drawn.  The routine returns
c   control when all of the data has been processed.
c
      include 'comconst'

      iendci=index(directory,' ')-1
      fullpath=directory(1:iendci)//'/'//filename
      iendfullpath=index(fullpath,' ')-1
      if (index(fullpath,'.bin').ne.0) then
         open (unit=iunit,file=fullpath,form='unformatted',
     &      status='old',err=33)
         ibin=1
      elseif (index(fullpath,'.ascii').ne.0) then
         open (unit=iunit,file=fullpath,form='formatted',
     &      status='old',err=33)
         ibin=0
      else
         write(iup,*)'   hires map data file name must end in'
         write(iup,*)'   either .ascii or .bin.'
      endif
c
 10   continue
      if (ibin.eq.1) then
         read(iunit,end=1000) numvals,xlatmin,xlatmax,xlonmin,xlonmax
         numpts=numvals/2
         read(iunit)(xlatmap(i),xlonmap(i),i=1,numpts)
      else
c         read(iunit,'(I4,14X,6F9.3)',end=1000) numvals,xlatmin,xlatmax,
c     &      xlonmin,xlonmax,xlatmap(1),xlonmap(1)
         read(iunit,*,end=1000) numvals,xlatmin,xlatmax,
     &      xlonmin,xlonmax,xlatmap(1),xlonmap(1)
         numpts=numvals/2
         read(iunit,'(8F9.3)')(xlatmap(i),xlonmap(i),i=2,numpts)
      endif
c
      call mapit(xlatmap(1),xlonmap(1),0)
      do i=2,numpts
         call mapit(xlatmap(i),xlonmap(i),2)
      enddo
      call mapiq
      goto 10
1000  continue
      close (iunit)
c
      return
c
 33   continue
c
      write(iup,*)'Could not find the custom map data file'
      write(iup,*)'specified by rip_root and outy.'
      write(iup,*)'Looked for file "',fullpath(1:iendfullpath),'"'
c
      return
c
      end
