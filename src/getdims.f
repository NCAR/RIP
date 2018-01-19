c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getdims(casename,iendc,xtimeavl,cxtimeavl,
     &   maxtavl,ncxc,miy,mjx,mkzh,mabpl,morpl,iup)
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*(*),fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   Create a file name for a file that is likely to exist,
c   just so we can read it and get the dimensions of the dataset.
c
      fname=casename(1:iendc)//'_'//cxtimeavl(1)(1:ncxc)//'_'//'ter'
      open(unit=25,file=fname,form='unformatted',status='old')
      read(25,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip
c
      miy=ihrip(4)
      mjx=ihrip(5)
      icd=ihrip(7)
      if (icd.ge.2) then
c
c     This is E-grid data.  The dimensions are for the compact E-grid
c     data arrays.  Change the dimensions to be the number of dot points
c     in a B-grid whose cross points are coincident with the E-grid H
c     and V points.  In the parlance of ripdp_wrfnmm.f, the values read
c     in here are miyec,mjxec, but we want to convert them to miyeb,
c     mjxeb.
c
         miy=miy+1
         mjx=2*mjx
      endif
      mkzh=ihrip(9)
      mmax=max(miy,mjx)
      mabpl=max(mmax,min(2*mmax-1,401))
      morpl=max(max(miy,mjx),mkzh)+1
c
      close (25)
      goto 200
c
 170  write(iup,*) 'The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  write(iup,*) 'Unexpected EOF reached when trying to read'
      write(iup,*) 'model data header.  Stopping.'
      stop
c
 200  continue
c
      return
      end
