c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,varname,miy,mjx,mkzh,maxtavl,ndim,iabort,arr,istat)
c
      dimension arr(miy,mjx,1+(ndim-2)*(mkzh-1))
      dimension xtimeavl(maxtavl)
      character fname*256,varname*10,cxtimeavl(maxtavl)*10,
     &   casename*(*)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      include 'comconst'
c
c     If you are attempting to re-code rip to read data directly from
c     native model output file(s), this would be the appropriate place
c     to first try to read a requested variable from that file(s).
c     However, if that fails, you should still have the routine attempt
c     to read the variable from a rip data file, in case it was
c     previously saved to such a file by rip.
c
      fname=casename(1:iendc)//'_'//
     &      cxtimeavl(nxt)(1:ncxc)//'_'//varname
      open(unit=25,err=30,file=fname,form='unformatted',
     &   status='old')
      read (25)
     &   vardesc,plchun,ihrip,rhrip,chrip
      ndimch=ihrip(6)
      if (ndim.ne.ndimch) then
         write(iup,*)'In attempting to read the file called ',fname
         write(iup,*)'you are trying to fill an array of dimension'
         write(iup,*)ndim,' with data of dimension ',ndimch
         stop
      endif
c
c   Check for B-grid/E-grid data
c
      icdch=ihrip(7)
      nlevs=1+(ndim-2)*(mkzh-1)
      if (icdch.le.1) then  ! B-grid data
         read (25) arr
         close (25)
         istat=1
c
c      Make sure rmsg values are EXACTLY rmsg.
c      (This is necessary for the Cray, if data file was IEEE)
c
         do k=1,nlevs
         do j=1,mjx-icdch
         do i=1,miy-icdch
            if (abs((arr(i,j,k)-rmsg)/rmsg).lt.1e-6) arr(i,j,k)=rmsg
         enddo
         enddo
         enddo
      else    ! E-grid data
c
c     In the parlance of ripdp_wrfnmm.f, the grid in RIP is a B-grid of
c     dimensions miyeb,mjxeb (i.e., miy,mjx in this routine are equal to
c     miyeb,mjxeb).  However, the data is stored in compact E-grid data
c     arrays of dimension miyec,mjxec.  See ripdp_wrfnmm.f comments for
c     explanation.
c
         miyeb=miy
         mjxeb=mjx
         miyef=miyeb-1
         mjxef=mjxeb-1
         miyec=miyef
         mjxec=(mjxef+1)/2
         read (25) (((arr(i,j,k),i=1,miyec),j=1,mjxec),k=1,nlevs)
         close (25)
         istat=1
c
c      Make sure rmsg values are EXACTLY rmsg.
c      (This is necessary for the Cray, if data file was IEEE)
c
         do k=1,nlevs
         do j=1,mjxec
         do i=1,miyec
            if (abs((arr(i,j,k)-rmsg)/rmsg).lt.1e-6) arr(i,j,k)=rmsg
         enddo
         enddo
         enddo
c
c      Spread/interpolate the data to the full E-grid (both H and V points)
c
         call egridfill(arr,icdch,miyef,mjxef,miy,mjx,nlevs)
c
c      If E-grid data was at V points, transfer to the B-grid
c      "dot" points.
c
         if (icdch.eq.2) then
            do k=1,nlevs
               call xtodot(arr(1,1,k),miy,mjx)
            enddo
         endif
c
      endif
c
      return
 30   istat=-1
      if (iabort.eq.1) then
         write(iup,*)'Fatal: couldn''t find file named'
         write(iup,*)'   ',fname
         stop
      endif
      return
      end
