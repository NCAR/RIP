c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine refmos3calc(fname,dbz,ght,miy,mjx,mkzh,
     &   nxmos,nymos,nzmos,ichoice)
c
      dimension dbz(miy,mjx,mkzh),ght(miy,mjx,mkzh)
      dimension reflz(miy,mjx,nzmos),np(miy,mjx,nzmos)
      integer*2 ireflllz(nxmos,nymos,nzmos),imaxr
      integer*2 ireflllz2(nxmos,nymos,nzmos)
      character fname*256,zlevchar(nzmos)*10
c
      include 'comconst'
c
      include "netcdf.inc"
c
      nf_status = nf_open (fname, nf_nowrite, ncid)
      nf_status = nf_inq_dimid (ncid, 'x', dimid)
      nf_status = nf_inq_dimlen (ncid, dimid, nxmosch)
      nf_status = nf_inq_dimid (ncid, 'y', dimid)
      nf_status = nf_inq_dimlen (ncid, dimid, nymosch)
      nf_status = nf_inq_dimid (ncid, 'mrefl_mosaic_levels', dimid)
      nf_status = nf_inq_dimlen (ncid, dimid, nzmosch)
      if (nxmosch.ne.nxmos.or.nymosch.ne.nymos.or.nzmosch.ne.nzmos) then
         stop
      endif
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'xMin', xmindeg)
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'yMin', ymindeg)
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'dx', dxdeg)
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'dy', dydeg)
c
c   Get reflectivity 3D array
c
      nf_status = nf_inq_varid (ncid, 'mrefl_mosaic', ivarid)
      nf_status = nf_get_var_int2 (ncid, ivarid, ireflllz)
c
c   Very special case: patch in a piece of a different mosaic time
c   because a key radar is missing at the time desired
c
      if (fname.eq.'data/20050213_1800_tile3_3d') then
         nf_status = nf_open ('data/20050213_1815_tile3_3d',
     &       nf_nowrite, ncid2)
         nf_status = nf_inq_varid (ncid2, 'mrefl_mosaic', ivarid)
         nf_status = nf_get_var_int2 (ncid2, ivarid, ireflllz2)
         print*,'Doing patch job.'
         do ly=1,nymos
            rlatmos=ymindeg+(ly-1)*dydeg
         do lx=1,nxmos
            rlonmos=xmindeg+(lx-1)*dxdeg
            dist=sqrt((rlatmos-44.85)**2+(.707*(rlonmos+93.57))**2)
            if (dist.le.1.4) then
               do lz=1,nzmos
                  ireflllz(lx,ly,lz)=
     &               max(ireflllz2(lx,ly,lz),ireflllz(lx,ly,lz))
               enddo
            endif
         enddo
         enddo
      endif
c
c     If ichoice=2, get the max value in each column, and assign it to
c     the lz=1 slab.  Also set nzmos=1.
c
      if (ichoice.eq.2) then
         do ly=1,nymos
         do lx=1,nxmos
            imaxr=-9991
            do lz=1,nzmos
               imaxr=max(imaxr,ireflllz(lx,ly,lz))
            enddo
            ireflllz(lx,ly,1)=imaxr
         enddo
         enddo
         nzmosdo=1
      else
         nzmosdo=nzmos
      endif
c
c   Get levels character array
c
      nf_status = nf_inq_varid (ncid, 'mrefl_mosaicLevels', ivarid)
      nf_status = nf_get_var_text (ncid, ivarid, zlevchar)
c
c     Go through the reflectivity array and add each value to a model
c     grid box.  Also keep track of the number of values added to each
c     grid box so that the average can be calculated for that grid box
c     after all is said and done.
c
      do lz=1,nzmosdo
      do j=1,mjx-1
      do i=1,miy-1
         reflz(i,j,lz)=0.
         np(i,j,lz)=0
      enddo
      enddo
      enddo
      do ly=1,nymos
         rlatmos=ymindeg+(ly-1)*dydeg
      do lx=1,nxmos
         rlonmos=xmindeg+(lx-1)*dxdeg
         call maptform(yival,xjval,rlatmos,rlonmos,-1)
         yival=1.+(yival-yicorn)*refrat
         xjval=1.+(xjval-xjcorn)*refrat
         iy=nint(yival)
         jx=nint(xjval)
         if (iy.ge.1.and.iy.le.miy-1.and.jx.ge.1.and.jx.le.mjx-1) then
            do lz=1,nzmosdo
c               if (lx.lt.10.and.ly.lt.10.and.lz.lt.10) then
c                  print*,'lx,ly,lz,ireflllz(lx,ly,lz)=',
c     &               lx,ly,lz,ireflllz(lx,ly,lz)
c               endif
c
c     Note, in 3D mosaic file, header info shows a "fill value" of
c     -9990, but in the data there are "missing data" values of both
c     -9990 and -990.  The former seem to appear at the bottoms of
c     columns (below-ground?)  and the latter above (missing data above
c     the ground?)
c
               if (ireflllz(lx,ly,lz).ne.-9990.and.
     &             ireflllz(lx,ly,lz).ne.-990.) then
                  refl=10.**(0.01*float(ireflllz(lx,ly,lz)))
               else
                  refl=.01   ! set to minimum value in this data set (equiv. to -20 dBZ)
               endif
               reflz(iy,jx,lz)=reflz(iy,jx,lz)+refl
               np(iy,jx,lz)=np(iy,jx,lz)+1
            enddo
         endif
      enddo
      enddo
      do lz=1,nzmosdo
      do j=1,mjx-1
      do i=1,miy-1
         if (np(i,j,lz).gt.0) then
            reflz(i,j,lz)=reflz(i,j,lz)/np(i,j,lz)
            reflz(i,j,lz)=10.*log10(reflz(i,j,lz))
         else
            reflz(i,j,lz)=-20.
         endif
      enddo
      enddo
      enddo
c
      if (ichoice.ne.2) then
c
c   Now interpolate reflectivity (in dBZ) to model height levels
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         if (ght(i,j,k).le.5000.) then
            rlev=1.+.002*(ght(i,j,k)-1000.)
         else
            rlev=9.+.001*(ght(i,j,k)-5000.)
         endif
         lzdn=int(rlev)
         lzup=lzdn+1
         if (lzdn.lt.1) then
            dbz(i,j,k)=reflz(i,j,1)
         elseif (lzup.gt.nzmos) then
            dbz(i,j,k)=reflz(i,j,nzmos)
         else
            facdn=lzup-rlev
            facup=rlev-lzdn
            dbz(i,j,k)=facdn*reflz(i,j,lzdn)+facup*reflz(i,j,lzup)
         endif
      enddo
      enddo
      enddo
c
      else
c
      do j=1,mjx-1
      do i=1,miy-1
         dbz(i,j,1)=reflz(i,j,1)
      enddo
      enddo
c
      endif
c
      return
      end
