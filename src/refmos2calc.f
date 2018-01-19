c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine refmos2calc(fname,dbz,miy,mjx,nxmos,nymos,ichoice)
c
      dimension dbz(miy,mjx), np(miy,mjx)
      integer*2 ireflll(nxmos,nymos)
      character fname*256
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
      if (nxmosch.ne.nxmos.or.nymosch.ne.nymos) then
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
c      print*,'nxmosch,nymosch,xmindeg,ymindeg,dxdeg,dydeg='
c      print*,nxmosch,nymosch,xmindeg,ymindeg,dxdeg,dydeg
c
c   Get reflectivity 2D array
c
      if (ichoice.eq.1) then
         nf_status = nf_inq_varid (ncid, 'cref_smooth', ivarid)
      elseif (ichoice.eq.2) then
         nf_status = nf_inq_varid (ncid, 'hsr', ivarid)
      endif
c      print*,'ivarid=',ivarid
      nf_status = nf_get_var_int2 (ncid, ivarid, ireflll)
c
c     Go through the reflectivity array and add each value to a model
c     grid box.  Also keep track of the number of values added to each
c     grid box so that the average can be calculated for that grid box
c     after all is said and done.
c
      do j=1,mjx-1
      do i=1,miy-1
         dbz(i,j)=0.
         np(i,j)=0
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
c
c     Note, in 2D mosaic file, header info shows a "fill value" of
c     -9990, but in the data the only apparent "missing data" value used
c     is -990.
c
            if (ireflll(lx,ly).ne.-9990.and.
     &          ireflll(lx,ly).ne.-990.) then
               refl=10.**(0.01*float(ireflll(lx,ly)))
            else
               refl=.01
            endif
            dbz(iy,jx)=dbz(iy,jx)+refl
            np(iy,jx)=np(iy,jx)+1
         endif
      enddo
      enddo
      do j=1,mjx-1
      do i=1,miy-1
         if (np(i,j).gt.0) then
            dbz(i,j)=dbz(i,j)/np(i,j)
            dbz(i,j)=10.*log10(dbz(i,j))
         else
            dbz(i,j)=-20.
         endif
      enddo
      enddo
c
      return
      end
