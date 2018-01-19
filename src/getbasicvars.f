c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getbasicvars(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,maxtavl,uuu,vvv,tmk,qvp,www,prs,ght,sfp,sfpsm,ter,
     &   dmap,xmap,cor,dskm,miy,mjx,mkzh)
c
c   All of the arguments from uuu through cor are basic arrays for
c   which RIP *must* be able to find a corresponding file in order to
c   work properly (exceptions are that water vapor (qvp) and vertical
c   velocity (www) are not abosolutely necessary--RIP will fill the
c   qvp and/or www arrays with zeros if corresponding files are
c   not found, and proceed happily along).
c
c   All other data fields (hydrometeor mixing ratios, rainfall,
c   surface radiative properties, or whatever else your model spits
c   out) are optional.  RIP can look for the files and plot the data
c   if requested (in subroutine fields), but does not require them.
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*(*)
c
      dimension uuu(miy,mjx,mkzh),vvv(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),www(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),sfp(miy,mjx),sfpsm(miy,mjx),ter(miy,mjx),
     &   dmap(miy,mjx),xmap(miy,mjx),cor(miy,mjx)
c
c   Read 2D variables
c
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'ter       ',miy,mjx,mkzh,maxtavl,2,1,ter,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'sfp       ',miy,mjx,mkzh,maxtavl,2,1,sfp,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'dmap      ',miy,mjx,mkzh,maxtavl,2,0,dmap,istat)
      if (istat.eq.-1) then
         do j=1,mjx-1
         do i=1,miy-1
            dmap(i,j)=xmap(i,j)
         enddo
         enddo
         call xtodot(dmap,miy,mjx)
      endif
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'xmap      ',miy,mjx,mkzh,maxtavl,2,1,xmap,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'cor       ',miy,mjx,mkzh,maxtavl,2,1,cor,istat)
c
c   Read 3D variables
c
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'uuu       ',miy,mjx,mkzh,maxtavl,3,1,uuu,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'vvv       ',miy,mjx,mkzh,maxtavl,3,1,vvv,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'tmk       ',miy,mjx,mkzh,maxtavl,3,1,tmk,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'prs       ',miy,mjx,mkzh,maxtavl,3,1,prs,istat)
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'ght       ',miy,mjx,mkzh,maxtavl,3,1,ght,istat)
c
c   Create sfpsm, a smoothed version of sfp (used for interpolating
c   fields to below-ground locations)
c
      do j=1,mjx-1
      do i=1,miy-1
         sfpsm(i,j)=sfp(i,j)
      enddo
      enddo
      nsm=nint(20./dskm)+400  ! causes smoothing if grid spacing < 40 km.
      call smooth(sfpsm,qvp(1,1,1),nsm,miy,miy-1,mjx-1) ! Temporarily use a slab
c                                                       ! of qvp for work spc
c
c   These 3D fields have associated arrays in RIP, but files
c   are not required.  Arrays are filled with zeros if no files found. 
c
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'qvp       ',miy,mjx,mkzh,maxtavl,3,0,qvp,istat)
      if (istat.eq.-1) then
         call fillarray(qvp,miy*mjx*mkzh,0.)
      else ! convert from g/kg to kg/kg
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            qvp(i,j,k)=.001*qvp(i,j,k)
         enddo
         enddo
         enddo
      endif
      call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,ncxc,
     &     'www       ',miy,mjx,mkzh,maxtavl,3,0,www,istat)
      if (istat.eq.-1) call fillarray(www,miy*mjx*mkzh,0.)
c
      return
      end
