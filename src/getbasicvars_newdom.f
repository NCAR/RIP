c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getbasicvars_newdom(uuu,vvv,tmk,qvp,www,prs,ght,sfp,
     &   sfpsm,ter,dmap,xmap,cor,miy,mjx,mkzh,yicorn,xjcorn,refrat,
     &   iunewdom)
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
      dimension uuu(miy,mjx,mkzh),vvv(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),www(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),sfp(miy,mjx),sfpsm(miy,mjx),ter(miy,mjx),
     &   dmap(miy,mjx),xmap(miy,mjx),cor(miy,mjx)
      character fname*256
c
c   RIP header variables (for source terrain RIP file)
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      rewind(iunewdom)
      do i=1,10
         read(iunewdom,*)
      enddo
      read(iunewdom,*)istate
c
      if (istate.eq.1) then
c
c     This is a dry quiescent atmosphere, flat terrain, real-earth
c     coriolis and map factor (based on actual domain location), uniform
c     surface pressure, isothermal (273 K), hydrostatically balanced
c
      read(iunewdom,*) ztop    ! (in km)
      dz=ztop/mkzh
c
c   Calculate ter, xmap, cor, sfp, ght, and prs
c
      twoomega=1.458426e-4
      pi=4.*atan(1.)
      rpd=pi/180.
c
      do j=1,mjx-1
      do i=1,miy-1
         ter(i,j)=0.
         riy=yicorn+(.5+i-1.)/refrat
         rjx=xjcorn+(.5+j-1.)/refrat
         xmap(i,j)=xmapcalc(riy,rjx)
         call maptform(riy,rjx,rlat,rlon,1)
         cor(i,j)=twoomega*sin(rlat*rpd)
         sfp(i,j)=1013.
         do k=1,mkzh
            ght(i,j,k)=(k-.5)*dz
            prs(i,j,k)=1013.*exp(-ght(i,j,k)/(287.*273./9.81))
         enddo
      enddo
      enddo
c
c   Calculate dmap
c
      do j=1,mjx
      do i=1,miy
         riy=yicorn+(i-1.)/refrat
         rjx=xjcorn+(j-1.)/refrat
         dmap(i,j)=xmapcalc(riy,rjx)
      enddo
      enddo
c
      call fillarray(uuu,miy*mjx*mkzh,10.)
      call fillarray(vvv,miy*mjx*mkzh,10.)
      call fillarray(tmk,miy*mjx*mkzh,273.)
      call fillarray(qvp,miy*mjx*mkzh,0.)
      call fillarray(www,miy*mjx*mkzh,0.)
c
c   Calculate sfpsm
c
      do j=1,mjx-1
      do i=1,miy-1
         sfpsm(i,j)=sfp(i,j)
      enddo
      enddo
c
      elseif (istate.eq.2) then
c
c     This is a 2-D midlatitude front, dry, hydrostatically and
c     geotrophically balanced, flat terrain, real-earth, constant
c     coriolis and no map factor, with a geostrophic deformation field.
c     Vertical coordinate is sigma-z [ (z-zsfc)/(ztop-zsfc) ], and levels
c     are evenly spaced in sigma-z.
c
      call contrive2(cor,dmap,xmap,ter,prs,ght,
     &   tmk,qvp,uuu,vvv,www,sfp,sfpsm,miy,mjx,mkzh)
c
c      elseif (istate.eq.3) then
c
c     This is a 2-D midlatitude front, based on work of Bannon (1984).
c
c      call contrive3(cor,dmap,xmap,ter,prs,ght,
c     &   tmk,qvp,uuu,vvv,www,sfp,sfpsm,miy,mjx,mkzh)
c
      endif
c
      return
      end
