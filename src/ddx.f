c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ddx(arr,icdin,deriv,xmap,dmap,icdout,miy,mjx)
c
      dimension arr(miy,mjx),deriv(miy,mjx),xmap(miy,mjx),dmap(miy,mjx)
c
      include 'comconst'
c
      dsi=1./ds
c
      if (icdin.eq.icdout) then
         do j=1,mjx-icdin
            jp=min(j+1,mjx-icdin)
            jm=max(j-1,1)
            dg=jp-jm
         do i=1,miy-icdin
            if (icdout.eq.1) then
               dssi=xmap(i,j)*dsi/dg
            else
               dssi=dmap(i,j)*dsi/dg
            endif
            deriv(i,j)=dssi*(arr(i,jp)-arr(i,jm))
         enddo
         enddo
      else
         hdsi=dsi*.5
         do j=2-icdout,mjx-1
            jp=j+icdout
            jm=jp-1
         do i=2-icdout,miy-1
            ip=i+icdout
            im=ip-1
            if (icdout.eq.1) then
               hdssi=xmap(i,j)*hdsi
            else
               hdssi=dmap(i,j)*hdsi
            endif
            deriv(i,j)=hdssi*
     &         (arr(ip,jp)+arr(im,jp)-arr(ip,jm)-arr(im,jm))
         enddo
         enddo
         if (icdout.eq.0) then
            do j=2,mjx-1
               jp=j
               jm=jp-1
               dssi=dmap(1,j)*dsi
               deriv(1,j)=dssi*(deriv(1,jp)-deriv(1,jm))
               dssi=dmap(miy,j)*dsi
               deriv(miy,j)=dssi*(deriv(miy-1,jp)-deriv(miy-1,jm))
            enddo
            do i=1,miy
               deriv(i,1)=deriv(i,2)
               deriv(i,mjx)=deriv(i,mjx-1)
            enddo
         endif
      endif
c
      return
      end
