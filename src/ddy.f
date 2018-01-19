c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine ddy(arr,icdin,deriv,xmap,dmap,icdout,miy,mjx)
c
      dimension arr(miy,mjx),deriv(miy,mjx),xmap(miy,mjx),dmap(miy,mjx)
c
      include 'comconst'
c
      dsi=1./ds
c
      if (icdin.eq.icdout) then
         do j=1,mjx-icdout
         do i=1,miy-icdout
            ip=min(i+1,miy-icdout)
            im=max(i-1,1)
            dg=ip-im
            if (icdout.eq.1.and.nproj.ne.4) then
               dssi=xmap(i,j)*dsi/dg
            elseif (icdout.eq.0.and.nproj.ne.4) then
               dssi=dmap(i,j)*dsi/dg
            else
               dssi=dsi/dg
            endif
            deriv(i,j)=dssi*(arr(ip,j)-arr(im,j))
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
            if (icdout.eq.1.and.nproj.ne.4) then
               hdssi=xmap(i,j)*hdsi
            elseif (icdout.eq.0.and.nproj.ne.4) then
               hdssi=dmap(i,j)*hdsi
            else
               hdssi=hdsi
            endif
            deriv(i,j)=hdssi*
     &         (arr(ip,jp)+arr(ip,jm)-arr(im,jp)-arr(im,jm))
         enddo
         enddo
         if (icdout.eq.0) then
            do i=2,miy-1
               ip=i
               im=ip-1
               if (nproj.ne.4) then
                  dssi1=dmap(i,1)*dsi
                  dssim=dmap(i,mjx)*dsi
               else
                  dssi1=dsi
                  dssim=dsi
               endif
               deriv(i,1)=dssi1*(deriv(ip,1)-deriv(im,1))
               deriv(i,mjx)=dssim*(deriv(ip,mjx-1)-deriv(im,mjx-1))
            enddo
            do j=1,mjx
               deriv(1,j)=deriv(2,j)
               deriv(miy,j)=deriv(miy-1,j)
            enddo
         endif
      endif
c
      return
      end
