c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine smoothrain(arr,miy,mjx,mkzh,pslab1,pslab2,
     &   mabpl,morpl,nsmth)
c
c   Smooths rain field in PBL.
c   This should "fix" ugly looking CFL (?) errors in rain field.
c
      dimension arr(miy,mjx,mkzh),pslab1(mabpl,morpl),
     &   pslab2(mabpl,morpl)
c
      do k=27,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pslab1(j,i)=arr(i,j,k)
         enddo
         enddo
         call smooth(pslab1,pslab2,nsmth,mabpl,mjx-1,miy-1)
         do j=1,mjx-1
         do i=1,miy-1
            arr(i,j,k)=pslab1(j,i)
         enddo
         enddo
      enddo
c
      return
      end
