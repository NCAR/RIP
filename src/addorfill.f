c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine addorfill(arr1,arr2,miy,mjx,mkzh,ndim,icd,fac1,fac2)
c
      dimension arr1(miy,mjx,1+(ndim-2)*(mkzh-1)),
     &   arr2(miy,mjx,1+(ndim-2)*(mkzh-1))
c
      include 'comconst'
c
      nslab=1+(ndim-2)*(mkzh-1)
      do k=1,nslab
      do j=1,mjx-icd
      do i=1,miy-icd
         arr2(i,j,k)=fac1*arr1(i,j,k)+fac2*arr2(i,j,k)
      enddo
      enddo
      enddo
      return
      end
