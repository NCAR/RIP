      subroutine egridfill(arr,icd,miyef,mjxef,miy,mjx,nlevs)
c
c     This routine takes an array that contains data only on E-grid
c     staggered points (H or V, depending on icd=3 or 2), and spreads
c     the data out to all E-grid points (H and V).  The input array is
c     at least large enough for the full grid (all H and V points) but
c     initially contains data that is "scrunched" in the x direction.
c     The input array may even be larger than necessary, which is why
c     there are two sets of horizontal dimensions: miyef,mjxef are the
c     dimensions of the full E-grid, and miy,mjx are the actual
c     dimensions of the array.
c     
      dimension arr(miy,mjx,nlevs)

         do k=1,nlevs
c
         if (icd.eq.2) then  ! E-grid velocity data
c
c         First spread out the data
c
            do i=1,miyef,2
            do j=mjxef-1,2,-2
               jec=j/2
               arr(i,j,k)=arr(i,jec,k)
            enddo
            enddo
            do i=2,miyef-1,2
            do j=mjxef,1,-2
               jec=(j+1)/2
               arr(i,j,k)=arr(i,jec,k)
            enddo
            enddo
c
c         Next, fill in the corner points
c
            arr(1,1,k)=.25*(3.*(arr(2,1,k)+arr(1,2,k))-
     &                         (arr(3,2,k)+arr(2,3,k)))
            arr(1,mjxef,k)=.25*(3.*(arr(2,mjxef,k)+arr(1,mjxef-1,k))-
     &                             (arr(3,mjxef-1,k)+arr(2,mjxef-2,k)))
            arr(miyef,1,k)=.25*(3.*(arr(miyef-1,1,k)+arr(miyef,2,k))-
     &                         (arr(miyef-2,2,k)+arr(miyef-1,3,k)))
            arr(miyef,mjxef,k)=
     &         .25*(3.*(arr(miyef-1,mjxef,k)+arr(miyef,mjxef-1,k))-
     &                 (arr(miyef-2,mjxef-1,k)+arr(miyef-1,mjxef-2,k)))
c
            do i=3,miyef-2,2  ! Fill in left and right side
               arr(i,1,k)=.5*(arr(i+1,1,k)+arr(i-1,1,k))
               arr(i,mjxef,k)=.5*(arr(i+1,mjxef,k)+arr(i-1,mjxef,k))
            enddo
c
            do j=3,mjxef-2,2  ! Fill in bototm and top sides.
               arr(1,j,k)=.5*(arr(1,j+1,k)+arr(1,j-1,k))
               arr(miyef,j,k)=.5*(arr(miyef,j+1,k)+arr(miyef,j-1,k))
            enddo
c
            do i=2,miyef-1,2  ! Fill in interior even rows
            do j=2,mjxef-1,2
               arr(i,j,k)=.25*(arr(i,j+1,k)+arr(i,j-1,k)+
     &                         arr(i+1,j,k)+arr(i-1,j,k))
            enddo
            enddo
c
            do i=3,miyef-2,2  ! Fill in interior odd rows
            do j=3,mjxef-2,2
               arr(i,j,k)=.25*(arr(i,j+1,k)+arr(i,j-1,k)+
     &                         arr(i+1,j,k)+arr(i-1,j,k))
            enddo
            enddo
c
         elseif (icd.eq.3) then  ! E-grid mass data
c
c         First spread out the data
c
            do i=2,miyef-1,2
            do j=mjxef-1,2,-2
               jec=j/2
               arr(i,j,k)=arr(i,jec,k)
            enddo
            enddo
            do i=1,miyef,2
            do j=mjxef,1,-2
               jec=(j+1)/2
               arr(i,j,k)=arr(i,jec,k)
            enddo
            enddo
c
            do i=2,miyef-1,2  ! Fill in left and right side
               arr(i,1,k)=.5*(arr(i+1,1,k)+arr(i-1,1,k))
               arr(i,mjxef,k)=.5*(arr(i+1,mjxef,k)+arr(i-1,mjxef,k))
            enddo
c
            do j=2,mjxef-1,2  ! Fill in bototm and top sides.
               arr(1,j,k)=.5*(arr(1,j+1,k)+arr(1,j-1,k))
               arr(miyef,j,k)=.5*(arr(miyef,j+1,k)+arr(miyef,j-1,k))
            enddo
c
            do i=2,miyef-1,2  ! Fill in interior even rows
            do j=3,mjxef-2,2
               arr(i,j,k)=.25*(arr(i,j+1,k)+arr(i,j-1,k)+
     &                         arr(i+1,j,k)+arr(i-1,j,k))
            enddo
            enddo
c
            do i=3,miyef-2,2  ! Fill in interior odd rows
            do j=2,mjxef-1,2
               arr(i,j,k)=.25*(arr(i,j+1,k)+arr(i,j-1,k)+
     &                         arr(i+1,j,k)+arr(i-1,j,k))
            enddo
            enddo
c
         endif
c
         enddo   ! end of k-loop
c
         return
         end
