c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine xtodot(slab,maxiy,maxjx)
c
c     This routine converts data that is on the B-grid mass grid (known
c     in MM5 lingo as "cross points") to the B-grid velocity staggered
c     grid (known in MM5 lingo as "dot points")
c
      dimension slab(maxiy,maxjx),bot(10000),rleft(10000)
c
c   Extrapolate out to top and bottom edges.
c
      do 200 j=2,maxjx-1
         bot(j)=(3.*(slab(1,j-1)+slab(1,j))-
     &      (slab(2,j-1)+slab(2,j)))/4.
         slab(maxiy,j)=(3.*(slab(maxiy-1,j-1)+slab(maxiy-1,j))-
     &      (slab(maxiy-2,j-1)+slab(maxiy-2,j)))/4.
  200 continue
c
c   Extrapolate out to left and right edges.
c
      do 300 i=2,maxiy-1
         rleft(i)=(3.*(slab(i-1,1)+slab(i,1))-
     &      (slab(i-1,2)+slab(i,2)))/4.
         slab(i,maxjx)=(3.*(slab(i-1,maxjx-1)+slab(i,maxjx-1))-
     &      (slab(i-1,maxjx-2)+slab(i,maxjx-2)))/4.
  300 continue
c
c   Extrapolate out to corners.
c
      rleft(1)=(3.*slab(1,1)-slab(2,2))/2.
      rleft(maxiy)=(3.*slab(maxiy-1,1)-slab(maxiy-2,2))/2.
      bot(maxjx)=(3.*slab(1,maxjx-1)-slab(2,maxjx-2))/2.
      slab(maxiy,maxjx)=(3.*slab(maxiy-1,maxjx-1)-
     &   slab(maxiy-2,maxjx-2))/2.
c
c   Interpolate in the interior.
c
      do 100 j=maxjx-1,2,-1
      do 100 i=maxiy-1,2,-1
         slab(i,j)=.25*(slab(i-1,j-1)+slab(i,j-1)+slab(i-1,j)+
     &      slab(i,j))
  100    continue
c
c   Put "bot" and "rleft" values into slab.
c
      do j=2,maxjx
         slab(1,j)=bot(j)
      enddo
      do i=1,maxiy
         slab(i,1)=rleft(i)
      enddo
c
      return
      end
