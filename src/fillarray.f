c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine fillarray(array,ndim,val)
      dimension array(ndim)
      do 10 i=1,ndim
         array(i)=val
   10 continue

      return
      end
