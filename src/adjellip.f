c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine adjellip(adj,thresh,ipttp,iall,nscrs,mkp,
     &   tot,totdiff)
c
      dimension adj(nscrs,mkp),thresh(nscrs,mkp),ipttp(nscrs,mkp)
c
      include 'comconst'
c
      totdiff=0.
      tot=0.
      npts=0
      do k=1,mkp
      do ls=1,nscrs
         if ((iall.eq.0.and.ipttp(ls,k).eq.1).or.
     &       (iall.eq.1.and.ipttp(ls,k).ge.1)) then
            npts=npts+1
            tot=tot+adj(ls,k)
            if (adj(ls,k).lt.thresh(ls,k)) then
               diff=thresh(ls,k)-adj(ls,k)
               totdiff=totdiff+diff
               adj(ls,k)=thresh(ls,k)
            endif
         endif
      enddo
      enddo
      totdiffe=totdiff
c
      do iter=1,3
         diffpp=totdiffe/npts
         do k=1,mkp
         do ls=1,nscrs
            if ((iall.eq.0.and.ipttp(ls,k).eq.1).or.
     &          (iall.eq.1.and.ipttp(ls,k).ge.1)) then
               diffppt=min(diffpp,adj(ls,k)-thresh(ls,k))
               adj(ls,k)=adj(ls,k)-diffppt
               totdiffe=totdiffe-diffppt
            endif
         enddo
         enddo
      enddo
c
      return
      end
