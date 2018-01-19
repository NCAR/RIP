c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cfad(work,ght,tmk,rcfad,miy,mjx,mkzh)
c
      dimension work(miy,mjx,mkzh),ght(miy,mjx,mkzh),tmk(miy,mjx,mkzh)
c


c
c   Do CFAD counting and print out
c
      if (ncfadbin.gt.0) then
         do ib=1,ncfadbin
            ncfadcount(ib)=0
         enddo
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg) then
               ibin=nint((pslab1(j,i)-rcfad(1,ipl))/rcfad(2,ipl))+1
               if (ibin.ge.1.and.ibin.le.ncfadbin) then
                  ncfadcount(ibin)=ncfadcount(ibin)+1
               endif
            endif
         enddo
         enddo
         do ib=1,ncfadbin
            write(58,*),rlevl(ilev,ipl),
     &                 rcfad(1,ipl)+(ib-1)*rcfad(2,ipl),
     &                 ncfadcount(ib)
         enddo
      endif
