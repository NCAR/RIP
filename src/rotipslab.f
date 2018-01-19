c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine rotipslab(ipslab,ipslabt,mabpl,morpl,njx,niy,irota)
c
      dimension ipslab(mabpl,morpl),ipslabt(mabpl,morpl)
c
      if (irota.eq.0) return
      if(irota.eq.90.or.irota.eq.-90) then
         nsv=niy
         niy=njx
         njx=nsv
      endif
      if(irota.eq.90) then
         do i=1,niy  ! new niy
         do j=1,njx  ! new njx
            ipslabt(j,i)=ipslab(i,njx-j+1)
         enddo
         enddo
      elseif(irota.eq.-90) then
         do i=1,niy  ! new niy
         do j=1,njx  ! new njx
            ipslabt(j,i)=ipslab(niy-i+1,j)
         enddo
         enddo
      elseif(irota.eq.180.or.irota.eq.-180) then
         do i=1,niy
         do j=1,njx
            ipslabt(j,i)=ipslab(njx-j+1,niy-i+1)
         enddo
         enddo
      endif
c
      do i=1,niy
      do j=1,njx
         ipslab(j,i)=ipslabt(j,i)
      enddo
      enddo
c
      return
      end
