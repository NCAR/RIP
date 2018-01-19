c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine rotpslab(pslab,pslabt,mabpl,morpl,njx,niy,irota)
c
      dimension pslab(mabpl,morpl),pslabt(mabpl,morpl)
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
            pslabt(j,i)=pslab(i,njx-j+1)
         enddo
         enddo
      elseif(irota.eq.-90) then
         do i=1,niy  ! new niy
         do j=1,njx  ! new njx
            pslabt(j,i)=pslab(niy-i+1,j)
         enddo
         enddo
      elseif(irota.eq.180.or.irota.eq.-180) then
         do i=1,niy
         do j=1,njx
            pslabt(j,i)=pslab(njx-j+1,niy-i+1)
         enddo
         enddo
      endif
c
      do i=1,niy
      do j=1,njx
         pslab(j,i)=pslabt(j,i)
      enddo
      enddo
c
      return
      end
