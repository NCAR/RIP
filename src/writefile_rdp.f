c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine writefile_rdp (var,varname,ndim,icd,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      dimension var(miy,mjx,1+(mkzh-1)*(ndim-2))
      character varname*10,fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c ensure data do not exceed ieee limits (but only do this for array
c elements for which there is actual data--set non-data elements to
c zero).
c
      if (icd.le.1) then
         icdilim=icd   ! B-grid data goes only to i=miy-icd
      else   ! E-grid data
         icdilim=0     ! E-grid data always goes all the way to i=miy
      endif
      kk = 1+(mkzh-1)*(ndim-2)
      do k=1,kk
         do i=1,miy-icdilim
            if (icd.le.1) then
c              B-grid data goes to j=mjx-icd
               icdjlim=icd
            elseif (icd.eq.2) then
c              E-grid vel. data goes to j=mjx-1 for odd-numbered rows
               icdjlim=mod(i,2)
            elseif (icd.eq.3) then
c              E-grid mass data goes to j=mjx-1 for even-numbered rows
               icdjlim=mod(i+1,2)
            endif
            do j=1,mjx-icdjlim
              if(var(i,j,k) .gt. 3.4e+38) var(i,j,k)= 3.4e+38
              if(var(i,j,k) .lt.-3.4e+38) var(i,j,k)=-3.4e+38
              if(abs(var(i,j,k)).lt. 1.2e-38) var(i,j,k)= 0.
            enddo
            if (icdjlim.eq.1) var(i,mjx,k)=0.
         enddo
         if (icdilim.eq.1) then
            do j=1,mjx
               var(miy,j,k)=0.
            enddo
         endif
      enddo
c
      iendv=index(varname,' ')-1
      if (iendv.eq.-1) iendc=LEN(varname)
      if (iendf1+iendv .gt. LEN(fname)) then
         print*, 'increase length of fname greater than ',iendf1+iendv
         STOP 'writefile_rdp: insufficient space for character var'
      endif
      fname(iendf1+1:)=varname
      iendt = LEN(fname)
      do n = iendt, 1, -1
         ich = ichar(fname(n:n))
         if (.not. ( (ich.ge.65 .and. ich.le.90) .or.     ! Letters A-Z
     &               (ich.ge.97 .and. ich.le.122) .or.    ! Letters a-z
     &               (ich.ge.45 .and. ich.le.58) .or.     ! digits 0-9, also [-./]
     &               (ich.eq.95) .or. ich.eq.0 ) ) then   ! underscore, null
            iendt = n
         endif
      enddo
      open(unit=65,file=fname(1:iendt),form='unformatted')
      ihrip(6)=ndim ! number of dimensions of this variable (2 or 3)
      ihrip(7)=icd  ! grid of this var. (1: cross point, 0: dot point)
      write(65) vardesc,plchun,ihrip,rhrip,chrip
      if (iexpanded.eq.1.and.iexpandedout.eq.0) then
         write(65) (((var(i,j,k),i=1+ioffexp,miy-ioffexp),
     &      j=1+joffexp,mjx-joffexp),k=1,1+(mkzh-1)*(ndim-2))
      else
         write(65) var
      endif
      close (65)
      return
      end
