c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine smoothcp(arr,icd,scr3a,prs,pslab1,pslab2,
     &   numpas,miy,mjx,mkzh,mabpl,morpl)
c
      dimension scr3a(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),arr(miy,mjx,mkzh),pslab1(mabpl,morpl),
     &   pslab2(mabpl,morpl)
c
      dimension plev(200),valc(200)
c
      include 'comconst'
c
      if (mod(numpas,100).eq.0) return
c
c   First, calculate pressure.
c
      pmax=-9e9
      pmin=9e9
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=prs(i,j,k)
         enddo
         enddo
         if (icd.eq.0) then
            do i=miy-1,2,-1
            do j=mjx-1,2,-1
               scr3a(i,j,k)=.25*(scr3a(i-1,j-1,k)+scr3a(i,j-1,k)+
     &                          scr3a(i-1,j,k)+scr3a(i,j,k))
            enddo
            enddo
            scr3a(1,mjx,k)=scr3a(2,mjx-1,k)
            scr3a(miy,mjx,k)=scr3a(miy-1,mjx-1,k)
            scr3a(1,1,k)=scr3a(2,2,k)
            scr3a(miy,1,k)=scr3a(miy-1,2,k)
            do i=2,miy-1
               scr3a(i,1,k)=scr3a(i,2,k)
               scr3a(i,mjx,k)=scr3a(i,mjx-1,k)
            enddo
            do j=2,mjx-1
               scr3a(1,j,k)=scr3a(2,j,k)
               scr3a(miy,j,k)=scr3a(miy-1,j,k)
            enddo
         endif
         do j=1,mjx-icd
         do i=1,miy-icd
            pmax=max(pmax,scr3a(i,j,k))
            pmin=min(pmin,scr3a(i,j,k))
         enddo
         enddo
      enddo
      dp=(pmax-pmin)/(mkzh-1.)
      dpi=1./dp
c
c   Define pressure levels for interpolation
c
      do k=1,mkzh
         plev(k)=pmin+(k-1.)*dp
      enddo
c
c   Do the interpolation to pressure
c
      do j = 1, mjx-icd
      do i = 1, miy-icd
c
      do kp=1,mkzh
         if (plev(kp).gt.scr3a(i,j,mkzh).or.
     &       plev(kp).lt.scr3a(i,j,1)) then
            valc(kp)=rmsg
         else
            do ks=1,mkzh-1
               if (plev(kp).le.scr3a(i,j,ks+1).and.
     &             plev(kp).ge.scr3a(i,j,ks)) then
                  if (arr(i,j,ks+1).ne.rmsg.and.
     &                arr(i,j,ks).ne.rmsg) then
                     valc(kp)=((plev(kp)-scr3a(i,j,ks))*arr(i,j,ks+1)+
     &                      (scr3a(i,j,ks+1)-plev(kp))*arr(i,j,ks))/
     &                     (scr3a(i,j,ks+1)-scr3a(i,j,ks))
                  else
                     valc(kp)=rmsg
                  endif
                  goto 30
               endif
            enddo
 30         continue
         endif
      enddo
c
      do k=1,mkzh
         arr(i,j,k)=valc(k)
      enddo
c
      enddo
      enddo
c
c   Do the smoothing
c
      njx=mjx-icd
      niy=miy-icd
      do k=1,mkzh
         do j=1,njx
         do i=1,niy
            pslab1(j,i)=arr(i,j,k)
         enddo
         enddo
         call smooth(pslab1,pslab2,numpas,mabpl,njx,niy)
         do j=1,njx
         do i=1,niy
            arr(i,j,k)=pslab1(j,i)
         enddo
         enddo
      enddo
c
c   Interpolate back to original model vertical levels.
c
      do j = 1, mjx-icd
      do i = 1, miy-icd
c
      kvmax=-10000
      kvmin=10000
      do kp=1,mkzh
         if (arr(i,j,kp).ne.rmsg) then
            kvmax=max(kvmax,kp)
            kvmin=min(kvmin,kp)
         endif
      enddo
      if (kvmax-kvmin.lt.1) then
         do ks=1,mkzh
            arr(i,j,ks)=rmsg
         enddo
      else
         do ks=1,mkzh
            if (scr3a(i,j,ks).gt.plev(kvmax)) then
               valc(ks)=((scr3a(i,j,ks)-plev(kvmax-1))*arr(i,j,kvmax)+
     &            (plev(kvmax)-scr3a(i,j,ks))*arr(i,j,kvmax-1))*dpi
            elseif (scr3a(i,j,ks).lt.plev(kvmin)) then
               valc(ks)=((scr3a(i,j,ks)-plev(kvmin))*arr(i,j,kvmin+1)+
     &            (plev(kvmin+1)-scr3a(i,j,ks))*arr(i,j,kvmin))*dpi
            else
               do kp=kvmin,kvmax-1
                  if (scr3a(i,j,ks).le.plev(kp+1).and.
     &                scr3a(i,j,ks).ge.plev(kp)) then
                     valc(ks)=((scr3a(i,j,ks)-plev(kp))*arr(i,j,kp+1)+
     &                   (plev(kp+1)-scr3a(i,j,ks))*arr(i,j,kp))*dpi
                     goto 40
                  endif
               enddo
 40            continue
            endif
         enddo
      endif
c
      do k=1,mkzh
         arr(i,j,k)=valc(k)
      enddo
c
      enddo
      enddo

      return
      end
