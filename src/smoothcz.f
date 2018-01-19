c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine smoothcz(arr,icd,scr3a,ght,pslab1,pslab2,
     &   numpas,miy,mjx,mkzh,mabpl,morpl)
c
      dimension scr3a(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),arr(miy,mjx,mkzh),pslab1(mabpl,morpl),
     &   pslab2(mabpl,morpl)
c
      dimension hlev(200),valc(200)
c
      include 'comconst'
c
      if (mod(numpas,100).eq.0) return
c
c   First, calculate height.
c
      hmax=-9e9
      hmin=9e9
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=ght(i,j,k)
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
            hmax=max(hmax,scr3a(i,j,k))
            hmin=min(hmin,scr3a(i,j,k))
         enddo
         enddo
      enddo
      dh=(hmax-hmin)/(mkzh-1.)
      dhi=1./dh
c
c   Define height levels for interpolation
c
      do k=1,mkzh
         hlev(k)=hmax-(k-1.)*dh
      enddo
c
c   Do the interpolation to height
c
      do j = 1, mjx-icd
      do i = 1, miy-icd
c
      do kh=1,mkzh
         if (hlev(kh).lt.scr3a(i,j,mkzh).or.
     &       hlev(kh).gt.scr3a(i,j,1)) then
            valc(kh)=rmsg
         else
            do ks=1,mkzh-1
               if (hlev(kh).ge.scr3a(i,j,ks+1).and.
     &             hlev(kh).le.scr3a(i,j,ks)) then
                  if (arr(i,j,ks+1).ne.rmsg.and.
     &                arr(i,j,ks).ne.rmsg) then
                     valc(kh)=((hlev(kh)-scr3a(i,j,ks))*arr(i,j,ks+1)+
     &                      (scr3a(i,j,ks+1)-hlev(kh))*arr(i,j,ks))/
     &                     (scr3a(i,j,ks+1)-scr3a(i,j,ks))
                  else
                     valc(kh)=rmsg
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
      do kh=1,mkzh
         if (arr(i,j,kh).ne.rmsg) then
            kvmax=max(kvmax,kh)
            kvmin=min(kvmin,kh)
         endif
      enddo
      if (kvmax-kvmin.lt.1) then
         do ks=1,mkzh
            arr(i,j,ks)=rmsg
         enddo
      else
         do ks=1,mkzh
            if (scr3a(i,j,ks).lt.hlev(kvmax)) then
               valc(ks)=((scr3a(i,j,ks)-hlev(kvmax-1))*arr(i,j,kvmax)+
     &            (hlev(kvmax)-scr3a(i,j,ks))*arr(i,j,kvmax-1))*dhi
            elseif (scr3a(i,j,ks).gt.hlev(kvmin)) then
               valc(ks)=((scr3a(i,j,ks)-hlev(kvmin))*arr(i,j,kvmin+1)+
     &            (hlev(kvmin+1)-scr3a(i,j,ks))*arr(i,j,kvmin))*dhi
            else
               do kh=kvmin,kvmax-1
                  if (scr3a(i,j,ks).ge.hlev(kh+1).and.
     &                scr3a(i,j,ks).le.hlev(kh)) then
                     valc(ks)=-((scr3a(i,j,ks)-hlev(kh))*arr(i,j,kh+1)+
     &                   (hlev(kh+1)-scr3a(i,j,ks))*arr(i,j,kh))*dhi
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
