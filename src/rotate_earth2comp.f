

      subroutine rotate_earth2comp(ue1,ve1,stdlon,xlon_u,xlon_v,nproj,
     &                              uu1,vv1)
c
c    This subroutine rotates the uue,vve winds from MPAS which are earth 
c    relative to grid relative (uuu,vvv).

      real diff
      real alpha

      include "commptf"

c      ROTATE U FIELD
      diff = stdlon - xlon_u
      if (diff > 180.) then
         diff = diff - 360.
      else if (diff < -180.) then
         diff = diff + 360.
      end if
c       Calculate the rotation angle, alpha, in radians
      if (nproj == 1) then
c       if using lambert conformal
         alpha = diff * cone_mptf * rpd_mptf * ihm_mptf
      else
         alpha = diff * rpd_mptf * ihm_mptf
      end if
c      Calculate U on computational grid
      uu1=ue1*cos(alpha)+ve1*sin(alpha)



c      ROTATE V FIELD
      diff = stdlon - xlon_v
      if (diff > 180.) then
         diff = diff - 360.
      else if (diff < -180.) then
         diff = diff + 360.
      end if
c       Calculate the rotation angle, alpha, in radians
      if (nproj == 1) then
c       if using lambert conformal
         alpha = diff * cone_mptf * rpd_mptf * ihm_mptf
      else
         alpha = diff *  rpd_mptf * ihm_mptf
      end if
c      Calculate V on computational grid
      vv1=ve1*cos(alpha)-ue1*sin(alpha)


      return
      end
