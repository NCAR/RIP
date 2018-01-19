c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine dbzcalc(qvp,qra,qsn,qgr,tmk,prs,dbz,miy,mjx,mkzh,
     &   in0r,in0s,in0g,iliqskin)
c
c     This routine computes equivalent reflectivity factor (in dBZ) at
c     each model grid point.  In calculating Ze, the RIP algorithm makes
c     assumptions consistent with those made in an early version
c     (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
c     model (i.e., the scheme known as "Resiner-2").  For each species:
c
c     1. Particles are assumed to be spheres of constant density.  The
c     densities of rain drops, snow particles, and graupel particles are
c     taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
c     rho_g = 400 kg m^-3, respectively. (l refers to the density of
c     liquid water.)
c
c     2. The size distribution (in terms of the actual diameter of the
c     particles, rather than the melted diameter or the equivalent solid
c     ice sphere diameter) is assumed to follow an exponential
c     distribution of the form N(D) = N_0 * exp( lambda*D ).
c
c     3. If in0X=0, the intercept parameter is assumed constant (as in
c     early Reisner-2), with values of 8x10^6, 2x10^7, and 4x10^6 m^-4,
c     for rain, snow, and graupel, respectively.  Various choices of
c     in0X are available (or can be added).  Currently, in0X=1 gives the
c     variable intercept for each species that is consistent with
c     Thompson, Rasmussen, and Manning (2004, Monthly Weather Review,
c     Vol. 132, No. 2, pp. 519-542.)
c
c     4. If iliqskin=1, frozen particles that are at a temperature above
c     freezing are assumed to scatter as a liquid particle.
c
c     More information on the derivation of simulated reflectivity in RIP
c     can be found in Stoelinga (2005, unpublished write-up).  Contact
c     Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.
c
      dimension qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   qsn(miy,mjx,mkzh),qgr(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),dbz(miy,mjx,mkzh)
      character meth*1
c
      include 'comconst'
c
c   Constant intercepts
c
      rn0_r = 8.e6    ! m^-4
      rn0_s = 2.e7    ! m^-4
      rn0_g = 4.e6    ! m^-4
c
c   Constants used to calculate variable intercepts
c
      r1=1.e-15
      ron=8.e6
      ron2=1.e10
      son=2.e7
      gon=5.e7
      ron_min = 8.e6
      ron_qr0 = 0.00010
      ron_delqr0 = 0.25*ron_qr0
      ron_const1r = (ron2-ron_min)*0.5
      ron_const2r = (ron2+ron_min)*0.5
c
c   Other constants:
c
      gamma_seven = 720.
      rho_r = rhowat ! 1000. kg m^-3
      rho_s = 100.   ! kg m^-3
      rho_g = 400.   ! kg m^-3
      alpha = 0.224
      factor_r = gamma_seven * 1.e18 * (1./(pi*rho_r))**1.75
      factor_s = gamma_seven * 1.e18 * (1./(pi*rho_s))**1.75
     &    * (rho_s/rhowat)**2 * alpha
      factor_g = gamma_seven * 1.e18 * (1./(pi*rho_g))**1.75
     &    * (rho_g/rhowat)**2 * alpha
c
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
c
         rhoair=prs(i,j,k)*100./
     &      (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))      ! air density
c
c      Adjust factor for brightband, where snow or graupel particle
c      scatters like liquid water (alpha=1.0) because it is assumed to
c      have a liquid skin.
c
         if (iliqskin.eq.1.and.tmk(i,j,k).gt.celkel) then
            factorb_s=factor_s/alpha
            factorb_g=factor_g/alpha
         else
            factorb_s=factor_s
            factorb_g=factor_g
         endif
c
c      Calculate variable intercept parameters if wanted
c
         if (in0s.eq.1) then  ! N_0s as in Thompson et al.
            temp_c = amin1(-0.001, tmk(i,j,k)-celkel)
            sonv = amin1(2.0e8, 2.0e6*exp(-0.12*temp_c))
         else
            sonv = rn0_s
         endif
c
         if (in0g.eq.1) then  ! N_0g as in Thompson et al.
            gonv = gon
            if (qgr(i,j,k).gt.r1) then
               gonv = 2.38*(pi*rho_g/
     +            (rhoair*qgr(i,j,k)))**0.92
               gonv = max(1.e4, min(gonv,gon))
            endif
         else
            gonv = rn0_g
         endif
c
         if (in0r.eq.1) then  ! N_0r as in Thompson et al.
            ronv = ron2
            if (qra(i,j,k).gt. r1) then
                ronv = ron_const1r*tanh((ron_qr0-qra(i,j,k))
     +             /ron_delqr0) + ron_const2r
            endif
c            if (i.eq.97.and.j.eq.100) then
c               print*,'k,ron2,qra(i,j,k),r1,ronv,rn0_r='
c               print*,'    ',k,ron2,qra(i,j,k),r1,ronv,rn0_r
c            endif
         else
            ronv = rn0_r
         endif
c
c      Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
c      the sum of z_e for each hydrometeor species:
c
         z_e =   factor_r  * (rhoair*qra(i,j,k))**1.75 /     ! rain
     &           ronv**.75
     &         + factorb_s * (rhoair*qsn(i,j,k))**1.75 /     ! snow
     &           sonv**.75
     &         + factorb_g * (rhoair*qgr(i,j,k))**1.75 /     ! graupel
     &           gonv**.75
c
c      Adjust small values of Z_e so that dBZ is no lower than -20
c
         z_e = max(z_e,.01)
c
c      Convert to dBZ
c
         dbz(i,j,k) = 10. * log10(z_e)
c
      enddo
      enddo
      enddo
c
      return
      end
