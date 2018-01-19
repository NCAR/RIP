c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine prcver(prmod,probs,mdate,miy,mjx)
c
c   This routine calculates the following precip verifications:
c      1) Categorical forecast score
c      2) Threat score
c      3) Bias score
c    All of these scores are calculated at the threshholds designated
c       in the data statement (in mm).
c
      parameter (nthresh=6)
c
      dimension prmod(miy,mjx),probs(miy,mjx)
c
      dimension thresh(nthresh),arobs(nthresh),armod(nthresh),
     &   arboth(nthresh),arneith(nthresh)
c
      include 'comconst'
c
      data thresh / 0.1000, 0.2154, 0.4642,
     &               1.000,  2.154,  4.642 /
c
      iu=iuprcver
      artot=0.
      do 100 kt=1,nthresh
         arobs(kt)=0.
         armod(kt)=0.
         arboth(kt)=0.
         arneith(kt)=0.
  100 continue
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
         if (probs(i,j).eq.rmsg) goto 1000
         artot=artot+1.
         do 900 kt = 1,nthresh
            if (probs(i,j).ge.thresh(kt).and.
     &          prmod(i,j).ge.thresh(kt)) arboth(kt)=arboth(kt)+1.
            if (probs(i,j).ge.thresh(kt)) arobs(kt)=arobs(kt)+1.
            if (prmod(i,j).ge.thresh(kt)) armod(kt)=armod(kt)+1.
            if (probs(i,j).lt.thresh(kt).and.
     &          prmod(i,j).lt.thresh(kt)) arneith(kt)=arneith(kt)+1.
  900    continue
 1000 continue
c
c   Calculate and print out the scores to the speicified unit number.
c
      write(iu,*)
      write(iu,*) '================================================='
      write(iu,*)
      write(iu,*) 'Precipitation verification valid at ',mdate,'.'
      write(iu,*) '    # of grid points =         ',(miy-1)*(mjx-1)
      write(iu,*) '    # of verification points = ',nint(artot)
      write(iu,*)
      write(iu,*) '  Thresh(mm)  CFS(%)  Threat(%)  Bias(%)'
      write(iu,*) '  ----------  ------  ---------  -------'
      do 2000 kt=1,nthresh
         cfs=999.
         threat=999.
         bias=999.
         denom=arobs(kt)+armod(kt)-arboth(kt)
         if (artot.gt..5) cfs=100.*(arboth(kt)+arneith(kt))/artot
         if (denom.gt..5) threat=100.*arboth(kt)/denom
         if (arobs(kt).gt..5) bias=100.*armod(kt)/arobs(kt)
         write(iu,'(3x,f10.4,2x,f5.1,2(5x,f5.1))')
     &      thresh(kt),cfs,threat,bias
 2000 continue
c
      return
      end
