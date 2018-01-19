c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine smooth(pslab,work,numpas,mabpl,njx,niy)
c
c   This is a smoothing routine, with several choices:
c
c   If numpas is between 1 and 99, then a 9-point weighted smoother is
c   applied numpas times.  The smoother follows equation 11-107 in
c   Haltiner and Williams. One pass completely removes 2-delta-x waves
c   on the interior.  On the outer row and column, and near missing
c   data points, smoothing is carried out in a manner that preserves
c   the domain average value of the field.
c
c   If numpas is between 101 and 199, then a smoother-desmoother is
c   applied (numpas-100) times.  One pass removes a large fraction
c   of the 2-delta-x component, but is not as harsh on longer
c   wavelengths as the 9-point smoother
c
c   If numpas is between 201 and 299, then the smoother-desmoother is
c   applied (numpas-200) times, and, after each pass, the data field
c   is forced to be non-negative.
c
c   If numpas is between 301 and 399, then a weighted
c   smoother is applied, in which the smoothed value
c   is given by a weighted average of values at
c   surrounding grid points.  The weighting function
c   is the Cressman weighting function:
c
c               w = ( D**2 - d**2 ) / ( D**2 + d**2 )
c
c   In the above, d is the distance (in grid increments)
c   of the neighboring point to the smoothing point, and
c   D is the radius of influence [in grid increments,
c   given by (numpas-300)].
c
c   If numpas is between 401 and 499, then the smoothing
c   is similar for numpas=301-399, except the weighting
c   function is the circular apperture diffraction function
c   (following a suggestion of Barnes et al. 1996):
c
c               w = bessel(3.8317*d/D)/(3.8317*d/D)
c
c   If numpas is between 501 and 599, then the smoothing
c   is similar for numpas=301-399, except the weighting
c   function is the product of the rectangular
c   apperture diffraction function in the x and y directions
c   (the function used in Barnes et al. 1996):
c
c               w = [sin(pi*x/D)/(pi*x/D)]*[sin(pi*y/D)/(pi*y/D)]
c
c   Note, the first index of pslab varies along the abcissa
c   (or x), and the second index varies along the ordinate (or y).
c
      parameter(beszero=3.8317)
c
      dimension pslab(mabpl,niy),work(mabpl,niy),xnu(2),fprint(150,150)
c
      include 'comconst'
c
      if (mod(numpas,100).eq.0) return
c
      if (numpas.le.99) then   ! 9-point smoother
c
      do ipas=1,numpas
c
      do i=1,niy
      do j=1,njx
         work(j,i)=0.
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            work(j,i)=rmsg
         else
            totgive=0.
            if (i.gt.1) then
               if (pslab(j,i-1).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j,i-1)=work(j,i-1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i-1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j-1,i-1)=work(j-1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i-1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j+1,i-1)=work(j+1,i-1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (i.lt.niy) then
               if (pslab(j,i+1).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j,i+1)=work(j,i+1)+give
                  totgive=totgive+give
               endif
               if (j.gt.1) then
                  if (pslab(j-1,i+1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j-1,i+1)=work(j-1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
               if (j.lt.njx) then
                  if (pslab(j+1,i+1).ne.rmsg) then 
                     give=.0625*pslab(j,i)
                     work(j+1,i+1)=work(j+1,i+1)+give
                     totgive=totgive+give
                  endif
               endif
            endif
            if (j.gt.1) then
               if (pslab(j-1,i).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j-1,i)=work(j-1,i)+give
                  totgive=totgive+give
               endif
            endif
            if (j.lt.njx) then
               if (pslab(j+1,i).ne.rmsg) then 
                  give=.125*pslab(j,i)
                  work(j+1,i)=work(j+1,i)+give
                  totgive=totgive+give
               endif
            endif
            work(j,i)=work(j,i)+pslab(j,i)-totgive
         endif
      enddo
      enddo
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
c
      enddo
c
      elseif (numpas.le.299) then   ! smoother-desmoother
c
      if (numpas.ge.200) then
         nump=numpas-200
         inn=1
      else
         nump=numpas-100
         inn=0
      endif
c
      if (nump.lt.1) return
c
c   Check if any data is missing.
c
      imsg=0
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).eq.rmsg) then
            imsg=1
            goto 15
         endif
      enddo
      enddo
 15   continue
c
      if (imsg.eq.1) then
c
c   Get average value of pslab.
c
      nval=0
      total=0.
      do 10 i=1,niy
      do 10 j=1,njx
         if (pslab(j,i).ne.rmsg) then
            total=total+pslab(j,i)
            nval=nval+1
         endif
   10 continue
      if (nval.eq.0) then
         write(*,*)'  All elements of this pslab are rmsg.'
         return
      endif
      avgval=total/nval
c
c   Set each element that is currently rmsg to avgval, and
c   keep track of them with work.
c
      do 20 i=1,niy
      do 20 j=1,njx
         if (pslab(j,i).eq.rmsg) then
            pslab(j,i)=avgval
            work(j,i)=1.
         else
            work(j,i)=0.
         endif
   20 continue
c
      endif
c
c     *** Do calculation and put into pslab array.
c
      xnu(1) = 0.50
      xnu(2) = -0.52
      je = njx - 1
      ie = niy - 1
      do 100 ipass = 1,nump*2
         kp=2-mod(ipass,2)          
c     
c        *** First, smooth in the njx direction.
c     
         do 60 j = 2,je
            asv = pslab(j,1)
            do 50 i = 2,ie
               aplus = pslab(j,i+1)
               cell = pslab(j,i)
               pslab(j,i)= pslab(j,i) + xnu(kp)*
     +            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   50       continue
   60    continue
c     
c        *** Now, smooth in the niy direction.
c     
         do 80 i = 2,ie
            asv = pslab(1,i)
            do 70 j = 2,je
               aplus = pslab(j+1,i)
               cell = pslab(j,i)
               pslab(j,i) = pslab(j,i) + xnu(kp)*
     +            ((asv + aplus)/2.0 - pslab(j,i))
               asv = cell
   70       continue
   80    continue
c
      if (inn.eq.1) then
c
c      Make non-negative.
c
         do i=1,niy
         do j=1,njx
            pslab(j,i)=max(0.,pslab(j,i))
         enddo
         enddo
      endif
c
  100 continue
c
      if (imsg.eq.1) then
c
c      Set rmsg elements back to rmsg
c
         do 200 i=1,niy
         do 200 j=1,njx
            pslab(j,i)=work(j,i)*rmsg + (1.-work(j,i))*pslab(j,i)
  200    continue
      endif
c
      elseif (numpas.le.599) then   ! weighted smoother
c
      idist=mod(numpas,100)
      if (idist.eq.0) return
      nfp=1+2*idist
      npsq=idist*idist
      if (numpas.le.399) then  ! Cressman function
         do i=1,nfp
         do j=1,nfp
            distsq=(i-idist-1.)**2+(j-idist-1.)**2
            fprint(j,i)=max((npsq-distsq)/(npsq+distsq),0.0)
         enddo
         enddo
      elseif (numpas.le.499) then   ! Circular diffraction function
         do i=1,nfp
         do j=1,nfp
            dist=beszero/idist*sqrt((i-idist-1.)**2+(j-idist-1.)**2)
            if (i.eq.idist+1.and.j.eq.idist+1) then
               fprint(j,i)=.5
            else
               fprint(j,i)=max(0.,bes(dist)/dist)
            endif
         enddo
         enddo
      elseif (numpas.le.599) then   ! Rect. diffraction function
         do i=1,nfp
         do j=1,nfp
            if (j.eq.idist+1) then
               xfac=1.
            else
               xdist=pi/idist*(j-idist-1.)
               xfac=sin(xdist)/xdist
            endif
            if (i.eq.idist+1) then
               yfac=1.
            else
               ydist=pi/idist*(i-idist-1.)
               yfac=sin(ydist)/ydist
            endif
            fprint(j,i)=xfac*yfac
         enddo
         enddo
      endif
c
      do i=1,niy
      do j=1,njx
         if (pslab(j,i).ne.rmsg) then
            tot=0.
            totwt=0.
            is=max(1,i-idist)
            ie=min(niy,i+idist)
            js=max(1,j-idist)
            je=min(njx,j+idist)
            do ireg=is,ie
               ifp=ireg-i+idist+1
            do jreg=js,je
               jfp=jreg-j+idist+1
               if (pslab(jreg,ireg).ne.rmsg) then
                  totwt=totwt+fprint(jfp,ifp)
                  tot=tot+fprint(jfp,ifp)*pslab(jreg,ireg)
               endif
            enddo
            enddo
            work(j,i)=tot/totwt
         else
            work(j,i)=rmsg
         endif
      enddo
      enddo
c
      do i=1,niy
      do j=1,njx
         pslab(j,i)=work(j,i)
      enddo
      enddo
c
      endif
c
      return
      end
