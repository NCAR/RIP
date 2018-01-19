c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine mconvert(mdate,mhour,idir,nsplityear)
c
c   mdate: an 8-digit integer specification for a date,
c      given as yymmddhh
c   mhour: an integer specificying the number of hours since
c      00 UTC 1 January 1 AD.
c
c   This routine converts an mdate to an mhour if idir=1, or vice versa
c   if idir=-1.
c
c   If idir=1, how do we know what century mdate refers to?  You
c   provide a year, called "nsplityear", denoted "aabb".  If mdate
c   is denoted "yymmddhh", then if yy >or= bb, the century is
c   assumed to be aa.  Otherwise it is assumed to be the century
c   after aa, or aa+1.
c
c   Leap year definition: every fourth year has a 29th day in February,
c      with the exception of century years not divisible by 400.
c
      dimension ndaypmo(12)
      integer yy,mm,dd,hh,aa,bb
      data ndaypmo /31,28,31,30,31,30,31,31,30,31,30,31/
c
      if (idir.eq.1) then
c
      yy=mdate/1000000
      bb=mod(nsplityear,100)
      aa=nsplityear-bb
      iyear=aa+yy
      if (yy.lt.bb) iyear=iyear+100
      iyearp=iyear-1
      idayp = iyearp*365 + iyearp/4 - iyearp/100 +iyearp/400
      mm=mod(mdate,1000000)/10000
      imonthp=mm-1
      if ((mod(iyear,4).eq.0.and.mod(iyear,100).ne.0).or.
     &    mod(iyear,400).eq.0)
     &   ndaypmo(2)=29
      do 5 i=1,imonthp
         idayp=idayp+ndaypmo(i)
 5    continue
      ndaypmo(2)=28
      dd=mod(mdate,10000)/100
      idayp=idayp+dd-1
      hh=mod(mdate,100)
      mhour=24*idayp+hh
c
      else
c
      nhour=mhour
c
c   Get an estimate of iyear that is guaranteed to be close to but
c   less than the current year
c
      iyear = max(0,nhour-48)*1.14079e-4
      ihour=24*(iyear*365+iyear/4-iyear/100+iyear/400)
 10   iyear=iyear+1
      ihourp=ihour
      ihour = 24*(iyear*365 + iyear/4 - iyear/100 +iyear/400)
      if (ihour.le.nhour) goto 10
      nhour=nhour-ihourp
      if ((mod(iyear,4).eq.0.and.mod(iyear,100).ne.0).or.
     &    mod(iyear,400).eq.0)
     &   ndaypmo(2)=29
      imo=0
      ihour=0
 20   imo=imo+1
      ihourp=ihour
      ihour=ihour+24*ndaypmo(imo)
      if (ihour.le.nhour) goto 20
      nhour=nhour-ihourp
      ndaypmo(2)=28
      iday = nhour/24 + 1
      ihour=mod(nhour,24)
      mdate=mod(iyear,100)*1000000+imo*10000+iday*100+ihour
c
      endif
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine mconvertccyy(mdate,nsplityear,iyear)
c
c   mdate: an 8-digit integer specification for a date,
c      given as yymmddhh
c   mhour: an integer specificying the number of hours since
c      00 UTC 1 January 1 AD.
c
c   How do we know what century mdate refers to?  You
c   provide a year, called "nsplityear", denoted "aabb".  If mdate
c   is denoted "yymmddhh", then if yy >or= bb, the century is
c   assumed to be aa.  Otherwise it is assumed to be the century
c   after aa, or aa+1.
c
      integer yy,aa,bb
c
      yy=mdate/1000000
      bb=mod(nsplityear,100)
      aa=nsplityear-bb
      iyear=aa+yy
      if (yy.lt.bb) iyear=iyear+100

      end
