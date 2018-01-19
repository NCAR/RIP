c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt,maxtavl,ncxc,nproj,miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn,xjcorn,rhourb,dsc,ds,refrat,iup)
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*(*),fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   Create a file name for a file that is likely to exist,
c   just so we can read it and get the header information from the dataset.
c
      fname=casename(1:iendc)//'_'//cxtimeavl(nxt)(1:ncxc)//'_'//'ter'
      open(unit=25,file=fname,form='unformatted',status='old')
      read(25,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip
c
      nproj=ihrip(1)
      miycors=ihrip(2)
      mjxcors=ihrip(3)
      icdch=ihrip(7)
c      inhyd=ihrip(8)
      mdateb=ihrip(10)
      call mconvert(mdateb,mhourb,1,1940)
      iice=ihrip(12)
c      iprog=ihrip(13)
c      if (iprog.lt.1.or.iprog.gt.50) iprog=6
c
c   iplevdata should be 3 for p-lev (though anything <4 will work),
c      and 6 for terrain-following (though anything >4 will work).
c
      iplevdata=ihrip(13)
      ilandset=ihrip(14)
      if (iplevdata.lt.1.or.iplevdata.gt.50) iplevdata=6
      true1=rhrip(1)
      true2=rhrip(2)
      xlatc=rhrip(3)
      xlonc=rhrip(4)
      dskmc=rhrip(5)
      dskm=rhrip(6)
      yicorn=rhrip(7)
      xjcorn=rhrip(8)
      if (icdch.ge.2) then
c
c      For E-grid data, yicorn and xjcorn were set to the location
c      of the lower-left H point in the E-grid domain, but we must
c      shift it down and left by half a nest grid space to give the location
c      of the lower-left dot point of the B-grid whose cross points
c      overlap the E-grid points.
c
         yicorn=yicorn-.5*(dskm/dskmc)
         xjcorn=xjcorn-.5*(dskm/dskmc)
      endif
c      refslp=rhrip(10)
c      refslt=rhrip(11)
c      reflaps=rhrip(12)
c      refstratt=rhrip(16)
c      if (refstratt.gt.500.or.refstratt.lt.0.) refstratt=0.1
      rhourb=rhrip(13)
      dsc=dskmc*1000.
      ds=dskm*1000.
      refrat=dsc/ds
c
      close (25)
      goto 200
c
 170  write(iup,*) 'The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  write(iup,*) 'Unexpected EOF reached when trying to read'
      write(iup,*) 'model data header.  Stopping.'
      stop
c
 200  continue
c
      return
      end
