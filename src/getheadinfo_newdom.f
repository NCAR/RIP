c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getheadinfo_newdom(nproj,miycors,mjxcors,mdateb,
     &   mhourb,iice,iplevdata,true1,true2,xlatc,xlonc,dskmc,dskm,
     &   yicorn,xjcorn,rhourb,dsc,ds,refrat,iunewdom,iup)
c
      rewind (iunewdom)
      read(iunewdom,*)
      read(iunewdom,*)
      read(iunewdom,*) nproj
      read(iunewdom,*) miycors,mjxcors
      read(iunewdom,*) mdateb,rhourb
      read(iunewdom,*) iice,iplevdata,ilandset
      read(iunewdom,*) true1,true2
      read(iunewdom,*) xlatc,xlonc
      read(iunewdom,*) dskmc,dskm
      read(iunewdom,*) yicorn,xjcorn
c
      call mconvert(mdateb,mhourb,1,1940)
      if (iplevdata.lt.1.or.iplevdata.gt.50) iplevdata=6
c
c   Set up map transformation stuff
c
      dsc=dskmc*1000.
      ds=dskm*1000.
      refrat=dsc/ds
c
      return
      end
