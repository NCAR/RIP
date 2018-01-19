c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine setripheader(ihrip,rhrip,chrip,miy,mjx,mkzh)
c
      include 'comconst'
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64
c
c   Set up stuff for RIP format data files.
c
      do i=1,32
         ihrip(i)=999999999
         rhrip(i)=9e9
         chrip(i)=' '
         chrip(i+32)=' '
      enddo
      chrip(1)=
     &'map projection (0: none/ideal, 1: LC, 2: PS, 3: ME, 4: SRCE)'
      ihrip(1)=nproj
      chrip(2)=
     &'number of dot points in the y-direction (coarse domain)'
      ihrip(2)=miycors
      chrip(3)=
     &'number of dot points in the x-direction (coarse domain)'
      ihrip(3)=mjxcors
      chrip(4)=
     &'number of dot points in the y-direction (this domain)'
      ihrip(4)=miy
      chrip(5)=
     &'number of dot points in the x-direction (this domain)'
      ihrip(5)=mjx
      chrip(6)=
     &'number of dimensions of this variable (2 or 3)'
      ihrip(6)=999    ! this is set separately for each variable
      chrip(7)=
     &'grid of this variable (1: cross point dom., 0: dot point dom.)'
      ihrip(7)=999    ! this is set separately for each variable
ccc      chrip(8)=
ccc     &'vertical coordinate (0: hydrostatic sigma, 1: nonhyd. sigma)'
ccc      ihrip(8)=inhyd
ccc      chrip(9)=
ccc     &'number of half sigma levels'
      chrip(9)=
     &'number of vertical levels in the data'
      ihrip(9)=mkzh
      chrip(10)=
     &'mdateb: YYMMDDHH (truncated hour) of hour-0 for this dataset'
      ihrip(10)=mdateb
      chrip(11)=
     &'mdate: YYMMDDHH (truncated hour) of this time'
      ihrip(11)=99999999    ! this is set separately for each time
      chrip(12)=
     &'ice physics (1: sep. arrays for ice fields, 0: no sep. arrays)'
      ihrip(12)=iice
ccc      chrip(13)=
ccc     &'Program #: 1:TER. 2:DG/RG. 3:RAW. 5:INT. 6:MOD. 11:MOD.(MM5V3)'
      chrip(13)=
     &'ver. coord. type: <or=3: hgt. or prs.; >or=4: terrain-following'
      ihrip(13)=iplevdata
      chrip(14)=
     &'landuse dataset (1: old, 13-cat; 2: USGS, 24-cat; 3: SiB, 16 )'
      ihrip(14)=ilandset
c
      ijmp=32
      chrip(ijmp+1)=
     &'first true latitude (deg.)'
      rhrip(1)=true1
      chrip(ijmp+2)=
     &'second true latitude (deg.)'
      rhrip(2)=true2
      chrip(ijmp+3)=
     &'central latitude of coarse domain (deg.)'
      rhrip(3)=xlatc
      chrip(ijmp+4)=
     &'central longitude of coarse domain(deg.)'
      rhrip(4)=xlonc
      chrip(ijmp+5)=
     &'grid distance of coarse domain (km)'
      rhrip(5)=dskmc
      chrip(ijmp+6)=
     &'grid distance of this domain (km)'
      rhrip(6)=dskm
      chrip(ijmp+7)=
     &'coarse dom. y-position of lower left corner of this domain'
      rhrip(7)=yicorn
      chrip(ijmp+8)=
     &'coarse dom. x-position of lower left corner of this domain'
      rhrip(8)=xjcorn
ccc      chrip(ijmp+9)=
ccc     &'pressure level (hPa) of the model top'
ccc      rhrip(9)=ptop
ccc      chrip(ijmp+10)=
ccc     &'reference sea-level pressure (Pa)'
ccc      rhrip(10)=refslp
ccc      chrip(ijmp+11)=
ccc     &'reference sea-level temperature (K)'
ccc      rhrip(11)=refslt
ccc      chrip(ijmp+12)=
ccc     &'reference lapse rate (dT/d(ln(p)), K)'
ccc      rhrip(12)=reflaps
      chrip(ijmp+13)=
     &'rhourb: diff (in h) between exact time and mdate of hour-0'
      rhrip(13)=rhourb
      chrip(ijmp+14)=
     &'rhour: diff (in h) between exact time and mdate of this data'
      rhrip(14)=9e9       ! this is set separately for each time
      chrip(ijmp+15)=
     &'xtime: exact time of this data relative to exact hour-0 (in h)'

      rhrip(15)=9e9       ! this is set separately for each time
ccc      chrip(ijmp+16)=
ccc     &'reference stratospheric constant temperature (K)'
ccc      rhrip(16)=refstratt
c
      return
      end
