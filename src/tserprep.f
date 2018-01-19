c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine tserprep(tseryi,tserxj,
     &   tserlocpoint,tserlocdesc,ntsers,ntsert,maxtsers,rip_root,
     &   miy,mjx)
c
c   This subroutine reads the time series station locations and
c      calculates the corresponding x/y locations on the grid.
c
      dimension tseryi(maxtsers),tserxj(maxtsers)
      character tserlocpoint(maxtsers)*58,
     &   tserlocdesc(maxtsers)*44,rip_root*256
c
      character tserlocsp(2)*20,string*48,icaoid*4,locdesc*44
c
      include 'comconst'
c
      i=0
      open (unit=iutserstn,file='tserstn.dat',form='formatted',
     &   status='old')
  300 read(iutserstn,'(a48)',end=350) string
      i=i+1
      if (i .gt. maxtsers) then
        write(iup,*) 'There are more than ',maxtsers,' stations in your
     &tserstn.dat file.'
        write(iup,*) 'Either reduce the number of stations or increase m
     &axtsers in driver.f and recompile.'
        stop '\tstopping in tserprep'
      endif
      ipos=1
      nterm=0
   37 nterm=nterm+1
      if (nterm.gt.2) then
         write(iup,*)'In tserprep: Too many terms for tserlocsp.'
         write(iup,*)'Location # = ',i
         stop
      endif
      call getchar(string,ipos,tserlocsp(nterm),0)
      if (string(ipos-1:ipos-1).eq.',') goto 37
      if (nterm.eq.1) tserlocsp(2)='missing             '
      call locinterp(tserlocsp,gridx,gridy,
     &   rlat,rlon,iwmo,icaoid,stelev,locdesc,rip_root)
      tserlocdesc(i)=locdesc
      tserlocpoint(i)=' '
      write(tserlocpoint(i),213)'x,y=      ,         lat,lon=',
     &   rlat,',',rlon
      itdone=0
      if (icaoid.ne.'XXXX') then
         write(tserlocpoint(i)(43:),214)'  stn=',icaoid
         itdone=1
      endif
      if (iwmo.ne.99999) then
         if (itdone.eq.0) then
            write(tserlocpoint(i)(43:),215)'   stn=',iwmo
         else
            write(tserlocpoint(i)(53:),216)',',iwmo
         endif
      endif
      sxgn=1.+(gridx-xjcorn)*refrat
      sygn=1.+(gridy-yicorn)*refrat
      if (sxgn.le..5.or.sxgn.ge.mjx-.5.or.
     &    sygn.le..5.or.sygn.ge.miy-.5) then
         write(iup,*)'The time series requested at the following'
         write(iup,*)'location is outside the cross-point domain:'
         write(iup,*) tserlocpoint(i)
         write(iup,*)'Please correct tserstn.dat and re-execute RIP.'
         stop
      endif
      tserxj(i)=sxgn
      tseryi(i)=sygn
      write(tserlocpoint(i)(5:10),'(f6.2)') sxgn
      write(tserlocpoint(i)(12:17),'(f6.2)') sygn
c      
  213 format(a28,f6.2,a1,f7.2)
  214 format(a6,a4)
  215 format(a7,i5)
  216 format(a1,i5)
c
      goto 300
c
  350 ntsers=i
      ntsert=0
      close (iutserstn)
      return
      end
