c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine gettimes(casename,iendc,xtimeavl,cxtimeavl,nxtavl,
     &   maxtavl,ncxc,iup)
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*(*),fname*256
c
c   Read in the available xtimes, in both character and
c   floating point arrays
c
      fname=casename(1:iendc)//'.xtimes'
      open(unit=25,file=fname,form='formatted',status='old')
      read(25,*)nxtavl
      if (nxtavl.gt.maxtavl) then
         write(iup,*)'There are ',nxtavl,' times in the ".xtime"'
         write(iup,*)'file, but maxtavl is only ',maxtavl,'.'
         write(iup,*)'Increase maxtavl in rip.f to be at least'
         write(iup,*)'as big as ',nxtavl,' and then recompile'
         write(iup,*)'and re-run rip.'
         stop
      endif
      do i=1,nxtavl
         read(25,'(a10)') cxtimeavl(i)
         if (i.eq.1) then
            ncxc=10
            if (cxtimeavl(i)(10:10).eq.' ') ncxc=9
         endif
         if (ncxc.eq.9) then
            read(cxtimeavl(i),'(f9.5)') xtimeavl(i)
         else
            read(cxtimeavl(i),'(f10.5)') xtimeavl(i)
         endif
      enddo
      close (25)
      do i=nxtavl+1,maxtavl
         cxtimeavl(i)=' '
         xtimeavl(i)=0.
      enddo
c
      return
      end
