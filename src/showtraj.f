      program showtraj
c
      parameter (maxtraj=7000,maxtrajtime=200)
c
      dimension stortr(maxtrajtime,maxtraj,3)
      character fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   Define constants.  These must be the same as in rip.
c
      rgas=287.04  !J/K/kg
      grav=9.81           ! m/s**2
      sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
c
      print*,'Enter the trajectory file name:'
      read(*,'(a256)') fname
      open (unit=10,file=fname,form='unformatted',status='old')
c      open (unit=11,file='trajprint',form='formatted',status='unknown')
      read(10) vardesc,plchun,ihrip,rhrip,chrip
      read (10) rtim,ctim,dttraj,ntraj
c     print*,'rtim,ctim,dttraj,ntraj=',rtim,ctim,dttraj,ntraj
      ntrajtime=nint(abs(rtim-ctim)/dttraj*3600) + 1
      if (rtim.lt.ctim) then
         trendtime=ctim
         trbegtime=rtim
         itm1=1
         itm2=ntrajtime
         itmi=1
      else
         trendtime=rtim
         trbegtime=ctim
         itm1=ntrajtime
         itm2=1
         itmi=-1
      endif
      do itm=itm1,itm2,itmi
         read(10) (stortr(itm,itr,1),itr=1,ntraj),
     &       (stortr(itm,itr,2),itr=1,ntraj),
     &       (stortr(itm,itr,3),itr=1,ntraj)
      enddo
      do itr=1,ntraj
         print*
         print*, 'Trajectory number ',itr,':'
         do itm=1,ntrajtime
            time=trbegtime+(itm-1)*dttraj/3600.
            xnest=1.+(stortr(itm,itr,2)-rhrip(8))*rhrip(5)/rhrip(6)
            ynest=1.+(stortr(itm,itr,1)-rhrip(7))*rhrip(5)/rhrip(6)
            write(6,'(a5,f7.3,a5,f7.3,a5,f7.3,a5,f8.5,a3)')'time=',time,
     &         '   x=',xnest,'   y=',ynest,
     &         '   z=',-.001*sclht*log(stortr(itm,itr,3)),' km'
         enddo
      enddo
      stop
      end
