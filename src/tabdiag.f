      program tabdiag
c
      parameter (maxtraj=7000,maxtrajtime=200,maxvar=20)
c
      dimension diag(maxtrajtime,maxtraj,maxvar),time(maxtrajtime)
      character colhead*128,fmt*128,sepline*128
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      character argum(16)*256,outfile*256
c
c   Get command line arguments.
c
      nargum=iargc()
      if (nargum.ne.2) then
         print*,'Usage: tabdiag diagnostic_file input_file'
         stop
      endif
      do i=1,nargum
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:8).eq.'tabdiag_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      open (unit=9,file=argum(2),form='formatted',status='old')
      read(9,*)colhead
      read(9,*)fmt
      iendcolh=lennonblank(colhead)
      do i=1,128
         sepline(i:i)='='
      enddo
c
      open (unit=10,file=argum(1),form='unformatted',status='old')
c
      read(10) vardesc,plchun,ihrip,rhrip,chrip
      read(10) rtim,ctim,dtfile,ntraj,nvar
      ndtm=nint(abs(rtim-ctim)*3600./dtfile)+1
      if (rtim.lt.ctim) then
         diagendtime=ctim
         diagbegtime=rtim
         idtm1=1
         idtm2=ndtm
         idtmi=1
      else
         diagendtime=rtim
         diagbegtime=ctim
         idtm1=ndtm
         idtm2=1
         idtmi=-1
      endif
c
      do idtm=1,ndtm
         time(idtm)=diagbegtime+(idtm-1)*dtfile/3600.
      enddo
c
c   Get the diagnostics
c      
      do idtm=idtm1,idtm2,idtmi
         do ivar=1,nvar
            read (10) (diag(idtm,itr,ivar),itr=1,ntraj)
         enddo
      enddo
c
      idd=index(argum(1),'.diag')
      if (idd.eq.0) then
         print*,'your diagnostic file should have a .diag suffix'
         stop
      else
         outfile=argum(1)(1:idd)//'tabdiag'
         open(unit=11,file=outfile,form='formatted',status='unknown')
      endif
c
      write(11,*)
      write(11,*) 'Trajectory diagnostics from file ',
     &   argum(1)(1:index(argum(1),' ')-1)
      write(11,*)
      write(11,*)
      do itr=1,ntraj
         write(11,*) 'Trajectory # ',itr,':'
         write(11,*)
         write(11,'(a)') colhead(1:iendcolh)
	 write(11,'(a)') sepline(1:iendcolh)
         do idtm=1,ndtm
            write(11,fmt) time(idtm),(diag(idtm,itr,ivar),ivar=1,nvar)
         enddo
         write(11,*)
      enddo
      stop
      end
