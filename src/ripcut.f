      program ripcut
c
c     This program reads in model output (in rip-format files) from a
c     source domain and creates a new file which has the data from the
c     source domain file cut to a subdomain.  The subdomain grid is
c     identical to the source domain grid in every way except it is a
c     subdomain.  A new case name is also given to the subdomain file,
c     which is the same as the original name, with "cut" added to the
c     casename.
c
      character argum(16)*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      character fname1*256,fname2*256
c
      dimension arr(1500,1500,60)
c
c   Get command line arguments.
c
      nargum=iargc()
      do i=1,nargum
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:7).eq.'ripcut_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if (nargum.ne.5) then
         print*,
     &      'Usage: ripcut source_file jx1 jx2 iy1 iy2'
         stop
      endif
c
c   Get information from header record of the source domain
c
      fname1=argum(1)
      print*,'fname1=',fname1
      open (unit=11,file=fname1,form='unformatted',status='old')
      read(11,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip
      goto 190
 170  print*,'Error in reading file .'
      stop
 180  print*,'Unexpected EOF encountered in reading file .'
      stop
 190  continue
      il=ihrip(4)
      jl=ihrip(5)
      kl=ihrip(9)
      ndim=ihrip(6)
      if (ndim.eq.3) then
         kldata=kl
      else
         kldata=1
      endif
      dskmc=rhrip(5)
      dskm=rhrip(6)
      yicorn=rhrip(7)
      xjcorn=rhrip(8)
c
c   Read the data array
c
      print*,'il,jl,kldata=',il,jl,kldata
      read(11) (((arr(i,j,k),i=1,il),j=1,jl),k=1,kldata)
c
c  Determine yicorn, xjcorn of target domain
c
      read(argum(2),'(i4)') jx1
      read(argum(3),'(i4)') jx2
      read(argum(4),'(i4)') iy1
      read(argum(5),'(i4)') iy2
      refrat=dskmc/dskm
      yicorn2=yicorn+(iy1-1.)/refrat
      xjcorn2=xjcorn+(jx1-1.)/refrat
c
c   Change RIP header values
c
      ihrip(4)=iy2-iy1+1
      ihrip(5)=jx2-jx1+1
      rhrip(7)=yicorn2
      rhrip(8)=xjcorn2
c
c   Create new file name
c      
      iendc=index(fname1,'_0')-1
      fname2=fname1(1:iendc)//'cut'//fname1(iendc+1:)
c
c   Write the new file
c
      open (unit=21,file=fname2,form='unformatted',status='unknown')
      write(21)
     &   vardesc,plchun,ihrip,rhrip,chrip
      write(21) (((arr(i,j,k),i=iy1,iy2),j=jx1,jx2),k=1,kldata)
c
      stop
      end
