      program ripcomp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c   This program reads in two rip data files and compares their        c
c   contents.                                                          c
c
      parameter (miymax=200,mjxmax=200,mkzhmax=35)
      dimension arr(miymax,mjxmax,mkzhmax),arr2(miymax,mjxmax,mkzhmax)
      character argum(16)*256,fname*256,fname2*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
      dimension ihrip2(32),rhrip2(32)
      character chrip2(64)*64,vardesc2*64,plchun2*24
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
      if (argum(1)(1:8).eq.'ripcomp_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if (nargum.ne.2) then
         print*,'Usage: ripcomp filename1 filename2'
         stop
      endif
      ishowd=0
      fname=argum(1)
      fname2=argum(2)
      iendc=index(fname,' ')-1
      iendc2=index(fname2,' ')-1
      if (iendc.eq.-1) iendc=256
      if (iendc2.eq.-1) iendc2=256
      open (unit=10,file=fname,form='unformatted',status='old')
      open (unit=11,file=fname2,form='unformatted',status='old')
c
      read(10,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip
      read(11,err=170,end=180)
     &   vardesc2,plchun2,ihrip2,rhrip2,chrip2
c      dimension ihrip(32),rhrip(32)
c      character chrip(64)*64,vardesc*64,plchun*24
      if (vardesc2.ne.vardesc) then
         print*,'vardesc not the same.'
         print*,'vardesc (1,2)=',vardesc,vardesc2
      endif
      if (plchun2.ne.plchun) then
         print*,'plchun not the same.'
         print*,'plchun (1,2)=',plchun,plchun2
      endif
      do i=1,32
         if (ihrip2(i).ne.ihrip(i)) then
            print*,'ihrip(',i,') not the same.'
            print*,'ihrip(',i,') (1,2)=',ihrip(i),ihrip2(i)
         endif
      enddo
      do i=1,32
         if (rhrip2(i).ne.rhrip(i)) then
            print*,'rhrip(',i,') not the same.'
            print*,'rhrip(',i,') (1,2)=',rhrip(i),rhrip2(i)
         endif
      enddo
      do i=1,64
         if (chrip2(i).ne.chrip(i)) then
            print*,'chrip(',i,') not the same.'
            print*,'chrip(',i,') (1,2)=',chrip(i),chrip2(i)
         endif
      enddo
      miy=ihrip(4)
      mjx=ihrip(5)
      mkzh=ihrip(9)
         if (ihrip(6).eq.3.and.ihrip2(6).eq.3) then
            read(10,err=190,end=195)
     &         (((arr(i,j,k),i=1,miy),j=1,mjx),k=1,mkzh)
            read(11,err=190,end=195)
     &         (((arr2(i,j,k),i=1,miy),j=1,mjx),k=1,mkzh)
            kmax=mkzh
         elseif (ihrip(6).eq.2.and.ihrip2(6).eq.2) then
            read(10,err=190,end=195)
     &         ((arr(i,j,1),i=1,miy),j=1,mjx)
            read(11,err=190,end=195)
     &         ((arr2(i,j,1),i=1,miy),j=1,mjx)
            kmax=1
         else
            print*,'files are not same dimension.'
            stop
         endif
      icd=ihrip(7)
      diffabmax=0.
      diffpmax=0.
      do k=1,kmax
      do j=1,mjx-icd
      do i=1,miy-icd
         if (arr(i,j,k).ne.0..or.arr2(i,j,k).ne.0.) then
            diffab=abs(arr2(i,j,k)-arr(i,j,k))
            diffp=200.*diffab/(arr2(i,j,k)+arr(i,j,k))
         else
            diffab=0.
            diffp=0.
         endif
         diffabmax=max(diffab,diffabmax)
         diffpmax=max(diffp,diffpmax)
      enddo
      enddo
      enddo
      print*,'max absolute difference=',diffabmax
      print*,'max percent difference=',diffpmax
      stop
c
 170  print*,'The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  print*,'Unexpected EOF reached when trying to read'
      print*,'model data header.  Stopping.'
      stop
c
 190  print*,'Error in reading the model data array. Stopping.'
      stop
c
 195  print*,'Unexpected EOF encountered while'
      print*,'reading the model data array. Stopping.'
      stop
      end
