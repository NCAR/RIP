      program ripshow
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c   This program reads in a rip data file and prints out the contents  c
c   of the header record.                                              c
c
      parameter (miymax=1500,mjxmax=3000,mkzhmax=100)
      dimension arr(miymax,mjxmax,mkzhmax),ii(5),jj(5),kk(5)
      character argum(16)*256,fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
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
      if (argum(1)(1:8).eq.'ripshow_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if (argum(1).eq.'-d') then
         if (nargum.ne.2) then
            print*,'Usage: ripshow [-d] filename'
            stop
         endif
         ishowd=1
         fname=argum(2)
      else
         if (nargum.ne.1) then
            print*,'Usage: ripshow [-d] filename'
            stop
         endif
         ishowd=0
         fname=argum(1)
      endif
      iendc=index(fname,' ')-1
      if (iendc.eq.-1) iendc=256
      open (unit=10,file=fname,form='unformatted',status='old')
c
      read(10,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip
      print*,'The variable and units are ',vardesc
      print*,'The units, in "PLOTCHAR" format, are ',plchun
      print*
      do i=1,32
         if (chrip(i)(1:5).eq.'mdate') then
             write(*,38) i,ihrip(i),chrip(i)
         elseif (chrip(i)(1:8).ne.'        '.or.
     &           ihrip(i).ne.999999999) then
             write(*,37) i,ihrip(i),chrip(i)
         endif
      enddo
      do i=1,32
         diff=rhrip(i)-9e9
         if (chrip(i+32)(1:8).ne.'        '.or.abs(diff).gt.10000.) then
             write(*,39) i,rhrip(i),chrip(i+32)
         endif
      enddo
 37   format('ihrip(',i2,')=',i9,5x,a)
 38   format('ihrip(',i2,')= ',i8.8,5x,a)
 39   format('rhrip(',i2,')=',f14.4,5x,a)
      miy=ihrip(4)
      mjx=ihrip(5)
      mkzh=ihrip(9)
      print*
c
c   Show data
c
      if (ishowd.eq.1) then
         if (miy.gt.miymax.or.mjx.gt.mjxmax.or.mkzh.gt.mkzhmax) then
            print*,'bump up parameter values to be'
            print*,'at least as large as:'
            print*,'miy,mjx,mkzh=',miy,mjx,mkzh
            stop
         endif
         if (ihrip(6).eq.3) then
            read(10,err=190,end=195)
     &         (((arr(i,j,k),i=1,miy),j=1,mjx),k=1,mkzh)
            kmax=5
         elseif (ihrip(6).eq.2) then
            read(10,err=190,end=195)
     &         ((arr(i,j,1),i=1,miy),j=1,mjx)
            kmax=1
         endif
         ii(1)=1
         ii(2)=(miy-ihrip(7))/4
         ii(3)=(miy-ihrip(7))/2
         ii(4)=(3*(miy-ihrip(7)))/4
         ii(5)=(miy-ihrip(7))
         jj(1)=1
         jj(2)=(mjx-ihrip(7))/4
         jj(3)=(mjx-ihrip(7))/2
         jj(4)=(3*(mjx-ihrip(7)))/4
         jj(5)=(mjx-ihrip(7))
         kk(1)=1
         kk(2)=1+mkzh/4
         kk(3)=1+mkzh/2
         kk(4)=1+(3*mkzh)/4
         kk(5)=mkzh
         do k=1,kmax
            if (ihrip(6).eq.3) then
               print*
               print*,'k=',kk(k)
               print*,'============================='
            elseif (ihrip(6).eq.2) then
               print*
            endif
            do i=5,1,-1
               write(*,'(5f16.7)')(arr(ii(i),jj(j),kk(k)),j=1,5)
            enddo
         enddo
      endif
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
