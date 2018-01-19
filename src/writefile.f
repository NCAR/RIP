c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine writefile (var,varname,numpas,ndim,icd,vardesc,
     &   plchun,ihrip,rhrip,chrip,casename,iendc,cxtimeavl,xtimeavl,
     &   nxt,ncxc,maxtavl,miy,mjx,mkzh)
c
      dimension var(miy,mjx,1+(mkzh-1)*(ndim-2))
      dimension xtimeavl(maxtavl)
      character fname*256,varname*10,cxtimeavl(maxtavl)*10,
     &   casename*(*),temp*4
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
      include 'comconst'
c
      iendv=index(varname,' ')-1
      if (iendv.eq.-1) iendv=10
      if (numpas.gt.0) then
         write(temp,'(i4)') 1000+numpas
         fname=casename(1:iendc)//'_'//
     &      cxtimeavl(nxt)(1:ncxc)//'_'//varname(1:iendv)//temp(2:4)
      else
         fname=casename(1:iendc)//'_'//
     &      cxtimeavl(nxt)(1:ncxc)//'_'//varname(1:iendv)
      endif
      open(unit=25,file=fname,form='unformatted',status='unknown')
      ihrip(6)=ndim ! number of dimensions of this variable (2 or 3)
      ihrip(7)=icd  ! grid of this var. (1:cross point, 0:dot point)
      write(25) vardesc,plchun,ihrip,rhrip,chrip
      write(25) var
      close (25)
      return
      end
