      subroutine mktrjpts(xjtraj,yitraj,zktraj,maxtraj,rmsg)
c
      dimension xjtraj(maxtraj),yitraj(maxtraj),zktraj(maxtraj)
c
      dimension zlev(100)
c
c   This routine makes a 3D rectangular grid of trajectory initial
c   points given: the x and y values of the lower left corner of
c   the traj. grid; the x and y values of another point defining the
c   positive x-axis of the traj. grid; the traj. grid
c   spacing; and the number of points in the traj. grid x
c   and y directions.  These values should be the first 7 values
c   in the xjtraj array.  The first value should be negative, indicating
c   that a grid is to be defined (rather than just points), but the
c   absolute value of that value will be used.  Any yitraj values given
c   are ignored.  The zktraj values specify the vertical levels of the
c   3D grid to be defined.
c
      xcorn=-xjtraj(1)
      ycorn=xjtraj(2)
      xpx=xjtraj(3)
      ypx=xjtraj(4)
      gsp=xjtraj(5)
      nx=xjtraj(6)
      ny=xjtraj(7)
      do i=1,maxtraj
         if (zktraj(i).ne.rmsg) then
            zlev(i)=zktraj(i)
         else
            nzlev=i-1
            goto 41
         endif
      enddo
 41   continue
c
      write(6,*) 'Rectangular grid of trajectories defined as follows:'
      write(6,*) 'xcorn = ',xcorn,' ycorn = ',ycorn
      write(6,*) 'xpx = ',xpx,' ypx = ',ypx
      write(6,*) 'gsp = ',gsp,' nx = ',nx,' ny = ',ny

      dx=xpx-xcorn
      dy=ypx-ycorn
      hypot=sqrt(dx*dx+dy*dy)
      cosang=dx/hypot
      sinang=dy/hypot
c
      np=nx*ny
      ip=0
      do i=1,ny
      do j=1,nx
         ip=ip+1
         xjtraj(ip)=xcorn+gsp*((j-1)*cosang-(i-1)*sinang)
         yitraj(ip)=ycorn+gsp*((i-1)*cosang+(j-1)*sinang)
         zktraj(ip)=zlev(1)
         do iz=2,nzlev
            ipa=(iz-1)*np+ip
	    if (ipa .gt. maxtraj) then
	      write(6,*) 'Number of trajectories exceeds maxtraj'
	      write(6,*) 'Estimated number of trajectories = ',np*nzlev
	      stop 'mktrjpts'
	    endif
c           write(0,*) 'ipa = ',ipa,' maxtraj = ',maxtraj
            xjtraj(ipa)=xjtraj(ip)
            yitraj(ipa)=yitraj(ip)
            zktraj(ipa)=zlev(iz)
         enddo
      enddo
      enddo
c
      return
      end
