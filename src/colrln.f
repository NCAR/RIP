c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine colrln (xcs,ycs,ncs,iai,iag,nai)
      dimension xcs(*),ycs(*),iai(*),iag(*)
      common /emap/ llcolor,lllinw,llndot,mpfillco(6),llmask,ioutype
c
c   This is a user-supplied subroutine that is used by NCAR Graphics
c   for making color-filled maps.
c
        call dashdb (llndot)
        call gqplci(ierr,ioldc)
        call gsplci(llcolor)
        call getusv('LW',ioldlw)
        call setusv('LW',lllinw*1000)
        itm=1
        do 101 i=1,nai
          if (iai(i).lt.0) itm=0
  101   continue
        if (itm.ne.0) then
          itm=0
          do 102 i=1,nai
            if (iag(i).eq.1) itm=iai(i)
  102     continue
          if (ioutype.eq.1) then
             iaci=mapaci(itm)
          else
             iaci=mpisci(itm)
          endif
          if ((llmask.eq. 1.and.iaci.eq.1).or.
     &        (llmask.eq.-1.and.iaci.ne.1).or.llmask.eq.0) then
            call plotit(0,0,0)
            call curved (xcs,ycs,ncs)
          endif
        endif
        call dashdb (65535)
        call gsplci(ioldc)
        call setusv('LW',ioldlw)
        return
      end
