c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getdash(idash,ndot)
c
      include 'comconst'
c
      isolid=int(abs(idash)/10.+.00001)
      iblank=abs(idash)-isolid*10
      itot=isolid+iblank
c
      if (itot.lt.1) then
         write(iup,*)'Bad idash - check pl.spec. file.'
         stop
      endif
c
      iss=isolid
      ibb=itot
      ndot=0
      do 100 i=1,16
   50    continue
         if (i.le.iss) then
            iswitch=1
         elseif (i.le.ibb) then
            iswitch=0
         else
            iss=iss+itot
            ibb=ibb+itot
            goto 50
         endif
         ndot=ndot+iswitch*2**(16-i)
  100 continue
      ndot=ndot*(idash/abs(idash))
      return
      end
