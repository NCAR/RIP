c                                                                     c
c*********************************************************************c
c                                                                     c
      function virtualhyd(temp,ratmix,ratmixhyd)
c
c   This function returns virtual temperature in K, given temperature
c      in K and vapor mixing ratio in kg/kg.
c   April 2002: Stoelinga added hydrometeor mixing ratio to virtual
c      temperature.  Presence of hydrometeors falling at terminal
c      velocity adds "weight" or density to the air, which can (and
c      should according to some) be accounted for in virtual temp.
c
      parameter (eps=0.622)
c
      virtualhyd = temp * ( (eps+ratmix)/(eps*(1.+ratmix)) - ratmixhyd )
c
      return
      end
