c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine rdcolt (rip_root,nco,conam,fred,fgreen,fblue)
c
c   This routine reads the color table defined by the user.
c
      dimension fred(0:255),fgreen(0:255),fblue(0:255)
      character*40  conam(0:255)*40, string*80, cocheck*11,
     &   rip_root*256, fname*256
      logical       red,green,blue
c
      include 'comconst'
c
      cocheck='COLOR TABLE'
c
c   First, check if a color table file exists in rip input file,
c   which should already be open and rewound to first record.
c
 150  read (iuinput,'(a80)',end=160) string
      if (index(string,cocheck).ne.0) then
         iun=iuinput
         goto 300
      endif
      goto 150
 160  continue
      rewind(iuinput)
c
c   No table found there, so check for the existence of the file
c   rip_root/color.tbl.  Borrow unit number for stationlist.
c
      iendrr=index(rip_root,' ')-1
      fname=rip_root(1:iendrr)//'/color.tbl'
      open(unit=iustnlist,file=fname,form='formatted',status='old',
     &   err=260)
 250  read (iustnlist,'(a80)',end=265) string
      if (index(string,cocheck).ne.0) then
         iun=iustnlist
         goto 300
      endif
      goto 250
 260  ierflag=1
      goto 270
 265  ierflag=2
      rewind(iuinput)
c
 270  write(iup,*)'RIP couldn''t find color table in rip input file'
      if (ierflag.eq.1) then
         write(iup,*)'nor could it find the file rip_root/color.tbl.'
      elseif (ierflag.eq.2) then
         write(iup,*)'or in the file rip_root/color.tbl.'
      endif
      stop
c
 300  continue
      read(iun,*)
      read(iun,*)
c
c  Read in an entire line from the color table
c
      ico=-1
10    continue
      read (iun,110) string(1:80)
      if (string(1:1) .eq. '-') goto 90
      ico=ico+1
      if (ico.eq.201) then
         write(iup,*) 'You''re trying to define more than 200 colors'
         write(iup,*) 'in your color table.  This will leave less than'
         write(iup,*) '55 color indices (per frame) for new color'
         write(iup,*) 'shades in color-filled plots, which may not'
         write(iup,*) 'be enough.  Consider yourself warned.'
      endif
c
c    Initialize flags
c
      red   = .false.
      green = .false.
      blue  = .false.
c
c  Parse the line for its components
c
      i = 1
      igotconam = 0
30    continue
      if ((string(i:i) .ne. '|') .and. (string(i:i) .ne. '!')) then
         i = i + 1
         goto 30
      endif
      if (igotconam.eq.0) then
         conam(ico)=string(1:i-1)
         igotconam=1
      endif
      i = i + 1
c
c    Parse to a non-blank character
c
40    continue
      if (string(i:i) .eq. ' ') then
         i = i + 1
         if (i.le.len(string)) goto 40
      endif
      if (blue)  goto 70
      if (green) goto 60
      if (red)   goto 50
c
c    Red component
c
      if (i.gt.(len(string)-3)) i = len(string)-3
      read (string(i:i+3),120) fred(ico)
      red = .true.
      i = i + 4
      goto 30
c
c    Green component
c
50    continue
      if (i.gt.(len(string)-3)) i = len(string)-3
      read (string(i:i+3),120) fgreen(ico)
      green = .true.
      i = i + 4
      goto 30
c
c    Blue component
c
60    continue
      if (i.gt.(len(string)-3)) i = len(string)-3
      read (string(i:i+3),120) fblue(ico)
      blue = .true.
      i = i + 4
      goto 30
c
c    Color number
c
70    continue
c      if (string(i+1:i+1) .ne. ' ') then
c         read (string(i:i+1),130) icotable
c      else
c         read (string(i:i),140) icotable
c      endif
c      if (icotable.ne.ico) then
c         write(iup,*)'The colors are not properly numbered',
c     &      ' in the table.'
c         stop
c      endif
c
      goto 10
90    write(iup,*)'Color Names and Fractions Defined'
      nco=ico
c
      if (iun.eq.iustnlist) then
         close (iustnlist)
      elseif (iun.eq.iuinput) then
         rewind (iuinput)
      endif
c
c  Format statements
c
110   format (A80)
120   format (F4.2)
130   format (I2)
140   format (I1)
c
      return
      end
