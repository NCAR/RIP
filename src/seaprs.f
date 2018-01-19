      SUBROUTINE SEAPRS(T,PRS,TER,SFP,IMX,JMX,KX,SLP,iup)
C
C     SECTION  DIAGNOSTIC
C     PURPOSE  COMPUTES SEA LEVEL PRESSURE FROM THE RULE
C              T1/T2=(P1/P2)**(GAMMA*R/G).
c        (this routine was swiped from the Graph program)
C
C     INPUT       T        TEMPERATURE                CROSS    3D
C                 PRS      PRESSURE                   CROSS    3D
C                 TER      TERRAIN                    CROSS    2D
C                 SFP      SURFACE PRESSURE           CROSS    2D
C                 IMX      DOT POINT DIMENSION N-S
C                 JMX      DOT POINT DIMENSION E-W
C                 KX       NUMBER OF VERTICAL LEVELS
C
C     OUTPUT      SLP      SEA LEVEL PRESSURE         CROSS    2D
C
      DIMENSION T(IMX,JMX,KX), PRS(IMX,JMX,KX), SFP(IMX,JMX) ,
     *          TER(IMX,JMX), SLP(IMX,JMX)
      PARAMETER (R=287.04,G=9.8,GAMMA=6.5E-3)
      PARAMETER (TC=273.16+17.5) ! T CRITICAL IN PSFC/PSLV
      PARAMETER (PCONST=100.)
C

      LOGICAL L1,L2,L3
C
C     ... SEA LEVEL PRESSURE
C
      XTERM=GAMMA*R/G
C
C     ... COMPUTE PRESSURE AT PCONST hPa ABOVE SURFACE (PL)
C
      KUPTO=KX/2
99    CONTINUE
      DO 600 J=1,JMX-1
      DO 600 I=1,IMX-1
         PL=SFP(I,J)-PCONST
         XKLEV=0.
C
C     ... FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J
C
         DO 125 K=KX-1,KUPTO,-1
            XK=FLOAT(K)
            XKHOLD=XKLEV
         if ( ((PRS(I,J,K).LT.PL) .AND. (PRS(I,J,K+1).GE.PL)) ) then
              xklev = xk
            else
              xklev = xkhold
            endif
c           XKLEV=CVMGT(XK,XKHOLD,
c    *         ((PRS(I,J,K).LT.PL) .AND. (PRS(I,J,K+1).GE.PL)))
125      CONTINUE
         IF(XKLEV.LT.1.) THEN
            WRITE(IUP,*)'ERROR FINDING PRESSURE LEVEL ',PCONST,' hPa ',
     *              'ABOVE THE SURFACE'
            WRITE(IUP,*)'LAST K LEVEL =',KUPTO
            IF(KUPTO.NE.1) THEN
               WRITE(IUP,*)'TRYING AGAIN WITH KUPTO=1'
               KUPTO=1
               GOTO 99
            ELSE
               WRITE(IUP,*)'I,J=',I,J
               WRITE(IUP,*)'PL=',PL
               WRITE(IUP,*)'PSFC=',SFP(I,J)
               STOP 'seaprs'
            END IF
         END IF
C
C     ... GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)
C         AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE
C
         KLO=NINT(XKLEV)+1
         KHI=NINT(XKLEV)
         PLO=PRS(I,J,KLO)
         PHI=PRS(I,J,KHI)
         TLO=T(I,J,KLO)
         THI=T(I,J,KHI)
         TL=THI-(THI-TLO)*ALOG(PL/PHI)/ALOG(PLO/PHI)
         TS=TL*(SFP(I,J)/PL)**XTERM
         TBAR=(TS+TL)*0.5
         HL=TER(I,J)-R/G*ALOG(PL/SFP(I,J))*TBAR
         T0=TL+GAMMA*HL
C
C     ... CORRECT SEA LEVEL TEMPERATURE IF TOO HOT
C
         L1=T0.LT.TC
         L2=TS.LE.TC
         L3=.NOT.L1
         T0HOLD=T0
          if (l1 .and. l2) then
            t0 = T0HOLD
          else
            if (l2 .and. l3) then
              t0 = tc
            else
              t0 = TC-0.005*(TS-TC)**2
            endif
          endif
c        T0=CVMGT(T0HOLD,
c    *      CVMGT(TC,TC-0.005*(TS-TC)**2,L2.AND.L3),
c    *      L1.AND.L2)
C
C     ... COMPUTE SEA LEVEL PRESSURE
C
         SLP(I,J)=SFP(I,J)*EXP(2.*G*TER(I,J)/(R*(TS+T0)))
600   CONTINUE
      RETURN
      END
