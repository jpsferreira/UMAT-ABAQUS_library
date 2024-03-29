      SUBROUTINE SDVWRITE(FRAC,RU0,DET,VARACT,DIRMAX,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER VV,POS1,POS2,POS3,I1
      DOUBLE PRECISION STATEV(NSDV),DET,VARACT
      DOUBLE PRECISION FRAC(4),RU0(NWP),DIRMAX(3)
C
        POS1=0
        DO I1=1,4
         POS2=POS1+I1
         STATEV(POS2)=FRAC(I1)
        ENDDO
C
        DO I1=1,NWP
          POS3=POS2+I1
          STATEV(POS3)=RU0(I1)
        ENDDO
        STATEV(POS3+1)=DET
        STATEV(POS3+2)=VARACT 
        STATEV(POS3+3)=DIRMAX(1)
        STATEV(POS3+4)=DIRMAX(2)
        STATEV(POS3+5)=DIRMAX(3)
      RETURN
C
      END SUBROUTINE SDVWRITE
