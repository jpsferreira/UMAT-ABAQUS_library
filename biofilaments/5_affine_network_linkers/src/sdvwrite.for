      SUBROUTINE SDVWRITE(DET,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER VV,POS1,POS2,POS3,I1
      DOUBLE PRECISION STATEV(NSDV),DET
C
        STATEV(1)=DET
      RETURN
C
      END SUBROUTINE SDVWRITE
