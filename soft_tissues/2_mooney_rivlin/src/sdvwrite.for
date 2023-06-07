      SUBROUTINE SDVWRITE(DET,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION STATEV(NSDV),DET
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
       STATEV(1)=DET

      RETURN
C
      END SUBROUTINE SDVWRITE
