      SUBROUTINE SDVWRITE(STATEV,DET,VV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION STATEV(NSDV),DET
      INTEGER VV,POS1
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
      POS1=9*VV 
      STATEV(POS1+1)=DET

      RETURN
C
      END SUBROUTINE SDVWRITE
