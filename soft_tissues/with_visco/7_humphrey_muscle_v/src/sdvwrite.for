      SUBROUTINE SDVWRITE(STATEV,DET,LAMBDA,ACT,VV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION STATEV(NSDV),DET,LAMBDA,ACT
      INTEGER VV,POS1
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
      POS1=9*VV 
      STATEV(POS1+1)=DET
      STATEV(POS1+2)=LAMBDA 
      STATEV(POS1+3)=ACT

      RETURN
C
      END SUBROUTINE SDVWRITE
