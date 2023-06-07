      SUBROUTINE SDVWRITE(DET,LAMBDA,ACT,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION STATEV(NSDV),DET,LAMBDA,ACT
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
       STATEV(1)=DET
       STATEV(2)=LAMBDA 
       STATEV(3)=ACT

      RETURN
C
      END SUBROUTINE SDVWRITE
