      SUBROUTINE SDVWRITE(DET,LAMBDA,SEF0,DMG,DMGNL,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION STATEV(NSDV),DET,LAMBDA,SEF0,DMG,DMGNL
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
       STATEV(1)=DET
       STATEV(2)=LAMBDA
       STATEV(3)=SEF0
       STATEV(4)=DMG
       STATEV(5)=DMGNL

      RETURN
C
      END SUBROUTINE SDVWRITE
