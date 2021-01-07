      SUBROUTINE SDVWRITE(DET,SEF0,DMG,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION STATEV(NSDV),DET,SEF0,DMG
C        write your sdvs here. they should be allocated 
C                after the viscous terms (check hvwrite)
       STATEV(1)=DET
       STATEV(2)=SEF0
       STATEV(3)=DMG

      RETURN
C
      END SUBROUTINE SDVWRITE
