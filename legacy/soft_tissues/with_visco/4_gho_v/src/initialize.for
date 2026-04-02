       SUBROUTINE INITIALIZE(STATEV, VV)
C
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C      
C      COMMON /KCOMMON/KBLOCK
C
C      DOUBLE PRECISION TIME(2),KSTEP
      INTEGER I1,POS,POS1,VV
      DOUBLE PRECISION STATEV(NSDV)
C        read your sdvs here. they should be allocated. 
C          after the viscous terms (only if you use viscosity check hvread)
        POS1=9*VV
C        DETERMINANT
        STATEV(POS1+1)=ONE    
C     
      RETURN
C
      END SUBROUTINE INITIALIZE
