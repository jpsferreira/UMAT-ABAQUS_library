       SUBROUTINE CONTRACTILE(FC,FFC,DW,DDW,FFCC,RU,RU0,LAMBDAC,LAMBDAP,
     1                 LAMBDA0,L,R0,MU0,BETA,B0,FFCMAX,FRIC,FRAC,DTIME)       
C
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION DTIME,FRAC(4)
      DOUBLE PRECISION FC,FFC,RU,RU0,DW,DDW,FFCC
      DOUBLE PRECISION LAMBDAP,LAMBDAC,L,R0,MU0,BETA,B0,FRIC,FFCMAX
      DOUBLE PRECISION LAMBDA0
C      CHECK STRETCH
        IF (RU0.GT.ZERO) THEN
        LAMBDAC=LAMBDAP+RU0  
        ELSE
        LAMBDAC=LAMBDAP
        ENDIF
C        IF(LAMBDAC.LT.LAMBDAP)THEN
C          LAMBDAC=LAMBDAP
C           RU=ZERO
C        ENDIF
C        FORCE WITH CONTRACTION
          CALL FIL(FC,FFC,DW,DDW,LAMBDAC,LAMBDA0,L,R0,MU0,BETA,B0)
C        PASSIVE FORCE          
C          CALL FIL(FC,FFC,DW,DDW,LAMBDAP,LAMBDA0,L,R0,MU0,BETA,B0)
C        RELATIVE SLIDING          
          CALL SLIDING(FFCC,RU,FFC,RU0,FFCMAX,FRIC,FRAC,DTIME)
C
      RETURN
C
      END SUBROUTINE CONTRACTILE

