      SUBROUTINE SETVOL(CVOL,PV,PPV,UNIT2,UNIT4S,NDI)
C>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),UNIT4S(NDI,NDI,NDI,NDI),
     1                 CVOL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION PV,PPV
C
      DO I1 = 1, NDI
        DO J1 = 1, NDI
         DO K1 = 1, NDI
           DO L1 = 1, NDI
             CVOL(I1,J1,K1,L1)=
     1                 PPV*UNIT2(I1,J1)*UNIT2(K1,L1)
     2                 -TWO*PV*UNIT4S(I1,J1,K1,L1)
           END DO
         END DO
        END DO
      END DO
C      
      RETURN
      END SUBROUTINE SETVOL
