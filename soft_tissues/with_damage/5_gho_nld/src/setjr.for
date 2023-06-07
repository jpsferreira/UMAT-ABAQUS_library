      SUBROUTINE SETJR(CJR,SIGMA,UNIT2,NDI)
C>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),
     1                 CJR(NDI,NDI,NDI,NDI),SIGMA(NDI,NDI)
C
      DO I1 = 1, NDI
        DO J1 = 1, NDI
         DO K1 = 1, NDI
           DO L1 = 1, NDI
           
              CJR(I1,J1,K1,L1)=
     1             (ONE/TWO)*(UNIT2(I1,K1)*SIGMA(J1,L1)
     2             +SIGMA(I1,K1)*UNIT2(J1,L1)+UNIT2(I1,L1)*SIGMA(J1,K1)
     3             +SIGMA(I1,L1)*UNIT2(J1,K1))
           END DO
         END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE SETJR
