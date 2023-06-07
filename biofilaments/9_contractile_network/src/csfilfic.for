      SUBROUTINE CSFILFIC(CFIC,RHO,LAMBDA,DW,DDW,M,RW,NDI)
C>    AFFINE NETWORK: 'FICTICIOUS' ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION CFIC(NDI,NDI,NDI,NDI),M(NDI)
      DOUBLE PRECISION RHO,AUX,DW,DDW,RW,LAMBDA,AUX0
C
      AUX0=DDW-(LAMBDA**(-ONE))*DW
      AUX=RHO*AUX0*RW*(LAMBDA**(-TWO))
      DO I1=1,NDI
       DO J1=1,NDI
        DO K1=1,NDI
         DO L1=1,NDI
          CFIC(I1,J1,K1,L1)=AUX*M(I1)*M(J1)*M(K1)*M(L1)
         END DO
        END DO
       END DO
      END DO
C
      RETURN
      END SUBROUTINE CSFILFIC
