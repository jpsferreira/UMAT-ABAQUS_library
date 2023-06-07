      SUBROUTINE PROJEUL(A,AA,PE,NDI)
C>    EULERIAN PROJECTION TENSOR
C      INPUTS:
C          IDENTITY TENSORS - A, AA
C      OUTPUTS:
C          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
C
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER I,J,K,L,NDI
C
      DOUBLE PRECISION A(NDI,NDI),AA(NDI,NDI,NDI,NDI),
     1                 PE(NDI,NDI,NDI,NDI)
C
      DO I=1,NDI
         DO J=1,NDI
          DO K=1,NDI
             DO L=1,NDI
              PE(I,J,K,L)=AA(I,J,K,L)-(ONE/THREE)*(A(I,J)*A(K,L))
           END DO
          END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PROJEUL
