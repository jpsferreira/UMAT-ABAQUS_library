      SUBROUTINE PUSH4(SPATIAL,MAT,F,DET,NDI)
C>        PIOLA TRANSFORMATION
C>      INPUT:
C>       MAT - MATERIAL ELASTICITY TENSOR
C>       F - DEFORMATION GRADIENT
C>       DET - DEFORMATION DETERMINANT
C>      OUTPUT:
C>       SPATIAL - SPATIAL ELASTICITY TENSOR
       IMPLICIT NONE
       INCLUDE 'param_umat.inc'
C
       INTEGER I1,J1,K1,L1,II1,JJ1,KK1,LL1,NDI
       DOUBLE PRECISION MAT(NDI,NDI,NDI,NDI),F(NDI,NDI)
       DOUBLE PRECISION SPATIAL(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION AUX,DET
C
       DO I1=1,NDI
        DO J1=1,NDI
         DO K1=1,NDI
          DO L1=1,NDI
           AUX=ZERO
           DO II1=1,NDI
            DO JJ1=1,NDI
             DO KK1=1,NDI
              DO LL1=1,NDI
              AUX=AUX+(DET**(-ONE))*
     +        F(I1,II1)*F(J1,JJ1)*
     +        F(K1,KK1)*F(L1,LL1)*MAT(II1,JJ1,KK1,LL1)
              END DO
             END DO
            END DO
           END DO
           SPATIAL(I1,J1,K1,L1)=AUX
          END DO
         END DO
        END DO
       END DO
C
       RETURN
      END SUBROUTINE PUSH4
