      SUBROUTINE CONTRACTION44(S,LT,RT,NDI)
C>       DOUBLE CONTRACTION BETWEEN 4TH ORDER TENSORS
C>      INPUT:
C>       LT - RIGHT 4TH ORDER TENSOR
C>       RT - LEFT  4TH ORDER TENSOR
C>      OUTPUT:
C>       S - DOUBLE CONTRACTED TENSOR (4TH ORDER)
       IMPLICIT NONE
       INCLUDE 'param_umat.inc'
C
       INTEGER I1,J1,K1,L1,M1,N1,NDI
C
       DOUBLE PRECISION LT(NDI,NDI,NDI,NDI),RT(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION S(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION AUX

C                     
       DO I1=1,NDI
        DO J1=1,NDI
         DO K1=1,NDI
          DO L1=1,NDI
           AUX=ZERO
           DO M1=1,NDI
            DO N1=1,NDI
                AUX=AUX+LT(I1,J1,M1,N1)*RT(M1,N1,K1,L1)
            END DO
           END DO
           S(I1,J1,K1,L1)=AUX
          END DO
         END DO
        END DO
       END DO
C
       RETURN
      END SUBROUTINE CONTRACTION44
