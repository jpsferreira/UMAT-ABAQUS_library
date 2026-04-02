SUBROUTINE contraction22(aux,lt,rt,ndi)
!>       DOUBLE CONTRACTION BETWEEN 2nd ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  2nd ODER TENSOR
!>      OUTPUT:
!>       aux - DOUBLE CONTRACTED TENSOR (scalar)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aux
INTEGER :: i1,j1


    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(i1,j1)*rt(j1,i1)
      END DO
    END DO
RETURN
END SUBROUTINE contraction22
