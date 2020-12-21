SUBROUTINE contraction24(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  4TH ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO k1=1,ndi
  DO l1=1,ndi
    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(k1,l1)*rt(i1,j1,k1,l1)
      END DO
    END DO
    s(k1,l1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction24
