SUBROUTINE projeul(a,aa,pe,ndi)



!>    EULERIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!      OUTPUTS:
!          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pe(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pe(i,j,k,l)=aa(i,j,k,l)-(one/three)*(a(i,j)*a(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projeul
