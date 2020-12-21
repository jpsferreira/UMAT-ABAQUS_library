SUBROUTINE onem(a,aa,aas,ndi)



!>      THIS SUBROUTINE GIVES:
!>          2ND ORDER IDENTITY TENSORS - A
!>          4TH ORDER IDENTITY TENSOR - AA
!>          4TH ORDER SYMMETRIC IDENTITY TENSOR - AAS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aas(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    IF (i == j) THEN
      a(i,j) = one
    ELSE
      a(i,j) = zero
    END IF
  END DO
END DO

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        aa(i,j,k,l)=a(i,k)*a(j,l)
        aas(i,j,k,l)=(one/two)*(a(i,k)*a(j,l)+a(i,l)*a(j,k))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE onem
