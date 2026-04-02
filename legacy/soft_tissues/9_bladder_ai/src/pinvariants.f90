SUBROUTINE pinvariants(a,inv4,ndi,st,lambda,barlambda,det)
!
use global
implicit none
!>    AND 4TH PSEUDO-INVARIANTS OF A TENSOR
INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: inv4
DOUBLE PRECISION, INTENT(IN)             :: st(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: lambda
DOUBLE PRECISION, INTENT(OUT)            :: barlambda
DOUBLE PRECISION, INTENT(IN)             :: det


INTEGER :: i,j
DOUBLE PRECISION :: scale1


inv4=zero
DO i=1,ndi
  DO j=1, ndi
    inv4=inv4+a(i,j)*st(i,j)
  END DO
END DO
!     STRETCH
scale1=det**(-one /three)
barlambda=DSQRT(inv4)
lambda=barlambda/scale1

RETURN
END SUBROUTINE pinvariants
