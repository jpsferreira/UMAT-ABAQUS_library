SUBROUTINE deffil(lambda,m,m0,f,ndi)



!>      SINGLE FILAMENT: STRETCH AND DEFORMED DIRECTION
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: lambda
DOUBLE PRECISION, INTENT(OUT)            :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: m0(ndi)
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)


INTEGER :: i1,j1

DOUBLE PRECISION :: aux

lambda=zero
DO i1=1,ndi
  aux=zero
  DO j1=1,ndi
    aux=aux+f(i1,j1)*m0(j1)
  END DO
  m(i1)=aux
END DO
lambda=dot_product(m,m)
lambda=SQRT(lambda)

RETURN
END SUBROUTINE deffil
