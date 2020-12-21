SUBROUTINE hfilfic(h,hh,pp,lambda,m,rw,ndi)



!>      NON-AFFINE NETWORK: STRUCTURE TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: h(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: hh(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pp
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN)             :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw



INTEGER :: i1,j1,k1,l1

DOUBLE PRECISION :: aux0,aux, pi,aux1

pi=four*ATAN(one)
aux0=four*pi
aux=(lambda**(pp-two))*rw
aux1=(pp-two)*(lambda**(pp-four))*rw

DO i1=1,ndi
  DO j1=1,ndi
    h(i1,j1)=aux*m(i1)*m(j1)
    DO k1=1,ndi
      DO l1=1,ndi
        hh(i1,j1,k1,l1)=aux1*m(i1)*m(j1)*m(k1)*m(l1)
      END DO
    END DO
  END DO
END DO

RETURN

END SUBROUTINE hfilfic
