SUBROUTINE csfilfic(cfic,rho,lambda,dw,ddw,m,rw,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rho
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN)             :: dw
DOUBLE PRECISION, INTENT(IN)             :: ddw
DOUBLE PRECISION, INTENT(IN)             :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw



INTEGER :: i1,j1,k1,l1

DOUBLE PRECISION :: aux, aux0

aux0=ddw-(lambda**(-one))*dw
aux=rho*aux0*rw*(lambda**(-two))
DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        cfic(i1,j1,k1,l1)=aux*m(i1)*m(j1)*m(k1)*m(l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE csfilfic
