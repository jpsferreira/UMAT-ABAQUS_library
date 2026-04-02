SUBROUTINE sigfilfic(sfic,rho,lambda,dw,m,rw,ndi)



!>    SINGLE FILAMENT:  'FICTICIOUS' CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi !number of dimensions
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi) !ficticious cauchy stress 
DOUBLE PRECISION, INTENT(IN)             :: rho  !angular density at m
DOUBLE PRECISION, INTENT(IN)             :: lambda !filament stretch
DOUBLE PRECISION, INTENT(IN)             :: dw !derivative of filament strain energy
DOUBLE PRECISION, INTENT(IN)             :: m(ndi) !direction vector
DOUBLE PRECISION, INTENT(IN)             :: rw ! integration weights


INTEGER :: i1,j1

DOUBLE PRECISION :: aux

aux=rho*lambda**(-one)*rw*dw
DO i1=1,ndi
  DO j1=1,ndi
    sfic(i1,j1)=aux*m(i1)*m(j1)
  END DO
END DO

RETURN
END SUBROUTINE sigfilfic
