SUBROUTINE sigfibfic(sfic,rho,dw,m,rw,ndi)
!>    SINGLE FIBER:  'FICTICIOUS' CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi !number of dimensions
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi) !ficticious cauchy stress 
DOUBLE PRECISION, INTENT(IN)             :: rho  !angular density at m
DOUBLE PRECISION, INTENT(IN)             :: dw !derivative of fiber strain energy
DOUBLE PRECISION, INTENT(IN)             :: m(ndi) !direction vector
DOUBLE PRECISION, INTENT(IN)             :: rw ! integration weight


INTEGER :: i1,j1

DOUBLE PRECISION :: aux

aux = rho*dw*rw
DO i1=1,ndi
  DO j1=1,ndi
    sfic(i1,j1)=aux*m(i1)*m(j1)
  END DO
END DO

RETURN
END SUBROUTINE sigfibfic
