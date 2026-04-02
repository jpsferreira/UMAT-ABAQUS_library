SUBROUTINE density(rho,ang,bb,erfi)



!>    SINGLE FILAMENT: DENSITY FUNCTION VALUE
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: rho
DOUBLE PRECISION, INTENT(IN OUT)         :: ang
DOUBLE PRECISION, INTENT(IN OUT)         :: bb
DOUBLE PRECISION, INTENT(IN OUT)         :: erfi



DOUBLE PRECISION :: pi,aux1,aux2

pi=four*ATAN(one)
aux1=SQRT(bb/(two*pi))
aux2=DEXP(bb*(COS(two*ang)+one))
rho=four*aux1*aux2*(erfi**(-one))
!      RHO=RHO*((FOUR*PI)**(-ONE))

RETURN
END SUBROUTINE density
