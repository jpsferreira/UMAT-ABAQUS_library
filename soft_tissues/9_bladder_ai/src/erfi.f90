SUBROUTINE erfi(erf,b,nterm)
!>    IMAGINARY ERROR FUNCTION OF SQRT(B); B IS THE DISPERSION PARAM
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: erf
DOUBLE PRECISION, INTENT(IN OUT)         :: b
INTEGER, INTENT(IN)                      :: nterm


DOUBLE PRECISION :: pi
DOUBLE PRECISION :: aux,aux1,aux2,aux3,aux4,fact
INTEGER :: i1,j1

pi=four*ATAN(one)
aux=SQRT(two*b)
aux1=two*aux
aux2=(two/three)*(aux**three)
aux4=zero
DO j1=3,nterm
  i1=j1-1
  CALL factorial (fact,i1)
  aux3=two*j1-one
  aux4=aux4+(aux**aux3)/(half*aux3*fact)
END DO

erf=pi**(-one/two)*(aux1+aux2+aux4)
RETURN
END SUBROUTINE erfi
