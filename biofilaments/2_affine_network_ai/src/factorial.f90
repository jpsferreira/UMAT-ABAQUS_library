SUBROUTINE factorial(fact,term)



!>    FACTORIAL
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: fact
INTEGER, INTENT(IN)                      :: term



INTEGER :: m

fact = 1

DO  m = 1, term
  fact = fact * m
END DO

RETURN
END SUBROUTINE factorial
