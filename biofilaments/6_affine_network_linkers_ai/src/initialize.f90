SUBROUTINE initialize(statev,phi0)


use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
DOUBLE PRECISION, INTENT(IN)            :: phi0

statev(1)=phi0

RETURN

END SUBROUTINE initialize
