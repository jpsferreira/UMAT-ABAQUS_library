SUBROUTINE sdvread(statev,phi_t,cr)



!>     READ STATE VARS
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: cr,phi_t
!DOUBLE PRECISION, INTENT(IN OUT)         :: dmudx(3,1)


phi_t=statev(1)
cr=statev(2)

RETURN

END SUBROUTINE sdvread
