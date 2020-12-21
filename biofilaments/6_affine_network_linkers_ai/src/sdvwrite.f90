SUBROUTINE sdvwrite(statev,phi_tau,cr,stress)



!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
IMPLICIT NONE

!DOUBLE PRECISION, INTENT(IN OUT)         :: det
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!DOUBLE PRECISION, INTENT(IN OUT)         :: dirmax(3)
DOUBLE PRECISION, INTENT(IN)             :: phi_tau,cr,stress
!DOUBLE PRECISION, INTENT(IN OUT)         :: dmudx(3,1)

statev(1)=phi_tau
statev(2)=cr
statev(3)=stress

RETURN

END SUBROUTINE sdvwrite
