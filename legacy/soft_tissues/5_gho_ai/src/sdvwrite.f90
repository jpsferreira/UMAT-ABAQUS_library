SUBROUTINE sdvwrite(det,lambda,statev)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1
!
DOUBLE PRECISION, INTENT(IN)             :: det,lambda
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=0
statev(pos1+1)=det
statev(pos1+2)=lambda

RETURN

END SUBROUTINE sdvwrite
