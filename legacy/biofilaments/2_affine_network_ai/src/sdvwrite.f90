SUBROUTINE sdvwrite(det,statev)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1
!
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=0
statev(pos1+1)=det

RETURN

END SUBROUTINE sdvwrite
