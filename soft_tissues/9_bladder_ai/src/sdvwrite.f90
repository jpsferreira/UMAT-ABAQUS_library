SUBROUTINE sdvwrite(det,lambda,frac_pld,prefang,statev)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1
!
DOUBLE PRECISION, INTENT(IN)             :: det,lambda,frac_pld,prefang
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=0
statev(pos1+1)=det
statev(pos1+2)=lambda
statev(pos1+3)=frac_pld
statev(pos1+4)=prefang


RETURN

END SUBROUTINE sdvwrite
