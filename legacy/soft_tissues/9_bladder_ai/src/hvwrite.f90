SUBROUTINE hvwrite(statev,hv,v1,ndi)



!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
DOUBLE PRECISION, INTENT(IN)             :: hv(ndi,ndi)
INTEGER, INTENT(IN)                      :: v1



INTEGER :: pos


pos=9*v1-9
statev(1+pos)=hv(1,1)
statev(2+pos)=hv(1,2)
statev(3+pos)=hv(1,3)
statev(4+pos)=hv(2,1)
statev(5+pos)=hv(2,2)
statev(6+pos)=hv(2,3)
statev(7+pos)=hv(3,1)
statev(8+pos)=hv(3,2)
statev(9+pos)=hv(3,3)

RETURN

END SUBROUTINE hvwrite
