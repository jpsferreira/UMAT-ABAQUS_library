SUBROUTINE hvread(hv,statev,v1,ndi)



!>    VISCOUS DISSIPATION: READ STATE VARS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi

DOUBLE PRECISION, INTENT(OUT)            :: hv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
INTEGER, INTENT(IN)                      :: v1



INTEGER :: pos


pos=9*v1-9
hv(1,1)=statev(1+pos)
hv(1,2)=statev(2+pos)
hv(1,3)=statev(3+pos)
hv(2,1)=statev(4+pos)
hv(2,2)=statev(5+pos)
hv(2,3)=statev(6+pos)
hv(3,1)=statev(7+pos)
hv(3,2)=statev(8+pos)
hv(3,3)=statev(9+pos)

RETURN

END SUBROUTINE hvread
