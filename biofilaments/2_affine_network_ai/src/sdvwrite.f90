SUBROUTINE sdvwrite(frac,ru0,det,varact,dirmax,statev)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: vv,pos1,pos2,pos3,i1
!
DOUBLE PRECISION, INTENT(IN)             :: frac(4)
DOUBLE PRECISION, INTENT(IN)             :: ru0(nwp)
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN)             :: varact
DOUBLE PRECISION, INTENT(IN)             :: dirmax(3)
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=0
DO i1=1,4
  pos2=pos1+i1
  statev(pos2)=frac(i1)
END DO

DO i1=1,nwp
  pos3=pos2+i1
  statev(pos3)=ru0(i1)
END DO
statev(pos3+1)=det
statev(pos3+2)=varact
statev(pos3+3)=dirmax(1)
statev(pos3+4)=dirmax(2)
statev(pos3+5)=dirmax(3)
RETURN

END SUBROUTINE sdvwrite
