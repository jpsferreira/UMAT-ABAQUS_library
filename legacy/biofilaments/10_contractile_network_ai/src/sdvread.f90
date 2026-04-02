SUBROUTINE sdvread(frac,ru0,statev)
use global
implicit none
!>    VISCOUS DISSIPATION: READ STATE VARS
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)


INTEGER :: vv,pos1,pos2,pos3,i1

DOUBLE PRECISION, INTENT(OUT)            :: frac(4)
DOUBLE PRECISION, INTENT(OUT)            :: ru0(nwp)


pos1=0
DO i1=1,4
  pos2=pos1+i1
  frac(i1)=statev(pos2)
END DO

DO i1=1,nwp
  pos3=pos2+i1
  ru0(i1)=statev(pos3)
END DO

RETURN

END SUBROUTINE sdvread
