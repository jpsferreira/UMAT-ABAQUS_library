SUBROUTINE initialize(statev)
use global
IMPLICIT NONE

!      DOUBLE PRECISION TIME(2),KSTEP
INTEGER :: i1,pos,pos1,pos2,pos3,vv


DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)


pos1=0
!       DETERMINANT
statev(pos1+1)=one
!        CONTRACTION VARIANCE
statev(pos1+2)=zero

RETURN

END SUBROUTINE initialize
