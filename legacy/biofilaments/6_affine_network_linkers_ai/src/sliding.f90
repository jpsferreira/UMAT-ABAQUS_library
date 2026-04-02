SUBROUTINE sliding(ffc,ru,ffc0,ru0,ffcmax,fric,frac,dtime)

use global
implicit none

DOUBLE PRECISION, INTENT(OUT)            :: ffc
DOUBLE PRECISION, INTENT(OUT)            :: ru
DOUBLE PRECISION, INTENT(IN)             :: ffc0
DOUBLE PRECISION, INTENT(IN OUT)         :: ru0
DOUBLE PRECISION, INTENT(IN)             :: ffcmax
DOUBLE PRECISION, INTENT(IN)         :: fric
DOUBLE PRECISION, INTENT(IN)             :: frac(4)
DOUBLE PRECISION, INTENT(IN)             :: dtime





DOUBLE PRECISION :: aux0,aux1,aux2, arg
!      INTEGER STATE

aux1=frac(3)
aux0=aux1+frac(4)
aux2=aux0
arg=ffc0/ffcmax

IF(arg < aux1) THEN
  ffc=aux1*ffcmax
ELSE IF (arg > aux2)THEN
  ffc=aux2*ffcmax
ELSE
  ffc=ffc0
END IF

ru=ru0+dtime*(fric**(-one))*(ffc-ffc0)
ru0=ru

RETURN

END SUBROUTINE sliding
