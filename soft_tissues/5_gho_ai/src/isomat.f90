SUBROUTINE isomat(sseiso,diso,c10,c01,cbari1,cbari2)



!>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: sseiso
DOUBLE PRECISION, INTENT(OUT)            :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: c10
DOUBLE PRECISION, INTENT(IN)             :: c01
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari1
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2


sseiso=c10*(cbari1-three)+c01*(cbari2-three)

diso(1)=c10
diso(2)=c01
diso(3)=zero
diso(4)=zero
diso(5)=zero

RETURN
END SUBROUTINE isomat
