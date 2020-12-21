SUBROUTINE vol(ssev,pv,ppv,k,det,dete,dets)



!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN OUT)         :: dete
DOUBLE PRECISION, INTENT(IN OUT)         :: dets



DOUBLE PRECISION :: g, aux

g=(one/four)*(det*det-one-two*log(det))
ssev=k*g

!no sweelling
pv=k*(one/two)*(det-one/det)
aux=k*(one/two)*(one+one/(det*det))
ppv=pv+det*aux
dete=det
dets=det

!sweelling
!pv=dets*k*(dlog(dete))/det
!aux=dets*k*(one/two)*(one+one/(dete*dete))
!ppv=pv+dete*aux

RETURN
END SUBROUTINE vol
