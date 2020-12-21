SUBROUTINE csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)



!>    ISOTROPIC MATRIX: SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: distgr(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



call push4(cisomatfic,cmisomatfic,distgr,det,ndi)

RETURN
END SUBROUTINE csisomatfic
