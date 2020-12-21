SUBROUTINE setvol(cvol,pv,ppv,unit2,unit4s,ndi)



!>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4s(ndi,ndi,ndi,ndi)


INTEGER :: i1,j1,k1,l1



DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= ppv*unit2(i1,j1)*unit2(k1,l1)  &
            -two*pv*unit4s(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setvol
