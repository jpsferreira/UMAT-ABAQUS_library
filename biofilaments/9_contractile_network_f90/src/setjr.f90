SUBROUTINE setjr(cjr,sigma,unit2,ndi)


use global
IMPLICIT NONE
!>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cjr(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sigma(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1


DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        
        cjr(i1,j1,k1,l1)= (one/two)*(unit2(i1,k1)*sigma(j1,l1)  &
            +sigma(i1,k1)*unit2(j1,l1)+unit2(i1,l1)*sigma(j1,k1)  &
            +sigma(i1,l1)*unit2(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setjr
