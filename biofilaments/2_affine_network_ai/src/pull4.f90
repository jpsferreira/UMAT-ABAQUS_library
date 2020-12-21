SUBROUTINE pull4(mat,spatial,finv,det,ndi)



!>        PULL-BACK TIMES DET OF 4TH ORDER TENSOR

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+det* finv(i1,ii1)*finv(j1,jj1)*  &
                    finv(k1,kk1)*finv(l1,ll1)*spatial(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        mat(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE pull4
