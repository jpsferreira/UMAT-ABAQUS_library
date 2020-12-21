SUBROUTINE pull2(pk,sig,finv,det,ndi)



!>       PULL-BACK TIMES DET OF A 2ND ORDER TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+det*finv(i1,ii1)*finv(j1,jj1)*sig(ii1,jj1)
      END DO
    END DO
    pk(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE pull2
