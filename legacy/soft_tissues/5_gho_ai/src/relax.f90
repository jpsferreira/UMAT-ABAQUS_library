SUBROUTINE relax(qv,hv,aux1,hv0,pkiso,dtime,tau,teta,ndi)



!>    VISCOUS DISSIPATION: STRESS RELAXATION TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: qv(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: hv(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aux1
DOUBLE PRECISION, INTENT(IN)             :: hv0(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: tau
DOUBLE PRECISION, INTENT(IN)             :: teta



INTEGER :: i1,j1

DOUBLE PRECISION :: aux

qv=zero
hv=zero

aux=DEXP(-dtime*((two*tau)**(-one)))
aux1=teta*aux
DO i1=1,ndi
  DO j1=1,ndi
    qv(i1,j1)=hv0(i1,j1)+aux1*pkiso(i1,j1)
    hv(i1,j1)=aux*(aux*qv(i1,j1)-teta*pkiso(i1,j1))
  END DO
END DO

RETURN
END SUBROUTINE relax
