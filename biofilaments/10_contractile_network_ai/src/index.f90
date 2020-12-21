SUBROUTINE indexx(stress,ddsdde,sig,tng,ntens,ndi)



!>    INDEXATION: FULL SIMMETRY  IN STRESSES AND ELASTICITY TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN)                      :: ntens
DOUBLE PRECISION, INTENT(OUT)            :: stress(ntens)
DOUBLE PRECISION, INTENT(OUT)            :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: tng(ndi,ndi,ndi,ndi)



INTEGER :: ii1(6),ii2(6), i1,j1


DOUBLE PRECISION :: pp1,pp2

ii1(1)=1
ii1(2)=2
ii1(3)=3
ii1(4)=1
ii1(5)=1
ii1(6)=2

ii2(1)=1
ii2(2)=2
ii2(3)=3
ii2(4)=2
ii2(5)=3
ii2(6)=3

DO i1=1,ntens
!       STRESS VECTOR
  stress(i1)=sig(ii1(i1),ii2(i1))
  DO j1=1,ntens
!       DDSDDE - FULLY SIMMETRY IMPOSED
    pp1=tng(ii1(i1),ii2(i1),ii1(j1),ii2(j1))
    pp2=tng(ii1(i1),ii2(i1),ii2(j1),ii1(j1))
    ddsdde(i1,j1)=(one/two)*(pp1+pp2)
  END DO
END DO

RETURN

END SUBROUTINE indexx
