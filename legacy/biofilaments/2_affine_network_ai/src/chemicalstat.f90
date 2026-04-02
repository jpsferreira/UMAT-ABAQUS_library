SUBROUTINE chemicalstat(frac,frac0,k,dtime)


use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: frac(4)
DOUBLE PRECISION, INTENT(IN)             :: frac0(4)
DOUBLE PRECISION, INTENT(IN)             :: k(7)
DOUBLE PRECISION, INTENT(IN)             :: dtime


DOUBLE PRECISION :: stiff(4,4)

DOUBLE PRECISION :: aux
INTEGER :: i1,j1

frac=zero
stiff=zero
stiff(1,1)=-k(1)
stiff(1,2)=k(2)
stiff(1,4)=k(7)
stiff(2,1)=k(1)
stiff(2,2)=-k(2)-k(3)
stiff(2,3)=k(4)
stiff(3,2)=k(3)
stiff(3,3)=-k(4)-k(5)
stiff(3,4)=k(6)
stiff(4,3)=k(5)
stiff(4,4)=-k(6)-k(7)

DO i1=1,4
  aux=zero
  DO j1=1,4
    aux=aux+stiff(i1,j1)*frac0(j1)
  END DO
  frac(i1)=aux*dtime+frac0(i1)
END DO

!      FRAC0=FRAC

RETURN
END SUBROUTINE chemicalstat
