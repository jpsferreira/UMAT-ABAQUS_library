SUBROUTINE setiso(ciso,cfic,pe,siso,sfic,unit2,ndi)


use global
IMPLICIT NONE

!>    ISOCHORIC SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ciso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pe(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz

cisoaux1=zero
cisoaux=zero

CALL contraction44(cisoaux1,pe,cfic,ndi)
CALL contraction44(cisoaux,cisoaux1,pe,ndi)

trfic=zero
DO i1=1,ndi
  trfic=trfic+sfic(i1,i1)
END DO

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=cisoaux(i1,j1,k1,l1)
        yy=trfic*pe(i1,j1,k1,l1)
        zz=siso(i1,j1)*unit2(k1,l1)+unit2(i1,j1)*siso(k1,l1)
        
        ciso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setiso
