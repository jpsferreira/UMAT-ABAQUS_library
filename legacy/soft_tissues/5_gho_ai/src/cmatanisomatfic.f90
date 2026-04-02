SUBROUTINE cmatanisomatfic(cmanisomatfic,m0,daniso,unit2,det,ndi)
!
use global
IMPLICIT NONE
!
!>    ANISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cmanisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: m0(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: daniso(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det
!
INTEGER :: i,j,k,l
DOUBLE PRECISION :: cinv4(ndi,ndi,ndi,ndi),cinv14(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: d2udi4,d2udi1di4
DOUBLE PRECISION :: imm(ndi,ndi,ndi,ndi),mmi(ndi,ndi,ndi,ndi),  &
    mm0(ndi,ndi,ndi,ndi)

!-----------------------------------------------------------------------------
!2ND DERIVATIVE OF SSEANISO IN ORDER TO I4
d2udi4=daniso(2)
!2ND DERIVATIVE OF SSEANISO IN ORDER TO I1 AND I4
d2udi1di4=daniso(3)

CALL tensorprod2(m0,m0,mm0,ndi)
CALL tensorprod2(unit2,m0,imm,ndi)
CALL tensorprod2(m0,unit2,mmi,ndi)

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        cinv4(i,j,k,l)=d2udi4*mm0(i,j,k,l)
        cinv14(i,j,k,l)=d2udi1di4*(imm(i,j,k,l)+mmi(i,j,k,l))
        cmanisomatfic(i,j,k,l)=four*(cinv4(i,j,k,l)+cinv14(i,j,k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE cmatanisomatfic
