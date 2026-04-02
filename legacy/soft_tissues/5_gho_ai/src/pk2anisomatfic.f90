SUBROUTINE pk2anisomatfic(afic,daniso,cbar,inv4,st0,ndi)
!
use global
implicit none
!>      ANISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
!       INPUT:
!       DANISO - ANISOTROPIC STRAIN-ENERGY DERIVATIVES
!       CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
!       INV1,INV4 -CBAR INVARIANTS
!       UNIT2 - 2ND ORDER IDENTITY TENSOR
!       OUTPUT:
!       AFIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: afic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: daniso(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: cbar(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: inv4
DOUBLE PRECISION, INTENT(IN)             :: st0(3,3)

DOUBLE PRECISION :: dudi4,di4dc(3,3)



!-----------------------------------------------------------------------------
!FIRST DERIVATIVE OF SSEANISO IN ORDER TO I4
dudi4=daniso(1)

di4dc=st0

afic=two*(dudi4*di4dc)

RETURN
END SUBROUTINE pk2anisomatfic
