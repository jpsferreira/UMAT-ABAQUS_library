SUBROUTINE anisomat(sseaniso,daniso,diso,k1,k2,kdisp,i4,i1)
!
use global
IMPLICIT NONE
!>     ANISOTROPIC PART : ISOCHORIC SEF AND DERIVATIVES

DOUBLE PRECISION, INTENT(OUT)            :: sseaniso
DOUBLE PRECISION, INTENT(OUT)            :: daniso(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: k1
DOUBLE PRECISION, INTENT(IN)             :: k2
DOUBLE PRECISION, INTENT(IN)             :: kdisp
DOUBLE PRECISION, INTENT(IN)             :: i4
DOUBLE PRECISION, INTENT(IN)             :: i1


DOUBLE PRECISION :: sseiso

DOUBLE PRECISION :: dudi1,d2ud2i1
DOUBLE PRECISION :: e1,ee2,ee3,dudi4,d2ud2i4,d2dudi1di4,d2dudi2di4

dudi1=diso(1)
d2ud2i1=diso(3)

e1=i4*(one-three*kdisp)+i1*kdisp-one

sseaniso=(k1/k2)*(DEXP(k2*e1*e1)-one)

IF(e1 > zero) THEN
  
  ee2=DEXP(k2*e1*e1)
  ee3=(one+two*k2*e1*e1)
  
  dudi1=dudi1+k1*kdisp*e1*ee2
  d2ud2i1=d2ud2i1+k1*kdisp*kdisp*ee3*ee2
  
  dudi4=k1*(one-three*kdisp)*e1*ee2
  
  d2ud2i4=k1*((one-three*kdisp)**two)*ee3*ee2
  
  d2dudi1di4=k1*(one-three*kdisp)*kdisp*ee3*ee2
  d2dudi2di4=zero
  
  
  
ELSE
  dudi4=zero
  d2ud2i4=zero
  d2dudi1di4=zero
  d2dudi2di4=zero
  
  d2ud2i1=zero
  
END IF
!FIRST DERIVATIVE OF SSEANISO IN ORDER TO I1
daniso(1)=dudi4
!FIRST DERIVATIVE OF SSEANISO IN ORDER TO I2
daniso(2)=d2ud2i4
!2ND DERIVATIVE OF SSEANISO IN ORDER TO I1
daniso(3)=d2dudi1di4
!2ND DERIVATIVE OF SSEANISO IN ORDER TO I2
daniso(4)=d2dudi2di4

diso(1)=dudi1
diso(3)=d2ud2i1

RETURN
END SUBROUTINE anisomat
