SUBROUTINE signetfic(sfic,f,mf0,rw,filprops,netprops,  &
        rho,lambda0,n,det,ndi)



!>    AFFINE NETWORK:  'FICTICIUOUS' CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(5)
DOUBLE PRECISION, INTENT(IN)             :: netprops(3)
DOUBLE PRECISION, INTENT(OUT)            :: rho
DOUBLE PRECISION, INTENT(OUT)            :: lambda0
DOUBLE PRECISION, INTENT(OUT)            :: n
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1
DOUBLE PRECISION :: sfilfic(ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta

!     FILAMENT
l      = filprops(1)
r0      = filprops(2)
mu0     = filprops(3)
beta    = filprops(4)
b0      = filprops(5)
!     NETWORK
n      = netprops(1)
lambda0 = netprops(2)
rho     = netprops(3)

pi=four*ATAN(one)
aux=n*(det**(-one))*four*pi
sfic=zero

DO i1=1,nwp
  
  mfi=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  CALL fil(fi,ffi,dwi,ddwi,lambdai,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdai,dwi,mfi,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
    END DO
  END DO
  
END DO

RETURN
END SUBROUTINE signetfic
