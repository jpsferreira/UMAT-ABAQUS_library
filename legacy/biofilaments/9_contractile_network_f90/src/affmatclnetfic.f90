SUBROUTINE affmatclnetfic(pkfic,cfic,f,mf0,rw,filprops,affprops,  &
        efi,noel,det,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' PK STRESS AND ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: sfilfic(ndi,ndi),   &
     cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,fi,ffi,r0
DOUBLE PRECISION :: r0c,etac,lambdaif,lambdaic
DOUBLE PRECISION :: l,r0f,mu0,b0,beta,lambda0,m,rho
DOUBLE PRECISION :: b,fric,ffmax,ang

!     FILAMENT
l       = filprops(1)
r0f     = filprops(2)
r0c     = filprops(3)
etac    = filprops(4)
mu0     = filprops(5)
beta    = filprops(6)
b0      = filprops(7)
lambda0 = filprops(8)
!     NETWORK
m       = affprops(1)
b       = affprops(2)

pi=four*ATAN(one)
aux=m*(det**(-one))*four*pi
cfic=zero
pkfic=zero
r0=r0f+r0c

rho=one

!        LOOP OVER THE INTEGRATION DIRECTIONS
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
!       CALL FIL(FI,FFI,DWI,DDWI,LAMBDAI,LAMBDA0,L,R0,MU0,BETA,B0)
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  rho = one
  
  IF((etac > zero).AND.(etac < one))THEN
    
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaic=zero
  END IF
  
  CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdai,dwi,mf0i,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdai,dwi,ddwi,mf0i,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      pkfic(j1,k1)=pkfic(j1,k1)+aux*sfilfic(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
END DO

RETURN
END SUBROUTINE affmatclnetfic
