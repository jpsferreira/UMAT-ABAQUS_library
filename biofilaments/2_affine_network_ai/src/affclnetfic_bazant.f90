SUBROUTINE affclnetfic_bazant(sfic,cfic,f,mf0,rw,filprops,affprops,  &
        efi,noel,det,ndi)

!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> BAZANT SPHERICAL INTEGRATION SCHEME
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,k1,l1,m1, im1
DOUBLE PRECISION :: sfilfic(ndi,ndi), cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,r0f,r0,mu0,b0,beta,lambda0,rho,n,fi,ffi,dtime
DOUBLE PRECISION :: r0c,etac,lambdaif
DOUBLE PRECISION :: b,fric,ffmax,ang, frac(4),ru0(nwp),ru
DOUBLE PRECISION :: vara,avga,maxa,aux0,ffic,suma,rho0,dirmax(ndi)
DOUBLE PRECISION :: aa,rr
INTEGER          :: nn,factor
real ( kind = 8 ), external :: polyterm_value_3d
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
n       = affprops(1)
b       = affprops(2)

pi=four*ATAN(one)
aux=n*(det**(-one))*four*pi
cfic=zero
sfic=zero

rho=one
r0=r0f+r0c

aa = zero
avga=zero
maxa=zero
suma=zero
dirmax=zero
!       CALL DENSITY(RHO0,ZERO,B,EFI)

!             OPEN (UNIT=20,FILE="projfil.out",action="write",
!     1 status="replace")

!        LOOP OVER THE INTEGRATION DIRECTIONS (BAZANT INTEGRATION SCHEME)
!bazant integration scheme
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  rho=one
  
  IF((etac > zero).AND.(etac < one))THEN
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaic=zero
  END IF
  
  CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,rwi,ndi)
  
  
  DO j1=1,ndi
    DO k1=1,ndi
      sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
  !bazant integration scheme
 ! aa = aa +rwi*dwi
END DO !end of discretization scheme

!!discrete angular integration test
!factor = 4
!nn = 0
!rr = zero
!call sphere01_quad_icos1c ( factor, polyterm_value_3d, nn, rr )
!write(*,*) 'bazant', four*pi*aa

RETURN
END SUBROUTINE affclnetfic_bazant
