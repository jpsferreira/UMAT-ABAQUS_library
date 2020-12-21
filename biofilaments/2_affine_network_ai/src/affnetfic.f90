SUBROUTINE affnetfic(sfic,cfic,f,mf0,rw,filprops,affprops,ru0,  &
        dtime,frac,efi,noel,vara,dirmax,det,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(5)
DOUBLE PRECISION, INTENT(IN OUT)         :: ru0(nwp)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: frac(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(OUT)            :: vara
DOUBLE PRECISION, INTENT(OUT)            :: dirmax(ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1, im1
DOUBLE PRECISION :: sfilfic(ndi,ndi),   &
     cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda0,rho,m,fi,ffi
DOUBLE PRECISION :: r0f,r0c,etac,lambdaif,lambdaicl
DOUBLE PRECISION :: b,fric,ffmax,ang, ru
DOUBLE PRECISION :: avga,maxa,aux0,ffic,suma,rho0

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
m       = affprops(2)
b       = affprops(3)
fric    = affprops(4)
ffmax   = affprops(5)

pi=four*ATAN(one)
aux=m*(det**(-one))*four*pi
cfic=zero
sfic=zero

rho=one
r0=r0f+r0c

avga=zero
maxa=zero
suma=zero
dirmax=zero
!       CALL DENSITY(RHO0,ZERO,B,EFI)

!             OPEN (UNIT=20,FILE="projfil.out",action="write",
!     1 status="replace")

!        LOOP OVER THE INTEGRATION DIRECTIONS
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  IF((etac > zero).AND.(etac < one))THEN
    
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaicl=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaicl=zero
  END IF
  
  ru=ru0(i1)
  CALL contractile(fi,ffi,dwi,ddwi,ffic,ru,ru,lambdaic,lambdaif,  &
      lambda0,l,r0,mu0,beta,b0,ffmax,fric,frac,dtime)
  ru0(i1)=ru
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  
!        AUX0=(FFIC/FFMAX)*(RHO)
  aux0=ru*rho
  
  IF (ru > zero) THEN
!        AVERAGE CONTRACTION LEVEL
    avga=avga+aux0
    suma=suma+one
    IF (aux0 > maxa) THEN
!        MAXIMUM CONTRACTION LEVEL
      maxa = aux0
      dirmax=mfi
      im1=i1
    END IF
  END IF
  
  CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,rwi,ndi)
  
  
  IF(ru > zero)THEN
    
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
    
  END IF
  
END DO
!      close(20)

IF (suma > zero) THEN
  avga=avga/nwp
END IF
vara=(maxa-avga)*((maxa)**(-one))
!        WRITE(*,*) VARA,MAXA,IM1,SUMA

RETURN
END SUBROUTINE affnetfic
