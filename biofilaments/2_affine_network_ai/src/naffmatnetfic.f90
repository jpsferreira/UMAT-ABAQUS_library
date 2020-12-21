SUBROUTINE naffmatnetfic(pkfic,cmfic,f,mf0,rw,filprops,naffprops,  &
        det,ndi)

!>    NON-AFFINE NETWORK:'FICTICIOUS' PK2 STRESS AND ELASTICITY TNSR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(6)
DOUBLE PRECISION, INTENT(IN)             :: naffprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: h(ndi,ndi),hh(ndi,ndi,ndi,ndi),  &
    hi(ndi,ndi),hhi(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dw,ddw,rwi,aux1,aux2,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda,lambda0,n,pp

!     FILAMENT
l      = filprops(1)
r0      = filprops(2)
mu0     = filprops(3)
beta    = filprops(4)
b0      = filprops(5)
lambda0 = filprops(6)
!     NETWORK
n       = naffprops(1)
pp      = naffprops(2)

lambda=zero
h=zero
hh=zero
hi=zero
hhi=zero
aux=n*(det**(-one))

DO i1=1,nwp
  
  mfi=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  lambda=lambda+(lambdai**pp)*rwi
  
  CALL hfilfic(hi,hhi,pp,lambdai,mf0i,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      h(j1,k1)=h(j1,k1)+hi(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          hh(j1,k1,l1,m1)=hh(j1,k1,l1,m1)+hhi(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
END DO

lambda=lambda**(pp**(-one))

CALL fil(fi,ffi,dw,ddw,lambda,lambda0,l,r0,mu0,beta,b0)

pkfic=zero
cmfic=zero
aux1=aux*dw*lambda**(one-pp)
aux2=aux*(ddw*(lambda**(two*(one-pp)))- (pp-one)*dw*(lambda**(one-two*pp)))

DO j1=1,ndi
  DO k1=1,ndi
    pkfic(j1,k1)=aux1*h(j1,k1)
    DO l1=1,ndi
      DO m1=1,ndi
        cmfic(j1,k1,l1,m1)=aux1*hh(j1,k1,l1,m1)+aux2*h(j1,k1)*h(l1,m1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE naffmatnetfic
