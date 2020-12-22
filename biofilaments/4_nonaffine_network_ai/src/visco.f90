SUBROUTINE visco(pk,cmat,vv,pkvol,pkiso,cmatvol,cmatiso,dtime,  &
        vscprops,statev,ndi)



!>    VISCOUS DISSIPATION: MAXWELL SPRINGS AND DASHPOTS SCHEME
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cmat(ndi,ndi,ndi,ndi)
INTEGER, INTENT(IN)                      :: vv
DOUBLE PRECISION, INTENT(IN)             :: pkvol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cmatvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmatiso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN)             :: vscprops(6)
DOUBLE PRECISION, INTENT(IN OUT)         :: statev(nsdv)



INTEGER :: i1,j1,k1,l1, v1

DOUBLE PRECISION :: q(ndi,ndi),qv(ndi,ndi),hv(ndi,ndi), hv0(ndi,ndi)
DOUBLE PRECISION :: teta,tau,aux,auxc

q=zero
qv=zero
hv=zero
auxc=zero

!     ( GENERAL MAXWELL DASHPOTS)
DO v1=1,vv
  
  tau=vscprops(2*v1-1)
  teta=vscprops(2*v1)
  
!      READ STATE VARIABLES
  CALL hvread(hv,statev,v1,ndi)
  hv0=hv
!        RALAXATION TENSORS
  CALL relax(qv,hv,aux,hv0,pkiso,dtime,tau,teta,ndi)
  auxc=auxc+aux
!        WRITE STATE VARIABLES
  CALL hvwrite(statev,hv,v1,ndi)
  
  q=q+qv
  
END DO

auxc=one+auxc
pk=pkvol+pkiso


DO i1=1,ndi
  DO j1=1,ndi
    pk(i1,j1)=pk(i1,j1)+q(i1,j1)
    DO k1=1,ndi
      DO l1=1,ndi
        cmat(i1,j1,k1,l1)= cmatvol(i1,j1,k1,l1)+ auxc*cmatiso(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO



RETURN
END SUBROUTINE visco
