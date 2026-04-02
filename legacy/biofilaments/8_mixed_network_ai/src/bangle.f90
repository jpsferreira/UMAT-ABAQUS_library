SUBROUTINE bangle(ang,f,mf,noel,pdir,ndi)

!>    ANGLE BETWEEN FILAMENT AND PREFERED DIRECTION

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ang
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf(ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pdir(ndi)
INTEGER, INTENT(IN OUT)                  :: noel
!
!
INTEGER :: inoel,i,j
DOUBLE PRECISION :: dnorm, mfa(ndi),aux
DOUBLE PRECISION :: c(ndi,ndi),egvc(ndi,ndi),egvl(ndi)
!
inoel=0
i=0
!DO i=1,nelem
!               ELEMENT IDENTIFICATION
!  IF(noel == INT(prefdir(i,1))) THEN
!    inoel=i
!  END IF
!END DO
!
!DO i=1,ndi
!  j=i+1
!       PREFERED ORIENTATION  ORIENTATION NORMALIZED
!  pdir(i)=prefdir(inoel,j)
!END DO
!        ALTERNATIVE APPROACH: BUNDLES FOLLOW PRINCIPAL DIRECTIONS
!c=matmul(transpose(f),f)
!CALL spectral(c,egvl,egvc)
!       WRITE(*,*) EGVC
!pdir(1)=egvc(1,1)
!pdir(2)=egvc(2,1)
!pdir(3)=egvc(3,1)
!        END OF ALTERNATIVE

!     PREFERED ORIENTATION
dnorm=dot_product(pdir,pdir)
dnorm=DSQRT(dnorm)
!     PREFERED ORIENTATION  NORMALIZED
pdir=pdir/dnorm

!       FILAMENT ORIENTATION
mfa=mf
dnorm=dot_product(mfa,mfa)
dnorm=dsqrt(dnorm)

!       FILAMENT ORIENTATION  NORMALIZED
mfa=mfa/dnorm
!        ANGLE BETWEEN PREFERED ORIENTATION AND FILAMENT - BANGLE
aux=dot_product(mfa,pdir)
!        if AUX.GT.ONE
!        endif
!        write(*,*) aux
ang=acos(aux)

RETURN
END SUBROUTINE bangle

