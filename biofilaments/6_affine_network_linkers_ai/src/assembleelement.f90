SUBROUTINE assembleelement(ndim,nnode,ndofel,  &
        ru,rc,kuu,kuc,kcu,kcc,  &
        rhs,amatrx)

!
! Subroutine to assemble the local elements residual and tangent
!
use global

INTEGER, INTENT(IN)                      :: ndim
INTEGER, INTENT(IN)                      :: nnode
INTEGER, INTENT(IN)                      :: ndofel
real(8), INTENT(IN)                       :: ru(ndim*nnode,1)
real(8), INTENT(IN)                       :: rc(nnode,1)
real(8), INTENT(IN)                       :: kuu(ndim*nnode,ndim*nnode)
real(8), INTENT(IN)                       :: kuc(ndim*nnode,nnode)
real(8), INTENT(IN)                       :: kcu(nnode,ndim*nnode)
real(8), INTENT(IN)                       :: kcc(nnode,nnode)
real(8), INTENT(OUT)                      :: rhs(ndofel,1)
real(8), INTENT(OUT)                      :: amatrx(ndofel,ndofel)


INTEGER :: i,j,k,l,m,n,a11,a12,b11,b12, ndofn

! Total number of degrees of freedom per node
!
ndofn = ndofel/nnode
! init
!
rhs = 0.d0
amatrx = 0.d0
IF(ndim == 2) THEN
!
! Assemble the element level residual
!
  DO i=1,nnode
    a11 = ndofn*(i-1)+1
    a12 = ndim*(i-1)+1
!
! displacement
!
    rhs(a11,1) = ru(a12,1)
    rhs(a11+1,1) = ru(a12+1,1)
!
! chemical potential
!
    rhs(a11+2,1) = rc(i,1)
!            ! temperature
!            !
!            rhs(A11+3,1) = Rt(i,1)
  END DO
!
! Assemble the element level tangent matrix
!
  DO i=1,nnode
    DO j=1,nnode
      a11 = ndofn*(i-1)+1
      a12 = ndim*(i-1)+1
      b11 = ndofn*(j-1)+1
      b12 = ndim*(j-1)+1
!
! displacement
!
      amatrx(a11,b11) = kuu(a12,b12)
      amatrx(a11,b11+1) = kuu(a12,b12+1)
      amatrx(a11+1,b11) = kuu(a12+1,b12)
      amatrx(a11+1,b11+1) = kuu(a12+1,b12+1)
!
! chemical potential
!
      amatrx(a11+2,b11+2) = kcc(i,j)
!               ! temperature
!               !
!               amatrx(A11+3,B11+3) = Ktt(i,j)
!
! displacement - chemical potential
!
      amatrx(a11,b11+2) = kuc(a12,j)
      amatrx(a11+1,b11+2) = kuc(a12+1,j)
!               ! displacement - temperature
!               !
!               amatrx(A11,B11+3) = Kut(A12,j)
!               amatrx(A11+1,B11+3) = Kut(A12+1,j)
!
! chemical potential - displacement
!
      amatrx(a11+2,b11) = kcu(i,b12)
      amatrx(a11+2,b11+1) = kcu(i,b12+1)
!
!               ! chemical potential - temperature
!               !
!               amatrx(A11+2,B11+3) = Kct(i,j)
!               !
!               ! temperature - displacement
!               !
!               amatrx(A11+3,B11) = Ktu(i,B12)
!               amatrx(A11+3,B11+1) = Ktu(i,B12+1)
!               !
!               ! temperature - chemical potential
!               !
!               amatrx(A11+3,B11+2) = Ktc(i,j)
    END DO
  END DO
!
ELSE IF(ndim == 3) THEN
!
! Assemble the element level residual
!
  DO i=1,nnode
    a11 = ndofn*(i-1)+1
    a12 = ndim*(i-1)+1
!
! displacement
!
    rhs(a11,1)   = ru(a12,1)
    rhs(a11+1,1) = ru(a12+1,1)
    rhs(a11+2,1) = ru(a12+2,1)
!
! chemical potential
!
    rhs(a11+3,1) = rc(i,1)
!
!             ! temperature
!            !
!            rhs(A11+4,1) = Rt(i,1)
  END DO
!
! Assembly the element level tangent matrix
!
  DO i=1,nnode
    DO j=1,nnode
      a11 = ndofn*(i-1)+1
      a12 = ndim*(i-1)+1
      b11 = ndofn*(j-1)+1
      b12 = ndim*(j-1)+1
!
! displacement
!
      amatrx(a11,b11)     = kuu(a12,b12)
      amatrx(a11,b11+1)   = kuu(a12,b12+1)
      amatrx(a11,b11+2)   = kuu(a12,b12+2)
      amatrx(a11+1,b11)   = kuu(a12+1,b12)
      amatrx(a11+1,b11+1) = kuu(a12+1,b12+1)
      amatrx(a11+1,b11+2) = kuu(a12+1,b12+2)
      amatrx(a11+2,b11)   = kuu(a12+2,b12)
      amatrx(a11+2,b11+1) = kuu(a12+2,b12+1)
      amatrx(a11+2,b11+2) = kuu(a12+2,b12+2)
!
! chemical potential
!
      amatrx(a11+3,b11+3) = kcc(i,j)
!                ! temperature
!               !
!               amatrx(A11+4,B11+4) = Ktt(i,j)
!
! displacement - chemical potential
!
      amatrx(a11,b11+3) = kuc(a12,j)
      amatrx(a11+1,b11+3) = kuc(a12+1,j)
      amatrx(a11+2,b11+3) = kuc(a12+2,j)
      
!               ! displacement - temperature
!               !
!               amatrx(A11,B11+4) = Kut(A12,j)
!               amatrx(A11+1,B11+4) = Kut(A12+1,j)
!               amatrx(A11+2,B11+4) = Kut(A12+2,j)
!
! chemical potential - displacement
!
      amatrx(a11+3,b11) = kcu(i,b12)
      amatrx(a11+3,b11+1) = kcu(i,b12+1)
      amatrx(a11+3,b11+2) = kcu(i,b12+2)
!
!
!               ! chemical potential - temperature
!               !
!               amatrx(A11+3,B11+4) = Kct(i,j)
!               !
!               ! temperature - displacement
!               !
!               amatrx(A11+4,B11) = Ktu(i,B12)
!               amatrx(A11+4,B11+1) = Ktu(i,B12+1)
!               amatrx(A11+4,B11+2) = Ktu(i,B12+2)
!               !
!               ! temperature - chemical potential
!               !
!               amatrx(A11+4,B11+3) = Ktc(i,j)
    END DO
  END DO
!
ELSE
  WRITE(*,*) 'How did you get nDim=',ndim
  CALL xit()
END IF

RETURN
END SUBROUTINE assembleelement
