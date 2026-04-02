SUBROUTINE tensorprod2(a,b,c,ndi)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: b(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: c(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        c(i,j,k,l)=a(i,j)*b(k,l)
      END DO
    END DO
  END DO
END DO

RETURN

END SUBROUTINE tensorprod2
