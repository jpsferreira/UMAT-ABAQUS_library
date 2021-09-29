!****************************************************************************



!     Utility subroutines
!****************************************************************************
!***************************************************************************

SUBROUTINE matinv3dd(a,a_inv,det_a,istat)
!
! Returns A_inv, the inverse and det_A, the determinant
! Note that the det is of the original matrix, not the
! inverse
!
use global

real(8), INTENT(IN)                       :: a(3,3)
real(8), INTENT(OUT)                      :: a_inv(3,3)
real(8), INTENT(OUT)                      :: det_a
INTEGER, INTENT(OUT)                     :: istat

!

!
real(8)  det_a_inv
!
istat = 1

det_a = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -  &
    a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) +  &
    a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3))

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: subroutine matInv3Dd:'
  WRITE(*,*) 'WARNING: det of mat=',det_a
  istat = 0
  RETURN
END IF
!
det_a_inv = 1.d0/det_a
!
a_inv(1,1) = det_a_inv*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
a_inv(1,2) = det_a_inv*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
a_inv(1,3) = det_a_inv*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
a_inv(2,1) = det_a_inv*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
a_inv(2,2) = det_a_inv*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
a_inv(2,3) = det_a_inv*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
a_inv(3,1) = det_a_inv*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
a_inv(3,2) = det_a_inv*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
a_inv(3,3) = det_a_inv*(a(1,1)*a(2,2)-a(2,1)*a(1,2))
!
RETURN
END SUBROUTINE matinv3dd
!!!

SUBROUTINE matinv2d(a,a_inv,det_a,istat)
!
! Returns A_inv, the inverse, and det_A, the determinant
! Note that the det is of the original matrix, not the
! inverse
!
use global

real(8), INTENT(IN)                       :: a(2,2)
real(8), INTENT(OUT)                      :: a_inv(2,2)
real(8), INTENT(OUT)                      :: det_a
INTEGER, INTENT(OUT)                     :: istat

!

!
real(8)  det_a_inv


istat = 1

det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: subroutine matInv2D:'
  WRITE(*,*) 'WARNING: det of mat=',det_a
  istat = 0
  RETURN
END IF

det_a_inv = 1.d0/det_a

a_inv(1,1) =  det_a_inv*a(2,2)
a_inv(1,2) = -det_a_inv*a(1,2)
a_inv(2,1) = -det_a_inv*a(2,1)
a_inv(2,2) =  det_a_inv*a(1,1)


RETURN
END SUBROUTINE matinv2d

!****************************************************************************

SUBROUTINE mdet(a,det)
!
! This subroutine calculates the determinant
! of a 3 by 3 matrix [A]
!
use global

real(8), INTENT(IN)                       :: a(3,3)
real(8), INTENT(OUT)                      :: det

!



det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)  &
    + a(1,3)*a(2,1)*a(3,2) - a(3,1)*a(2,2)*a(1,3)  &
    - a(3,2)*a(2,3)*a(1,1) - a(3,3)*a(2,1)*a(1,2)


RETURN
END SUBROUTINE mdet

!****************************************************************************

SUBROUTINE onem0(a)
!
! This subroutine stores the identity matrix in the
! 3 by 3 matrix [A]
!
use global

real(8), INTENT(OUT)                      :: a(3,3)

!
INTEGER :: i,j
!



DO i=1,3
  DO j=1,3
    IF (i == j) THEN
      a(i,j) = 1.0
    ELSE
      a(i,j) = 0.0
    END IF
  END DO
END DO


RETURN
END SUBROUTINE onem0
!***************************************************************************
