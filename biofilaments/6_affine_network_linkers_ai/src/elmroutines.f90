!****************************************************************************



!     Element subroutines
!****************************************************************************

SUBROUTINE xint2d1pt(xi,w,ninttpt)
!
! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 1 gauss point for integration
!
!  xi(nInttPt,2): xi,eta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
!
use global

real(8), INTENT(OUT)                      :: xi(1,2)
real(8), INTENT(OUT)                      :: w(1)
INTEGER, INTENT(OUT)                     :: ninttpt

!
INTEGER :: ndim
!


! Initialize
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 1


! Gauss weights
!
w = 4.d0


! Gauss pt location in master element
!
xi(1,1) = 0.d0
xi(1,2) = 0.d0


RETURN
END SUBROUTINE xint2d1pt

!************************************************************************

SUBROUTINE xint2d4pt(xi,w,ninttpt)
!
! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 4 gauss points for integration
!
!  xi(nInttPt,2): xi,eta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
!
use global

real(8), INTENT(OUT)                      :: xi(4,2)
real(8), INTENT(OUT)                      :: w(4)
INTEGER, INTENT(OUT)                     :: ninttpt

!
INTEGER :: ndim
!



! Initialize
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 4


! Gauss weights
!
w(1) = 1.d0
w(2) = 1.d0
w(3) = 1.d0
w(4) = 1.d0


! Gauss pt locations in master element
!
xi(1,1) = -DSQRT(1.d0/3.d0)
xi(1,2) = -DSQRT(1.d0/3.d0)
xi(2,1) = DSQRT(1.d0/3.d0)
xi(2,2) = -DSQRT(1.d0/3.d0)
xi(3,1) = -DSQRT(1.d0/3.d0)
xi(3,2) = DSQRT(1.d0/3.d0)
xi(4,1) = DSQRT(1.d0/3.d0)
xi(4,2) = DSQRT(1.d0/3.d0)


RETURN
END SUBROUTINE xint2d4pt

!***********************************************************************

SUBROUTINE xint3d1pt(xi,w,ninttpt)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using a 2 gauss points for integration
!
!  xi(nInttPt,3): xi,eta,zeta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights

use global

real(8), INTENT(OUT)                      :: xi(1,3)
real(8), INTENT(OUT)                      :: w(1)
INTEGER, INTENT(OUT)                     :: ninttpt


INTEGER :: ndim




! Init
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 1


! Gauss weights
!
w(1) = 8.d0


! Gauss pt locations in master element
!
xi(1,1) = 0.d0
xi(1,2) = 0.d0
xi(1,3) = 0.d0

RETURN
END SUBROUTINE xint3d1pt

!***********************************************************************

SUBROUTINE xint3d8pt(xi,w,ninttpt)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 8 gauss points for integration
!
!  xi(nInttPt,3): xi,eta,zeta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

real(8), INTENT(OUT)                      :: xi(8,3)
real(8), INTENT(OUT)                      :: w(8)
INTEGER, INTENT(OUT)                     :: ninttpt


INTEGER :: ndim




! Init
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 8


! Gauss weights
!
w(1) = 1.d0
w(2) = 1.d0
w(3) = 1.d0
w(4) = 1.d0
w(5) = 1.d0
w(6) = 1.d0
w(7) = 1.d0
w(8) = 1.d0


! Gauss pt locations in master element
!
xi(1,1) = -DSQRT(1.d0/3.d0)
xi(1,2) = -DSQRT(1.d0/3.d0)
xi(1,3) = -DSQRT(1.d0/3.d0)
xi(2,1) = DSQRT(1.d0/3.d0)
xi(2,2) = -DSQRT(1.d0/3.d0)
xi(2,3) = -DSQRT(1.d0/3.d0)
xi(3,1) = -DSQRT(1.d0/3.d0)
xi(3,2) = DSQRT(1.d0/3.d0)
xi(3,3) = -DSQRT(1.d0/3.d0)
xi(4,1) = DSQRT(1.d0/3.d0)
xi(4,2) = DSQRT(1.d0/3.d0)
xi(4,3) = -DSQRT(1.d0/3.d0)
xi(5,1) = -DSQRT(1.d0/3.d0)
xi(5,2) = -DSQRT(1.d0/3.d0)
xi(5,3) = DSQRT(1.d0/3.d0)
xi(6,1) = DSQRT(1.d0/3.d0)
xi(6,2) = -DSQRT(1.d0/3.d0)
xi(6,3) = DSQRT(1.d0/3.d0)
xi(7,1) = -DSQRT(1.d0/3.d0)
xi(7,2) = DSQRT(1.d0/3.d0)
xi(7,3) = DSQRT(1.d0/3.d0)
xi(8,1) = DSQRT(1.d0/3.d0)
xi(8,2) = DSQRT(1.d0/3.d0)
xi(8,3) = DSQRT(1.d0/3.d0)


RETURN
END SUBROUTINE xint3d8pt

!***********************************************************************

SUBROUTINE xintsurf2d1pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(1)
real(8), INTENT(OUT)                      :: ylocal(1)
real(8), INTENT(OUT)                      :: w(1)



! Gauss weights
!
w(1) = two


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = zero
  ylocal(1) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = zero
ELSE IF(face == 3) THEN
  xlocal(1) = zero
  ylocal(1) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = zero
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d1pt

!***********************************************************************

SUBROUTINE xintsurf2d2pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(2)
real(8), INTENT(OUT)                      :: ylocal(2)
real(8), INTENT(OUT)                      :: w(2)



! Gauss weights
!
w(1) = one
w(2) = one


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(one/three)
  xlocal(2) = one
  ylocal(2) = DSQRT(one/three)
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = DSQRT(one/three)
  xlocal(2) = -one
  ylocal(2) = -DSQRT(one/three)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d2pt

!***********************************************************************

SUBROUTINE xintsurf2d3pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(3)
real(8), INTENT(OUT)                      :: ylocal(3)
real(8), INTENT(OUT)                      :: w(3)



real(8), PARAMETER :: five=5.d0
real(8), PARAMETER :: eight=8.d0
real(8), PARAMETER :: nine=9.d0


! Gauss weights
!
w(1) = five/nine
w(2) = eight/nine
w(3) = five/nine


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(three/five)
  ylocal(1) = -one
  xlocal(2) = zero
  ylocal(2) = -one
  xlocal(2) = DSQRT(three/five)
  ylocal(2) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(three/five)
  xlocal(2) = one
  ylocal(2) = zero
  xlocal(3) = one
  ylocal(3) = DSQRT(three/five)
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(three/five)
  ylocal(1) = one
  xlocal(2) = zero
  ylocal(2) = one
  xlocal(3) = DSQRT(three/five)
  ylocal(3) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = DSQRT(three/five)
  xlocal(2) = -one
  ylocal(2) = zero
  xlocal(3) = -one
  ylocal(3) = -DSQRT(three/five)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d3pt

!***********************************************************************

SUBROUTINE xintsurf3d1pt(face,xlocal,ylocal,zlocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 1 gauss point for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  zLocal(nInttPt): z coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(1)
real(8), INTENT(OUT)                      :: ylocal(1)
real(8), INTENT(OUT)                      :: zlocal(1)
real(8), INTENT(OUT)                      :: w(1)




! Gauss weights
!
w(1) = four


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = zero
  ylocal(1) = zero
  zlocal(1) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = zero
  ylocal(1) = zero
  zlocal(1) = one
ELSE IF(face == 3) THEN
  xlocal(1) = zero
  ylocal(1) = -one
  zlocal(1) = zero
ELSE IF(face == 4) THEN
  xlocal(1) = one
  ylocal(1) = zero
  zlocal(1) = zero
ELSE IF(face == 5) THEN
  xlocal(1) = zero
  ylocal(1) = one
  zlocal(1) = zero
ELSE IF(face == 6) THEN
  xlocal(1) = -one
  ylocal(1) = zero
  zlocal(1) = zero
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4,5,6'
! !write(80,*) 'face.ne.1,2,3,4,5,6'
  CALL xit()
END IF

END SUBROUTINE xintsurf3d1pt

!***********************************************************************

SUBROUTINE xintsurf3d4pt(face,xlocal,ylocal,zlocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 4 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  yLocal(nInttPt): z coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(4)
real(8), INTENT(OUT)                      :: ylocal(4)
real(8), INTENT(OUT)                      :: zlocal(4)
real(8), INTENT(OUT)                      :: w(4)


! Gauss weights
!
w(1) = one
w(2) = one
w(3) = one
w(4) = one


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -DSQRT(one/three)
  zlocal(2) = -one
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = -one
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = DSQRT(one/three)
  zlocal(4) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -DSQRT(one/three)
  zlocal(2) = one
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = one
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = DSQRT(one/three)
  zlocal(4) = one
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -one
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -one
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = -one
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = -one
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 4) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = one
  ylocal(2) = DSQRT(one/three)
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = one
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = one
  ylocal(4) = -DSQRT(one/three)
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 5) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = one
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = one
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = one
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = one
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 6) THEN
  xlocal(1) = -one
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = -one
  ylocal(2) = DSQRT(one/three)
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = -one
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -one
  ylocal(4) = -DSQRT(one/three)
  zlocal(4) = DSQRT(one/three)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4,5,6'
! !write(80,*) 'face.ne.1,2,3,4,5,6'
  CALL xit()
END IF

END SUBROUTINE xintsurf3d4pt

!************************************************************************

SUBROUTINE calcshape2dlinear(ninttpt,xi_int,intpt,sh,dshxi)
!
! Calculate the shape functions and their derivatives at the
! given integration point in the master element


! Calculate the shape functions and their derivatives at the
! given integration point in the master element
!
!                          eta
!   4-----------3          |
!   |           |          |
!   |           |          |
!   |           |          |
!   |           |          |
!   |           |          O--------- xi
!   1-----------2        origin at center
!
!
! sh(i) = shape function of node i at the intpt.
! dshxi(i,j) = derivative wrt j direction of shape fn of node i
!
use global
INTEGER, INTENT(IN OUT)                  :: ninttpt
real(8), INTENT(IN)                       :: xi_int(ninttpt,2)
INTEGER, INTENT(IN OUT)                  :: intpt
real(8), INTENT(OUT)                      :: sh(4)
real(8), INTENT(OUT)                      :: dshxi(4,2)

!
INTEGER :: ndim
!
real(8)  xi,eta
!
real(8), PARAMETER :: fourth=1.d0/4.d0


! Location in the master element
!
xi = xi_int(intpt,1)
eta = xi_int(intpt,2)


! The shape functions
!
sh(1) = fourth*(one - xi)*(one - eta)
sh(2) = fourth*(one + xi)*(one - eta)
sh(3) = fourth*(one + xi)*(one + eta)
sh(4) = fourth*(one - xi)*(one + eta)


! The first derivatives
!
dshxi(1,1) = -fourth*(one - eta)
dshxi(1,2) = -fourth*(one - xi)
dshxi(2,1) = fourth*(one - eta)
dshxi(2,2) = -fourth*(one + xi)
dshxi(3,1) = fourth*(one + eta)
dshxi(3,2) = fourth*(one + xi)
dshxi(4,1) = -fourth*(one + eta)
dshxi(4,2) = fourth*(one - xi)


RETURN
END SUBROUTINE calcshape2dlinear

!***********************************************************************

SUBROUTINE calcshape3dlinear(ninttpt,xi_int,intpt,sh,dshxi)
!
!
! Calculate the shape functions and their derivatives at the
! given integration point in the master element
!
! This subroutine uses a 8-node linear 3D element as shown
!
!      8-----------7
!     /|          /|       zeta
!    / |         / |
!   5-----------6  |       |     eta
!   |  |        |  |       |   /
!   |  |        |  |       |  /
!   |  4--------|--3       | /
!   | /         | /        |/
!   |/          |/         O--------- xi
!   1-----------2        origin at cube center
!
!
! sh(i) = shape function of node i at the intpt.
! dshxi(i,j) = derivative wrt j direction of shape fn of node i
! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i
use global

INTEGER, INTENT(IN OUT)                  :: ninttpt
real(8), INTENT(IN)                       :: xi_int(ninttpt,3)
INTEGER, INTENT(IN OUT)                  :: intpt
real(8), INTENT(OUT)                      :: sh(8)
real(8), INTENT(OUT)                      :: dshxi(8,3)


INTEGER :: ndim, i,j


real(8) d2shxi(8,3,3),xi,eta,zeta

real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Location in the master element
!
xi = xi_int(intpt,1)
eta = xi_int(intpt,2)
zeta = xi_int(intpt,3)
!
! The shape functions
!
sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
!
! The first derivatives
!
dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
dshxi(1,3) = -eighth*(one - xi)*(one - eta)
dshxi(2,1) = eighth*(one - eta)*(one - zeta)
dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
dshxi(2,3) = -eighth*(one + xi)*(one - eta)
dshxi(3,1) = eighth*(one + eta)*(one - zeta)
dshxi(3,2) = eighth*(one + xi)*(one - zeta)
dshxi(3,3) = -eighth*(one + xi)*(one + eta)
dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
dshxi(4,2) = eighth*(one - xi)*(one - zeta)
dshxi(4,3) = -eighth*(one - xi)*(one + eta)
dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
dshxi(5,3) = eighth*(one - xi)*(one - eta)
dshxi(6,1) = eighth*(one - eta)*(one + zeta)
dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
dshxi(6,3) = eighth*(one + xi)*(one - eta)
dshxi(7,1) = eighth*(one + eta)*(one + zeta)
dshxi(7,2) = eighth*(one + xi)*(one + zeta)
dshxi(7,3) = eighth*(one + xi)*(one + eta)
dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
dshxi(8,2) = eighth*(one - xi)*(one + zeta)
dshxi(8,3) = eighth*(one - xi)*(one + eta)
!
! The second derivatives
!
d2shxi = zero
d2shxi(1,1,2) = eighth*(one - zeta)
d2shxi(1,2,1) = d2shxi(1,1,2)
d2shxi(1,1,3) = eighth*(one - eta)
d2shxi(1,3,1) = d2shxi(1,1,3)
d2shxi(1,2,3) = eighth*(one - xi)
d2shxi(1,3,2) = d2shxi(1,2,3)
d2shxi(2,1,2) = -eighth*(one - zeta)
d2shxi(2,2,1) = d2shxi(2,1,2)
d2shxi(2,1,3) = -eighth*(one - eta)
d2shxi(2,3,1) = d2shxi(2,1,3)
d2shxi(2,2,3) = eighth*(one + xi)
d2shxi(2,3,2) = d2shxi(2,2,3)
d2shxi(3,1,2) = eighth*(one - zeta)
d2shxi(3,2,1) = d2shxi(2,1,2)
d2shxi(3,1,3) = -eighth*(one + eta)
d2shxi(3,3,1) = d2shxi(2,1,3)
d2shxi(3,2,3) = -eighth*(one + xi)
d2shxi(3,3,2) = d2shxi(2,2,3)
d2shxi(4,1,2) = -eighth*(one - zeta)
d2shxi(4,2,1) = d2shxi(2,1,2)
d2shxi(4,1,3) = eighth*(one + eta)
d2shxi(4,3,1) = d2shxi(2,1,3)
d2shxi(4,2,3) = -eighth*(one - xi)
d2shxi(4,3,2) = d2shxi(2,2,3)
d2shxi(5,1,2) = eighth*(one + zeta)
d2shxi(5,2,1) = d2shxi(2,1,2)
d2shxi(5,1,3) = -eighth*(one - eta)
d2shxi(5,3,1) = d2shxi(2,1,3)
d2shxi(5,2,3) = -eighth*(one - xi)
d2shxi(5,3,2) = d2shxi(2,2,3)
d2shxi(6,1,2) = eighth*(one + zeta)
d2shxi(6,2,1) = d2shxi(2,1,2)
d2shxi(6,1,3) = eighth*(one - eta)
d2shxi(6,3,1) = d2shxi(2,1,3)
d2shxi(6,2,3) = -eighth*(one + xi)
d2shxi(6,3,2) = d2shxi(2,2,3)
d2shxi(7,1,2) = eighth*(one + zeta)
d2shxi(7,2,1) = d2shxi(2,1,2)
d2shxi(7,1,3) = eighth*(one + eta)
d2shxi(7,3,1) = d2shxi(2,1,3)
d2shxi(7,2,3) = eighth*(one + xi)
d2shxi(7,3,2) = d2shxi(2,2,3)
d2shxi(8,1,2) = -eighth*(one + zeta)
d2shxi(8,2,1) = d2shxi(2,1,2)
d2shxi(8,1,3) = -eighth*(one + eta)
d2shxi(8,3,1) = d2shxi(2,1,3)
d2shxi(8,2,3) = eighth*(one - xi)
d2shxi(8,3,2) = d2shxi(2,2,3)

RETURN
END SUBROUTINE calcshape3dlinear

!************************************************************************


SUBROUTINE computesurf(xlocal,ylocal,face,coords,sh,ds)

! This subroutine computes the shape functions, derivatives
!  of shape functions, and the length ds, so that one can
!  do the numerical integration on the boundary for fluxes
!  on the 4-node quadrilateral elements

use global

real(8), INTENT(IN OUT)                   :: xlocal
real(8), INTENT(IN OUT)                   :: ylocal
INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(IN)                       :: coords(2,4)
real(8), INTENT(OUT)                      :: sh(4)
real(8), INTENT(OUT)                      :: ds



real(8)  dshxi(4,2), dxdxi,dxdeta,dydxi
real(8) dydeta, shape,normal(2,1)
real(8), PARAMETER :: fourth=1.d0/4.d0

sh(1) = fourth*(one - xlocal)*(one - ylocal)
sh(2) = fourth*(one + xlocal)*(one - ylocal)
sh(3) = fourth*(one + xlocal)*(one + ylocal)
sh(4) = fourth*(one - xlocal)*(one + ylocal)

dshxi(1,1) = -fourth*(one - ylocal)
dshxi(1,2) = -fourth*(one - xlocal)
dshxi(2,1) = fourth*(one - ylocal)
dshxi(2,2) = -fourth*(one + xlocal)
dshxi(3,1) = fourth*(one + ylocal)
dshxi(3,2) = fourth*(one + xlocal)
dshxi(4,1) = -fourth*(one + ylocal)
dshxi(4,2) = fourth*(one - xlocal)

dxdxi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)  &
    + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
dxdeta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)  &
    + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
dydxi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)  &
    + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
dydeta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)  &
    + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


! Jacobian of the mapping
!
IF((face == 2).OR.(face == 4)) THEN
  ds = DSQRT(dxdeta*dxdeta + dydeta*dydeta)
ELSE IF((face == 1).OR.(face == 3)) THEN
  ds = DSQRT(dxdxi*dxdxi + dydxi*dydxi)
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF


! Surface normal, outward pointing in this case. Useful for
!  ``follower'' type loads. The normal is referential or spatial
!  depending on which coords were supplied to this subroutine
!  (NOT fully tested)
!
IF((face == 2).OR.(face == 4)) THEN
  normal(1,1) = dydeta/DSQRT(dxdeta*dxdeta + dydeta*dydeta)
  normal(2,1) = -dxdeta/DSQRT(dxdeta*dxdeta + dydeta*dydeta)
  IF(face == 4) normal = -normal
ELSE IF((face == 1).OR.(face == 3)) THEN
  normal(1,1) = dydxi/DSQRT(dxdxi*dxdxi + dydxi*dydxi)
  normal(2,1) = -dxdxi/DSQRT(dxdxi*dxdxi + dydxi*dydxi)
  IF(face == 3) normal = -normal
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF

RETURN
END SUBROUTINE computesurf

!***********************************************************************

SUBROUTINE computesurf3d(xlocal,ylocal,zlocal,face,coords,sh,da)

! This subroutine computes the shape functions, derivatives
!  of shape functions, and the area dA, so that one can
!  do the numerical integration on the boundary for fluxes
!  on the 8-node brick elements

use global

real(8), INTENT(IN OUT)                   :: xlocal
real(8), INTENT(IN OUT)                   :: ylocal
real(8), INTENT(IN OUT)                   :: zlocal
INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(IN)                       :: coords(3,8)
real(8), INTENT(OUT)                      :: sh(8)
real(8), INTENT(OUT)                      :: da

INTEGER :: stat,i,j,k

real(8)  dshxi(8,3), dsh(8,3)
real(8)  mapj(3,3),mag,normal(3,1)

real(8) dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta
real(8) dzdzeta


real(8), PARAMETER :: eighth=1.d0/8.d0


! The shape functions
!
sh(1) = eighth*(one - xlocal)*(one - ylocal)*(one - zlocal)
sh(2) = eighth*(one + xlocal)*(one - ylocal)*(one - zlocal)
sh(3) = eighth*(one + xlocal)*(one + ylocal)*(one - zlocal)
sh(4) = eighth*(one - xlocal)*(one + ylocal)*(one - zlocal)
sh(5) = eighth*(one - xlocal)*(one - ylocal)*(one + zlocal)
sh(6) = eighth*(one + xlocal)*(one - ylocal)*(one + zlocal)
sh(7) = eighth*(one + xlocal)*(one + ylocal)*(one + zlocal)
sh(8) = eighth*(one - xlocal)*(one + ylocal)*(one + zlocal)


! Shape function derivatives
!
dshxi(1,1) = -eighth*(one - ylocal)*(one - zlocal)
dshxi(1,2) = -eighth*(one - xlocal)*(one - zlocal)
dshxi(1,3) = -eighth*(one - xlocal)*(one - ylocal)
dshxi(2,1) = eighth*(one - ylocal)*(one - zlocal)
dshxi(2,2) = -eighth*(one + xlocal)*(one - zlocal)
dshxi(2,3) = -eighth*(one + xlocal)*(one - ylocal)
dshxi(3,1) = eighth*(one + ylocal)*(one - zlocal)
dshxi(3,2) = eighth*(one + xlocal)*(one - zlocal)
dshxi(3,3) = -eighth*(one + xlocal)*(one + ylocal)
dshxi(4,1) = -eighth*(one + ylocal)*(one - zlocal)
dshxi(4,2) = eighth*(one - xlocal)*(one - zlocal)
dshxi(4,3) = -eighth*(one - xlocal)*(one + ylocal)
dshxi(5,1) = -eighth*(one - ylocal)*(one + zlocal)
dshxi(5,2) = -eighth*(one - xlocal)*(one + zlocal)
dshxi(5,3) = eighth*(one - xlocal)*(one - ylocal)
dshxi(6,1) = eighth*(one - ylocal)*(one + zlocal)
dshxi(6,2) = -eighth*(one + xlocal)*(one + zlocal)
dshxi(6,3) = eighth*(one + xlocal)*(one - ylocal)
dshxi(7,1) = eighth*(one + ylocal)*(one + zlocal)
dshxi(7,2) = eighth*(one + xlocal)*(one + zlocal)
dshxi(7,3) = eighth*(one + xlocal)*(one + ylocal)
dshxi(8,1) = -eighth*(one + ylocal)*(one + zlocal)
dshxi(8,2) = eighth*(one - xlocal)*(one + zlocal)
dshxi(8,3) = eighth*(one - xlocal)*(one + ylocal)


dxdxi = zero
dxdeta = zero
dxdzeta = zero
dydxi = zero
dydeta = zero
dydzeta = zero
dzdxi = zero
dzdeta = zero
dzdzeta = zero
DO k=1,8
  dxdxi = dxdxi + dshxi(k,1)*coords(1,k)
  dxdeta = dxdeta + dshxi(k,2)*coords(1,k)
  dxdzeta = dxdzeta + dshxi(k,3)*coords(1,k)
  dydxi = dydxi + dshxi(k,1)*coords(2,k)
  dydeta = dydeta + dshxi(k,2)*coords(2,k)
  dydzeta = dydzeta + dshxi(k,3)*coords(2,k)
  dzdxi = dzdxi + dshxi(k,1)*coords(3,k)
  dzdeta = dzdeta + dshxi(k,2)*coords(3,k)
  dzdzeta = dzdzeta + dshxi(k,3)*coords(3,k)
END DO


! Jacobian of the mapping
!
IF((face == 1).OR.(face == 2)) THEN
! zeta = constant on this face
  da = DSQRT( (dydxi*dzdeta - dydeta*dzdxi)**two  &
      + (dxdxi*dzdeta - dxdeta*dzdxi)**two + (dxdxi*dydeta - dxdeta*dydxi)**two  &
      )
ELSE IF((face == 3).OR.(face == 5)) THEN
! eta = constant on this face
  da = DSQRT( (dydxi*dzdzeta - dydzeta*dzdxi)**two  &
      + (dxdxi*dzdzeta - dxdzeta*dzdxi)**two  &
      + (dxdxi*dydzeta - dxdzeta*dydxi)**two )
ELSE IF((face == 4).OR.(face == 6)) THEN
! xi = constant on this face
  da = DSQRT( (dydeta*dzdzeta - dydzeta*dzdeta)**two  &
      + (dxdeta*dzdzeta - dxdzeta*dzdeta)**two  &
      + (dxdeta*dydzeta - dxdzeta*dydeta)**two )
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF


! Surface normal, outward pointing in this case. Useful for
!  ``follower'' type loads. The normal is referential or spatial
!  depending on which coords were supplied to this subroutine
!  (NOT fully tested)
!
IF((face == 1).OR.(face == 2)) THEN
! zeta = constant on this face
  normal(1,1) = dydxi*dzdeta - dydeta*dzdxi
  normal(2,1) = dxdxi*dzdeta - dxdeta*dzdxi
  normal(3,1) = dxdxi*dydeta - dxdeta*dydxi
  IF(face == 1) normal = -normal
ELSE IF((face == 3).OR.(face == 5)) THEN
! eta = constant on this face
  normal(1,1) = dydxi*dzdzeta - dydzeta*dzdxi
  normal(2,1) = dxdxi*dzdzeta - dxdzeta*dzdxi
  normal(3,1) = dxdxi*dydzeta - dxdzeta*dydxi
  IF(face == 5) normal = -normal
ELSE IF((face == 4).OR.(face == 6)) THEN
! xi = constant on this face
  normal(1,1) = dydeta*dzdzeta - dydzeta*dzdeta
  normal(2,1) = dxdeta*dzdzeta - dxdzeta*dzdeta
  normal(3,1) = dxdeta*dydzeta - dxdzeta*dydeta
  IF(face == 6) normal = -normal
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF
mag = DSQRT(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
normal(1,1) = normal(1,1)/mag
normal(2,1) = normal(2,1)/mag
normal(3,1) = normal(3,1)/mag

END SUBROUTINE computesurf3d

!***********************************************************************

SUBROUTINE mapshape2d(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.
!
use global
INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,2)
real(8), INTENT(IN)                       :: coords(3,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,2)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat

!
INTEGER :: i,j,k, ieror
!
real(8)  mapj(2,2), mapj_inv(2,2)
!
real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,2
  DO j=1,2
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv2d(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape2D'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


RETURN
END SUBROUTINE mapshape2d

!*************************************************************************

SUBROUTINE mapshape2da(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.
!
! This subroutine is exactly the same as the regular mapShape2D
!with the exception that coords(2,nNode) here and coords(3,nNode)
!in the regular.  I have noticed that a "heat transfer" and
!"static" step uses MCRD=2,
!but for "coupled-temperature-displacement"
!you will get MCRD=3, even for a plane analysis.
use global

INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,2)
real(8), INTENT(IN)                       :: coords(2,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,2)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat

!
INTEGER :: i,j,k, ieror
!
real(8)  mapj(2,2), mapj_inv(2,2)
!
real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,2
  DO j=1,2
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv2d(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape2Da'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


RETURN
END SUBROUTINE mapshape2da

!***********************************************************************

SUBROUTINE mapshape3d(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.  This subroutine works for both 8-node
!  linear and 20-node quadratic 3D elements.
!
use global
INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,3)
real(8), INTENT(IN)                       :: coords(3,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,3)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat


INTEGER :: i,j,k, ieror


real(8) mapj(3,3),mapj_inv(3,3)

real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,3
  DO j=1,3
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv3dd(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape3D'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


! The second derivatives may be calculated.
!

RETURN
END SUBROUTINE mapshape3d
