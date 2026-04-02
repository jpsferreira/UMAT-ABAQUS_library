SUBROUTINE pullforce(zero0, a, b, machep, t,  &
        lambda,lambda0,l,r0,mu0,beta,b0)



!>    SINGLE FILAMENT: COMPUTES PULLING FORCE FOR A GIVEN STRETCH
!*********************************************************************72

!     ZERO SEEKS THE ROOT OF A FUNCTION F(X) IN AN INTERVAL [A,B].

!     DISCUSSION:

!     THE INTERVAL [A,B] MUST BE A CHANGE OF SIGN INTERVAL FOR F.
!     THAT IS, F(A) AND F(B) MUST BE OF OPPOSITE SIGNS.  THEN
!     ASSUMING THAT F IS CONTINUOUS IMPLIES THE EXISTENCE OF AT LEAST
!     ONE VALUE C BETWEEN A AND B FOR WHICH F(C) = 0.

!     THE LOCATION OF THE ZERO IS DETERMINED TO WITHIN AN ACCURACY
!     OF 6 * MACHEPS * ABS ( C ) + 2 * T.


!     LICENSING:

!     THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.

!     MODIFIED:

!     11 FEBRUARY 2013

!     AUTHOR:

!     RICHARD BRENT
!     MODIFICATIONS BY JOHN BURKARDT

!     REFERENCE:

!     RICHARD BRENT,
!     ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES,
!     DOVER, 2002,
!     ISBN: 0-486-41998-3,
!     LC: QA402.5.B74.

!     PARAMETERS:

!     INPUT, DOUBLE PRECISION A, B, THE ENDPOINTS OF THE CHANGE OF SIGN
!     INTERVAL.
!     INPUT, DOUBLE PRECISION MACHEP, AN ESTIMATE FOR THE RELATIVE
!     MACHINE PRECISION.

!     INPUT, DOUBLE PRECISION T, A POSITIVE ERROR TOLERANCE.

!     INPUT, EXTERNAL DOUBLE PRECISION F, THE NAME OF A USER-SUPPLIED
!     FUNCTION, OF THE FORM "FUNCTION G ( F )", WHICH EVALUATES THE
!     FUNCTION WHOSE ZERO IS BEING SOUGHT.

!     OUTPUT, DOUBLE PRECISION ZERO, THE ESTIMATED VALUE OF A ZERO OF
!     THE FUNCTION G.
use global

DOUBLE PRECISION, INTENT(OUT)            :: zero0
DOUBLE PRECISION, INTENT(IN)             :: a
DOUBLE PRECISION, INTENT(IN)             :: b
DOUBLE PRECISION, INTENT(IN)             :: machep
DOUBLE PRECISION, INTENT(IN)             :: t
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: l
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION :: c
DOUBLE PRECISION :: d
DOUBLE PRECISION :: e
DOUBLE PRECISION :: fa
DOUBLE PRECISION :: fb
DOUBLE PRECISION :: fc
DOUBLE PRECISION :: m

DOUBLE PRECISION :: p
DOUBLE PRECISION :: q
DOUBLE PRECISION :: r
DOUBLE PRECISION :: s
DOUBLE PRECISION :: sa
DOUBLE PRECISION :: sb

DOUBLE PRECISION :: tol




!     MAKE LOCAL COPIES OF A AND B.

sa = a
sb = b
CALL evalg(fa,sa,lambda,lambda0,l,r0,mu0,beta,b0)
CALL evalg(fb,sb,lambda,lambda0,l,r0,mu0,beta,b0)
!      FA = F ( SA )
!      FB = F ( SB )

10    CONTINUE

c = sa
fc = fa
e = sb - sa
d = e

20    CONTINUE

IF ( ABS ( fc ) < ABS ( fb ) ) THEN
  sa = sb
  sb = c
  c = sa
  fa = fb
  fb = fc
  fc = fa
END IF

30    CONTINUE

tol = 2.0D+00 * machep * ABS ( sb ) + t
m = 0.5D+00 * ( c - sb )
IF ( ABS ( m ) <= tol .OR. fb == 0.0D+00 ) GO TO 140
IF ( ABS ( e ) >= tol .AND. ABS ( fa ) > ABS ( fb ) ) GO TO 40

e = m
d = e
GO TO 100

40    CONTINUE

s = fb / fa
IF ( sa /= c ) GO TO 50

p = 2.0D+00 * m * s
q = 1.0D+00 - s
GO TO 60

50    CONTINUE

q = fa / fc
r = fb / fc
p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

60    CONTINUE

IF ( p <= 0.0D+00 ) GO TO 70

q = - q
GO TO 80

70    CONTINUE

p = - p

80    CONTINUE

s = e
e = d
IF ( 2.0D+00 * p >= 3.0D+00 * m * q - ABS ( tol * q ) .OR.  &
    p >= ABS ( 0.5D+00 * s * q ) ) GO TO 90

d = p / q
GO TO 100

90    CONTINUE

e = m
d = e

100   CONTINUE

sa = sb
fa = fb
IF ( ABS ( d ) <= tol ) GO TO 110
sb = sb + d
GO TO 130

110   CONTINUE

IF ( m <= 0.0D+00 ) GO TO 120
sb = sb + tol
GO TO 130

120   CONTINUE

sb = sb - tol

130   CONTINUE

!      FB = F ( SB )
CALL evalg(fb,sb,lambda,lambda0,l,r0,mu0,beta,b0)
IF ( fb > 0.0D+00 .AND. fc > 0.0D+00 ) GO TO 10
IF ( fb <= 0.0D+00 .AND. fc <= 0.0D+00 ) GO TO 10
GO TO 20

140   CONTINUE

zero0 = sb

RETURN
END SUBROUTINE pullforce

!*********************************************************************72
