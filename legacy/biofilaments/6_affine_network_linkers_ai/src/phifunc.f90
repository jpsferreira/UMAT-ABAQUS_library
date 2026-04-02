SUBROUTINE phifunc(phi,f,df,args,nargs)



! This subroutine serves as the function we would like to solve
! for the polymer volume fraction by finding phi such that ``f=0''
use global

real(8), INTENT(IN)                       :: phi
real(8), INTENT(OUT)                      :: f
real(8), INTENT(OUT)                      :: df
real(8), INTENT(IN)                       :: args(nargs)
INTEGER, INTENT(IN OUT)                  :: nargs


INTEGER :: material
INTEGER, PARAMETER :: neohookean=1
INTEGER, PARAMETER :: langevin=2

real(8)  mu,mu0,rgas,theta,chi,vmol,gshear,kbulk
real(8) detf, rt


! Obtain relevant quantities
!
mu     = args(1)
mu0    = args(2)
rgas   = args(3)
theta  = args(4)
chi    = args(5)
vmol   = args(6)
kbulk  = args(7)
detf   = args(8)


! Compute the useful quantity
!
rt = rgas*theta


! Compute the residual
!
f = (mu0 - mu)/rt + DLOG(one - phi) + phi + chi*phi*phi  &
    - ((kbulk*vmol)/rt)*DLOG(detf*phi)  &
    + ((kbulk*vmol)/(two*rt))*(DLOG(detf*phi)**two)


! Compute the tangent
!
IF(phi > 0.999D0) THEN
  df = zero
ELSE
  df = one - (one/(one - phi)) + two*chi*phi - (kbulk*vmol)/(rt*phi)  &
      + ((kbulk*vmol)/(rt*phi))*DLOG(detf*phi)
END IF


RETURN
END SUBROUTINE phifunc
