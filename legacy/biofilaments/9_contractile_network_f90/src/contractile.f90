SUBROUTINE contractile(fc,ffc,dw,ddw,ffcc,ru,ru0,lambdac,lambdap,  &
        lambda0,l,r0,mu0,beta,b0,ffcmax,fric,frac,dtime)


use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(IN OUT)         :: fc
DOUBLE PRECISION, INTENT(IN OUT)         :: ffc
DOUBLE PRECISION, INTENT(IN OUT)         :: dw
DOUBLE PRECISION, INTENT(IN OUT)         :: ddw
DOUBLE PRECISION, INTENT(IN OUT)         :: ffcc
DOUBLE PRECISION, INTENT(IN OUT)         :: ru
DOUBLE PRECISION, INTENT(IN OUT)             :: ru0
DOUBLE PRECISION, INTENT(OUT)            :: lambdac
DOUBLE PRECISION, INTENT(IN OUT)             :: lambdap
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: l
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION, INTENT(IN OUT)         :: ffcmax
DOUBLE PRECISION, INTENT(IN OUT)         :: fric
DOUBLE PRECISION, INTENT(IN OUT)         :: frac(nch)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime


!      CHECK STRETCH
IF (ru0 > zero) THEN
  lambdac=lambdap+ru0
ELSE
  lambdac=lambdap
END IF
!        IF(LAMBDAC.LT.LAMBDAP)THEN
!          LAMBDAC=LAMBDAP
!           RU=ZERO
!        ENDIF
!        FORCE WITH CONTRACTION
CALL fil(fc,ffc,dw,ddw,lambdac,lambda0,l,r0,mu0,beta,b0)
!        PASSIVE FORCE
!          CALL FIL(FC,FFC,DW,DDW,LAMBDAP,LAMBDA0,L,R0,MU0,BETA,B0)
!        RELATIVE SLIDING
CALL sliding(ffcc,ru,ffc,ru0,ffcmax,fric,frac,dtime)

RETURN

END SUBROUTINE contractile

