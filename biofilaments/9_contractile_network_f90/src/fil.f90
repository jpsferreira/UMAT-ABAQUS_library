SUBROUTINE fil(f,ff,dw,ddw,lambda,lambda0,ll,r0,mu0,beta,b0)



!>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: f
DOUBLE PRECISION, INTENT(OUT)            :: ff
DOUBLE PRECISION, INTENT(OUT)            :: dw
DOUBLE PRECISION, INTENT(OUT)            :: ddw
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: ll
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0




DOUBLE PRECISION :: a,b,machep,t
DOUBLE PRECISION :: aux, pi,alpha
DOUBLE PRECISION :: aux0,aux1,aux2,aux3,aux4,aux5,aux6,y

a=zero
b=1.0E09
machep=2.2204E-16
t=1.0E-6
f=zero

CALL pullforce(f, a, b, machep, t, lambda,lambda0,ll,r0,mu0,beta,b0)

pi=four*ATAN(one)
ff=f*ll*(pi*pi*b0)**(-one)
alpha=pi*pi*b0*(ll*ll*mu0)**(-one)

aux0=beta/alpha
aux=alpha*ff
aux1=one+ff+aux*ff
aux2=one+two*aux
aux3=one+aux
aux4=lambda0*r0*r0*mu0*(ll**(-one))
aux5=((one+aux)*(aux1**(-one)))**beta
aux6=one-r0*((ll)**(-one))

y=aux0*(aux2*aux2*(aux1**(-one)))-beta*(aux2*(aux3**(-one)))-two

dw=lambda0*r0*f
ddw=aux4*((one+y*aux5*aux6)**(-one))

RETURN
END SUBROUTINE fil
