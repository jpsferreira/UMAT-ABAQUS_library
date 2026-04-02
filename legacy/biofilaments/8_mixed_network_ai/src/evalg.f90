SUBROUTINE evalg(g,f,lambda,lambda0,l,r0,mu0,beta,b0)



!>     ESTABLISHMENT OF G(F)=LHS-RHS(F) THAT RELATES
!>       STRETCHFORCE RELATIONSHIP OF A SINGLE EXNTESIBLE FILAMENT
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: g
DOUBLE PRECISION, INTENT(IN)             :: f
DOUBLE PRECISION, INTENT(IN)             :: lambda
DOUBLE PRECISION, INTENT(IN)             :: lambda0
DOUBLE PRECISION, INTENT(IN)             :: l
DOUBLE PRECISION, INTENT(IN)             :: r0
DOUBLE PRECISION, INTENT(IN)             :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0


DOUBLE PRECISION :: lhs,rhs


DOUBLE PRECISION :: aux0,aux1,aux,aux2,aux3,aux4,pi

pi=four*ATAN(one)
aux0=one-r0/l
aux1=l*l*((pi*pi*b0)**(-one))
aux=f/mu0
aux2=one+aux
aux3=one+two*aux
aux4=one+f*aux1+f*aux*aux1

rhs=one+aux-aux0*(aux2**beta)*aux3*(aux4**(-beta))
lhs=lambda*lambda0*r0*(l**(-one))

g=lhs-rhs

RETURN
END SUBROUTINE evalg
