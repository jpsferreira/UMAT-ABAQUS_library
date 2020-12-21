SUBROUTINE chempot(phi_t,phi_tau,dpdt,dphidmu,dphidotdmu,mfluid,  &
        dmdmu,dmdj,spucmod,spcumodfac,vmol,mu_tau,theta,  &
        detfe,detfs,detf,cr,kbulk,dffprops,iden,dtime,ndi)



!>    ASDASDASD
use global
IMPLICIT NONE

INTEGER :: nargs,ndi
PARAMETER(nargs=8)
DOUBLE PRECISION :: phi_tau,phi_t,phi_per,phi_m,dtime
DOUBLE PRECISION :: d,mu_tau,mu0,rgas,theta,chi,vmol,kbulk,detf,cr
DOUBLE PRECISION :: dpdt,dpdt_per,dpdt_m,dphidmu,dphidj,dphidotdmu
DOUBLE PRECISION :: detfe,detfs,deltamu,mfluid,dmdmu,dmdj,dmdphi
DOUBLE PRECISION :: spcumodfac(ndi,ndi),spucmod(ndi,ndi),  &
    iden(ndi,ndi),args(nargs)
DOUBLE PRECISION :: dffprops(5),cbari1

!      ! COMPUTE THE POLYMER VOLUME FRACTION
!      !
chi      = dffprops(1)
d        = dffprops(2)
mu0      = dffprops(3)
vmol     = dffprops(4)
rgas     = dffprops(5)
!
args(1)  = mu_tau
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
!
CALL solvephi(phi_tau,args,nargs,phi_t)
! COMPUTE THE ELASTIC VOLUME RATIO, DETFE
!
cr=(one/vmol)*(one*(phi_tau**(-one))-one)
detfe = detf*phi_tau
!
! COMPUTE THE SWELLING VOLUME RATIO, DETFS
!
detfs = one/phi_tau
! COMPUTE THE TIME RATE OF THE POLYMER VOLUME FRACTION USING
!  A FINITE DIFFERENCE IN TIME
!
dpdt = (phi_tau - phi_t)/dtime


! COMPUTE THE DERIVATIVE OF THE POLYMER VOLUME FRACTION WITH
!  RESPECT TO THE CHEMICAL POTENTIAL.  COMPUTED VIA IMPLICIT
!  DIFFERENTIATION ON THE CHEMICAL POTENTIAL EQUATION.
!
dphidmu =  (one/(rgas*theta))/ (  &
    (one/(phi_tau - one)) + one + two*chi*phi_tau  &
    - ((vmol*kbulk)/(rgas*theta*phi_tau))  &
    + ((vmol*kbulk)/(rgas*theta*phi_tau))*DLOG(detf*phi_tau) )
! COMPUTE THE DERIVATIVE OF THE POLYMER VOLUME FRACTION WITH
!  RESPECT TO THE CHEMICAL POTENTIAL.  COMPUTED VIA IMPLICIT
!  DIFFERENTIATION ON THE CHEMICAL POTENTIAL EQUATION.
!
dphidj = ( (vmol*kbulk)/(rgas*theta*detf)  &
    - ((vmol*kbulk)/(rgas*theta*detf))*DLOG(detf*phi_tau) )/  &
    ( (one/(phi_tau - one)) + one + two*chi*phi_tau  &
    - ((vmol*kbulk)/(rgas*theta*phi_tau))  &
    + ((vmol*kbulk)/(rgas*theta*phi_tau))*DLOG(detf*phi_tau) )
! COMPUTE THE PERTURBATION ON THE CHEMICAL POTENTIAL
!
IF(DABS(mu_tau) > one) THEN
  deltamu = DABS(mu_tau)*1.d-8
ELSE
  deltamu = 1.d-8
END IF

! COMPUTE A PERTURBED POLYMER VOLUME FRACTION
!
args(1)  = mu_tau + deltamu
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
CALL solvephi(phi_per,args,nargs,phi_t)

! COMPUTE A PERTURBED POLYMER VOLUME FRACTION
!
args(1)  = mu_tau - deltamu
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
CALL solvephi(phi_m,args,nargs,phi_t)


! COMPUTE THE DERIVATIVE OF THE TIME RATE OF CHANGE OF THE
!  POLYMER VOLUME FRACTION WITH RESPECT TO THE CHEMICAL POTENTIAL
!
dpdt_per = (phi_per - phi_t)/dtime
dpdt_m   = (phi_m - phi_t)/dtime
dphidotdmu = (dpdt_per - dpdt_m)/(two*deltamu)

! COMPUTE THE FLUID PERMEABILITY AT THIS INTEG. POINT
!
!MFLUID = D/(VMOL*RGAS*THETA)
!
! TO DO M = (D*CR)/(R*T), USE THE FOLLOWING LINE
!MFLUID = (D*(ONE/PHI_TAU - ONE))/(VMOL*RGAS*THETA)
!
! TO DO M = (D*C)/(R*T), USE THE FOLLOWING LINE
!   MFLUID = (cbari1**(cbari1)/3.d0)*
! 1                (D*(ONE/PHI_TAU - ONE))/(DETF*VMOL*RGAS*THETA)
mfluid = (d*(one/phi_tau - one))/(detf*vmol*rgas*theta)
!write(*,*) cbari1
! COMPUTE THE TANGENTS OF THE FLUID MOBILITY
!
!DMDPHI = ZERO
!
! TO DO M = (D*CR)/(R*T), USE THE FOLLOWING LINE
!DMDPHI = -(D/(VMOL*PHI_TAU*PHI_TAU*RGAS*THETA))
!
! TO DO M = (D*C)/(R*T), USE THE FOLLOWING LINE
dmdphi = -(d/(detf*vmol*phi_tau*phi_tau*rgas*theta))
!
dmdmu = dmdphi*dphidmu
dmdj  = dmdphi*dphidj
! COMPUTE THE DISPLACEMENT - CHEMICAL POTENTIAL MODULUS
!
spucmod = (kbulk/(detfe*phi_tau))*iden*dphidmu
! COMPUTE THE CHEMICAL POTENTIAL - DISPLACEMENT MODULUS
!
spcumodfac = mfluid*iden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PHIT=ONE
!      PHITAU=PHIT
!C     TIME DERIVATIVE OF THE FRACTION CONCENTRATION (FINITE DIF)
!      DPDT= (PHITAU - PHIT)/DTIME
!C     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
!      DPHIDMU=ZERO
!C     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
!      DPHIDMU=ZERO
!      DDUDDC=ONE
!      DDCDDU=ONE

!      CHI    = PROPS(3)
!      D      = PROPS(4)
!      MU0    = PROPS(5)
!      VMOL   = PROPS(6)
!      RGAS   = PROPS(7)
!
!
!        DDUDD
!        DDCDDU
!         MUTAU
!        THETA0
!        PHIT
!        PHITAU
!        DPDT
!        DPHIDMU
!        DPHITDMU,
!         MOB
!        DMDMU
!        DMDJ
!        VMOL = ONE



END SUBROUTINE
