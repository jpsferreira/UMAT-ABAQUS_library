SUBROUTINE integ(props,nprops,dtime,  &
        f_tau,mu_tau,phi_t,theta,  &
        t_tau,sptanmod,  &
        phi_tau)



! This subroutine computes everything required for the time integration
! of the problem.
!
! Inputs:
!  1) material parameters, props(nprops)
!  2) time increment, dtime
!  3) deformation gradient, F_tau(3,3)
!  4) chemical potential, mu_tau
!  5) old polymer volume fraction, phi_t
!  6) temperature, theta
!
! Outputs:
!  1) Cauchy stress, T_tau(3,3)
!  2) spatial tangent modulus, SpTanMod(3,3,3,3)
!  3) polymer volume fraction, phi_tau
!  4) time rate of polymer volume fraction, dPdt
!  5) derivative of the phi with mu, DphiDmu
!  6) derivative of the time rate of phi with mu, DphidotDmu
!  7) scalar fluid permeability, Mfluid
!  8) derivative of permeability with chemical potential, DmDmu
!  9) volume of a mole of fluid, Vmol
! 10) displacement - chemical potential modulus terms
! 11) chemical potential - displacement modulus terms
use global


real(8), INTENT(IN OUT)                   :: props(nprops)
INTEGER, INTENT(IN OUT)                  :: nprops
real(8), INTENT(IN OUT)                   :: dtime
real(8), INTENT(IN OUT)                   :: f_tau(3,3)
real(8), INTENT(IN OUT)                   :: mu_tau
real(8), INTENT(IN OUT)                   :: phi_t
real(8), INTENT(IN OUT)                   :: theta
real(8), INTENT(OUT)                      :: t_tau(3,3)
real(8), INTENT(OUT)                      :: sptanmod(3,3,3,3)
real(8), INTENT(OUT)                      :: phi_tau

INTEGER :: i,j,k,l,m,n, stat
INTEGER, PARAMETER :: nargs=8

real(8) iden(3,3)
real(8)  tr_tau(3,3), dtrdf(3,3,3,3),gshear,kbulk
real(8)  chi,d,mu0,vmol,rgas,detf,finvt(3,3),dpdt
real(8) b_tau(3,3),trb_tau,c_tau(3,3),trc_tau,args(nargs),detfe
real(8) deltamu,dphidmu,dpdt_per,dpdt_m,dphidotdmu,mfluid,finv(3,3)
real(8) phi_per,phi_m, dmdmu,dphidj,spucmod(3,3),dmdphi,dmdj
real(8) spcumodfac(3,3),detfs




! Identity tensor
!
CALL onem0(iden)



! Compute the inverse of F, its determinant, and its transpose
!
CALL matinv3dd(f_tau,finv,detf,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero'
  CALL xit()
END IF
finvt = transpose(finv)


! Compute the left Cauchy-Green tensor and its trace
!
b_tau = matmul(f_tau,transpose(f_tau))
trb_tau = b_tau(1,1) + b_tau(2,2) + b_tau(3,3)


! Compute the right Cauchy-Green tensor and its trace
!
c_tau = matmul(transpose(f_tau),f_tau)
trc_tau = c_tau(1,1) + c_tau(2,2) + c_tau(3,3)


! Compute the Cauchy stress
!
t_tau = (gshear*(b_tau-iden) + kbulk*DLOG(detf)*iden)/detf

! Compute the 1st Piola stress
!
tr_tau = gshear*(f_tau - finvt) + kbulk*DLOG(detf)*finvt


! Compute dTRdF, the so-called material tangent modulus
!
dtrdf = zero
DO i=1,3
  DO j = 1,3
    DO k = 1,3
      DO l = 1,3
        dtrdf(i,j,k,l) = dtrdf(i,j,k,l) + gshear*iden(i,k)*iden(j,l)  &
            + gshear*finv(l,i)*finv(j,k) + kbulk*finv(j,i)*finv(l,k)  &
            - kbulk*DLOG(detf)*finv(l,i)*finv(j,k)
      END DO
    END DO
  END DO
END DO
!
! Calculate the so-called spatial tangent modulus, based
!  on the push forward of the material tangent modulus
!
sptanmod = zero
DO i=1,3
  DO j=1,3
    DO k=1,3
      DO l=1,3
        DO m=1,3
          DO n=1,3
            sptanmod(i,j,k,l) = sptanmod(i,j,k,l) +  &
                (dtrdf(i,m,k,n)*f_tau(j,m)*f_tau(l,n))/detf
          END DO
        END DO
      END DO
    END DO
  END DO
END DO
phi_tau=detf
RETURN
END SUBROUTINE integ
