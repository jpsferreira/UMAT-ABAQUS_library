!> @brief Isotropic hyperelastic strain energy functions.
!>
!> Each model computes the isochoric strain energy and fills the DISO(5) array:
!>   DISO(1) = dW/dI1
!>   DISO(2) = dW/dI2
!>   DISO(3) = d2W/dI1^2
!>   DISO(4) = d2W/dI2^2
!>   DISO(5) = d2W/(dI1 dI2)
module mod_hyperelastic
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, FOUR, HALF
  use mod_tensor, only: spectral
  use mod_kinematics, only: invariants
  implicit none
  private
  public :: sef_neo_hooke, sef_mooney_rivlin, sef_humphrey, sef_ogden
  public :: n_params_iso

contains

  !> Return the number of material parameters expected for a given isotropic model type.
  pure function n_params_iso(iso_type) result(n)
    integer, intent(in) :: iso_type
    integer :: n
    select case (iso_type)
    case (1) ! Neo-Hooke
      n = 1
    case (2) ! Mooney-Rivlin
      n = 2
    case (3) ! Ogden: first param is N_terms, then 2*N pairs
      n = -1  ! variable, handled by caller
    case (4) ! Humphrey
      n = 2
    case default
      n = 0
    end select
  end function n_params_iso

  !> Neo-Hookean: W = C10 * (I1bar - 3)
  subroutine sef_neo_hooke(sseiso, diso, c10, cbari1, cbari2)
    real(dp), intent(out) :: sseiso, diso(5)
    real(dp), intent(in)  :: c10, cbari1, cbari2

    sseiso  = c10 * (cbari1 - THREE)
    diso(1) = c10
    diso(2) = ZERO
    diso(3) = ZERO
    diso(4) = ZERO
    diso(5) = ZERO
  end subroutine sef_neo_hooke

  !> Mooney-Rivlin: W = C10*(I1bar - 3) + C01*(I2bar - 3)
  subroutine sef_mooney_rivlin(sseiso, diso, c10, c01, cbari1, cbari2)
    real(dp), intent(out) :: sseiso, diso(5)
    real(dp), intent(in)  :: c10, c01, cbari1, cbari2

    sseiso  = c10*(cbari1 - THREE) + c01*(cbari2 - THREE)
    diso(1) = c10
    diso(2) = c01
    diso(3) = ZERO
    diso(4) = ZERO
    diso(5) = ZERO
  end subroutine sef_mooney_rivlin

  !> Exponential (Humphrey): W = C10 * (exp[C01*(I1bar - 3)] - 1)
  subroutine sef_humphrey(sseiso, diso, c10, c01, cbari1, cbari2)
    real(dp), intent(out) :: sseiso, diso(5)
    real(dp), intent(in)  :: c10, c01, cbari1, cbari2
    real(dp) :: expterm

    expterm = exp(c01 * (cbari1 - THREE))
    sseiso  = c10 * (expterm - ONE)
    diso(1) = c10 * c01 * expterm
    diso(2) = ZERO
    diso(3) = c01 * c01 * c10 * expterm
    diso(4) = ZERO
    diso(5) = ZERO
  end subroutine sef_humphrey

  !> Ogden model: W = sum_i (2*mu_i/alpha_i^2) * (lam1^alpha + lam2^alpha + lam3^alpha - 3)
  !> Uses spectral decomposition; params = [mu_1, alpha_1, mu_2, alpha_2, ...]
  subroutine sef_ogden(sseiso, diso, c, cbar, params, nterms)
    real(dp), intent(out) :: sseiso, diso(5)
    real(dp), intent(in)  :: c(3,3), cbar(3,3)
    real(dp), intent(in)  :: params(:)
    integer,  intent(in)  :: nterms

    real(dp) :: ps(3), an(3,3), ci1, ci2, cbari1, cbari2
    real(dp) :: check12, check13, check23, tol, mincheck
    real(dp) :: alpha, gamma, coef, aa, bb
    real(dp) :: dudi1, dudi2, d2ud2i1, d2dudi1di2, d2ud2i2
    integer  :: neig, i1, k1

    ! Ogden two-eigenvalue case variables
    real(dp) :: di1da, di2da, di1db, di2db, d
    real(dp) :: dadi1, dadi2, dbdi1, dbdi2
    real(dp) :: d2di1da2, d2di1dadb, d2di1db2
    real(dp) :: d2di2da2, d2di2dadb, d2di2db2
    real(dp) :: d2dadi1di2, d2dbdi1di2, d2dad2i1, d2dbd2i1
    real(dp) :: d2dad2i2, d2dbd2i2
    real(dp) :: duda, dudb, d2duda, d2dudadb, d2dudb
    real(dp) :: v1, v2, v3, v4, v5, v6, v7
    real(dp) :: gamma1, gamma2, gamma3
    real(dp) :: gamma01, gamma02, gamma03, gamma04, gamma05
    real(dp) :: dd(3)

    ! Ogden three-equal case variables
    real(dp) :: c10, c01, c11, c20, c02

    call invariants(c, ci1, ci2)
    call invariants(cbar, cbari1, cbari2)
    call spectral(cbar, ps, an)

    ! Classify eigenvalue degeneracy
    tol = 1.0e-5_dp
    check12 = abs(ps(1) - ps(2))
    check13 = abs(ps(1) - ps(3))
    check23 = abs(ps(2) - ps(3))

    if (check12 < tol .and. check13 < tol .and. check23 < tol) then
      neig = 3
    else if (check12 < tol .or. check13 < tol .or. check23 < tol) then
      neig = 2
      mincheck = check12
      aa = (mincheck / TWO)**2
      bb = (ps(1) + ps(2)) / TWO
      if (check13 < mincheck) then
        mincheck = check13
        aa = (mincheck / TWO)**2
        bb = (ps(1) + ps(3)) / TWO
      else if (check23 < mincheck) then
        mincheck = check23
        aa = (mincheck / TWO)**2
        bb = (ps(2) + ps(3)) / TWO
      end if
    else
      neig = 1
    end if

    ! Initialize accumulators
    sseiso     = ZERO
    dudi1      = ZERO
    dudi2      = ZERO
    d2ud2i1    = ZERO
    d2dudi1di2 = ZERO
    d2ud2i2    = ZERO

    if (neig == 1) then
      ! --- Three distinct eigenvalues ---
      dd(1) = (ps(1) - ps(2)) * (ps(1) - ps(3))
      dd(2) = (ps(2) - ps(1)) * (ps(2) - ps(3))
      dd(3) = (ps(3) - ps(1)) * (ps(3) - ps(2))

      do i1 = 1, nterms
        alpha = params(2*i1)
        gamma = HALF * alpha
        coef  = TWO * params(2*i1-1) / (FOUR * gamma**2)
        ! Strain energy sums over all 3 eigenvalues per Ogden term
        do k1 = 1, 3
          sseiso = sseiso + coef * (ps(k1)**gamma - ONE)
        end do
        do k1 = 1, 3
          dudi1 = dudi1 + coef*gamma*(ps(k1)**(gamma+ONE)) / dd(k1)
          dudi2 = dudi2 - coef*gamma*(ps(k1)**gamma) / dd(k1)
          d2ud2i1 = d2ud2i1 &
            + coef*gamma*(gamma+ONE)*(ps(k1)**(gamma+TWO)) / (dd(k1)**2) &
            + coef*TWO*gamma*(ps(k1)**(gamma+ONE)) &
              * (ONE - ps(k1)**3) / (dd(k1)**3)
          d2dudi1di2 = d2dudi1di2 &
            - coef*gamma*gamma*(ps(k1)**(gamma+ONE)) / (dd(k1)**2) &
            - coef*TWO*gamma*(ps(k1)**gamma) &
              * (ONE - ps(k1)**3) / (dd(k1)**3)
          d2ud2i2 = d2ud2i2 &
            + coef*gamma*gamma*(ps(k1)**gamma) / (dd(k1)**2) &
            + coef*gamma*(ps(k1)**gamma) &
              * (cbari2 - THREE*ps(k1)**2) / (dd(k1)**3)
        end do
      end do

    else if (neig == 2) then
      ! --- Two equal eigenvalues ---
      di1da = bb**(-4) + TWO * bb**(-6.0_dp) * aa
      di2da = TWO * bb**(-3) - ONE + FOUR * bb**(-5.0_dp) * aa
      di1db = TWO - TWO*bb**(-3) - FOUR*bb**(-5.0_dp)*aa &
            - 6.0_dp*bb**(-7.0_dp)*aa*aa
      di2db = TWO*bb - TWO*bb**(-2) - 6.0_dp*bb**(-4.0_dp)*aa &
            - 10.0_dp*bb**(-6.0_dp)*aa*aa
      d = di1da*di2db - di1db*di2da

      dadi1  =  di2db / d
      dadi2  = -di1db / d
      dbdi1  = -di2da / d
      dbdi2  =  di1da / d

      d2di1da2  = TWO * bb**(-6.0_dp)
      d2di1dadb = -FOUR*bb**(-5.0_dp) - 12.0_dp*bb**(-7.0_dp)*aa
      d2di1db2  = 6.0_dp*bb**(-4) + 20.0_dp*bb**(-6.0_dp)*aa &
                + 42.0_dp*bb**(-8.0_dp)*aa*aa
      d2di2da2  = FOUR * bb**(-5.0_dp)
      d2di2dadb = -6.0_dp*bb**(-4.0_dp) - 20.0_dp*bb**(-6.0_dp)*aa
      d2di2db2  = TWO + 4.0_dp*bb**(-3) + 24.0_dp*bb**(-5.0_dp)*aa &
                + 60.0_dp*bb**(-7.0_dp)*aa*aa

      ! Second derivatives of a, b w.r.t. I1, I2
      v1 = dadi1*dadi2;  v2 = dbdi1*dbdi2
      v3 = dadi1*dbdi2 + dbdi1*dadi2
      v4 = d2di1da2*v1 + d2di1dadb*v3 + d2di1db2*v2
      v5 = d2di2da2*v1 + d2di2dadb*v3 + d2di2db2*v2
      d2dadi1di2 = -dadi1*v4 - dadi2*v5
      d2dbdi1di2 = -dbdi1*v4 - dbdi2*v5

      v1 = dadi1*dadi1;  v2 = dbdi1*dbdi1
      v6 = dadi1*dbdi1 + dbdi1*dadi1
      v4 = d2di1da2*v1 + d2di1dadb*v6 + d2di1db2*v2
      v5 = d2di2da2*v1 + d2di2dadb*v6 + d2di2db2*v2
      d2dad2i1 = -dadi1*v4 - dadi2*v5
      d2dbd2i1 = -dbdi1*v4 - dbdi2*v5

      v1 = dadi2*dadi2;  v2 = dbdi2*dbdi2
      v7 = dadi2*dbdi2 + dbdi2*dadi2
      v4 = d2di1da2*v1 + d2di1dadb*v7 + d2di1db2*v2
      v5 = d2di2da2*v1 + d2di2dadb*v7 + d2di2db2*v2
      d2dad2i2 = -dadi1*v4 - dadi2*v5
      d2dbd2i2 = -dbdi1*v4 - dbdi2*v5

      do i1 = 1, nterms
        alpha = params(2*i1)
        gamma = HALF * alpha
        coef  = TWO * params(2*i1-1) / (FOUR * gamma**2)
        gamma1 = gamma + ONE;  gamma2 = gamma1 + ONE
        gamma3 = gamma2 + ONE
        gamma01 = gamma - ONE;  gamma02 = gamma01 - ONE
        gamma03 = gamma02 - ONE; gamma04 = gamma03 - ONE
        gamma05 = gamma04 - ONE

        do k1 = 1, 3
          sseiso = sseiso + coef * (ps(k1)**gamma - ONE)
        end do

        duda = gamma*bb**(-TWO*gamma1) &
             + gamma*gamma01*bb**gamma02 &
             + (ONE/6.0_dp)*gamma*gamma01*gamma02*gamma03*bb**gamma04*aa &
             + gamma*gamma1*bb**(-TWO*gamma2)*aa
        dudb = TWO*gamma*bb**gamma01 - TWO*gamma*bb**(-gamma1-gamma) &
             - TWO*gamma*gamma1*bb**(-TWO*gamma-THREE)*aa &
             + gamma*gamma01*gamma02*bb**gamma03*aa &
             + (ONE/12.0_dp)*gamma*gamma01*gamma02*gamma03 &
               *gamma04*bb**gamma05*aa*aa &
             - gamma*gamma1*gamma2*bb**(-TWO*gamma-5.0_dp)*aa*aa

        d2duda = (ONE/6.0_dp)*gamma*gamma01*gamma02*gamma03*bb**gamma04 &
               + gamma*gamma1*bb**(-TWO*gamma2)
        d2dudadb = -TWO*gamma*gamma1*bb**(-TWO*gamma-THREE) &
                 + gamma*gamma01*gamma02*bb**gamma03 &
                 + (ONE/6.0_dp)*gamma*gamma01*gamma02*gamma03 &
                   *gamma04*bb**gamma05*aa &
                 - TWO*gamma*gamma1*gamma2*bb**(-TWO*gamma2-ONE)*aa
        d2dudb = TWO*gamma*gamma01*bb**gamma02 &
               + TWO*gamma*(TWO*gamma+ONE)*bb**(-TWO*gamma1) &
               + TWO*gamma*gamma1*(TWO*gamma+THREE)*bb**(-TWO*gamma2)*aa &
               + gamma*gamma01*gamma02*gamma03*bb**gamma04*aa &
               + gamma*gamma1*gamma2*(TWO*gamma+5.0_dp) &
                 *bb**(-TWO*gamma3)*aa*aa &
               + (ONE/12.0_dp)*gamma*gamma01*gamma02*gamma03 &
                 *gamma04*gamma05*bb**(gamma-6.0_dp)*aa*aa

        dudi1 = dudi1 + coef*(duda*dadi1 + dudb*dbdi1)
        dudi2 = dudi2 + coef*(duda*dadi2 + dudb*dbdi2)
        d2dudi1di2 = d2dudi1di2 + coef*( &
          d2duda*dadi1*dadi2 + d2dudadb*v3 &
        + d2dudb*dbdi1*dbdi2 + duda*d2dadi1di2 + dudb*d2dbdi1di2)
        d2ud2i1 = d2ud2i1 + coef*( &
          d2duda*dadi1*dadi1 + d2dudadb*v6 &
        + d2dudb*dbdi1*dbdi1 + duda*d2dad2i1 + dudb*d2dbd2i1)
        d2ud2i2 = d2ud2i2 + coef*( &
          d2duda*dadi2*dadi2 + d2dudadb*v7 &
        + d2dudb*dbdi2*dbdi2 + duda*d2dad2i2 + dudb*d2dbd2i2)
      end do

    else
      ! --- Three equal eigenvalues ---
      do i1 = 1, nterms
        alpha = params(2*i1)
        gamma = HALF * alpha
        coef  = TWO * params(2*i1-1) / (alpha**2)
        c10 = HALF*(ONE + gamma)
        c01 = HALF*(ONE - gamma)
        c11 = (-ONE/120.0_dp)*(gamma*gamma - ONE)*(gamma*gamma - FOUR)
        c20 = (ONE/240.0_dp)*(gamma*gamma - ONE)*(gamma*gamma + 5.0_dp*gamma + 6.0_dp)
        c02 = (ONE/240.0_dp)*(gamma*gamma - ONE)*(gamma*gamma - 5.0_dp*gamma + 6.0_dp)

        do k1 = 1, 3
          sseiso = sseiso + coef*(ps(k1)**gamma - ONE)
        end do
        dudi1 = dudi1 + coef*gamma**2*( &
          c10 + TWO*c20*(ci1 - THREE) &
        + c11*((ci2 - THREE) - 0.125_dp*(ci1+ci2-6.0_dp)**2))
        dudi2 = dudi2 + coef*gamma**2*( &
          c01 + TWO*c02*(ci2 - THREE) &
        + c11*((ci1 - THREE) - 0.125_dp*(ci1+ci2-6.0_dp)**2))
        d2ud2i1 = d2ud2i1 + coef*gamma**2*( &
          TWO*c20 - 0.25_dp*c11*(ci1+ci2-6.0_dp))
        d2ud2i2 = d2ud2i2 + coef*gamma**2*( &
          TWO*c02 - 0.25_dp*c11*(ci1+ci2-6.0_dp))
        d2dudi1di2 = d2dudi1di2 + coef*gamma**2*( &
          c11*(ONE - 0.25_dp*(ci1+ci2-6.0_dp)))
      end do
    end if

    diso(1) = dudi1
    diso(2) = dudi2
    diso(3) = d2ud2i1
    diso(4) = d2ud2i2
    diso(5) = d2dudi1di2
  end subroutine sef_ogden

end module mod_hyperelastic
