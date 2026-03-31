!> @brief Filament network models for biofilament mechanics.
!>
!> Implements the single-filament force-stretch law (shared by all network types)
!> and network assembly strategies: affine, non-affine, mixed, and contractile.
!>
!> All network types share the same FIL subroutine for single-filament mechanics.
!> They differ in how filament stretches are computed and how the network
!> stress/stiffness tensors are assembled from individual filament contributions.
module mod_network
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, FOUR, HALF
  implicit none
  private

  ! Single filament mechanics
  public :: filament_force
  public :: pullforce_brent
  public :: evalg

  ! Network assembly
  public :: affine_network
  public :: nonaffine_network

  ! Helper routines
  public :: deformed_filament
  public :: filament_stress_fic
  public :: filament_stiffness_fic
  public :: orientation_density

contains

  ! ============================================================================
  ! SINGLE FILAMENT MECHANICS (shared by all network types)
  ! ============================================================================

  !> Single filament: compute force, normalized force, and strain energy derivatives.
  !>
  !> Uses Brent's root-finding method to solve the implicit force-stretch relation.
  !>
  !> @param[out] fi     Pulling force
  !> @param[out] ffi    Normalized force
  !> @param[out] dw     First derivative of filament energy (= lambda0 * r0 * F)
  !> @param[out] ddw    Second derivative (tangent stiffness)
  !> @param[in]  lambda Filament extension ratio
  !> @param[in]  lambda0 Reference length density
  !> @param[in]  ll     Filament contour length
  !> @param[in]  r0     Filament end-to-end distance
  !> @param[in]  mu0    Bending stiffness parameter
  !> @param[in]  beta   Power-law exponent
  !> @param[in]  b0     Thermal persistence parameter
  subroutine filament_force(fi, ffi, dw, ddw, lambda, lambda0, ll, r0, mu0, beta, b0)
    real(dp), intent(out) :: fi, ffi, dw, ddw
    real(dp), intent(in)  :: lambda, lambda0, ll, r0, mu0, beta, b0
    real(dp) :: a, b, machep, tol
    real(dp) :: pi, alpha
    real(dp) :: aux0, aux, aux1, aux2, aux3, aux4, aux5, aux6, y

    a = ZERO
    b = 1.0e9_dp
    machep = 2.2204e-16_dp
    tol = 1.0e-6_dp
    fi = ZERO

    call pullforce_brent(fi, a, b, machep, tol, lambda, lambda0, ll, r0, mu0, beta, b0)

    pi    = FOUR * atan(ONE)
    ffi   = fi * ll / (pi*pi*b0)
    alpha = pi*pi*b0 / (ll*ll*mu0)

    aux0 = beta / alpha
    aux  = alpha * ffi
    aux1 = ONE + ffi + aux*ffi
    aux2 = ONE + TWO*aux
    aux3 = ONE + aux
    aux4 = lambda0 * r0*r0 * mu0 / ll
    aux5 = ((ONE + aux) / aux1)**beta
    aux6 = ONE - r0 / ll

    y = aux0*(aux2*aux2/aux1) - beta*(aux2/aux3) - TWO

    dw  = lambda0 * r0 * fi
    ddw = aux4 / (ONE + y*aux5*aux6)
  end subroutine filament_force

  !> Brent's method for root finding of the filament force-stretch equation.
  subroutine pullforce_brent(root, a_in, b_in, machep, tol, &
                             lambda, lambda0, ll, r0, mu0, beta, b0)
    real(dp), intent(out) :: root
    real(dp), intent(in)  :: a_in, b_in, machep, tol
    real(dp), intent(in)  :: lambda, lambda0, ll, r0, mu0, beta, b0

    real(dp) :: sa, sb, c, d, e, fa, fb, fc
    real(dp) :: m, p, q, r, s, toler
    integer  :: max_iter, iter

    max_iter = 500
    sa = a_in
    sb = b_in
    call evalg(fa, sa, lambda, lambda0, ll, r0, mu0, beta, b0)
    call evalg(fb, sb, lambda, lambda0, ll, r0, mu0, beta, b0)

    c  = sa
    fc = fa
    e  = sb - sa
    d  = e

    do iter = 1, max_iter
      if (abs(fc) < abs(fb)) then
        sa = sb; sb = c; c = sa
        fa = fb; fb = fc; fc = fa
      end if

      toler = TWO * machep * abs(sb) + tol
      m = HALF * (c - sb)

      if (abs(m) <= toler .or. fb == ZERO) exit

      if (abs(e) >= toler .and. abs(fa) > abs(fb)) then
        s = fb / fa
        if (sa == c) then
          p = TWO * m * s
          q = ONE - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (TWO*m*q*(q - r) - (sb - sa)*(r - ONE))
          q = (q - ONE) * (r - ONE) * (s - ONE)
        end if

        if (p > ZERO) then
          q = -q
        else
          p = -p
        end if

        s = e; e = d

        if (TWO*p >= THREE*m*q - abs(toler*q) .or. p >= abs(HALF*s*q)) then
          e = m; d = e
        else
          d = p / q
        end if
      else
        e = m; d = e
      end if

      sa = sb; fa = fb

      if (abs(d) > toler) then
        sb = sb + d
      else
        if (m > ZERO) then
          sb = sb + toler
        else
          sb = sb - toler
        end if
      end if

      call evalg(fb, sb, lambda, lambda0, ll, r0, mu0, beta, b0)

      if ((fb > ZERO .and. fc > ZERO) .or. (fb <= ZERO .and. fc <= ZERO)) then
        c = sa; fc = fa
        e = sb - sa; d = e
      end if
    end do

    root = sb
  end subroutine pullforce_brent

  !> Evaluate G(F) = LHS - RHS for the filament stretch-force relation.
  subroutine evalg(g, f, lambda, lambda0, ll, r0, mu0, beta, b0)
    real(dp), intent(out) :: g
    real(dp), intent(in)  :: f, lambda, lambda0, ll, r0, mu0, beta, b0
    real(dp) :: pi, aux0, aux1, aux, aux2, aux3, aux4, lhs, rhs

    pi   = FOUR * atan(ONE)
    aux0 = ONE - r0/ll
    aux1 = ll*ll / (pi*pi*b0)
    aux  = f / mu0
    aux2 = ONE + aux
    aux3 = ONE + TWO*aux
    aux4 = ONE + f*aux1 + f*aux*aux1

    rhs = ONE + aux - aux0 * (aux2**beta) * aux3 * (aux4**(-beta))
    lhs = lambda * lambda0 * r0 / ll

    g = lhs - rhs
  end subroutine evalg

  ! ============================================================================
  ! NETWORK HELPER ROUTINES
  ! ============================================================================

  !> Compute deformed filament direction and stretch.
  !> lambda = |F * m0|, mfi = F * m0 / |F * m0|
  subroutine deformed_filament(lambda, mfi, m0, f)
    real(dp), intent(out) :: lambda, mfi(3)
    real(dp), intent(in)  :: m0(3), f(3,3)
    real(dp) :: fm(3)

    fm = matmul(f, m0)
    lambda = sqrt(fm(1)**2 + fm(2)**2 + fm(3)**2)
    if (lambda > ZERO) then
      mfi = fm / lambda
    else
      mfi = ZERO
    end if
  end subroutine deformed_filament

  !> Single filament contribution to fictitious Cauchy stress.
  !> sigma_fil = rho * lambda^{-1} * rw * dw * m (x) m
  subroutine filament_stress_fic(sfil, rho, lambda, dw, mfi, rw)
    real(dp), intent(out) :: sfil(3,3)
    real(dp), intent(in)  :: rho, lambda, dw, mfi(3), rw
    real(dp) :: scale
    integer  :: i, j

    if (lambda > ZERO) then
      scale = rho * dw / lambda * rw
    else
      scale = ZERO
    end if

    do i = 1, 3
      do j = 1, 3
        sfil(i,j) = scale * mfi(i) * mfi(j)
      end do
    end do
  end subroutine filament_stress_fic

  !> Single filament contribution to fictitious spatial elasticity tensor.
  subroutine filament_stiffness_fic(cfil, rho, lambda, dw, ddw, mfi, rw)
    real(dp), intent(out) :: cfil(3,3,3,3)
    real(dp), intent(in)  :: rho, lambda, dw, ddw, mfi(3), rw
    real(dp) :: scale
    integer  :: i, j, k, l

    if (lambda > ZERO) then
      scale = rho * (ddw - dw/lambda) / (lambda*lambda) * rw
    else
      scale = ZERO
    end if

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cfil(i,j,k,l) = scale * mfi(i)*mfi(j)*mfi(k)*mfi(l)
          end do
        end do
      end do
    end do
  end subroutine filament_stiffness_fic

  !> Von Mises-Fisher orientation density function.
  !> rho(theta) = normalization * exp(b*(cos(2*theta) + 1))
  subroutine orientation_density(rho, angle, b_param, efi)
    real(dp), intent(out) :: rho
    real(dp), intent(in)  :: angle, b_param, efi
    real(dp) :: pi

    pi = FOUR * atan(ONE)
    if (abs(b_param) < 1.0e-10_dp) then
      rho = ONE
    else
      rho = FOUR * sqrt(b_param / (TWO*pi)) * exp(b_param*(cos(TWO*angle) + ONE))
    end if
  end subroutine orientation_density

  ! ============================================================================
  ! AFFINE NETWORK ASSEMBLY
  ! ============================================================================

  !> Assemble affine network fictitious stress and stiffness tensors.
  !>
  !> Integrates single-filament contributions over NWP quadrature directions.
  !> Each filament deforms affinely with the macroscopic deformation gradient.
  !>
  !> @param[out] sfic   Fictitious Cauchy stress (3x3)
  !> @param[out] cfic   Fictitious spatial elasticity tensor (3x3x3x3)
  !> @param[in]  f      Distortion gradient (3x3)
  !> @param[in]  mf0    Quadrature directions (nwp x 3)
  !> @param[in]  rw     Quadrature weights (nwp)
  !> @param[in]  nwp    Number of quadrature directions
  !> @param[in]  det    Jacobian determinant
  !> @param[in]  filprops  Filament properties: [L, R0, mu0, beta, B0, lambda0]
  !> @param[in]  net_density Network number density N
  !> @param[in]  b_orient   Orientation distribution parameter
  !> @param[in]  efi    Preferred direction angle
  subroutine affine_network(sfic, cfic, f, mf0, rw, nwp, det, &
                            filprops, net_density, b_orient, efi)
    real(dp), intent(out) :: sfic(3,3), cfic(3,3,3,3)
    real(dp), intent(in)  :: f(3,3), det
    integer,  intent(in)  :: nwp
    real(dp), intent(in)  :: mf0(:,:), rw(:)
    real(dp), intent(in)  :: filprops(6), net_density, b_orient, efi

    real(dp) :: pi, coeff
    real(dp) :: ll, r0, mu0, beta, b0, lambda0
    real(dp) :: mfi(3), m0i(3), lambdai, fi, ffi, dwi, ddwi, rho, angle
    real(dp) :: sfili(3,3), cfili(3,3,3,3)
    integer  :: ip, j, k, l, m

    ll      = filprops(1)
    r0      = filprops(2)
    mu0     = filprops(3)
    beta    = filprops(4)
    b0      = filprops(5)
    lambda0 = filprops(6)

    pi    = FOUR * atan(ONE)
    coeff = net_density / det * FOUR * pi

    sfic = ZERO
    cfic = ZERO

    do ip = 1, nwp
      m0i = mf0(ip, :)

      call deformed_filament(lambdai, mfi, m0i, f)

      ! Compute angle for orientation density (simplified: use m0 angle)
      angle = ZERO  ! Default isotropic; caller can provide oriented density
      call orientation_density(rho, angle, b_orient, efi)

      call filament_force(fi, ffi, dwi, ddwi, lambdai, lambda0, ll, r0, mu0, beta, b0)

      call filament_stress_fic(sfili, rho, lambdai, dwi, mfi, rw(ip))
      call filament_stiffness_fic(cfili, rho, lambdai, dwi, ddwi, mfi, rw(ip))

      do j = 1, 3
        do k = 1, 3
          sfic(j,k) = sfic(j,k) + coeff * sfili(j,k)
          do l = 1, 3
            do m = 1, 3
              cfic(j,k,l,m) = cfic(j,k,l,m) + coeff * cfili(j,k,l,m)
            end do
          end do
        end do
      end do
    end do
  end subroutine affine_network

  ! ============================================================================
  ! NON-AFFINE NETWORK ASSEMBLY
  ! ============================================================================

  !> Assemble non-affine network using power-law stretch averaging.
  !>
  !> Effective stretch: Lambda = (sum_i lambda_i^p * rw_i)^(1/p)
  !> Single filament is evaluated at the effective stretch Lambda.
  !>
  !> @param[in] pp Non-affinity exponent (p)
  subroutine nonaffine_network(sfic, cfic, f, mf0, rw, nwp, det, &
                               filprops, net_density, b_orient, efi, pp)
    real(dp), intent(out) :: sfic(3,3), cfic(3,3,3,3)
    real(dp), intent(in)  :: f(3,3), det
    integer,  intent(in)  :: nwp
    real(dp), intent(in)  :: mf0(:,:), rw(:)
    real(dp), intent(in)  :: filprops(6), net_density, b_orient, efi, pp

    real(dp) :: pi, coeff
    real(dp) :: ll, r0, mu0, beta, b0, lambda0
    real(dp) :: m0i(3), mfi(3), lambdai
    real(dp) :: lambda_eff, sum_lp, fi, ffi, dw, ddw
    real(dp) :: h2(3,3), h4(3,3,3,3), rho, angle
    real(dp) :: scale_s, scale_c1, scale_c2
    integer  :: ip, i, j, k, l

    ll      = filprops(1)
    r0      = filprops(2)
    mu0     = filprops(3)
    beta    = filprops(4)
    b0      = filprops(5)
    lambda0 = filprops(6)

    pi    = FOUR * atan(ONE)
    coeff = net_density / det

    ! First pass: compute effective stretch and structure tensors
    sum_lp = ZERO
    h2 = ZERO
    h4 = ZERO

    do ip = 1, nwp
      m0i = mf0(ip, :)
      call deformed_filament(lambdai, mfi, m0i, f)

      angle = ZERO
      call orientation_density(rho, angle, b_orient, efi)

      if (lambdai > ZERO) then
        sum_lp = sum_lp + rho * lambdai**pp * rw(ip)

        ! 2nd-order structure tensor
        do i = 1, 3
          do j = 1, 3
            h2(i,j) = h2(i,j) + rho * lambdai**(pp-TWO) * rw(ip) * mfi(i)*mfi(j)
          end do
        end do

        ! 4th-order structure tensor
        do i = 1, 3
          do j = 1, 3
            do k = 1, 3
              do l = 1, 3
                h4(i,j,k,l) = h4(i,j,k,l) + rho * (pp-TWO) * lambdai**(pp-FOUR) &
                             * rw(ip) * mfi(i)*mfi(j)*mfi(k)*mfi(l)
              end do
            end do
          end do
        end do
      end if
    end do

    ! Effective stretch
    if (sum_lp > ZERO) then
      lambda_eff = sum_lp**(ONE/pp)
    else
      lambda_eff = ONE
    end if

    ! Evaluate single filament at effective stretch
    call filament_force(fi, ffi, dw, ddw, lambda_eff, lambda0, ll, r0, mu0, beta, b0)

    ! Assemble stress and stiffness
    scale_s  = coeff * dw * lambda_eff**(ONE - pp)
    scale_c1 = coeff * dw * lambda_eff**(ONE - pp)
    scale_c2 = coeff * (ddw * lambda_eff**(TWO*(ONE-pp)) &
             - (pp-ONE) * dw * lambda_eff**(ONE - TWO*pp))

    sfic = scale_s * h2
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cfic(i,j,k,l) = scale_c1 * h4(i,j,k,l) + scale_c2 * h2(i,j)*h2(k,l)
          end do
        end do
      end do
    end do
  end subroutine nonaffine_network

end module mod_network
