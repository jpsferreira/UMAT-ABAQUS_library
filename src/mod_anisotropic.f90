!> @brief Anisotropic fiber-reinforcement models.
!>
!> Implements fiber direction handling and anisotropic strain energy contributions.
!> Each model fills DANISO(4):
!>   DANISO(1) = dW/dI4
!>   DANISO(2) = d2W/dI4^2
!>   DANISO(3) = d2W/(dI1 dI4)
!>   DANISO(4) = d2W/(dI2 dI4)
module mod_anisotropic
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, FOUR, HALF
  implicit none
  private
  public :: fiber_direction, sef_aniso_hgo, sef_aniso_humphrey
  public :: sef_aniso_humphrey_act
  public :: sef_aniso_hgo_ai
  public :: n_params_aniso

contains

  !> Return number of parameters per fiber family for a given aniso model.
  pure function n_params_aniso(aniso_type) result(n)
    integer, intent(in) :: aniso_type
    integer :: n
    select case (aniso_type)
    case (1) ! HGO: K1, K2, kappa
      n = 3
    case (2) ! Humphrey fiber: K1, K2
      n = 2
    case (3) ! HGO AI: K1, K2, bdisp, factor
      n = 4
    case (4) ! Humphrey AI: K1, K2, bdisp, factor
      n = 4
    case (5) ! Humphrey + activation: K1, K2, T0M
      n = 3
    case default
      n = 0
    end select
  end function n_params_aniso

  !> Compute fiber structural tensor from orientation vector.
  !> vorif = undeformed fiber direction (3-vector)
  !> Returns st0 = v0 (x) v0 (structural tensor in reference config)
  !> and deformed direction vd = F * v0 / |F * v0|
  subroutine fiber_direction(vorif, f, st0, vd)
    real(dp), intent(in)  :: vorif(3), f(3,3)
    real(dp), intent(out) :: st0(3,3), vd(3)
    real(dp) :: v0(3), norm0, fv(3), normfv
    integer  :: i, j

    ! Normalize
    norm0 = sqrt(vorif(1)**2 + vorif(2)**2 + vorif(3)**2)
    if (norm0 > ZERO) then
      v0 = vorif / norm0
    else
      v0 = ZERO
    end if

    ! Structural tensor
    do i = 1, 3
      do j = 1, 3
        st0(i,j) = v0(i) * v0(j)
      end do
    end do

    ! Deformed fiber direction
    fv = matmul(f, v0)
    normfv = sqrt(fv(1)**2 + fv(2)**2 + fv(3)**2)
    if (normfv > ZERO) then
      vd = fv / normfv
    else
      vd = ZERO
    end if
  end subroutine fiber_direction

  !> HGO (Holzapfel-Gasser-Ogden) anisotropic model with dispersion.
  !> W_aniso = (K1/K2) * (exp[K2*E1^2] - 1)
  !> with E1 = I4*(1 - 3*kdisp) + I1*kdisp - 1
  !> This version also modifies DISO(1) and DISO(3) to account for the
  !> coupled I1-I4 contribution when kdisp > 0.
  subroutine sef_aniso_hgo(sseaniso, daniso, diso, k1, k2, kdisp, inv4, cbari1)
    real(dp), intent(out)   :: sseaniso, daniso(4)
    real(dp), intent(inout) :: diso(5)
    real(dp), intent(in)    :: k1, k2, kdisp, inv4, cbari1
    real(dp) :: e1, ee2, ee3, dudi4, d2ud2i4, d2dudi1di4, d2dudi2di4
    real(dp) :: dudi1, d2ud2i1

    dudi1   = diso(1)
    d2ud2i1 = diso(3)

    e1 = inv4*(ONE - THREE*kdisp) + cbari1*kdisp - ONE
    sseaniso = (k1/k2) * (exp(k2*e1*e1) - ONE)

    if (e1 > ZERO) then
      ee2 = exp(k2 * e1*e1)
      ee3 = ONE + TWO*k2*e1*e1

      ! Update isotropic derivatives for coupled contribution
      dudi1   = dudi1 + k1*kdisp*e1*ee2
      d2ud2i1 = d2ud2i1 + k1*kdisp*kdisp*ee3*ee2

      dudi4       = k1*(ONE - THREE*kdisp)*e1*ee2
      d2ud2i4     = k1*((ONE - THREE*kdisp)**2)*ee3*ee2
      d2dudi1di4  = k1*(ONE - THREE*kdisp)*kdisp*ee3*ee2
      d2dudi2di4  = ZERO
    else
      dudi4      = ZERO
      d2ud2i4    = ZERO
      d2dudi1di4 = ZERO
      d2dudi2di4 = ZERO
    end if

    daniso(1) = dudi4
    daniso(2) = d2ud2i4
    daniso(3) = d2dudi1di4
    daniso(4) = d2dudi2di4

    ! Feed back coupled terms
    diso(1) = dudi1
    diso(3) = d2ud2i1
  end subroutine sef_aniso_hgo

  !> Humphrey-type exponential fiber model (stretch-based).
  !> W_aniso = K1 * (exp[K2*(lambda-1)^2] - 1)
  !> where lambda = sqrt(I4). Only active for lambda > 1 (tension).
  subroutine sef_aniso_humphrey(sseaniso, daniso, diso, k1, k2, inv4, cbari1)
    real(dp), intent(out)   :: sseaniso, daniso(4)
    real(dp), intent(inout) :: diso(5)
    real(dp), intent(in)    :: k1, k2, inv4, cbari1
    real(dp) :: lambda, e1, e2, e3, e4, e5

    lambda = sqrt(inv4)
    sseaniso = k1 * (exp(k2*(lambda - ONE)**2) - ONE)

    e1 = lambda - ONE
    if (e1 > ZERO) then
      e3 = exp(k2 * e1*e1)
      e2 = TWO * k1 * k2
      e4 = ONE + TWO*k2*e1*e1*lambda
      e5 = FOUR * inv4**(THREE/TWO)

      daniso(1) = e2*e1*e3 / (TWO*lambda)       ! dW/dI4
      daniso(2) = e2*e3*e4 / e5                  ! d2W/dI4^2
    else
      daniso(1) = ZERO
      daniso(2) = ZERO
    end if
    daniso(3) = ZERO  ! no I1-I4 coupling
    daniso(4) = ZERO  ! no I2-I4 coupling
  end subroutine sef_aniso_humphrey

  !> Humphrey fiber model with muscle activation.
  !> Adds active stress contribution: U_act = ACT * T0M * (-(4/3)*E1^3 + lambda_bar)
  !> where E1 = sqrt(I4) - 1, ACT from external field (PREDEF).
  !> Only active for lambda_bar in [0.5, 1.5] (physiological range).
  !> Modifies sseaniso and daniso IN PLACE (call AFTER sef_aniso_humphrey).
  subroutine sef_aniso_humphrey_act(sseaniso, daniso, cbari4, predef, dpred, t0m)
    real(dp), intent(inout) :: sseaniso, daniso(4)
    real(dp), intent(in)    :: cbari4, predef(1), dpred(1), t0m
    real(dp) :: barlambda, e1, act, duact, d2uact

    barlambda = sqrt(cbari4)
    e1 = barlambda - ONE
    act = predef(1) + dpred(1)

    if (barlambda > HALF .and. barlambda < 1.5_dp) then
      duact  = act * t0m * (-FOUR * e1*e1 + ONE)
      d2uact = act * t0m * (-8.0_dp * e1)
      sseaniso = sseaniso + act * t0m * (-(FOUR/THREE) * e1**3 + barlambda)
    else
      duact  = ZERO
      d2uact = ZERO
    end if

    ! Chain rule: d/dI4 = d/d(lambda_bar) * d(lambda_bar)/dI4
    ! d(lambda_bar)/dI4 = (1/2) * I4^(-1/2)
    daniso(1) = daniso(1) + duact * HALF * cbari4**(-HALF)
    daniso(2) = daniso(2) + d2uact * (-0.25_dp) * cbari4**(-1.5_dp) &
                           + d2uact * (HALF * cbari4**(-HALF))**2
  end subroutine sef_aniso_humphrey_act

  ! ============================================================================
  ! DISCRETE ANGULAR INTEGRATION (AI) FOR ANISOTROPIC MODELS
  ! ============================================================================

  !> Compute imaginary error function of sqrt(b) for density normalization.
  subroutine erfi_func(erf_val, b, nterm)
    real(dp), intent(out) :: erf_val
    real(dp), intent(in)  :: b
    integer,  intent(in)  :: nterm
    real(dp) :: pi_val, x, term_sum, fact, aux3
    integer  :: j, i

    pi_val = FOUR * atan(ONE)
    x = sqrt(TWO * b)

    term_sum = TWO * x + (TWO / THREE) * x**3
    do j = 3, nterm
      i = j - 1
      fact = ONE
      do i = 1, j - 1
        fact = fact * real(i, dp)
      end do
      aux3 = real(2*j - 1, dp)
      term_sum = term_sum + x**aux3 / (HALF * aux3 * fact)
    end do

    erf_val = pi_val**(-HALF) * term_sum
  end subroutine erfi_func

  !> Von Mises-Fisher density function for fiber orientation distribution.
  pure subroutine von_mises_density(rho, angle, bdisp, efi)
    real(dp), intent(out) :: rho
    real(dp), intent(in)  :: angle, bdisp, efi
    real(dp) :: pi_val, aux1, aux2

    pi_val = FOUR * atan(ONE)
    aux1 = sqrt(bdisp / (TWO * pi_val))
    aux2 = exp(bdisp * (cos(TWO * angle) + ONE))
    rho = FOUR * aux1 * aux2 / efi
  end subroutine von_mises_density

  !> HGO model with discrete angular integration (icosahedron-based).
  !> Returns fictitious Cauchy stress and spatial elasticity tensor directly
  !> (not invariant derivatives). Must be added to sfic/cfic in the caller.
  subroutine sef_aniso_hgo_ai(sfic_aniso, cfic_aniso, w_aniso, &
                               distgr, det, k1, k2, bdisp, factor_ai, vorif)
    use mod_icosahedron, only: icos_shape, sphere01_triangle_project, &
                               sphere01_triangle_vertices_to_area, &
                               ICOS_POINT_NUM, ICOS_EDGE_NUM, ICOS_FACE_NUM, &
                               ICOS_FACE_ORDER_MAX
    real(dp), intent(out) :: sfic_aniso(3,3), cfic_aniso(3,3,3,3), w_aniso
    real(dp), intent(in)  :: distgr(3,3), det, k1, k2, bdisp, vorif(3)
    integer,  intent(in)  :: factor_ai

    real(dp) :: point_coord(3, ICOS_POINT_NUM)
    integer  :: edge_point(2, ICOS_EDGE_NUM)
    integer  :: face_order(ICOS_FACE_NUM)
    integer  :: face_point(ICOS_FACE_ORDER_MAX, ICOS_FACE_NUM)

    real(dp) :: a_xyz(3), b_xyz(3), c_xyz(3)
    real(dp) :: a2_xyz(3), b2_xyz(3), c2_xyz(3)
    real(dp) :: node_xyz(3), ai
    real(dp) :: mfi(3), m0i(3), lambdai, lambda2, ei
    real(dp) :: dwi, ddwi, wi, rho, angle
    real(dp) :: pd(3), pd_norm, efi_val
    real(dp) :: pi_val, aux_s, aux_c, aux_rho
    integer  :: face, ia, ib, ic, f1, f2, f3
    integer  :: j, k, l, m

    pi_val = FOUR * atan(ONE)

    ! Compute erfi normalization factor
    call erfi_func(efi_val, bdisp, 60)

    ! Compute deformed preferred direction
    pd = matmul(distgr, vorif)
    pd_norm = sqrt(dot_product(pd, pd))
    if (pd_norm > ZERO) pd = pd / pd_norm

    sfic_aniso = ZERO
    cfic_aniso = ZERO
    w_aniso = ZERO

    aux_s = TWO / det
    aux_c = FOUR * det**(-FOUR/THREE)

    ! Set up icosahedron
    call icos_shape(point_coord, edge_point, face_order, face_point)

    ! Loop over icosahedron faces
    do face = 1, ICOS_FACE_NUM
      ia = face_point(1, face)
      ib = face_point(2, face)
      ic = face_point(3, face)

      a_xyz = point_coord(:, ia)
      b_xyz = point_coord(:, ib)
      c_xyz = point_coord(:, ic)

      ! Same-orientation subtriangles
      do f3 = 1, 3*factor_ai - 2, 3
        do f2 = 1, 3*factor_ai - f3 - 1, 3
          f1 = 3*factor_ai - f3 - f2

          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3, node_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1+2, f2-1, f3-1, a2_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1-1, f2+2, f3-1, b2_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1-1, f2-1, f3+2, c2_xyz)
          call sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz, ai)

          m0i = node_xyz
          ! Deform fiber direction
          mfi = matmul(distgr, m0i)
          lambdai = sqrt(dot_product(mfi, mfi))
          if (lambdai > ZERO) mfi = mfi / lambdai

          ! Angle to preferred direction and density
          aux_rho = dot_product(mfi, pd)
          aux_rho = max(-ONE, min(ONE, aux_rho))
          angle = acos(aux_rho)
          call von_mises_density(rho, angle, bdisp, efi_val)

          ! Scaled weight
          ai = ai / (TWO * pi_val)

          ! Fiber strain-like quantity: E = lambda^2 - 1
          lambda2 = lambdai * lambdai
          ei = lambda2 - ONE

          if (ei >= ZERO) then
            ! HGO fiber strain energy and derivatives
            wi   = (k1 / (TWO * k2)) * (exp(k2 * ei * ei) - ONE)
            dwi  = k1 * ei * exp(k2 * ei * ei)
            ddwi = k1 * exp(k2 * ei * ei) * (TWO * k2 * ei * ei + ONE)

            ! Accumulate stress: sfic += aux_s * rho * dw * ai * m (x) m
            do j = 1, 3
              do k = 1, 3
                sfic_aniso(j,k) = sfic_aniso(j,k) + aux_s * rho * dwi * ai * mfi(j) * mfi(k)
                do l = 1, 3
                  do m = 1, 3
                    cfic_aniso(j,k,l,m) = cfic_aniso(j,k,l,m) &
                      + aux_c * rho * ddwi * ai * mfi(j) * mfi(k) * mfi(l) * mfi(m)
                  end do
                end do
              end do
            end do

            w_aniso = w_aniso + rho * ai * wi
          end if
        end do
      end do

      ! Opposite-orientation subtriangles
      do f3 = 2, 3*factor_ai - 4, 3
        do f2 = 2, 3*factor_ai - f3 - 2, 3
          f1 = 3*factor_ai - f3 - f2

          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3, node_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1-2, f2+1, f3+1, a2_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1+1, f2-2, f3+1, b2_xyz)
          call sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1+1, f2+1, f3-2, c2_xyz)
          call sphere01_triangle_vertices_to_area(a2_xyz, b2_xyz, c2_xyz, ai)

          m0i = node_xyz
          mfi = matmul(distgr, m0i)
          lambdai = sqrt(dot_product(mfi, mfi))
          if (lambdai > ZERO) mfi = mfi / lambdai

          aux_rho = dot_product(mfi, pd)
          aux_rho = max(-ONE, min(ONE, aux_rho))
          angle = acos(aux_rho)
          call von_mises_density(rho, angle, bdisp, efi_val)

          ai = ai / (TWO * pi_val)
          lambda2 = lambdai * lambdai
          ei = lambda2 - ONE

          if (ei >= ZERO) then
            wi   = (k1 / (TWO * k2)) * (exp(k2 * ei * ei) - ONE)
            dwi  = k1 * ei * exp(k2 * ei * ei)
            ddwi = k1 * exp(k2 * ei * ei) * (TWO * k2 * ei * ei + ONE)

            do j = 1, 3
              do k = 1, 3
                sfic_aniso(j,k) = sfic_aniso(j,k) + aux_s * rho * dwi * ai * mfi(j) * mfi(k)
                do l = 1, 3
                  do m = 1, 3
                    cfic_aniso(j,k,l,m) = cfic_aniso(j,k,l,m) &
                      + aux_c * rho * ddwi * ai * mfi(j) * mfi(k) * mfi(l) * mfi(m)
                  end do
                end do
              end do
            end do

            w_aniso = w_aniso + rho * ai * wi
          end if
        end do
      end do
    end do
  end subroutine sef_aniso_hgo_ai

end module mod_anisotropic
