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
  public :: n_params_aniso

contains

  !> Return number of parameters per fiber family for a given aniso model.
  pure function n_params_aniso(aniso_type) result(n)
    integer, intent(in) :: aniso_type
    integer :: n
    select case (aniso_type)
    case (1) ! HGO
      n = 3  ! K1, K2, KDISP
    case (2) ! Humphrey fiber
      n = 2  ! K1, K2
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

end module mod_anisotropic
