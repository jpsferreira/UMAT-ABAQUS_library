!> @brief Kinematics module for finite-strain continuum mechanics.
!>
!> Computes deformation measures: distortion gradient, Cauchy-Green tensors,
!> invariants, pseudo-invariants, stretch/rotation tensors, and projection tensors.
module mod_kinematics
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, HALF
  use mod_tensor, only: matinv3d, onem, spectral, contraction22
  implicit none
  private
  public :: compute_determinant
  public :: distortion_gradient
  public :: cauchy_green
  public :: invariants
  public :: pseudo_invariants
  public :: stretch_tensors
  public :: rotation_tensor
  public :: projection_euler
  public :: projection_lagrange

contains

  !> Compute determinant of a 3x3 matrix via cofactor expansion.
  function compute_determinant(f) result(det)
    real(dp), intent(in) :: f(3,3)
    real(dp) :: det
    det = f(1,1)*(f(2,2)*f(3,3) - f(2,3)*f(3,2)) &
        - f(1,2)*(f(2,1)*f(3,3) - f(2,3)*f(3,1)) &
        + f(1,3)*(f(2,1)*f(3,2) - f(2,2)*f(3,1))
  end function compute_determinant

  !> Compute distortion (volume-preserving) gradient: Fbar = det^(-1/3) * F
  subroutine distortion_gradient(f, fbar, det)
    real(dp), intent(in)  :: f(3,3)
    real(dp), intent(out) :: fbar(3,3)
    real(dp), intent(out) :: det
    real(dp) :: scale

    det = compute_determinant(f)
    scale = det**(-ONE/THREE)
    fbar = scale * f
  end subroutine distortion_gradient

  !> Compute right (C) and left (B) Cauchy-Green deformation tensors from F.
  subroutine cauchy_green(f, c, b)
    real(dp), intent(in)  :: f(3,3)
    real(dp), intent(out) :: c(3,3), b(3,3)

    c = matmul(transpose(f), f)
    b = matmul(f, transpose(f))
  end subroutine cauchy_green

  !> Compute 1st and 2nd invariants of a 3x3 symmetric tensor.
  subroutine invariants(a, i1, i2)
    real(dp), intent(in)  :: a(3,3)
    real(dp), intent(out) :: i1, i2
    real(dp) :: trA2

    i1 = a(1,1) + a(2,2) + a(3,3)
    call contraction22(trA2, a, a)
    i2 = HALF * (i1*i1 - trA2)
  end subroutine invariants

  !> Compute pseudo-invariant I4 = C:M0 (fiber stretch squared) and fiber stretches.
  subroutine pseudo_invariants(cbar, st0, det, inv4, lambda, barlambda)
    real(dp), intent(in)  :: cbar(3,3), st0(3,3), det
    real(dp), intent(out) :: inv4, lambda, barlambda

    call contraction22(inv4, cbar, st0)
    barlambda = sqrt(inv4)
    lambda = barlambda / (det**(-ONE/THREE))
  end subroutine pseudo_invariants

  !> Compute right (U) and left (V) stretch tensors via spectral decomposition.
  subroutine stretch_tensors(c, b, u, v)
    real(dp), intent(in)  :: c(3,3), b(3,3)
    real(dp), intent(out) :: u(3,3), v(3,3)
    real(dp) :: omega(3), eigvec(3,3), diag(3,3)
    integer :: i

    ! Right stretch U from C = U^2
    call spectral(c, omega, eigvec)
    diag = ZERO
    do i = 1, 3
      diag(i,i) = sqrt(omega(i))
    end do
    u = matmul(matmul(eigvec, diag), transpose(eigvec))

    ! Left stretch V from B = V^2
    call spectral(b, omega, eigvec)
    diag = ZERO
    do i = 1, 3
      diag(i,i) = sqrt(omega(i))
    end do
    v = matmul(matmul(eigvec, diag), transpose(eigvec))
  end subroutine stretch_tensors

  !> Compute rotation tensor: R = F * U^(-1)
  subroutine rotation_tensor(f, u, rot)
    real(dp), intent(in)  :: f(3,3), u(3,3)
    real(dp), intent(out) :: rot(3,3)
    real(dp) :: uinv(3,3)

    call matinv3d(u, uinv)
    rot = matmul(f, uinv)
  end subroutine rotation_tensor

  !> Eulerian (spatial) deviatoric projection tensor: PE = Isym - (1/3) I⊗I
  subroutine projection_euler(unit2, unit4s, pe)
    real(dp), intent(in)  :: unit2(3,3), unit4s(3,3,3,3)
    real(dp), intent(out) :: pe(3,3,3,3)
    integer :: i, j, k, l

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            pe(i,j,k,l) = unit4s(i,j,k,l) - (ONE/THREE)*unit2(i,j)*unit2(k,l)
          end do
        end do
      end do
    end do
  end subroutine projection_euler

  !> Lagrangian (material) deviatoric projection tensor: PL = I - (1/3) C^{-1}⊗C
  subroutine projection_lagrange(c, unit4, pl)
    real(dp), intent(in)  :: c(3,3), unit4(3,3,3,3)
    real(dp), intent(out) :: pl(3,3,3,3)
    real(dp) :: cinv(3,3)
    integer :: i, j, k, l

    call matinv3d(c, cinv)
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            pl(i,j,k,l) = unit4(i,j,k,l) - (ONE/THREE)*cinv(i,j)*c(k,l)
          end do
        end do
      end do
    end do
  end subroutine projection_lagrange

end module mod_kinematics
