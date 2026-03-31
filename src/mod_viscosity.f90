!> @brief Viscoelastic models based on generalized Maxwell elements.
!>
!> Implements stress relaxation via internal variables (hidden stress tensors)
!> using an exponential integration scheme. Supports multiple Maxwell branches.
module mod_viscosity
  use mod_constants, only: dp, ZERO, ONE, TWO, HALF
  implicit none
  private
  public :: visco_maxwell
  public :: relax_branch
  public :: hv_read, hv_write

  ! Maximum number of viscous branches
  integer, parameter, public :: MAX_VISCO_BRANCHES = 3

contains

  !> Apply generalized Maxwell viscoelastic contribution.
  !>
  !> Modifies PK2 stress and material elasticity tensor by adding
  !> viscous overstress from each Maxwell branch.
  !>
  !> @param[out]    pk2      Total PK2 stress (pkvol + pkiso + viscous)
  !> @param[out]    cmat     Total material elasticity tensor
  !> @param[in]     n_branches Number of Maxwell branches
  !> @param[in]     pkvol    Volumetric PK2 stress
  !> @param[in]     pkiso    Isochoric PK2 stress
  !> @param[in]     cmatvol  Volumetric material elasticity
  !> @param[in]     cmatiso  Isochoric material elasticity
  !> @param[in]     dtime    Time increment
  !> @param[in]     vscprops Viscous properties: [tau_1, theta_1, tau_2, theta_2, ...]
  !> @param[in,out] statev   State variables (hidden stress tensors stored here)
  !> @param[in]     sdv_offset Starting index in statev for viscous state variables
  subroutine visco_maxwell(pk2, cmat, n_branches, pkvol, pkiso, &
                           cmatvol, cmatiso, dtime, vscprops, statev, sdv_offset)
    real(dp), intent(out)   :: pk2(3,3), cmat(3,3,3,3)
    integer,  intent(in)    :: n_branches, sdv_offset
    real(dp), intent(in)    :: pkvol(3,3), pkiso(3,3)
    real(dp), intent(in)    :: cmatvol(3,3,3,3), cmatiso(3,3,3,3)
    real(dp), intent(in)    :: dtime, vscprops(2*MAX_VISCO_BRANCHES)
    real(dp), intent(inout) :: statev(*)

    real(dp) :: q(3,3), qv(3,3), hv(3,3), hv0(3,3)
    real(dp) :: tau, theta, aux_c, branch_aux
    integer  :: v1, i, j, k, l, pos

    q     = ZERO
    aux_c = ZERO

    do v1 = 1, n_branches
      tau   = vscprops(2*v1 - 1)
      theta = vscprops(2*v1)

      ! Read hidden stress from state variables
      pos = sdv_offset + 9*(v1-1)
      call hv_read(hv, statev, pos)
      hv0 = hv

      ! Compute relaxation
      call relax_branch(qv, hv, branch_aux, hv0, pkiso, dtime, tau, theta)
      aux_c = aux_c + branch_aux

      ! Write updated hidden stress
      call hv_write(statev, hv, pos)

      q = q + qv
    end do

    ! Assemble total stress and stiffness
    aux_c = ONE + aux_c
    pk2 = pkvol + pkiso

    do i = 1, 3
      do j = 1, 3
        pk2(i,j) = pk2(i,j) + q(i,j)
        do k = 1, 3
          do l = 1, 3
            cmat(i,j,k,l) = cmatvol(i,j,k,l) + aux_c * cmatiso(i,j,k,l)
          end do
        end do
      end do
    end do
  end subroutine visco_maxwell

  !> Single Maxwell branch relaxation (exponential integration).
  !>
  !> @param[out] qv        Viscous stress contribution from this branch
  !> @param[out] hv        Updated hidden stress
  !> @param[out] aux_stiff Stiffness scaling contribution from this branch
  !> @param[in]  hv0       Previous hidden stress
  !> @param[in]  pkiso     Current isochoric PK2 stress
  !> @param[in]  dtime     Time increment
  !> @param[in]  tau       Relaxation time
  !> @param[in]  theta     Viscous fraction
  subroutine relax_branch(qv, hv, aux_stiff, hv0, pkiso, dtime, tau, theta)
    real(dp), intent(out) :: qv(3,3), hv(3,3), aux_stiff
    real(dp), intent(in)  :: hv0(3,3), pkiso(3,3), dtime, tau, theta
    real(dp) :: decay
    integer :: i, j

    decay     = exp(-dtime / (TWO*tau))
    aux_stiff = theta * decay

    do i = 1, 3
      do j = 1, 3
        qv(i,j) = hv0(i,j) + aux_stiff * pkiso(i,j)
        hv(i,j) = decay * (decay*qv(i,j) - theta*pkiso(i,j))
      end do
    end do
  end subroutine relax_branch

  !> Read a 3x3 hidden stress tensor from the state variable array.
  subroutine hv_read(hv, statev, pos)
    real(dp), intent(out) :: hv(3,3)
    real(dp), intent(in)  :: statev(*)
    integer,  intent(in)  :: pos
    integer :: i, j, idx

    idx = pos
    do i = 1, 3
      do j = 1, 3
        idx = idx + 1
        hv(i,j) = statev(idx)
      end do
    end do
  end subroutine hv_read

  !> Write a 3x3 hidden stress tensor to the state variable array.
  subroutine hv_write(statev, hv, pos)
    real(dp), intent(inout) :: statev(*)
    real(dp), intent(in)    :: hv(3,3)
    integer,  intent(in)    :: pos
    integer :: i, j, idx

    idx = pos
    do i = 1, 3
      do j = 1, 3
        idx = idx + 1
        statev(idx) = hv(i,j)
      end do
    end do
  end subroutine hv_write

end module mod_viscosity
