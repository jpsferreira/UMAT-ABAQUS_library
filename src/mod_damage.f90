!> @brief Continuum damage mechanics models.
!>
!> Implements damage evolution laws that multiplicatively reduce
!> stress and stiffness based on strain energy history.
module mod_damage
  use mod_constants, only: dp, ZERO, ONE
  implicit none
  private
  public :: damage_sigmoid

contains

  !> Sigmoid damage evolution: d = 1 - 1/(1 + exp[beta*(W - psi_half)])
  !>
  !> Damage only increases — if current strain energy W < stored max W0,
  !> damage remains at its previous value.
  !>
  !> @param[in]     sef       Current strain energy
  !> @param[in,out] sef0      Maximum historical strain energy (updated if sef > sef0)
  !> @param[in]     dmg       Current damage variable (0 = undamaged, 1 = fully damaged)
  !> @param[out]    dmg_red   Damage reduction factor (1 - d)
  !> @param[out]    dmg_diff  Derivative d(dmg_red)/d(sef) for consistent tangent
  !> @param[in]     beta_d    Controls sharpness of damage activation
  !> @param[in]     psi_half  Energy threshold for 50% damage
  subroutine damage_sigmoid(sef, sef0, dmg, dmg_red, dmg_diff, beta_d, psi_half)
    real(dp), intent(in)    :: sef, beta_d, psi_half
    real(dp), intent(inout) :: sef0, dmg
    real(dp), intent(out)   :: dmg_red, dmg_diff

    dmg_diff = ZERO
    dmg_red  = ONE - dmg

    if (sef >= sef0) then
      dmg_red  = ONE / (ONE + exp(beta_d * (sef - psi_half)))
      dmg_diff = -beta_d * dmg_red * (ONE - dmg_red)
      sef0     = sef
      dmg      = ONE - dmg_red
    end if
  end subroutine damage_sigmoid

end module mod_damage
