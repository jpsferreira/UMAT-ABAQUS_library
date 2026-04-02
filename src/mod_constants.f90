!> @brief Global constants and precision parameters for the UMAT builder.
module mod_constants
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)

  real(dp), parameter :: ZERO  = 0.0_dp
  real(dp), parameter :: HALF  = 0.5_dp
  real(dp), parameter :: ONE   = 1.0_dp
  real(dp), parameter :: TWO   = 2.0_dp
  real(dp), parameter :: THREE = 3.0_dp
  real(dp), parameter :: FOUR  = 4.0_dp
  real(dp), parameter :: SIX   = 6.0_dp
  real(dp), parameter :: THIRD = ONE / THREE

  ! Model type IDs for PROPS-based dispatch
  ! Isotropic models
  integer, parameter :: ISO_NONE        = 0
  integer, parameter :: ISO_NEO_HOOKE   = 1
  integer, parameter :: ISO_MOONEY_RIVLIN = 2
  integer, parameter :: ISO_OGDEN       = 3
  integer, parameter :: ISO_HUMPHREY    = 4

  ! Anisotropic models
  integer, parameter :: ANISO_NONE          = 0
  integer, parameter :: ANISO_HGO           = 1
  integer, parameter :: ANISO_HUMPHREY      = 2
  integer, parameter :: ANISO_HGO_AI        = 3
  integer, parameter :: ANISO_HUMPHREY_AI   = 4
  integer, parameter :: ANISO_HUMPHREY_ACT  = 5

  ! Network models
  integer, parameter :: NET_NONE           = 0
  integer, parameter :: NET_AFFINE         = 1
  integer, parameter :: NET_NONAFFINE      = 2
  integer, parameter :: NET_MIXED          = 3
  integer, parameter :: NET_CONTRACTILE    = 4
  integer, parameter :: NET_AFFINE_AI      = 5
  integer, parameter :: NET_NONAFFINE_AI   = 6

  ! Damage models
  integer, parameter :: DMG_NONE        = 0
  integer, parameter :: DMG_SIGMOID     = 1

end module mod_constants
