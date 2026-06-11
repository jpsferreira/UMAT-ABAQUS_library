!> @brief Mesh-dependent UEL configuration and the global SDV transfer array.
!>
!> globalSdv carries state variables from the UEL to UVARM so they can be
!> visualized on a dummy mesh that overlays the real (UEL) mesh.
module uel_config
  implicit none

  ! Mesh-dependent settings. generate.py rewrites the parameter VALUES below
  ! when emitting uel.f90 for a model directory — keep each on its own line.
  integer, parameter :: numElem    = 1      ! number of UEL elements in the real mesh
  integer, parameter :: ElemOffset = 1000   ! element-number offset of the dummy mesh
  logical, parameter :: use_fbar   = .true. ! F-bar (de Souza Neto 1996) on the 8pt brick
  integer, parameter :: nIntPt     = 8      ! volume integration points (1 or 8)
  real(8), allocatable :: globalSdv(:,:,:)  ! (element, integ point, sdv)

end module uel_config
