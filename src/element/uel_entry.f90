! ==============================================================================
! ABAQUS-callable entry points for the user element layer.
!
!   uel    dispatches on jtype to the element kernels (module uel_element)
!   uvarm  transfers globalSdv entries onto the dummy mesh for visualization
!
! These must be EXTERNAL (non-module) procedures so ABAQUS can resolve them.
! ==============================================================================

subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
               props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, &
               time, dtime, kstep, kinc, jelem, params, ndload, jdltyp, &
               adlmag, predef, npredf, lflags, mlvarx, ddlmag, mdload, &
               pnewdt, jprops, njprop, period)

  use mod_constants, only: dp, ZERO
  use uel_config, only: nIntPt
  use uel_element, only: u3d8

  implicit none

  ! Variables passed into the UEL
  integer, intent(in)  :: ndofel, nrhs, nsvars, nprops, mcrd, nnode, jtype
  integer, intent(in)  :: kstep, kinc, jelem, ndload, npredf, mlvarx
  integer, intent(in)  :: mdload, njprop
  integer, intent(in)  :: jdltyp(mdload,1), lflags(4), jprops(njprop)

  ! Variables defined here, passed back to ABAQUS
  real(dp), intent(inout) :: rhs(mlvarx,1), amatrx(ndofel,ndofel)
  real(dp), intent(inout) :: svars(nsvars), energy(8), pnewdt

  real(dp), intent(in) :: props(nprops), coords(mcrd,nnode), u(ndofel)
  real(dp), intent(in) :: du(mlvarx,1), v(ndofel), a(ndofel)
  real(dp), intent(in) :: time(2), dtime, params(1), adlmag(mdload,1)
  real(dp), intent(in) :: predef(2,npredf,nnode), ddlmag(mdload,1), period

  integer :: nDim

  !----------------------------------------------------------------------
  ! Perform initial checks
  !
  ! Check the procedure type: static general (1,2) or coupled
  ! displacement procedures (64,65,72,73) are supported, the element
  ! itself being purely mechanical
  select case (lflags(1))
  case (1, 2, 64, 65, 72, 73)
    ! all is good
  case default
    write(*,*) 'Abaqus does not have the right procedure'
    write(*,*) 'go back and check the procedure type'
    write(*,*) 'lflags(1)=', lflags(1)
    call exit(1)
  end select

  ! Make sure Abaqus knows you are doing a large deformation problem
  if (lflags(2) == 0) then
    !
    ! lflags(2)=0 -> small disp.
    ! lflags(2)=1 -> large disp.
    !
    write(*,*) 'Abaqus thinks you are doing'
    write(*,*) 'a small displacement analysis'
    write(*,*) 'go in and set nlgeom=yes'
    call exit(1)
  end if

  ! Check if it is a general step or a linear perturbation step
  if (lflags(4) == 1) then
    !
    ! lflags(4)=0 -> general step
    ! lflags(4)=1 -> linear perturbation step
    !
    write(*,*) 'Abaqus thinks you are doing'
    write(*,*) 'a linear perturbation step'
    call exit(1)
  end if

  ! Do nothing if a "dummy" step
  if (dtime == ZERO) return
  !
  ! Done with initial checks
  !----------------------------------------------------------------------

  if (jtype == 3) then
    !
    ! This is a 3D analysis
    !
    nDim = 3
    call u3d8(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
              props, nprops, coords, mcrd, nnode, u, du, v, a, jtype, &
              time, dtime, kstep, kinc, jelem, params, ndload, jdltyp, &
              adlmag, predef, npredf, lflags, mlvarx, ddlmag, mdload, &
              pnewdt, jprops, njprop, period, nDim, nIntPt)
  else
    !
    ! We have a problem...
    !
    write(*,*) 'Element type not supported, jtype=', jtype
    call exit(1)
  end if
  !
  ! Done with this element; RHS and AMATRX are returned as output from
  ! the element routine called above

end subroutine uel


subroutine uvarm(uvar, direct, t, time, dtime, cmname, orname, &
                 nuvarm, noel, npt, layer, kspt, kstep, kinc, ndi, nshr, &
                 coord, jmac, jmatyp, matlayo, laccfla)

  ! This subroutine is used to transfer SDV's from the UEL onto the dummy
  ! mesh for viewing. Note that an offset of ElemOffset is used between
  ! the real mesh and the dummy mesh.

  use uel_config, only: globalSdv, ElemOffset

  implicit none

  character*80 :: cmname, orname
  integer :: nuvarm, noel, npt, layer, kspt, kstep, kinc, ndi, nshr
  integer :: matlayo, laccfla
  integer :: jmac(*), jmatyp(*)
  double precision :: uvar(nuvarm), direct(3,3), t(3,3), time(2), dtime
  double precision :: coord(*)

  integer :: i1

  do i1 = 1, nuvarm
    uvar(i1) = globalSdv(noel-ElemOffset, npt, i1)
  end do

end subroutine uvarm
