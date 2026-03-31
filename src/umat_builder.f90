!> @brief Composable UMAT builder — a single UMAT entry point that assembles
!>        material behaviour from building blocks selected at runtime via PROPS.
!>
!> Instead of maintaining separate UMAT files for every combination of
!> isotropic + anisotropic + network + damage + viscosity, the user specifies
!> which components to combine through the ABAQUS material property array.
!>
!> ============================================================================
!> PROPS LAYOUT
!> ============================================================================
!>
!>  PROPS(1)  = KBULK         Bulk modulus (required)
!>  PROPS(2)  = ISO_TYPE       Isotropic model: 0=none, 1=NH, 2=MR, 3=Ogden, 4=Humphrey
!>  PROPS(3)  = ANISO_TYPE     Anisotropic model: 0=none, 1=HGO, 2=Humphrey-fiber
!>  PROPS(4)  = N_FIBER_FAM    Number of fiber families (0, 1, or 2)
!>  PROPS(5)  = NETWORK_TYPE   Network model: 0=none, 1=affine, 2=non-affine
!>  PROPS(6)  = DAMAGE_TYPE    Damage model: 0=none, 1=sigmoid
!>  PROPS(7)  = N_VISCO        Number of Maxwell branches (0 = no viscosity)
!>
!>  PROPS(8:) = Material parameters, packed sequentially:
!>
!>    [isotropic params] [aniso params per family] [network params] [damage params] [visco params]
!>
!>  --- Isotropic parameters ---
!>    NH:       C10
!>    MR:       C10, C01
!>    Ogden:    N_TERMS, mu_1, alpha_1, mu_2, alpha_2, ...
!>    Humphrey: C10, C01
!>
!>  --- Anisotropic parameters (repeated per fiber family) ---
!>    HGO:            K1, K2, KDISP,  fiber_x, fiber_y, fiber_z
!>    Humphrey-fiber: K1, K2,         fiber_x, fiber_y, fiber_z
!>
!>  --- Network parameters ---
!>    Affine:     PHI, N, B_orient, EFI, L, R0, mu0, beta, B0, lambda0
!>    Non-affine: PHI, N, B_orient, EFI, PP, L, R0, mu0, beta, B0, lambda0
!>    (PHI = network volume fraction; remaining = filament + network params)
!>
!>  --- Damage parameters ---
!>    Sigmoid: beta_d, psi_half
!>
!>  --- Viscous parameters (per branch) ---
!>    tau_1, theta_1, tau_2, theta_2, ...
!>
!> ============================================================================
!> STATE VARIABLES (NSTATEV)
!> ============================================================================
!>   STATEV(1)         = Jacobian determinant
!>   STATEV(2)         = Damage variable (if damage active)
!>   STATEV(3)         = Max historical strain energy (if damage active)
!>   STATEV(4:4+9*VV)  = Hidden stress tensors for viscous branches
!> ============================================================================

! NOTE: This module uses ABAQUS double precision via aba_param.inc.
! When aba_param.inc is included, ABAQUS maps "real(8)" to its internal
! double-precision type. Our dp parameter matches this.

subroutine umat(stress, statev, ddsdde, sse, spd, scd, &
     rpl, ddsddt, drplde, drpldt, &
     stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
     ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
     celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

  use mod_constants
  use mod_tensor,       only: onem, push2, push4
  use mod_kinematics,   only: distortion_gradient, cauchy_green, invariants, &
                               pseudo_invariants, stretch_tensors, rotation_tensor, &
                               projection_euler, projection_lagrange
  use mod_continuum,    only: vol_energy, pk2_vol, sig_vol, met_vol, set_vol, &
                               pk2_isomatfic, sig_isomatfic, cmat_isomatfic, cs_isomatfic, &
                               pk2_iso, sig_iso, met_iso, set_iso, &
                               pk2_anisofic, cmat_anisofic, set_jr, indexx
  use mod_hyperelastic, only: sef_neo_hooke, sef_mooney_rivlin, sef_humphrey, sef_ogden
  use mod_anisotropic,  only: fiber_direction, sef_aniso_hgo, sef_aniso_humphrey
  use mod_damage,       only: damage_sigmoid
  use mod_viscosity,    only: visco_maxwell, MAX_VISCO_BRANCHES

  implicit none

  ! ABAQUS-required interface
  character(len=8), intent(in) :: cmname
  integer, intent(in) :: ndi, nshr, ntens, nstatev, nprops, noel, npt, &
                          layer, kspt, kstep, kinc
  real(dp), intent(inout) :: stress(ntens), statev(nstatev), ddsdde(ntens,ntens)
  real(dp), intent(inout) :: sse, spd, scd, rpl, drpldt, pnewdt
  real(dp), intent(in)    :: ddsddt(ntens), drplde(ntens)
  real(dp), intent(in)    :: stran(ntens), dstran(ntens), time(2), predef(1), dpred(1)
  real(dp), intent(in)    :: props(nprops), coords(3), drot(3,3)
  real(dp), intent(in)    :: dfgrd0(3,3), dfgrd1(3,3)
  real(dp), intent(in)    :: dtime, temp, dtemp, celent

  ! --- Configuration from PROPS header ---
  real(dp) :: kbulk
  integer  :: iso_type, aniso_type, n_fiber_fam, network_type, damage_type, n_visco
  integer  :: ip  ! parameter index pointer

  ! --- Identity and projection tensors ---
  real(dp) :: unit2(3,3), unit4(3,3,3,3), unit4s(3,3,3,3)
  real(dp) :: proje(3,3,3,3), projl(3,3,3,3)

  ! --- Kinematics ---
  real(dp) :: distgr(3,3), c(3,3), b(3,3), cbar(3,3), bbar(3,3)
  real(dp) :: distgrinv(3,3), dfgrd1inv(3,3)
  real(dp) :: ubar(3,3), vbar(3,3), rot(3,3)
  real(dp) :: det, cbari1, cbari2

  ! --- Volumetric ---
  real(dp) :: ssev, pv, ppv
  real(dp) :: pkvol(3,3), svol(3,3), cmvol(3,3,3,3), cvol(3,3,3,3)

  ! --- Isochoric isotropic ---
  real(dp) :: sseiso, diso(5)
  real(dp) :: pkmatfic(3,3), sisomatfic(3,3)
  real(dp) :: cmisomatfic(3,3,3,3), cisomatfic(3,3,3,3)
  real(dp) :: iso_params(20)  ! enough for Ogden with many terms
  integer  :: n_ogden_terms

  ! --- Isochoric anisotropic ---
  real(dp) :: sseaniso, daniso(4)
  real(dp) :: pkmatficaniso(3,3), sanisomatfic(3,3)
  real(dp) :: cmanisomatfic(3,3,3,3), canisomatfic(3,3,3,3)
  real(dp) :: vorif(3), vd(3), st0(3,3)
  real(dp) :: cbari4, lambda_fib, barlambda
  real(dp) :: k1_fib, k2_fib, kdisp_fib
  integer  :: ifam

  ! --- Combined fictitious tensors ---
  real(dp) :: pkfic(3,3), sfic(3,3), cmfic(3,3,3,3), cfic(3,3,3,3)

  ! --- Isochoric projected ---
  real(dp) :: pkiso(3,3), siso(3,3), cmiso(3,3,3,3), ciso(3,3,3,3)

  ! --- Total ---
  real(dp) :: pk2_total(3,3), sigma(3,3)
  real(dp) :: ddsigdde(3,3,3,3), ddpkdde(3,3,3,3)
  real(dp) :: cjr(3,3,3,3)

  ! --- Damage ---
  real(dp) :: dmg, dmg_red, dmg_diff, sef0_hist, beta_d, psi_half

  ! --- Viscosity ---
  real(dp) :: vscprops(6)
  real(dp) :: pk2_visc(3,3), cmat_visc(3,3,3,3)
  integer  :: sdv_offset_visco

  ! --- Temporaries ---
  integer :: i, j

  ! ============================================================================
  ! READ CONFIGURATION FROM PROPS
  ! ============================================================================
  kbulk       = props(1)
  iso_type    = nint(props(2))
  aniso_type  = nint(props(3))
  n_fiber_fam = nint(props(4))
  network_type= nint(props(5))
  damage_type = nint(props(6))
  n_visco     = nint(props(7))
  ip = 8  ! pointer to start of material parameters

  ! ============================================================================
  ! INITIALIZATIONS
  ! ============================================================================
  call onem(unit2, unit4, unit4s)

  ! Zero all working arrays
  pkmatfic = ZERO; sisomatfic = ZERO
  cmisomatfic = ZERO; cisomatfic = ZERO
  pkmatficaniso = ZERO; sanisomatfic = ZERO
  cmanisomatfic = ZERO; canisomatfic = ZERO
  pkfic = ZERO; sfic = ZERO; cmfic = ZERO; cfic = ZERO
  sigma = ZERO; ddsigdde = ZERO
  sseiso = ZERO; diso = ZERO
  sseaniso = ZERO; daniso = ZERO

  ! Initialize state variables on first increment
  if (time(1) == ZERO .and. kstep == 1) then
    statev(1) = ONE
    if (damage_type > 0 .and. nstatev >= 3) then
      statev(2) = ZERO  ! damage
      statev(3) = ZERO  ! max SEF
    end if
  end if

  ! ============================================================================
  ! KINEMATICS
  ! ============================================================================
  call distortion_gradient(dfgrd1, distgr, det)
  call cauchy_green(dfgrd1, c, b)
  call cauchy_green(distgr, cbar, bbar)
  call invariants(cbar, cbari1, cbari2)
  call stretch_tensors(cbar, bbar, ubar, vbar)
  call rotation_tensor(distgr, ubar, rot)
  call projection_euler(unit2, unit4s, proje)
  call projection_lagrange(c, unit4, projl)

  ! ============================================================================
  ! CONSTITUTIVE RELATIONS
  ! ============================================================================

  ! --- Volumetric contribution (always present) ---
  call vol_energy(kbulk, det, ssev, pv, ppv)

  ! --- Isochoric isotropic contribution ---
  select case (iso_type)
  case (ISO_NEO_HOOKE)
    sseiso = ZERO; diso = ZERO
    call sef_neo_hooke(sseiso, diso, props(ip), cbari1, cbari2)
    ip = ip + 1

  case (ISO_MOONEY_RIVLIN)
    call sef_mooney_rivlin(sseiso, diso, props(ip), props(ip+1), cbari1, cbari2)
    ip = ip + 2

  case (ISO_OGDEN)
    n_ogden_terms = nint(props(ip))
    ip = ip + 1
    call sef_ogden(sseiso, diso, c, cbar, props(ip:ip+2*n_ogden_terms-1), n_ogden_terms)
    ip = ip + 2*n_ogden_terms

  case (ISO_HUMPHREY)
    call sef_humphrey(sseiso, diso, props(ip), props(ip+1), cbari1, cbari2)
    ip = ip + 2

  case default
    sseiso = ZERO; diso = ZERO
  end select

  ! Compute isotropic fictitious stress and stiffness
  if (iso_type > 0) then
    call pk2_isomatfic(pkmatfic, diso, cbar, cbari1, unit2)
    call sig_isomatfic(sisomatfic, pkmatfic, distgr, det)
    call cmat_isomatfic(cmisomatfic, cbar, cbari1, cbari2, diso, unit2, unit4, det)
    call cs_isomatfic(cisomatfic, cmisomatfic, distgr, det)
  end if

  ! Start accumulating fictitious tensors
  pkfic  = pkmatfic
  sfic   = sisomatfic
  cmfic  = cmisomatfic
  cfic   = cisomatfic

  ! --- Anisotropic contribution (fiber families) ---
  do ifam = 1, n_fiber_fam
    select case (aniso_type)
    case (ANISO_HGO)
      k1_fib   = props(ip)
      k2_fib   = props(ip+1)
      kdisp_fib= props(ip+2)
      vorif    = props(ip+3:ip+5)
      ip = ip + 6

      ! Fiber direction and pseudo-invariant
      call fiber_direction(vorif, distgr, st0, vd)
      call pseudo_invariants(cbar, st0, det, cbari4, lambda_fib, barlambda)

      ! Anisotropic SEF (modifies diso for coupled terms)
      call sef_aniso_hgo(sseaniso, daniso, diso, k1_fib, k2_fib, kdisp_fib, cbari4, cbari1)

    case (ANISO_HUMPHREY)
      k1_fib = props(ip)
      k2_fib = props(ip+1)
      vorif  = props(ip+2:ip+4)
      ip = ip + 5

      call fiber_direction(vorif, distgr, st0, vd)
      call pseudo_invariants(cbar, st0, det, cbari4, lambda_fib, barlambda)

      call sef_aniso_humphrey(sseaniso, daniso, diso, k1_fib, k2_fib, cbari4, cbari1)

    case default
      cycle
    end select

    ! Anisotropic fictitious stress and stiffness
    call pk2_anisofic(pkmatficaniso, daniso, st0)
    call push2(sanisomatfic, pkmatficaniso, distgr, det)
    call cmat_anisofic(cmanisomatfic, st0, daniso, unit2, det)
    call push4(canisomatfic, cmanisomatfic, distgr, det)

    ! Accumulate
    pkfic = pkfic + pkmatficaniso
    sfic  = sfic  + sanisomatfic
    cmfic = cmfic + cmanisomatfic
    cfic  = cfic  + canisomatfic
  end do

  ! --- Network contribution ---
  ! (Network models operate in the spatial frame and produce their own
  !  fictitious stress/stiffness. They are combined with the matrix
  !  contribution using a volume fraction PHI.)
  ! For network models, the user provides PHI in the network params,
  ! and the matrix contribution is scaled by (1-PHI).
  ! This section can be extended by calling affine_network / nonaffine_network
  ! from mod_network and blending with the matrix contribution.
  ! Currently, the builder focuses on the soft-tissue models;
  ! network integration requires quadrature data loaded via UEXTERNALDB.

  ! --- Damage modification ---
  if (damage_type == DMG_SIGMOID) then
    beta_d   = props(ip)
    psi_half = props(ip+1)
    ip = ip + 2

    dmg      = statev(2)
    sef0_hist= statev(3)

    call damage_sigmoid(sseiso + sseaniso, sef0_hist, dmg, dmg_red, dmg_diff, beta_d, psi_half)

    ! Update state
    statev(2) = dmg
    statev(3) = sef0_hist

    ! Apply damage to fictitious tensors
    pkfic = dmg_red * pkfic
    sfic  = dmg_red * sfic
    cmfic = dmg_red * cmfic + dmg_diff * cmfic  ! consistent tangent correction
    cfic  = dmg_red * cfic  + dmg_diff * cfic
  end if

  ! ============================================================================
  ! STRESS MEASURES
  ! ============================================================================

  ! --- Volumetric ---
  call pk2_vol(pkvol, pv, c)
  call sig_vol(svol, pv, unit2)

  ! --- Isochoric (projected) ---
  call pk2_iso(pkiso, pkfic, projl, det)
  call sig_iso(siso, sfic, proje)

  ! --- Total ---
  pk2_total = pkvol + pkiso
  sigma     = svol  + siso

  ! ============================================================================
  ! ELASTICITY TENSORS
  ! ============================================================================

  ! --- Material elasticity ---
  call met_vol(cmvol, c, pv, ppv, det)
  call met_iso(cmiso, cmfic, projl, pkiso, pkfic, c, unit2, det)
  ddpkdde = cmvol + cmiso

  ! --- Spatial elasticity ---
  call set_vol(cvol, pv, ppv, unit2, unit4s)
  call set_iso(ciso, cfic, proje, siso, sfic, unit2)
  call set_jr(cjr, sigma, unit2)
  ddsigdde = cvol + ciso + cjr

  ! ============================================================================
  ! VISCOUS CONTRIBUTION (optional, operates on material frame)
  ! ============================================================================
  if (n_visco > 0) then
    vscprops = ZERO
    do i = 1, min(n_visco, MAX_VISCO_BRANCHES)
      vscprops(2*i-1) = props(ip)      ! tau
      vscprops(2*i)   = props(ip+1)    ! theta
      ip = ip + 2
    end do

    sdv_offset_visco = 3  ! after det, damage, sef0
    if (damage_type == 0) sdv_offset_visco = 1  ! after det only

    call visco_maxwell(pk2_visc, cmat_visc, n_visco, pkvol, pkiso, &
                       cmvol, cmiso, dtime, vscprops, statev, sdv_offset_visco)

    ! Push viscous material tangent to spatial frame for ABAQUS
    call push4(cvol, cmat_visc, dfgrd1, det)  ! reuse cvol as temp
    call set_jr(cjr, sigma, unit2)  ! recompute if stress changed

    ! For viscous case, override total with viscous-augmented version
    pk2_total = pk2_visc
    ! Spatial stiffness: push-forward the total material tangent + JR
    call push4(ddsigdde, cmat_visc, dfgrd1, det)
    call push2(sigma, pk2_total, dfgrd1, det)
    call set_jr(cjr, sigma, unit2)
    ddsigdde = ddsigdde + cjr
  end if

  ! ============================================================================
  ! STRAIN ENERGY
  ! ============================================================================
  sse = ssev + sseiso + sseaniso

  ! ============================================================================
  ! VOIGT INDEXING (ABAQUS interface)
  ! ============================================================================
  call indexx(stress, ddsdde, sigma, ddsigdde, ntens)

  ! ============================================================================
  ! STATE VARIABLES
  ! ============================================================================
  statev(1) = det

end subroutine umat
