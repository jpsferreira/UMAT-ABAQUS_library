!> @brief Composable UMAT builder — a single UMAT entry point that assembles
!>        material behaviour from building blocks selected at runtime via PROPS.
!>
!> ============================================================================
!> PROPS LAYOUT  (see PROPS_REFERENCE.md for full documentation)
!> ============================================================================
!>
!>  PROPS(1)  = KBULK         Bulk modulus
!>  PROPS(2)  = ISO_TYPE       0=none, 1=NH, 2=MR, 3=Ogden, 4=Humphrey
!>  PROPS(3)  = ANISO_TYPE     0=none, 1=HGO, 2=Humphrey-fiber
!>  PROPS(4)  = N_FIBER_FAM    Number of fiber families (0, 1, or 2)
!>  PROPS(5)  = NETWORK_TYPE   0=none, 1=affine, 2=non-affine
!>  PROPS(6)  = DAMAGE_TYPE    0=none, 1=sigmoid
!>  PROPS(7)  = N_VISCO        Number of Maxwell branches (0 = none)
!>  PROPS(8:) = [iso] [aniso x N_FAM] [network] [damage] [visco x N_VISCO]
!>
!> ============================================================================
!> STATE VARIABLES
!> ============================================================================
!>   STATEV(1)                     = det(F)
!>   STATEV(2)                     = damage variable d       (if damage)
!>   STATEV(3)                     = max historical SEF      (if damage)
!>   STATEV(sdv_visco : +9*N_V-1) = hidden stress tensors    (if visco)
!> ============================================================================

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

  ! --------------------------------------------------------------------------
  ! ABAQUS UMAT interface — DOUBLE PRECISION to match ABAQUS calling convention.
  ! INTENT(IN OUT) on most arguments matches the existing F90 UMAT convention.
  ! --------------------------------------------------------------------------
  integer, intent(in out) :: ndi, nshr, ntens, nstatev, nprops
  integer, intent(in out) :: noel, npt, layer, kspt, kstep, kinc

  double precision, intent(in out) :: stress(ntens)
  double precision, intent(in out) :: statev(nstatev)
  double precision, intent(in out) :: ddsdde(ntens, ntens)
  double precision, intent(in out) :: sse, spd, scd
  double precision, intent(in out) :: rpl, drpldt
  double precision, intent(in out) :: ddsddt(ntens), drplde(ntens)
  double precision, intent(in out) :: stran(ntens), dstran(ntens)
  double precision, intent(in out) :: time(2)
  double precision, intent(in out) :: dtime, temp, dtemp
  double precision, intent(in out) :: predef(1), dpred(1)
  character(len=8), intent(in out) :: cmname
  double precision, intent(in)     :: props(nprops)
  double precision, intent(in out) :: coords(3)
  double precision, intent(in out) :: drot(3,3)
  double precision, intent(in out) :: pnewdt, celent
  double precision, intent(in out) :: dfgrd0(3,3), dfgrd1(3,3)

  ! --- Configuration from PROPS header ---
  real(dp) :: kbulk
  integer  :: iso_type, aniso_type, n_fiber_fam, network_type, damage_type, n_visco
  integer  :: ip  ! parameter index pointer

  ! --- Identity and projection tensors ---
  real(dp) :: unit2(3,3), unit4(3,3,3,3), unit4s(3,3,3,3)
  real(dp) :: proje(3,3,3,3), projl(3,3,3,3)

  ! --- Kinematics ---
  real(dp) :: distgr(3,3), c(3,3), b(3,3), cbar(3,3), bbar(3,3)
  real(dp) :: ubar(3,3), vbar(3,3)
  real(dp) :: det, cbari1, cbari2

  ! --- Volumetric ---
  real(dp) :: ssev, pv, ppv
  real(dp) :: pkvol(3,3), svol(3,3), cmvol(3,3,3,3), cvol(3,3,3,3)

  ! --- Isochoric isotropic ---
  real(dp) :: sseiso, diso(5)
  real(dp) :: pkmatfic(3,3), sisomatfic(3,3)
  real(dp) :: cmisomatfic(3,3,3,3), cisomatfic(3,3,3,3)
  integer  :: n_ogden_terms

  ! --- Isochoric anisotropic ---
  real(dp) :: sseaniso, daniso(4)
  real(dp) :: pkmatficaniso(3,3), sanisomatfic(3,3)
  real(dp) :: cmanisomatfic(3,3,3,3), canisomatfic(3,3,3,3)
  real(dp) :: vorif(3), vd(3), st0(3,3)
  real(dp) :: cbari4, lambda_fib, barlambda
  real(dp) :: k1_fib, k2_fib, kdisp_fib
  integer  :: ifam

  ! --- Network ---
  real(dp) :: phi_net, net_density, b_orient, efi, pp_naff
  real(dp) :: filprops(6)
  real(dp) :: snet(3,3), cnet(3,3,3,3)

  ! --- Combined fictitious tensors ---
  real(dp) :: pkfic(3,3), sfic(3,3), cmfic(3,3,3,3), cfic(3,3,3,3)

  ! --- Isochoric projected ---
  real(dp) :: pkiso(3,3), siso(3,3), cmiso(3,3,3,3), ciso(3,3,3,3)

  ! --- Total ---
  real(dp) :: sigma(3,3), ddsigdde(3,3,3,3)
  real(dp) :: cjr(3,3,3,3)

  ! --- Damage ---
  real(dp) :: dmg, dmg_red, dmg_diff, sef0_hist, beta_d, psi_half

  ! --- Viscosity ---
  real(dp) :: vscprops(2*MAX_VISCO_BRANCHES)
  real(dp) :: pk2_visc(3,3), cmat_visc(3,3,3,3)
  real(dp) :: sigma_visc(3,3), cspatial_visc(3,3,3,3)
  integer  :: sdv_offset_visco

  ! --- Temporaries ---
  integer :: i

  ! --- Interface for external network subroutine ---
  interface
    subroutine network_contribution(network_type, props, ip, distgr, det, &
                                    phi_net, snet, cnet)
      use mod_constants, only: dp
      implicit none
      integer,  intent(in)    :: network_type
      double precision, intent(in) :: props(*)
      integer,  intent(inout) :: ip
      real(dp), intent(in)    :: distgr(3,3), det
      real(dp), intent(out)   :: phi_net, snet(3,3), cnet(3,3,3,3)
    end subroutine network_contribution
  end interface

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
  snet = ZERO; cnet = ZERO

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
    ! no isotropic contribution
  end select

  ! Isotropic fictitious stress and stiffness (if any iso model active)
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
      k1_fib    = props(ip)
      k2_fib    = props(ip+1)
      kdisp_fib = props(ip+2)
      vorif     = props(ip+3:ip+5)
      ip = ip + 6

      call fiber_direction(vorif, distgr, st0, vd)
      call pseudo_invariants(cbar, st0, det, cbari4, lambda_fib, barlambda)
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

    ! Accumulate into total fictitious tensors
    pkfic = pkfic + pkmatficaniso
    sfic  = sfic  + sanisomatfic
    cmfic = cmfic + cmanisomatfic
    cfic  = cfic  + canisomatfic
  end do

  ! --- Network contribution ---
  ! Network models require quadrature data loaded via UEXTERNALDB.
  ! They produce their own fictitious stress/stiffness in the spatial frame,
  ! blended with the matrix contribution using volume fraction PHI.
  if (network_type > 0) then
    call network_contribution(network_type, props, ip, distgr, det, &
                              phi_net, snet, cnet)
    ! Blend: total = (1-PHI)*matrix + PHI*network
    ! Network stress/stiffness are already in fictitious spatial frame
    sfic  = (ONE - phi_net) * sfic  + snet
    cfic  = (ONE - phi_net) * cfic  + cnet
    pkfic = (ONE - phi_net) * pkfic  ! PK2 only has matrix part
  end if

  ! --- Damage modification ---
  if (damage_type == DMG_SIGMOID) then
    beta_d   = props(ip)
    psi_half = props(ip+1)
    ip = ip + 2

    dmg       = statev(2)
    sef0_hist = statev(3)

    call damage_sigmoid(sseiso + sseaniso, sef0_hist, dmg, dmg_red, dmg_diff, beta_d, psi_half)

    statev(2) = dmg
    statev(3) = sef0_hist

    ! Apply damage reduction to fictitious tensors
    ! Consistent tangent: C_damaged = (1-d)*C + d'*S (x) dW/dC
    ! Simplified: scale both stress and stiffness by damage reduction
    pkfic = dmg_red * pkfic
    sfic  = dmg_red * sfic
    cmfic = dmg_red * cmfic
    cfic  = dmg_red * cfic
  end if

  ! ============================================================================
  ! STRESS MEASURES (elastic, no viscosity yet)
  ! ============================================================================

  ! --- Volumetric ---
  call pk2_vol(pkvol, pv, c)
  call sig_vol(svol, pv, unit2)

  ! --- Isochoric (projected) ---
  call pk2_iso(pkiso, pkfic, projl, det)
  call sig_iso(siso, sfic, proje)

  ! --- Total elastic Cauchy stress ---
  sigma = svol + siso

  ! ============================================================================
  ! ELASTICITY TENSORS (elastic)
  ! ============================================================================

  ! --- Material elasticity ---
  call met_vol(cmvol, c, pv, ppv, det)
  call met_iso(cmiso, cmfic, projl, pkiso, pkfic, c, unit2, det)

  ! --- Spatial elasticity ---
  call set_vol(cvol, pv, ppv, unit2, unit4s)
  call set_iso(ciso, cfic, proje, siso, sfic, unit2)
  call set_jr(cjr, sigma, unit2)
  ddsigdde = cvol + ciso + cjr

  ! ============================================================================
  ! VISCOUS CONTRIBUTION (optional)
  ! ============================================================================
  ! The viscous model operates on the material frame (PK2 + material tangent).
  ! It modifies PK2 by adding viscous overstress and scales the isochoric
  ! material tangent. We then push-forward the result to the spatial frame.
  if (n_visco > 0) then
    vscprops = ZERO
    do i = 1, min(n_visco, MAX_VISCO_BRANCHES)
      vscprops(2*i-1) = props(ip)      ! tau
      vscprops(2*i)   = props(ip+1)    ! theta
      ip = ip + 2
    end do

    ! Offset: after det (1) + damage vars (2 if active)
    sdv_offset_visco = 1
    if (damage_type > 0) sdv_offset_visco = 3

    ! visco_maxwell returns modified PK2 (pk2_visc) and material tangent (cmat_visc)
    call visco_maxwell(pk2_visc, cmat_visc, n_visco, pkvol, pkiso, &
                       cmvol, cmiso, dtime, vscprops, statev, sdv_offset_visco)

    ! Push-forward to get spatial quantities
    call push2(sigma_visc, pk2_visc, dfgrd1, det)
    call push4(cspatial_visc, cmat_visc, dfgrd1, det)

    ! Jaumann rate correction with viscous Cauchy stress
    call set_jr(cjr, sigma_visc, unit2)

    ! Override elastic results with viscous-augmented versions
    sigma    = sigma_visc
    ddsigdde = cspatial_visc + cjr
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


! ============================================================================
! INTERNAL: Network dispatch (reads PROPS, calls the appropriate network model)
! ============================================================================
subroutine network_contribution(network_type, props, ip, distgr, det, &
                                phi_net, snet, cnet)
  use mod_constants, only: dp, ZERO, ONE
  use mod_network,   only: affine_network, nonaffine_network

  implicit none
  integer,  intent(in)    :: network_type
  double precision, intent(in) :: props(*)
  integer,  intent(inout) :: ip
  real(dp), intent(in)    :: distgr(3,3), det
  real(dp), intent(out)   :: phi_net, snet(3,3), cnet(3,3,3,3)

  real(dp) :: net_density, b_orient, efi, pp_naff
  real(dp) :: filprops(6)

  ! Quadrature data loaded via UEXTERNALDB into COMMON blocks
  integer, parameter :: MAX_NWP = 720
  real(dp) :: mf0(MAX_NWP, 3), rw(MAX_NWP)
  integer  :: nwp_active
  common /kfil/  mf0
  common /kfilr/ rw
  common /knwp/  nwp_active

  snet = ZERO
  cnet = ZERO

  select case (network_type)
  case (1) ! Affine
    phi_net     = props(ip)
    net_density = props(ip+1)
    b_orient    = props(ip+2)
    efi         = props(ip+3)
    filprops(1) = props(ip+4)   ! L
    filprops(2) = props(ip+5)   ! R0
    filprops(3) = props(ip+6)   ! mu0
    filprops(4) = props(ip+7)   ! beta
    filprops(5) = props(ip+8)   ! B0
    filprops(6) = props(ip+9)   ! lambda0
    ip = ip + 10

    call affine_network(snet, cnet, distgr, mf0, rw, nwp_active, det, &
                        filprops, net_density, b_orient, efi)

  case (2) ! Non-affine
    phi_net     = props(ip)
    net_density = props(ip+1)
    b_orient    = props(ip+2)
    efi         = props(ip+3)
    pp_naff     = props(ip+4)   ! non-affinity exponent
    filprops(1) = props(ip+5)
    filprops(2) = props(ip+6)
    filprops(3) = props(ip+7)
    filprops(4) = props(ip+8)
    filprops(5) = props(ip+9)
    filprops(6) = props(ip+10)
    ip = ip + 11

    call nonaffine_network(snet, cnet, distgr, mf0, rw, nwp_active, det, &
                           filprops, net_density, b_orient, efi, pp_naff)

  case default
    phi_net = ZERO
  end select
end subroutine network_contribution
