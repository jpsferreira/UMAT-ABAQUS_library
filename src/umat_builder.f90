!> @brief Composable material builder — one constitutive core that assembles
!>        material behaviour from building blocks selected at runtime via PROPS,
!>        exposed through two thin interface wrappers:
!>
!>          umat(...)          ABAQUS *USER MATERIAL entry point. Adds the
!>                             Jaumann-rate term and collapses to Voigt
!>                             (STRESS/DDSDDE), as ABAQUS built-in elements
!>                             require.
!>          material_uel(...)  UEL-facing entry point (module mod_material).
!>                             Returns the full Cauchy stress sigma(3,3) and
!>                             the assembly tangent
!>                               a_ijkl = c_ijkl + delta_ik sigma_jl
!>                             (push-forward of dP/dF) that a displacement
!>                             UEL needs for Kuu = int G^T A G dv. Note this
!>                             geometric term is NOT the Jaumann correction.
!>
!> ============================================================================
!> PROPS LAYOUT  (see PROPS_REFERENCE.md for full documentation)
!> ============================================================================
!>
!>  PROPS(1)  = KBULK         Bulk modulus
!>  PROPS(2)  = ISO_TYPE       0=none, 1=NH, 2=MR, 3=Ogden, 4=Humphrey
!>  PROPS(3)  = ANISO_TYPE     0=none, 1=HGO, 2=Humphrey-fiber,
!>                             3=HGO-AI, 4=Humphrey-AI, 5=Humphrey-activation
!>  PROPS(4)  = N_FIBER_FAM    Number of fiber families (0, 1, or 2)
!>  PROPS(5)  = NETWORK_TYPE   0=none, 1=affine(RW), 2=non-affine(RW), 3=mixed(AI),
!>                             4=contractile(AI), 5=affine(AI), 6=non-affine(AI)
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

module mod_material
  implicit none

contains

  !> Constitutive core shared by the UMAT and UEL wrappers.
  !>
  !> Returns the total Cauchy stress and the BARE spatial elasticity tensor
  !> (volumetric + isochoric, damage-softened, visco-augmented) WITHOUT any
  !> rate-convention term — each wrapper adds its own (set_jr or set_geometric).
  !>
  !> istat: 0 = ok; 1 = inverted/degenerate element (det(F) <= 0), pnewdt has
  !> been set to HALF and sigma/cspatial are zero — the caller should return
  !> without touching its outputs so ABAQUS cuts the increment back.
  subroutine material_core(sigma, cspatial, sse, istat, &
                           statev, nstatev, props, nprops, dfgrd1, &
                           time, dtime, predef, dpred, pnewdt, &
                           noel, npt, kstep, kinc)

    use mod_constants
    use mod_tensor,       only: onem, push2, push4, pull2, pull4, matinv3d
    use mod_kinematics,   only: distortion_gradient, cauchy_green, invariants, &
                                 pseudo_invariants, stretch_tensors, rotation_tensor, &
                                 projection_euler, projection_lagrange
    use mod_continuum,    only: vol_energy, pk2_vol, sig_vol, met_vol, set_vol, &
                                 pk2_isomatfic, sig_isomatfic, cmat_isomatfic, cs_isomatfic, &
                                 pk2_iso, sig_iso, met_iso, set_iso, &
                                 pk2_anisofic, cmat_anisofic
    use mod_hyperelastic, only: sef_neo_hooke, sef_mooney_rivlin, sef_humphrey, sef_ogden
    use mod_anisotropic,  only: fiber_direction, sef_aniso_hgo, sef_aniso_humphrey, &
                                 sef_aniso_humphrey_act, sef_aniso_hgo_ai, &
                                 sef_aniso_humphrey_ai
    use mod_damage,       only: damage_sigmoid
    use mod_viscosity,    only: visco_maxwell, MAX_VISCO_BRANCHES

    implicit none

    real(dp), intent(out)   :: sigma(3,3), cspatial(3,3,3,3)
    double precision, intent(out)   :: sse
    integer,  intent(out)   :: istat
    integer,  intent(in)    :: nstatev, nprops
    double precision, intent(in out) :: statev(nstatev)
    double precision, intent(in)     :: props(nprops)
    double precision, intent(in)     :: dfgrd1(3,3)
    double precision, intent(in)     :: time(2), dtime
    double precision, intent(in)     :: predef(1), dpred(1)
    double precision, intent(in out) :: pnewdt
    integer,  intent(in)    :: noel, npt, kstep, kinc

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
    logical  :: iso_fic_active

    ! --- Isochoric anisotropic ---
    real(dp) :: sseaniso, daniso(4)
    real(dp) :: pkmatficaniso(3,3), sanisomatfic(3,3)
    real(dp) :: cmanisomatfic(3,3,3,3), canisomatfic(3,3,3,3)
    real(dp) :: vorif(3), vd(3), st0(3,3)
    real(dp) :: cbari4, lambda_fib, barlambda
    real(dp) :: k1_fib, k2_fib, kdisp_fib, t0m_act, bdisp_fib
    real(dp) :: sfic_aniso_ai(3,3), cfic_aniso_ai(3,3,3,3), w_aniso_ai
    real(dp) :: finv(3,3)
    integer  :: ifam, factor_aniso

    ! --- Network ---
    real(dp) :: phi_net
    real(dp) :: snet(3,3), cnet(3,3,3,3)

    ! --- Combined fictitious tensors ---
    real(dp) :: pkfic(3,3), sfic(3,3), cmfic(3,3,3,3), cfic(3,3,3,3)

    ! --- Isochoric projected ---
    real(dp) :: pkiso(3,3), siso(3,3), cmiso(3,3,3,3), ciso(3,3,3,3)

    ! --- Damage ---
    real(dp) :: dmg, dmg_red, dmg_diff, sef0_hist, beta_d, psi_half
    real(dp) :: pkiso0(3,3), siso0(3,3)   ! undamaged isochoric stress (for consistent tangent)
    integer  :: di, dj, dk, dl

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
                                      phi_net, snet, cnet, &
                                      statev, nstatev, dtime, time, kstep, noel, &
                                      damage_type, n_visco)
        use mod_constants, only: dp
        implicit none
        integer,  intent(in)    :: network_type
        double precision, intent(in) :: props(*)
        integer,  intent(inout) :: ip
        real(dp), intent(in)    :: distgr(3,3), det
        real(dp), intent(out)   :: phi_net, snet(3,3), cnet(3,3,3,3)
        double precision, intent(inout) :: statev(*)
        integer,  intent(in)    :: nstatev, kstep, noel
        double precision, intent(in) :: dtime, time(2)
        integer,  intent(in)    :: damage_type, n_visco
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
    istat = 0
    sse   = ZERO
    call onem(unit2, unit4, unit4s)

    ! Zero all working arrays
    pkmatfic = ZERO; sisomatfic = ZERO
    cmisomatfic = ZERO; cisomatfic = ZERO
    pkmatficaniso = ZERO; sanisomatfic = ZERO
    cmanisomatfic = ZERO; canisomatfic = ZERO
    pkfic = ZERO; sfic = ZERO; cmfic = ZERO; cfic = ZERO
    sigma = ZERO; cspatial = ZERO
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
      ! Explicitly zero the viscous hidden-stress tensors (9 per branch). ABAQUS
      ! zero-inits STATEV, but being explicit guards restarts / re-used arrays.
      if (n_visco > 0) then
        sdv_offset_visco = 1
        if (damage_type > 0) sdv_offset_visco = 3
        do i = sdv_offset_visco + 1, min(sdv_offset_visco + 9*n_visco, nstatev)
          statev(i) = ZERO
        end do
      end if
    end if

    ! ============================================================================
    ! KINEMATICS
    ! ============================================================================
    call distortion_gradient(dfgrd1, distgr, det)

    ! Guard against an inverted/degenerate element (det(F) <= 0). det**(-1/3) and
    ! C^{-1} are undefined there (NaN); rather than poison the Newton iteration,
    ! ask ABAQUS to cut the increment back and return.
    if (det <= ZERO) then
      pnewdt = HALF
      istat = 1
      return
    end if

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

    ! Whether the isotropic fictitious tensors must be built. An HGO fiber family
    ! with dispersion (kdisp > 0) feeds the I1 channel back into diso, so this is
    ! also switched on inside the fiber loop below.
    iso_fic_active = (iso_type > 0)

    ! Fictitious totals start at ZERO (already zeroed during initialization). The
    ! isotropic contribution is added AFTER the fiber loop so that diso reflects
    ! any HGO I1-I4 dispersion coupling fed back by sef_aniso_hgo.

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
        iso_fic_active = .true.  ! HGO dispersion couples into the I1 (diso) channel

      case (ANISO_HUMPHREY)
        k1_fib = props(ip)
        k2_fib = props(ip+1)
        vorif  = props(ip+2:ip+4)
        ip = ip + 5

        call fiber_direction(vorif, distgr, st0, vd)
        call pseudo_invariants(cbar, st0, det, cbari4, lambda_fib, barlambda)
        call sef_aniso_humphrey(sseaniso, daniso, diso, k1_fib, k2_fib, cbari4, cbari1)

      case (ANISO_HUMPHREY_ACT)
        k1_fib = props(ip)
        k2_fib = props(ip+1)
        vorif  = props(ip+2:ip+4)
        t0m_act = props(ip+5)
        ip = ip + 6

        call fiber_direction(vorif, distgr, st0, vd)
        call pseudo_invariants(cbar, st0, det, cbari4, lambda_fib, barlambda)
        call sef_aniso_humphrey(sseaniso, daniso, diso, k1_fib, k2_fib, cbari4, cbari1)
        call sef_aniso_humphrey_act(sseaniso, daniso, cbari4, predef, dpred, t0m_act)

      case (ANISO_HGO_AI)
        k1_fib    = props(ip)
        k2_fib    = props(ip+1)
        bdisp_fib = props(ip+2)
        factor_aniso = nint(props(ip+3))
        vorif     = props(ip+4:ip+6)
        ip = ip + 7

        call sef_aniso_hgo_ai(sfic_aniso_ai, cfic_aniso_ai, w_aniso_ai, &
                               distgr, det, k1_fib, k2_fib, bdisp_fib, factor_aniso, vorif)

        ! AI aniso returns spatial sfic/cfic. Pull back to the material frame so the
        ! fiber also enters pkfic/cmfic — needed for the viscous path (PK2/material
        ! tangent) and to mirror the non-AI fiber (D1 fix). Then skip the daniso path.
        call matinv3d(distgr, finv)
        call pull2(pkmatficaniso, sfic_aniso_ai, finv, det)
        call pull4(cmanisomatfic, cfic_aniso_ai, finv, det)
        pkfic = pkfic + pkmatficaniso
        sfic  = sfic  + sfic_aniso_ai
        cmfic = cmfic + cmanisomatfic
        cfic  = cfic  + cfic_aniso_ai
        cycle

      case (ANISO_HUMPHREY_AI)
        k1_fib    = props(ip)
        k2_fib    = props(ip+1)
        bdisp_fib = props(ip+2)
        factor_aniso = nint(props(ip+3))
        vorif     = props(ip+4:ip+6)
        ip = ip + 7

        call sef_aniso_humphrey_ai(sfic_aniso_ai, cfic_aniso_ai, w_aniso_ai, &
                                    distgr, det, k1_fib, k2_fib, bdisp_fib, factor_aniso, vorif)

        ! Pull back to the material frame (D1 fix — see ANISO_HGO_AI above).
        call matinv3d(distgr, finv)
        call pull2(pkmatficaniso, sfic_aniso_ai, finv, det)
        call pull4(cmanisomatfic, cfic_aniso_ai, finv, det)
        pkfic = pkfic + pkmatficaniso
        sfic  = sfic  + sfic_aniso_ai
        cmfic = cmfic + cmanisomatfic
        cfic  = cfic  + cfic_aniso_ai
        cycle

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

    ! --- Isochoric isotropic fictitious tensors ---
    ! Built AFTER the fiber loop so diso includes any HGO I1-I4 dispersion coupling
    ! fed back by sef_aniso_hgo, then added to the (so far anisotropic-only) totals.
    if (iso_fic_active) then
      call pk2_isomatfic(pkmatfic, diso, cbar, cbari1, unit2)
      call sig_isomatfic(sisomatfic, pkmatfic, distgr, det)
      call cmat_isomatfic(cmisomatfic, cbar, cbari1, cbari2, diso, unit2, unit4, det)
      call cs_isomatfic(cisomatfic, cmisomatfic, distgr, det)
      pkfic = pkfic + pkmatfic
      sfic  = sfic  + sisomatfic
      cmfic = cmfic + cmisomatfic
      cfic  = cfic  + cisomatfic
    end if

    ! --- Network contribution ---
    ! Compute the network now (this advances the PROPS pointer past the network
    ! params), but DEFER the matrix/network blend until AFTER damage, so damage
    ! scales only the matrix — consistent with the damage criterion, which uses
    ! the matrix energy and excludes the network (D3). Network models require
    ! quadrature data loaded via UEXTERNALDB; snet/cnet are fictitious spatial.
    if (network_type > 0) then
      call network_contribution(network_type, props, ip, distgr, det, &
                                phi_net, snet, cnet, &
                                statev, nstatev, dtime, time, kstep, noel, &
                                damage_type, n_visco)
    else
      phi_net = ZERO
    end if

    ! --- Damage modification (matrix only; applied BEFORE the network blend) ---
    if (damage_type == DMG_SIGMOID) then
      beta_d   = props(ip)
      psi_half = props(ip+1)
      ip = ip + 2

      dmg       = statev(2)
      sef0_hist = statev(3)

      call damage_sigmoid(sseiso + sseaniso, sef0_hist, dmg, dmg_red, dmg_diff, beta_d, psi_half)

      statev(2) = dmg
      statev(3) = sef0_hist

      ! Capture the UNDAMAGED matrix isochoric stresses before scaling. They drive
      ! the consistent-tangent softening term d'(W)*S0 (x) S0, added after assembly.
      call pk2_iso(pkiso0, pkfic, projl, det)
      call sig_iso(siso0, sfic, proje)

      ! Secant scaling of the MATRIX (the network, blended below, is left undamaged).
      pkfic = dmg_red * pkfic
      sfic  = dmg_red * sfic
      cmfic = dmg_red * cmfic
      cfic  = dmg_red * cfic

      ! Scale strain energy for consistency with damaged stress
      sseiso   = dmg_red * sseiso
      sseaniso = dmg_red * sseaniso
    end if

    ! --- Blend the (now damaged) matrix with the pristine network ---
    ! total = (1-PHI)*matrix + PHI*network; the network is NOT damaged here and
    ! (in the viscous branch below) is NOT relaxed.
    if (network_type > 0) then
      sfic  = (ONE - phi_net) * sfic  + phi_net * snet
      cfic  = (ONE - phi_net) * cfic  + phi_net * cnet
      pkfic = (ONE - phi_net) * pkfic   ! network has no PK2 part
    end if

    ! ============================================================================
    ! STRESS MEASURES (elastic, no viscosity yet)
    ! ============================================================================

    ! --- Volumetric ---
    call pk2_vol(pkvol, pv, c, det)
    call sig_vol(svol, pv, unit2)

    ! --- Isochoric (projected) ---
    call pk2_iso(pkiso, pkfic, projl, det)
    call sig_iso(siso, sfic, proje)

    ! --- Total elastic Cauchy stress ---
    sigma = svol + siso

    ! ============================================================================
    ! ELASTICITY TENSORS (elastic, no rate-convention term)
    ! ============================================================================

    ! --- Material elasticity ---
    call met_vol(cmvol, c, pv, ppv, det)
    call met_iso(cmiso, cmfic, projl, pkiso, pkfic, c, unit2, det)

    ! --- Spatial elasticity ---
    call set_vol(cvol, pv, ppv, unit2, unit4s)
    call set_iso(ciso, cfic, proje, siso, sfic, unit2)

    ! --- Consistent damage softening term ---
    ! sigma = sigma_vol + d_red(W)*sigma_iso0, so the tangent gains
    !   material: C_m += d'(W) * S_iso0 (x) S_iso0
    !   spatial:  c   += d'(W) * J * sigma_iso0 (x) sigma_iso0
    ! (dmg_diff = d(d_red)/dW < 0; it is ZERO unless damage is actively growing).
    if (damage_type == DMG_SIGMOID .and. dmg_diff /= ZERO) then
      do di = 1, 3
        do dj = 1, 3
          do dk = 1, 3
            do dl = 1, 3
              cmiso(di,dj,dk,dl) = cmiso(di,dj,dk,dl) &
                                 + dmg_diff * pkiso0(di,dj) * pkiso0(dk,dl)
              ciso(di,dj,dk,dl)  = ciso(di,dj,dk,dl) &
                                 + dmg_diff * det * siso0(di,dj) * siso0(dk,dl)
            end do
          end do
        end do
      end do
    end if

    cspatial = cvol + ciso

    ! ============================================================================
    ! VISCOUS CONTRIBUTION (optional)
    ! ============================================================================
    ! The viscous model operates on the material frame (PK2 + material tangent).
    ! It modifies PK2 by adding viscous overstress and scales the isochoric
    ! material tangent. We then push-forward the result to the spatial frame.
    ! Note: viscosity applies to the continuum (matrix) part only. If a network
    ! model is active, the network spatial contribution is added back afterwards.
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
      call visco_maxwell(pk2_visc, cmat_visc, &
                         min(n_visco, MAX_VISCO_BRANCHES), pkvol, pkiso, &
                         cmvol, cmiso, dtime, vscprops, statev, sdv_offset_visco)

      ! Push-forward to get spatial quantities
      call push2(sigma_visc, pk2_visc, dfgrd1, det)
      call push4(cspatial_visc, cmat_visc, dfgrd1, det)

      ! Override elastic results with viscous-augmented versions
      sigma    = sigma_visc
      cspatial = cspatial_visc

      ! Re-add the network (viscosity relaxes the matrix only; the network is
      ! pristine — NOT relaxed). Use snet/cnet directly, NOT the blended sfic/cfic,
      ! so the (1-phi)*matrix + phi*network blend is not double-counted (D2).
      if (network_type > 0) then
        call sig_iso(siso, snet, proje)
        call set_iso(ciso, cnet, proje, siso, snet, unit2)
        sigma    = sigma    + phi_net * siso
        cspatial = cspatial + phi_net * ciso
      end if
    end if

    ! ============================================================================
    ! STRAIN ENERGY
    ! ============================================================================
    sse = ssev + sseiso + sseaniso

    ! ============================================================================
    ! STATE VARIABLES
    ! ============================================================================
    statev(1) = det

  end subroutine material_core


  !> UEL-facing entry point: full Cauchy stress and the displacement-UEL
  !> assembly tangent  a_ijkl = c_ijkl + delta_ik sigma_jl  (push-forward of
  !> dP/dF). Use this from user elements that build Kuu = int G^T A G dv;
  !> a has no minor symmetry, so declare the element Unsymm.
  subroutine material_uel(sigma, atangent, sse, &
                          statev, nstatev, props, nprops, dfgrd1, &
                          time, dtime, predef, dpred, pnewdt, &
                          noel, npt, kstep, kinc)

    use mod_constants, only: dp, ZERO
    use mod_tensor,    only: onem
    use mod_continuum, only: set_geometric

    implicit none

    real(dp), intent(out)   :: sigma(3,3), atangent(3,3,3,3)
    double precision, intent(out)   :: sse
    integer,  intent(in)    :: nstatev, nprops
    double precision, intent(in out) :: statev(nstatev)
    double precision, intent(in)     :: props(nprops)
    double precision, intent(in)     :: dfgrd1(3,3)
    double precision, intent(in)     :: time(2), dtime
    double precision, intent(in)     :: predef(1), dpred(1)
    double precision, intent(in out) :: pnewdt
    integer,  intent(in)    :: noel, npt, kstep, kinc

    real(dp) :: cspatial(3,3,3,3), cgeo(3,3,3,3)
    real(dp) :: unit2(3,3), unit4(3,3,3,3), unit4s(3,3,3,3)
    integer  :: istat

    call material_core(sigma, cspatial, sse, istat, &
                       statev, nstatev, props, nprops, dfgrd1, &
                       time, dtime, predef, dpred, pnewdt, &
                       noel, npt, kstep, kinc)

    if (istat /= 0) then
      ! Inverted element: pnewdt already cut; hand back zeros so the element
      ! assembles cleanly while ABAQUS abandons the increment.
      atangent = ZERO
      return
    end if

    call onem(unit2, unit4, unit4s)
    call set_geometric(cgeo, sigma, unit2)
    atangent = cspatial + cgeo

  end subroutine material_uel

end module mod_material


! ==============================================================================
! ABAQUS *USER MATERIAL entry point (built-in elements).
! Thin wrapper: constitutive core + Jaumann rate term + Voigt indexing.
! ==============================================================================
subroutine umat(stress, statev, ddsdde, sse, spd, scd, &
     rpl, ddsddt, drplde, drpldt, &
     stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
     ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
     celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

  use mod_constants, only: dp
  use mod_tensor,    only: onem
  use mod_continuum, only: set_jr, indexx
  use mod_material,  only: material_core

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

  real(dp) :: sigma(3,3), cspatial(3,3,3,3), ddsigdde(3,3,3,3), cjr(3,3,3,3)
  real(dp) :: unit2(3,3), unit4(3,3,3,3), unit4s(3,3,3,3)
  integer  :: istat

  call material_core(sigma, cspatial, sse, istat, &
                     statev, nstatev, props, nprops, dfgrd1, &
                     time, dtime, predef, dpred, pnewdt, &
                     noel, npt, kstep, kinc)

  ! Inverted element: pnewdt already cut; leave stress/ddsdde untouched so
  ! ABAQUS abandons the increment (matches the pre-refactor behavior).
  if (istat /= 0) return

  ! Jaumann rate correction with the TOTAL Cauchy stress (ABAQUS UMAT
  ! material Jacobian convention), then collapse to Voigt.
  call onem(unit2, unit4, unit4s)
  call set_jr(cjr, sigma, unit2)
  ddsigdde = cspatial + cjr

  call indexx(stress, ddsdde, sigma, ddsigdde, ntens)

end subroutine umat


! ============================================================================
! INTERNAL: Network dispatch (reads PROPS, calls the appropriate network model)
! ============================================================================
subroutine network_contribution(network_type, props, ip, distgr, det, &
                                phi_net, snet, cnet, &
                                statev, nstatev, dtime, time, kstep, noel, &
                                damage_type, n_visco)
  use mod_constants, only: dp, ZERO, ONE
  use mod_network,   only: affine_network, nonaffine_network, &
                           affine_network_ai, nonaffine_network_ai, &
                           mixed_network_ai, &
                           contractile_network_ai, chemical_state

  implicit none
  integer,  intent(in)    :: network_type
  double precision, intent(in) :: props(*)
  integer,  intent(inout) :: ip
  real(dp), intent(in)    :: distgr(3,3), det
  real(dp), intent(out)   :: phi_net, snet(3,3), cnet(3,3,3,3)
  double precision, intent(inout) :: statev(*)
  integer,  intent(in)    :: nstatev, kstep, noel
  double precision, intent(in) :: dtime, time(2)
  integer,  intent(in)    :: damage_type, n_visco

  real(dp) :: net_density, b_orient, efi, pp_naff
  real(dp) :: density_aff, density_naff
  real(dp) :: filprops(8)
  real(dp) :: prefdir_net(3), prefdir_ai(3)
  integer  :: factor_ai

  ! Contractile-specific variables
  real(dp) :: fric_val, ffmax_val
  real(dp) :: kch(7), frac(4), frac_new(4)
  real(dp) :: ru0(2000)  ! sliding displacement workspace (must be >= nwp)
  integer  :: nwp, sdv_off

  ! Quadrature data loaded via UEXTERNALDB into COMMON blocks (for RW types)
  integer, parameter :: MAX_NWP  = 720
  integer, parameter :: MAX_ELEM = 100000
  real(dp) :: mf0(MAX_NWP, 3), rw(MAX_NWP)
  real(dp) :: prefdir_db(MAX_ELEM, 4)
  integer  :: nwp_active
  common /kfil/  mf0
  common /kfilr/ rw
  common /kfilp/ prefdir_db
  common /knwp/  nwp_active

  snet = ZERO
  cnet = ZERO

  select case (network_type)
  case (1) ! Affine (quadrature weights)
    phi_net        = props(ip)
    net_density    = props(ip+1)
    b_orient       = props(ip+2)
    efi            = props(ip+3)
    prefdir_net(1) = props(ip+4)   ! preferred direction
    prefdir_net(2) = props(ip+5)
    prefdir_net(3) = props(ip+6)
    filprops(1)    = props(ip+7)   ! L
    filprops(2)    = props(ip+8)   ! R0F
    filprops(3)    = props(ip+9)   ! mu0
    filprops(4)    = props(ip+10)  ! beta
    filprops(5)    = props(ip+11)  ! B0
    filprops(6)    = props(ip+12)  ! lambda0
    filprops(7)    = props(ip+13)  ! R0C
    filprops(8)    = props(ip+14)  ! ETAC
    ip = ip + 15

    call affine_network(snet, cnet, distgr, mf0, rw, nwp_active, det, &
                        filprops, net_density, b_orient, efi, prefdir_net)

  case (2) ! Non-affine (quadrature weights)
    phi_net        = props(ip)
    net_density    = props(ip+1)
    b_orient       = props(ip+2)
    efi            = props(ip+3)
    pp_naff        = props(ip+4)   ! non-affinity exponent
    prefdir_net(1) = props(ip+5)   ! preferred direction
    prefdir_net(2) = props(ip+6)
    prefdir_net(3) = props(ip+7)
    filprops(1)    = props(ip+8)
    filprops(2)    = props(ip+9)
    filprops(3)    = props(ip+10)
    filprops(4)    = props(ip+11)
    filprops(5)    = props(ip+12)
    filprops(6)    = props(ip+13)
    filprops(7)    = props(ip+14)
    filprops(8)    = props(ip+15)
    ip = ip + 16

    call nonaffine_network(snet, cnet, distgr, mf0, rw, nwp_active, det, &
                           filprops, net_density, b_orient, efi, pp_naff, prefdir_net)

  case (5) ! Affine — angular integration (icosahedron)
    phi_net        = props(ip)
    net_density    = props(ip+1)
    b_orient       = props(ip+2)
    efi            = props(ip+3)
    factor_ai      = nint(props(ip+4))  ! refinement factor
    prefdir_ai(1)  = props(ip+5)        ! preferred direction
    prefdir_ai(2)  = props(ip+6)
    prefdir_ai(3)  = props(ip+7)
    filprops(1)    = props(ip+8)   ! L
    filprops(2)    = props(ip+9)   ! R0F
    filprops(3)    = props(ip+10)  ! mu0
    filprops(4)    = props(ip+11)  ! beta
    filprops(5)    = props(ip+12)  ! B0
    filprops(6)    = props(ip+13)  ! lambda0
    filprops(7)    = props(ip+14)  ! R0C
    filprops(8)    = props(ip+15)  ! ETAC
    ip = ip + 16

    call affine_network_ai(snet, cnet, distgr, det, &
                           filprops, net_density, b_orient, efi, factor_ai, prefdir_ai)

  case (6) ! Non-affine — angular integration (icosahedron)
    phi_net        = props(ip)
    net_density    = props(ip+1)
    pp_naff        = props(ip+2)   ! non-affinity exponent
    factor_ai      = nint(props(ip+3))  ! refinement factor
    filprops(1)    = props(ip+4)   ! L
    filprops(2)    = props(ip+5)   ! R0F
    filprops(3)    = props(ip+6)   ! mu0
    filprops(4)    = props(ip+7)   ! beta
    filprops(5)    = props(ip+8)   ! B0
    filprops(6)    = props(ip+9)   ! lambda0
    filprops(7)    = props(ip+10)  ! R0C
    filprops(8)    = props(ip+11)  ! ETAC
    ip = ip + 12

    call nonaffine_network_ai(snet, cnet, distgr, det, &
                              filprops, net_density, pp_naff, factor_ai)

  case (3) ! Mixed — affine + non-affine (AI)
    phi_net        = props(ip)
    density_naff   = props(ip+1)
    pp_naff        = props(ip+2)
    density_aff    = props(ip+3)
    b_orient       = props(ip+4)
    efi            = props(ip+5)
    factor_ai      = nint(props(ip+6))
    prefdir_ai(1)  = props(ip+7)
    prefdir_ai(2)  = props(ip+8)
    prefdir_ai(3)  = props(ip+9)
    filprops(1)    = props(ip+10)  ! L
    filprops(2)    = props(ip+11)  ! R0F
    filprops(3)    = props(ip+12)  ! mu0
    filprops(4)    = props(ip+13)  ! beta
    filprops(5)    = props(ip+14)  ! B0
    filprops(6)    = props(ip+15)  ! lambda0
    filprops(7)    = props(ip+16)  ! R0C
    filprops(8)    = props(ip+17)  ! ETAC
    ip = ip + 18

    call mixed_network_ai(snet, cnet, distgr, det, &
                          filprops, density_naff, pp_naff, density_aff, &
                          b_orient, efi, factor_ai, prefdir_ai)

  case (4) ! Contractile (AI)
    phi_net        = props(ip)
    net_density    = props(ip+1)
    b_orient       = props(ip+2)
    efi            = props(ip+3)
    fric_val       = props(ip+4)
    ffmax_val      = props(ip+5)
    factor_ai      = nint(props(ip+6))
    prefdir_ai(1)  = props(ip+7)
    prefdir_ai(2)  = props(ip+8)
    prefdir_ai(3)  = props(ip+9)
    filprops(1:8)  = props(ip+10:ip+17)
    kch(1:7)       = props(ip+18:ip+24)
    ip = ip + 25

    ! State variable offset for contractile data
    ! After: det(1) + damage(2 if active) + visco(9*N_VISCO)
    sdv_off = 1
    if (damage_type > 0) sdv_off = sdv_off + 2
    sdv_off = sdv_off + 9 * n_visco

    ! Update chemical state (forward Euler)
    frac(1:4) = statev(sdv_off+1 : sdv_off+4)
    ! Initialize chemical state on first increment
    if (time(1) == ZERO .and. kstep == 1) then
      frac(1) = ONE
      frac(2:4) = ZERO
    end if
    call chemical_state(frac_new, frac, kch, dtime)
    statev(sdv_off+1 : sdv_off+4) = frac_new

    ! Derive nwp from icosahedron refinement factor (20 * factor^2)
    nwp = 20 * factor_ai * factor_ai
    if (nwp > 2000) nwp = 2000  ! workspace limit
    if (sdv_off + 4 + nwp > nstatev) then
      nwp = nstatev - sdv_off - 4
    end if
    if (nwp > 0) then
      ru0(1:nwp) = statev(sdv_off+5 : sdv_off+4+nwp)
    end if

    call contractile_network_ai(snet, cnet, distgr, det, &
                                filprops, net_density, b_orient, efi, &
                                fric_val, ffmax_val, factor_ai, prefdir_ai, &
                                kch, dtime, ru0, frac_new, nwp)

    ! Write back updated RU0 state variables
    if (nwp > 0) then
      statev(sdv_off+5 : sdv_off+4+nwp) = ru0(1:nwp)
    end if

  case default
    phi_net = ZERO
  end select
end subroutine network_contribution
