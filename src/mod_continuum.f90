!> @brief Stress measures, elasticity tensors, and Voigt indexing.
!>
!> Implements volumetric/isochoric stress split, material/spatial elasticity tensors,
!> the Jaumann rate correction, and 3D-to-Voigt index conversion for ABAQUS.
module mod_continuum
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, FOUR, HALF
  use mod_tensor, only: matinv3d, contraction42, contraction44, contraction22, &
                         push2, push4
  implicit none
  private

  ! Volumetric
  public :: vol_energy, pk2_vol, sig_vol, met_vol, set_vol

  ! Isochoric stress & elasticity
  public :: pk2_isomatfic, sig_isomatfic, cmat_isomatfic, cs_isomatfic
  public :: pk2_iso, sig_iso, met_iso, set_iso

  ! Anisotropic stress & elasticity
  public :: pk2_anisofic, cmat_anisofic

  ! Jaumann rate
  public :: set_jr

  ! Voigt indexing
  public :: indexx

contains

  ! ============================================================================
  ! VOLUMETRIC CONTRIBUTION
  ! ============================================================================

  !> Volumetric strain energy and its 1st/2nd derivatives.
  !> G = (1/4)(J^2 - 1 - 2*ln(J)), W_vol = K*G
  subroutine vol_energy(kbulk, det, ssev, pv, ppv)
    real(dp), intent(in)  :: kbulk, det
    real(dp), intent(out) :: ssev, pv, ppv
    real(dp) :: g, aux

    g    = 0.25_dp * (det*det - ONE - TWO*log(det))
    ssev = kbulk * g
    pv   = kbulk * HALF * (det - ONE/det)
    aux  = kbulk * HALF * (ONE + ONE/(det*det))
    ppv  = pv + det * aux
  end subroutine vol_energy

  !> Volumetric 2nd Piola-Kirchhoff stress: S_vol = p * J * C^{-1}
  subroutine pk2_vol(pkvol, pv, c)
    real(dp), intent(out) :: pkvol(3,3)
    real(dp), intent(in)  :: pv, c(3,3)
    real(dp) :: cinv(3,3)

    call matinv3d(c, cinv)
    pkvol = pv * cinv
  end subroutine pk2_vol

  !> Volumetric Cauchy stress: sigma_vol = p * I
  subroutine sig_vol(svol, pv, unit2)
    real(dp), intent(out) :: svol(3,3)
    real(dp), intent(in)  :: pv, unit2(3,3)

    svol = pv * unit2
  end subroutine sig_vol

  !> Volumetric material elasticity tensor.
  subroutine met_vol(cmvol, c, pv, ppv, det)
    real(dp), intent(out) :: cmvol(3,3,3,3)
    real(dp), intent(in)  :: c(3,3), pv, ppv, det
    real(dp) :: cinv(3,3)
    integer :: i, j, k, l

    call matinv3d(c, cinv)
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cmvol(i,j,k,l) = det*ppv * cinv(i,j)*cinv(k,l) &
                            - det*pv * (cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k))
          end do
        end do
      end do
    end do
  end subroutine met_vol

  !> Volumetric spatial elasticity tensor.
  subroutine set_vol(cvol, pv, ppv, unit2, unit4s)
    real(dp), intent(out) :: cvol(3,3,3,3)
    real(dp), intent(in)  :: pv, ppv, unit2(3,3), unit4s(3,3,3,3)
    integer :: i, j, k, l

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cvol(i,j,k,l) = ppv * unit2(i,j)*unit2(k,l) &
                           - TWO*pv * unit4s(i,j,k,l)
          end do
        end do
      end do
    end do
  end subroutine set_vol

  ! ============================================================================
  ! ISOCHORIC ISOTROPIC CONTRIBUTION
  ! ============================================================================

  !> Isotropic "fictitious" 2PK stress from invariant derivatives.
  !> Sfic = 2*(dW/dI1 + I1*dW/dI2)*I - 2*dW/dI2*Cbar
  subroutine pk2_isomatfic(fic, diso, cbar, cbari1, unit2)
    real(dp), intent(out) :: fic(3,3)
    real(dp), intent(in)  :: diso(5), cbar(3,3), unit2(3,3), cbari1
    real(dp) :: aux1, aux2
    integer  :: i, j

    aux1 = TWO * (diso(1) + cbari1*diso(2))
    aux2 = -TWO * diso(2)
    do i = 1, 3
      do j = 1, 3
        fic(i,j) = aux1*unit2(i,j) + aux2*cbar(i,j)
      end do
    end do
  end subroutine pk2_isomatfic

  !> Push-forward of fictitious PK2 to get fictitious Cauchy stress.
  subroutine sig_isomatfic(sfic, pkfic, distgr, det)
    real(dp), intent(out) :: sfic(3,3)
    real(dp), intent(in)  :: pkfic(3,3), distgr(3,3), det
    call push2(sfic, pkfic, distgr, det)
  end subroutine sig_isomatfic

  !> Isotropic "fictitious" material elasticity tensor (invariant-based).
  subroutine cmat_isomatfic(cmfic, cbar, cbari1, cbari2, diso, unit2, unit4, det)
    real(dp), intent(out) :: cmfic(3,3,3,3)
    real(dp), intent(in)  :: cbar(3,3), diso(5), unit2(3,3), unit4(3,3,3,3)
    real(dp), intent(in)  :: cbari1, cbari2, det
    real(dp) :: aux1, aux2, aux3, aux4
    integer  :: i, j, k, l

    aux1 = FOUR*(diso(3) + TWO*cbari1*diso(5) + diso(2) + cbari1*cbari1*diso(4))
    aux2 = -FOUR*(diso(5) + cbari1*diso(4))
    aux3 = FOUR*diso(4)
    aux4 = -FOUR*diso(2)

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cmfic(i,j,k,l) = aux1*unit2(i,j)*unit2(k,l) &
                            + aux2*(unit2(i,j)*cbar(k,l) + cbar(i,j)*unit2(k,l)) &
                            + aux3*cbar(i,j)*cbar(k,l) &
                            + aux4*unit4(i,j,k,l)
          end do
        end do
      end do
    end do
  end subroutine cmat_isomatfic

  !> Push-forward of fictitious material elasticity to spatial frame.
  subroutine cs_isomatfic(csfic, cmfic, distgr, det)
    real(dp), intent(out) :: csfic(3,3,3,3)
    real(dp), intent(in)  :: cmfic(3,3,3,3), distgr(3,3), det
    call push4(csfic, cmfic, distgr, det)
  end subroutine cs_isomatfic

  !> Isochoric PK2 stress: S_iso = J^{-2/3} * PL : S_fic
  subroutine pk2_iso(pkiso, pkfic, pl, det)
    real(dp), intent(out) :: pkiso(3,3)
    real(dp), intent(in)  :: pkfic(3,3), pl(3,3,3,3), det
    real(dp) :: scale2
    integer  :: i, j

    call contraction42(pkiso, pl, pkfic)
    scale2 = det**(-TWO/THREE)
    do i = 1, 3
      do j = 1, 3
        pkiso(i,j) = scale2 * pkiso(i,j)
      end do
    end do
  end subroutine pk2_iso

  !> Isochoric Cauchy stress: sigma_iso = PE : sigma_fic
  subroutine sig_iso(siso, sfic, pe)
    real(dp), intent(out) :: siso(3,3)
    real(dp), intent(in)  :: sfic(3,3), pe(3,3,3,3)
    call contraction42(siso, pe, sfic)
  end subroutine sig_iso

  !> Isochoric material elasticity tensor.
  subroutine met_iso(cmiso, cmfic, pl, pkiso, pkfic, c, unit2, det)
    real(dp), intent(out) :: cmiso(3,3,3,3)
    real(dp), intent(in)  :: cmfic(3,3,3,3), pl(3,3,3,3)
    real(dp), intent(in)  :: pkiso(3,3), pkfic(3,3), c(3,3), unit2(3,3), det
    real(dp) :: cinv(3,3), cisoaux(3,3,3,3), cisoaux1(3,3,3,3)
    real(dp) :: plt(3,3,3,3), pll(3,3,3,3)
    real(dp) :: trfic, xx, yy, zz, scl, scl2
    integer  :: i, j, k, l
    real(dp) :: pkfic_scaled(3,3)

    call matinv3d(c, cinv)

    ! PL : Cfic
    call contraction44(cisoaux1, pl, cmfic)

    ! Transpose of PL
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            plt(i,j,k,l) = pl(k,l,i,j)
          end do
        end do
      end do
    end do

    ! (PL : Cfic) : PL^T
    call contraction44(cisoaux, cisoaux1, plt)

    ! Trace of fictitious stress with C
    scl  = det**(-TWO/THREE)
    scl2 = scl * scl
    pkfic_scaled = scl * pkfic
    call contraction22(trfic, pkfic_scaled, c)

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            xx = scl2 * cisoaux(i,j,k,l)
            pll(i,j,k,l) = HALF*(cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k)) &
                          - (ONE/THREE)*cinv(i,j)*cinv(k,l)
            yy = trfic * pll(i,j,k,l)
            zz = pkiso(i,j)*cinv(k,l) + cinv(i,j)*pkiso(k,l)
            cmiso(i,j,k,l) = xx + (TWO/THREE)*yy - (TWO/THREE)*zz
          end do
        end do
      end do
    end do
  end subroutine met_iso

  !> Isochoric spatial elasticity tensor.
  subroutine set_iso(ciso, cfic, pe, siso, sfic, unit2)
    real(dp), intent(out) :: ciso(3,3,3,3)
    real(dp), intent(in)  :: cfic(3,3,3,3), pe(3,3,3,3)
    real(dp), intent(in)  :: siso(3,3), sfic(3,3), unit2(3,3)
    real(dp) :: cisoaux(3,3,3,3), cisoaux1(3,3,3,3), trfic, xx, yy, zz
    integer :: i, j, k, l

    cisoaux1 = ZERO
    cisoaux  = ZERO

    call contraction44(cisoaux1, pe, cfic)
    call contraction44(cisoaux, cisoaux1, pe)

    call contraction22(trfic, sfic, unit2)

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            xx = cisoaux(i,j,k,l)
            yy = trfic * (pe(i,j,k,l))
            zz = siso(i,j)*unit2(k,l) + unit2(i,j)*siso(k,l)
            ciso(i,j,k,l) = xx + (TWO/THREE)*yy - (TWO/THREE)*zz
          end do
        end do
      end do
    end do
  end subroutine set_iso

  ! ============================================================================
  ! ANISOTROPIC CONTRIBUTION (FICTITIOUS FRAME)
  ! ============================================================================

  !> Anisotropic "fictitious" 2PK stress: Afic = 2 * dW/dI4 * M0
  subroutine pk2_anisofic(afic, daniso, st0)
    real(dp), intent(out) :: afic(3,3)
    real(dp), intent(in)  :: daniso(4), st0(3,3)

    afic = TWO * daniso(1) * st0
  end subroutine pk2_anisofic

  !> Anisotropic "fictitious" material elasticity tensor.
  subroutine cmat_anisofic(cmaniso, st0, daniso, unit2, det)
    real(dp), intent(out) :: cmaniso(3,3,3,3)
    real(dp), intent(in)  :: st0(3,3), daniso(4), unit2(3,3), det
    real(dp) :: d2udi4, d2udi1di4
    integer :: i, j, k, l

    d2udi4    = daniso(2)
    d2udi1di4 = daniso(3)

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cmaniso(i,j,k,l) = FOUR * ( &
              d2udi4 * st0(i,j)*st0(k,l) &
            + d2udi1di4 * (unit2(i,j)*st0(k,l) + st0(i,j)*unit2(k,l)) )
          end do
        end do
      end do
    end do
  end subroutine cmat_anisofic

  ! ============================================================================
  ! JAUMANN RATE CORRECTION
  ! ============================================================================

  !> Jaumann rate contribution for the spatial elasticity tensor.
  subroutine set_jr(cjr, sigma, unit2)
    real(dp), intent(out) :: cjr(3,3,3,3)
    real(dp), intent(in)  :: sigma(3,3), unit2(3,3)
    integer :: i, j, k, l

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            cjr(i,j,k,l) = HALF * ( &
              unit2(i,k)*sigma(j,l) + sigma(i,k)*unit2(j,l) &
            + unit2(i,l)*sigma(j,k) + sigma(i,l)*unit2(j,k) )
          end do
        end do
      end do
    end do
  end subroutine set_jr

  ! ============================================================================
  ! VOIGT INDEXING (ABAQUS INTERFACE)
  ! ============================================================================

  !> Convert 3x3 stress and 3x3x3x3 stiffness to ABAQUS Voigt notation.
  subroutine indexx(stress, ddsdde, sigma, ddsigdde, ntens)
    integer,  intent(in)  :: ntens
    real(dp), intent(out) :: stress(ntens), ddsdde(ntens, ntens)
    real(dp), intent(in)  :: sigma(3,3), ddsigdde(3,3,3,3)

    ! Stress: [11, 22, 33, 12, 13, 23]
    stress(1) = sigma(1,1)
    stress(2) = sigma(2,2)
    stress(3) = sigma(3,3)
    if (ntens > 3) then
      stress(4) = sigma(1,2)
    end if
    if (ntens > 4) then
      stress(5) = sigma(1,3)
      stress(6) = sigma(2,3)
    end if

    ! Stiffness with full symmetry
    ddsdde(1,1) = ddsigdde(1,1,1,1)
    ddsdde(1,2) = ddsigdde(1,1,2,2)
    ddsdde(1,3) = ddsigdde(1,1,3,3)
    ddsdde(2,1) = ddsigdde(2,2,1,1)
    ddsdde(2,2) = ddsigdde(2,2,2,2)
    ddsdde(2,3) = ddsigdde(2,2,3,3)
    ddsdde(3,1) = ddsigdde(3,3,1,1)
    ddsdde(3,2) = ddsigdde(3,3,2,2)
    ddsdde(3,3) = ddsigdde(3,3,3,3)

    if (ntens > 3) then
      ddsdde(1,4) = ddsigdde(1,1,1,2)
      ddsdde(2,4) = ddsigdde(2,2,1,2)
      ddsdde(3,4) = ddsigdde(3,3,1,2)
      ddsdde(4,1) = ddsigdde(1,2,1,1)
      ddsdde(4,2) = ddsigdde(1,2,2,2)
      ddsdde(4,3) = ddsigdde(1,2,3,3)
      ddsdde(4,4) = ddsigdde(1,2,1,2)
    end if

    if (ntens > 4) then
      ddsdde(1,5) = ddsigdde(1,1,1,3)
      ddsdde(2,5) = ddsigdde(2,2,1,3)
      ddsdde(3,5) = ddsigdde(3,3,1,3)
      ddsdde(4,5) = ddsigdde(1,2,1,3)
      ddsdde(1,6) = ddsigdde(1,1,2,3)
      ddsdde(2,6) = ddsigdde(2,2,2,3)
      ddsdde(3,6) = ddsigdde(3,3,2,3)
      ddsdde(4,6) = ddsigdde(1,2,2,3)

      ddsdde(5,1) = ddsigdde(1,3,1,1)
      ddsdde(5,2) = ddsigdde(1,3,2,2)
      ddsdde(5,3) = ddsigdde(1,3,3,3)
      ddsdde(5,4) = ddsigdde(1,3,1,2)
      ddsdde(5,5) = ddsigdde(1,3,1,3)
      ddsdde(5,6) = ddsigdde(1,3,2,3)

      ddsdde(6,1) = ddsigdde(2,3,1,1)
      ddsdde(6,2) = ddsigdde(2,3,2,2)
      ddsdde(6,3) = ddsigdde(2,3,3,3)
      ddsdde(6,4) = ddsigdde(2,3,1,2)
      ddsdde(6,5) = ddsigdde(2,3,1,3)
      ddsdde(6,6) = ddsigdde(2,3,2,3)
    end if
  end subroutine indexx

end module mod_continuum
