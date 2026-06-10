!> @brief 8-node brick user element (U3D8) for large-deformation mechanics.
!>
!> Locking is avoided for the fully-integrated element with the F-bar method:
!>   de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
!>   Design of simple low order finite elements for large strain analysis of
!>   nearly incompressible solids. Int. J. Solids Structures, 33, 3277-3296.
!>
!> The constitutive response is provided by material_uel (mod_material).
module uel_element
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, HALF, THIRD
  use uel_config, only: numElem, globalSdv, use_fbar
  use uel_shape, only: calcShape3DLinear, mapShape3D, xint3D1pt, xint3D8pt, &
                       mdet, identity3
  use mod_material, only: material_uel
  implicit none
  private
  public :: u3d8, AssembleElement

contains

  subroutine u3d8(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
                  props, nprops, coords, mcrd, nNode, Uall, DUall, Vel, Accn, &
                  jtype, time, dtime, kstep, kinc, jelem, params, ndload, &
                  jdltyp, adlmag, predef, npredf, lflags, mlvarx, ddlmag, &
                  mdload, pnewdt, jprops, njprop, period, nDim, nIntt)

    ! Variables passed in
    integer, intent(in)  :: ndofel, nrhs, nsvars, nprops, mcrd, nNode, jtype
    integer, intent(in)  :: kstep, kinc, jelem, ndload, npredf, mlvarx
    integer, intent(in)  :: mdload, njprop, nDim, nIntt
    integer, intent(in)  :: jdltyp(mdload,1), lflags(4), jprops(njprop)

    ! Variables defined here, passed back to ABAQUS
    real(dp), intent(inout) :: rhs(mlvarx,1), amatrx(ndofel,ndofel)
    real(dp), intent(inout) :: svars(nsvars), energy(8), pnewdt

    real(dp), intent(in) :: props(nprops), coords(mcrd,nNode), Uall(ndofel)
    real(dp), intent(in) :: DUall(mlvarx,1), Vel(ndofel), Accn(ndofel)
    real(dp), intent(in) :: time(2), dtime, params(1), adlmag(mdload,1)
    real(dp), intent(in) :: predef(2,npredf,nNode), ddlmag(mdload,1), period

    ! Locals
    integer :: i, j, k, kk, jj, intpt, nInttPt, stat, nlSdv, ngSdv
    real(dp), allocatable :: statev(:)
    real(dp) :: u(nNode,3), du(nNode,ndofel), uOld(nNode,ndofel)
    real(dp) :: coordsC(mcrd,nNode)
    real(dp) :: Iden(3,3), Le, Ru(3*nNode,1), body(3)
    real(dp) :: Kuu(3*nNode,3*nNode), sh0(nNode), detMapJ0
    real(dp) :: dshxi(nNode,3), dsh0(nNode,3), dshC0(nNode,3), detMapJ0C
    real(dp) :: Fc_tau(3,3), Fc_t(3,3), detFc_tau, detFc_t
    real(dp) :: xi(nIntt,3), xi0(nIntt,3), w(nIntt)
    real(dp) :: sh(nNode), detMapJ, dsh(nNode,3), detMapJC, dshC(nNode,3)
    real(dp) :: F_tau(3,3), F_t(3,3), detF_tau, detF_t
    real(dp) :: T_tau(3,3), SpTanMod(3,3,3,3), sse_ip
    real(dp) :: predloc(1), dpredloc(1)
    real(dp) :: Smat(6,1), Bmat(6,3*nNode), BodyForceRes(3*nNode,1)
    real(dp) :: Gmat(9,3*nNode), G0mat(9,3*nNode), Amat(9,9), Qmat(9,9)

    ! Get element parameters
    nlSdv = jprops(1) ! number of local sdv's per integ point
    ngSdv = jprops(2) ! number of global sdv's per integ point
    allocate(statev(nlSdv))
    xi0 = ZERO

    ! Allocate memory for the globalSdv's used for viewing results
    ! on the dummy mesh
    if (.not. allocated(globalSdv)) then
      allocate(globalSdv(numElem,nIntt,ngSdv), stat=stat)
      if (stat /= 0) then
        write(*,*) 'U3D8: error when allocating globalSdv, stat=', stat
        write(*,*) '  numElem=', numElem, ' nIntt=', nIntt, ' ngSdv=', ngSdv
        call exit(1)
      end if
      write(*,*) 'U3D8: globalSdv allocated with numElem=', numElem, &
                 ' nIntt=', nIntt, ' ngSdv=', ngSdv
    end if

    ! Identity tensor
    call identity3(Iden)

    ! Initialize the residual and tangent matrices to zero
    Ru = ZERO
    Kuu = ZERO

    ! Body forces
    body(1:3) = ZERO

    ! Obtain nodal displacements
    k = 0
    do i = 1, nNode
      do j = 1, nDim
        k = k + 1
        u(i,j) = Uall(k)
        du(i,j) = DUall(k,1)
        uOld(i,j) = u(i,j) - du(i,j)
      end do
    end do

    ! Obtain current nodal coordinates
    do i = 1, nNode
      do j = 1, nDim
        coordsC(j,i) = coords(j,i) + u(i,j)
      end do
    end do

    ! Impose any time-stepping changes on the increments of displacement;
    ! displacement increment check based on element diagonal
    Le = sqrt(((coordsC(1,1) - coordsC(1,7))**TWO) &
            + ((coordsC(2,1) - coordsC(2,7))**TWO) &
            + ((coordsC(3,1) - coordsC(3,7))**TWO))
    !
    do i = 1, nNode
      do j = 1, nDim
        if (abs(du(i,j)) > 10.0_dp*Le) then
          pnewdt = HALF
          return
        end if
      end do
    end do

    !--------------------------------------------------------------------
    ! Take this opportunity to perform calculations at the element
    ! centroid. Get the deformation gradient for use in the F-bar method.
    !
    ! Obtain shape functions and their local gradients at the element
    ! centroid, that means xi=eta=zeta=0.0, and nInttPt=1
    if (nNode == 8) then
      call calcShape3DLinear(1, xi0, 1, sh0, dshxi)
    else
      write(*,*) 'Incorrect number of nodes: nNode.ne.8'
      call exit(1)
    end if

    ! Map shape functions from local to global reference coordinate system
    call mapShape3D(nNode, dshxi, coords, dsh0, detMapJ0, stat)
    if (stat == 0) then
      pnewdt = HALF
      return
    end if

    ! Map shape functions from local to global current coordinate system
    call mapShape3D(nNode, dshxi, coordsC, dshC0, detMapJ0C, stat)
    if (stat == 0) then
      pnewdt = HALF
      return
    end if

    ! Calculate the deformation gradient at the element centroid at the
    ! beginning and end of the increment for use in the F-bar method
    Fc_tau = Iden
    Fc_t = Iden
    do i = 1, nDim
      do j = 1, nDim
        do k = 1, nNode
          ! F at the end of the increment
          Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
          ! F at the beginning of the increment
          Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
        end do
      end do
    end do
    !
    call mdet(Fc_tau, detFc_tau)
    call mdet(Fc_t, detFc_t)
    !
    ! With the deformation gradient known at the element centroid
    ! we are now able to implement the F-bar method later
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Begin the loop over integration points
    !
    ! Obtain integration point local coordinates and weights
    if (nIntt == 1) then
      call xint3D1pt(xi, w, nInttPt) ! 1-pt integration
    else if (nIntt == 8) then
      call xint3D8pt(xi, w, nInttPt) ! 8-pt integration
    else
      write(*,*) 'Invalid number of int points, nIntt=', nIntt
      call exit(1)
    end if

    ! Loop over integration points
    jj = 0 ! jj is used for tracking the state variables
    do intpt = 1, nInttPt

      ! Obtain state variables from previous increment
      if ((kinc <= 1) .and. (kstep == 1)) then
        ! this is the first increment of the first step:
        ! give initial conditions
        statev = ZERO
      else
        ! this is not the first increment, read old values
        statev = svars(1+jj:nlSdv+jj)
      end if

      ! Obtain shape functions and their local gradients
      if (nNode == 8) then
        call calcShape3DLinear(nInttPt, xi, intpt, sh, dshxi)
      else
        write(*,*) 'Incorrect number of nodes: nNode.ne.8'
        call exit(1)
      end if

      ! Map shape functions from local to global reference coordinate system
      call mapShape3D(nNode, dshxi, coords, dsh, detMapJ, stat)
      if (stat == 0) then
        pnewdt = HALF
        return
      end if

      ! Map shape functions from local to global current coordinate system
      call mapShape3D(nNode, dshxi, coordsC, dshC, detMapJC, stat)
      if (stat == 0) then
        pnewdt = HALF
        return
      end if

      ! Obtain the deformation gradient at this integration point
      F_tau = Iden
      F_t = Iden
      do i = 1, nDim
        do j = 1, nDim
          do k = 1, nNode
            F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
            F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
          end do
        end do
      end do
      !
      ! Modify the deformation gradient for the F-bar method, only when
      ! using the 8-node fully-integrated linear element
      if (use_fbar .and. nNode == 8 .and. nIntt == 8) then
        call mdet(F_tau, detF_tau)
        call mdet(F_t, detF_t)
        F_tau = ((detFc_tau/detF_tau)**THIRD)*F_tau
        F_t = ((detFc_tau/detF_tau)**THIRD)*F_t
      end if

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !
      ! Constitutive update at this integration point: Cauchy stress,
      ! UEL assembly tangent, and state variables. pnewdt is governed by
      ! material_uel (element-inversion cutback).
      !
      ! Field-variable interpolation to the integration point is not
      ! implemented; the first nodal value is forwarded as-is.
      predloc(1) = predef(1,1,1)
      dpredloc(1) = ZERO
      !
      call material_uel(T_tau, SpTanMod, sse_ip, &
                        statev, nlSdv, props, nprops, F_tau, &
                        time, dtime, predloc, dpredloc, pnewdt, &
                        jelem, intpt, kstep, kinc)
      !
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      ! Save the state variables at this integ point at the end
      ! of the increment
      svars(1+jj:nlSdv+jj) = statev
      jj = jj + nlSdv ! setup for the next intpt

      ! Save the state variables at this integ point in the global
      ! array used for plotting field output
      globalSdv(jelem,intpt,1:ngSdv) = statev(1:ngSdv)

      ! Compute/update the displacement residual vector
      Smat(1,1) = T_tau(1,1)
      Smat(2,1) = T_tau(2,2)
      Smat(3,1) = T_tau(3,3)
      Smat(4,1) = T_tau(1,2)
      Smat(5,1) = T_tau(2,3)
      Smat(6,1) = T_tau(1,3)
      !
      Bmat = ZERO
      do kk = 1, nNode
        Bmat(1,1+nDim*(kk-1)) = dshC(kk,1)
        Bmat(2,2+nDim*(kk-1)) = dshC(kk,2)
        Bmat(3,3+nDim*(kk-1)) = dshC(kk,3)
        Bmat(4,1+nDim*(kk-1)) = dshC(kk,2)
        Bmat(4,2+nDim*(kk-1)) = dshC(kk,1)
        Bmat(5,2+nDim*(kk-1)) = dshC(kk,3)
        Bmat(5,3+nDim*(kk-1)) = dshC(kk,2)
        Bmat(6,1+nDim*(kk-1)) = dshC(kk,3)
        Bmat(6,3+nDim*(kk-1)) = dshC(kk,1)
      end do
      !
      BodyForceRes = ZERO
      do kk = 1, nNode
        BodyForceRes(1+nDim*(kk-1),1) = sh(kk)*body(1)
        BodyForceRes(2+nDim*(kk-1),1) = sh(kk)*body(2)
        BodyForceRes(3+nDim*(kk-1),1) = sh(kk)*body(3)
      end do
      !
      Ru = Ru + detMapJC*w(intpt)* &
           ( &
           -matmul(transpose(Bmat),Smat) &
           + BodyForceRes &
           )

      ! Compute/update the displacement tangent matrix
      Gmat = ZERO
      do kk = 1, nNode
        Gmat(1,1+nDim*(kk-1)) = dshC(kk,1)
        Gmat(2,2+nDim*(kk-1)) = dshC(kk,1)
        Gmat(3,3+nDim*(kk-1)) = dshC(kk,1)
        Gmat(4,1+nDim*(kk-1)) = dshC(kk,2)
        Gmat(5,2+nDim*(kk-1)) = dshC(kk,2)
        Gmat(6,3+nDim*(kk-1)) = dshC(kk,2)
        Gmat(7,1+nDim*(kk-1)) = dshC(kk,3)
        Gmat(8,2+nDim*(kk-1)) = dshC(kk,3)
        Gmat(9,3+nDim*(kk-1)) = dshC(kk,3)
      end do

      G0mat = ZERO
      do kk = 1, nNode
        G0mat(1,1+nDim*(kk-1)) = dshC0(kk,1)
        G0mat(2,2+nDim*(kk-1)) = dshC0(kk,1)
        G0mat(3,3+nDim*(kk-1)) = dshC0(kk,1)
        G0mat(4,1+nDim*(kk-1)) = dshC0(kk,2)
        G0mat(5,2+nDim*(kk-1)) = dshC0(kk,2)
        G0mat(6,3+nDim*(kk-1)) = dshC0(kk,2)
        G0mat(7,1+nDim*(kk-1)) = dshC0(kk,3)
        G0mat(8,2+nDim*(kk-1)) = dshC0(kk,3)
        G0mat(9,3+nDim*(kk-1)) = dshC0(kk,3)
      end do

      Amat = ZERO
      Amat(1,1) = SpTanMod(1,1,1,1)
      Amat(1,2) = SpTanMod(1,1,2,1)
      Amat(1,3) = SpTanMod(1,1,3,1)
      Amat(1,4) = SpTanMod(1,1,1,2)
      Amat(1,5) = SpTanMod(1,1,2,2)
      Amat(1,6) = SpTanMod(1,1,3,2)
      Amat(1,7) = SpTanMod(1,1,1,3)
      Amat(1,8) = SpTanMod(1,1,2,3)
      Amat(1,9) = SpTanMod(1,1,3,3)
      Amat(2,1) = SpTanMod(2,1,1,1)
      Amat(2,2) = SpTanMod(2,1,2,1)
      Amat(2,3) = SpTanMod(2,1,3,1)
      Amat(2,4) = SpTanMod(2,1,1,2)
      Amat(2,5) = SpTanMod(2,1,2,2)
      Amat(2,6) = SpTanMod(2,1,3,2)
      Amat(2,7) = SpTanMod(2,1,1,3)
      Amat(2,8) = SpTanMod(2,1,2,3)
      Amat(2,9) = SpTanMod(2,1,3,3)
      Amat(3,1) = SpTanMod(3,1,1,1)
      Amat(3,2) = SpTanMod(3,1,2,1)
      Amat(3,3) = SpTanMod(3,1,3,1)
      Amat(3,4) = SpTanMod(3,1,1,2)
      Amat(3,5) = SpTanMod(3,1,2,2)
      Amat(3,6) = SpTanMod(3,1,3,2)
      Amat(3,7) = SpTanMod(3,1,1,3)
      Amat(3,8) = SpTanMod(3,1,2,3)
      Amat(3,9) = SpTanMod(3,1,3,3)
      Amat(4,1) = SpTanMod(1,2,1,1)
      Amat(4,2) = SpTanMod(1,2,2,1)
      Amat(4,3) = SpTanMod(1,2,3,1)
      Amat(4,4) = SpTanMod(1,2,1,2)
      Amat(4,5) = SpTanMod(1,2,2,2)
      Amat(4,6) = SpTanMod(1,2,3,2)
      Amat(4,7) = SpTanMod(1,2,1,3)
      Amat(4,8) = SpTanMod(1,2,2,3)
      Amat(4,9) = SpTanMod(1,2,3,3)
      Amat(5,1) = SpTanMod(2,2,1,1)
      Amat(5,2) = SpTanMod(2,2,2,1)
      Amat(5,3) = SpTanMod(2,2,3,1)
      Amat(5,4) = SpTanMod(2,2,1,2)
      Amat(5,5) = SpTanMod(2,2,2,2)
      Amat(5,6) = SpTanMod(2,2,3,2)
      Amat(5,7) = SpTanMod(2,2,1,3)
      Amat(5,8) = SpTanMod(2,2,2,3)
      Amat(5,9) = SpTanMod(2,2,3,3)
      Amat(6,1) = SpTanMod(3,2,1,1)
      Amat(6,2) = SpTanMod(3,2,2,1)
      Amat(6,3) = SpTanMod(3,2,3,1)
      Amat(6,4) = SpTanMod(3,2,1,2)
      Amat(6,5) = SpTanMod(3,2,2,2)
      Amat(6,6) = SpTanMod(3,2,3,2)
      Amat(6,7) = SpTanMod(3,2,1,3)
      Amat(6,8) = SpTanMod(3,2,2,3)
      Amat(6,9) = SpTanMod(3,2,3,3)
      Amat(7,1) = SpTanMod(1,3,1,1)
      Amat(7,2) = SpTanMod(1,3,2,1)
      Amat(7,3) = SpTanMod(1,3,3,1)
      Amat(7,4) = SpTanMod(1,3,1,2)
      Amat(7,5) = SpTanMod(1,3,2,2)
      Amat(7,6) = SpTanMod(1,3,3,2)
      Amat(7,7) = SpTanMod(1,3,1,3)
      Amat(7,8) = SpTanMod(1,3,2,3)
      Amat(7,9) = SpTanMod(1,3,3,3)
      Amat(8,1) = SpTanMod(2,3,1,1)
      Amat(8,2) = SpTanMod(2,3,2,1)
      Amat(8,3) = SpTanMod(2,3,3,1)
      Amat(8,4) = SpTanMod(2,3,1,2)
      Amat(8,5) = SpTanMod(2,3,2,2)
      Amat(8,6) = SpTanMod(2,3,3,2)
      Amat(8,7) = SpTanMod(2,3,1,3)
      Amat(8,8) = SpTanMod(2,3,2,3)
      Amat(8,9) = SpTanMod(2,3,3,3)
      Amat(9,1) = SpTanMod(3,3,1,1)
      Amat(9,2) = SpTanMod(3,3,2,1)
      Amat(9,3) = SpTanMod(3,3,3,1)
      Amat(9,4) = SpTanMod(3,3,1,2)
      Amat(9,5) = SpTanMod(3,3,2,2)
      Amat(9,6) = SpTanMod(3,3,3,2)
      Amat(9,7) = SpTanMod(3,3,1,3)
      Amat(9,8) = SpTanMod(3,3,2,3)
      Amat(9,9) = SpTanMod(3,3,3,3)

      Qmat = ZERO
      Qmat(1,1) = THIRD*(Amat(1,1)+Amat(1,5)+Amat(1,9)) &
                  - (TWO/THREE)*T_tau(1,1)
      Qmat(2,1) = THIRD*(Amat(2,1)+Amat(2,5)+Amat(2,9)) &
                  - (TWO/THREE)*T_tau(2,1)
      Qmat(3,1) = THIRD*(Amat(3,1)+Amat(3,5)+Amat(3,9)) &
                  - (TWO/THREE)*T_tau(3,1)
      Qmat(4,1) = THIRD*(Amat(4,1)+Amat(4,5)+Amat(4,9)) &
                  - (TWO/THREE)*T_tau(1,2)
      Qmat(5,1) = THIRD*(Amat(5,1)+Amat(5,5)+Amat(5,9)) &
                  - (TWO/THREE)*T_tau(2,2)
      Qmat(6,1) = THIRD*(Amat(6,1)+Amat(6,5)+Amat(6,9)) &
                  - (TWO/THREE)*T_tau(3,2)
      Qmat(7,1) = THIRD*(Amat(7,1)+Amat(7,5)+Amat(7,9)) &
                  - (TWO/THREE)*T_tau(1,3)
      Qmat(8,1) = THIRD*(Amat(8,1)+Amat(8,5)+Amat(8,9)) &
                  - (TWO/THREE)*T_tau(2,3)
      Qmat(9,1) = THIRD*(Amat(9,1)+Amat(9,5)+Amat(9,9)) &
                  - (TWO/THREE)*T_tau(3,3)
      Qmat(1,5) = Qmat(1,1)
      Qmat(2,5) = Qmat(2,1)
      Qmat(3,5) = Qmat(3,1)
      Qmat(4,5) = Qmat(4,1)
      Qmat(5,5) = Qmat(5,1)
      Qmat(6,5) = Qmat(6,1)
      Qmat(7,5) = Qmat(7,1)
      Qmat(8,5) = Qmat(8,1)
      Qmat(9,5) = Qmat(9,1)
      Qmat(1,9) = Qmat(1,1)
      Qmat(2,9) = Qmat(2,1)
      Qmat(3,9) = Qmat(3,1)
      Qmat(4,9) = Qmat(4,1)
      Qmat(5,9) = Qmat(5,1)
      Qmat(6,9) = Qmat(6,1)
      Qmat(7,9) = Qmat(7,1)
      Qmat(8,9) = Qmat(8,1)
      Qmat(9,9) = Qmat(9,1)

      if (use_fbar .and. nNode == 8 .and. nIntt == 8) then
        !
        ! This is the tangent using the F-bar method with the
        ! 8-node fully-integrated linear element
        !
        Kuu = Kuu + detMapJC*w(intpt)* &
              ( &
              matmul(matmul(transpose(Gmat),Amat),Gmat) &
              + matmul(transpose(Gmat),matmul(Qmat,(G0mat-Gmat))) &
              )
      else
        !
        ! This is the tangent NOT using the F-bar method with all
        ! other elements
        !
        Kuu = Kuu + detMapJC*w(intpt)* &
              ( &
              matmul(matmul(transpose(Gmat),Amat),Gmat) &
              )
      end if

    end do
    !
    ! End the loop over integration points
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Return to Abaqus the RHS vector and the stiffness matrix
    !
    call AssembleElement(nDim, nNode, ndofel, Ru, Kuu, rhs, amatrx)
    !
    ! End return of RHS and AMATRX
    !--------------------------------------------------------------------

  end subroutine u3d8

  !> Assemble the local element residual and tangent (mechanical only):
  !> maps Ru -> rhs and Kuu -> amatrx.
  subroutine AssembleElement(nDim, nNode, ndofel, Ru, Kuu, rhs, amatrx)

    integer, intent(in)   :: nDim, nNode, ndofel
    real(dp), intent(in)  :: Ru(nDim*nNode,1), Kuu(nDim*nNode,nDim*nNode)
    real(dp), intent(out) :: rhs(ndofel,1), amatrx(ndofel,ndofel)

    integer :: i, j, A11, A12, B11, B12, nDofN

    ! Total number of degrees of freedom per node
    nDofN = ndofel/nNode

    ! init
    rhs = ZERO
    amatrx = ZERO

    if (nDim == 3) then
      !
      ! Assemble the element level residual
      !
      do i = 1, nNode
        A11 = nDofN*(i-1) + 1
        A12 = nDim*(i-1) + 1
        !
        ! displacement
        !
        rhs(A11,1)   = Ru(A12,1)
        rhs(A11+1,1) = Ru(A12+1,1)
        rhs(A11+2,1) = Ru(A12+2,1)
      end do
      !
      ! Assemble the element level tangent matrix
      !
      do i = 1, nNode
        do j = 1, nNode
          A11 = nDofN*(i-1) + 1
          A12 = nDim*(i-1) + 1
          B11 = nDofN*(j-1) + 1
          B12 = nDim*(j-1) + 1
          !
          ! displacement
          !
          amatrx(A11,B11)     = Kuu(A12,B12)
          amatrx(A11,B11+1)   = Kuu(A12,B12+1)
          amatrx(A11,B11+2)   = Kuu(A12,B12+2)
          amatrx(A11+1,B11)   = Kuu(A12+1,B12)
          amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
          amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
          amatrx(A11+2,B11)   = Kuu(A12+2,B12)
          amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
          amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
        end do
      end do
      !
    else
      write(*,*) 'AssembleElement: unsupported nDim=', nDim
      call exit(1)
    end if

  end subroutine AssembleElement

end module uel_element
