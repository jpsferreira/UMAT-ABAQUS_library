!> @brief Shape functions, Gauss schemes, and small matrix helpers for the
!>        8-node brick user element (volume integration path only).
module uel_shape
  use mod_constants, only: dp, ZERO, ONE, THREE
  use mod_tensor, only: matinv3d
  implicit none
  private
  public :: calcShape3DLinear, mapShape3D, xint3D1pt, xint3D8pt, mdet, identity3

  real(dp), parameter :: EIGHTH = ONE/8.0_dp

contains

  !> Integration point location and weight for 1-point Gauss integration
  !> on the 8-node brick.
  subroutine xint3D1pt(xi, w, nInttPt)
    real(dp), intent(out) :: xi(1,3), w(1)
    integer, intent(out)  :: nInttPt

    ! Init
    w = ZERO
    xi = ZERO

    ! Number of Gauss points
    nInttPt = 1

    ! Gauss weights
    w(1) = 8.0_dp

    ! Gauss pt locations in master element
    xi(1,1) = ZERO
    xi(1,2) = ZERO
    xi(1,3) = ZERO

  end subroutine xint3D1pt

  !> Integration point locations and weights for 8-point Gauss integration
  !> on the 8-node brick.
  subroutine xint3D8pt(xi, w, nInttPt)
    real(dp), intent(out) :: xi(8,3), w(8)
    integer, intent(out)  :: nInttPt

    ! Init
    w = ZERO
    xi = ZERO

    ! Number of Gauss points
    nInttPt = 8

    ! Gauss weights
    w(1) = ONE
    w(2) = ONE
    w(3) = ONE
    w(4) = ONE
    w(5) = ONE
    w(6) = ONE
    w(7) = ONE
    w(8) = ONE

    ! Gauss pt locations in master element
    xi(1,1) = -sqrt(ONE/THREE)
    xi(1,2) = -sqrt(ONE/THREE)
    xi(1,3) = -sqrt(ONE/THREE)
    xi(2,1) = sqrt(ONE/THREE)
    xi(2,2) = -sqrt(ONE/THREE)
    xi(2,3) = -sqrt(ONE/THREE)
    xi(3,1) = -sqrt(ONE/THREE)
    xi(3,2) = sqrt(ONE/THREE)
    xi(3,3) = -sqrt(ONE/THREE)
    xi(4,1) = sqrt(ONE/THREE)
    xi(4,2) = sqrt(ONE/THREE)
    xi(4,3) = -sqrt(ONE/THREE)
    xi(5,1) = -sqrt(ONE/THREE)
    xi(5,2) = -sqrt(ONE/THREE)
    xi(5,3) = sqrt(ONE/THREE)
    xi(6,1) = sqrt(ONE/THREE)
    xi(6,2) = -sqrt(ONE/THREE)
    xi(6,3) = sqrt(ONE/THREE)
    xi(7,1) = -sqrt(ONE/THREE)
    xi(7,2) = sqrt(ONE/THREE)
    xi(7,3) = sqrt(ONE/THREE)
    xi(8,1) = sqrt(ONE/THREE)
    xi(8,2) = sqrt(ONE/THREE)
    xi(8,3) = sqrt(ONE/THREE)

  end subroutine xint3D8pt

  !> Shape functions and their local derivatives at a given integration
  !> point of the 8-node linear brick master element.
  !>
  !>      8-----------7
  !>     /|          /|       zeta
  !>    / |         / |
  !>   5-----------6  |       |     eta
  !>   |  |        |  |       |   /
  !>   |  |        |  |       |  /
  !>   |  4--------|--3       | /
  !>   | /         | /        |/
  !>   |/          |/         O--------- xi
  !>   1-----------2        origin at cube center
  !>
  !> sh(i) = shape function of node i at the intpt
  !> dshxi(i,j) = derivative wrt j direction of shape fn of node i
  subroutine calcShape3DLinear(nInttPt, xi_int, intpt, sh, dshxi)
    integer, intent(in)   :: nInttPt, intpt
    real(dp), intent(in)  :: xi_int(nInttPt,3)
    real(dp), intent(out) :: sh(8), dshxi(8,3)

    real(dp) :: xi, eta, zeta

    ! Location in the master element
    xi = xi_int(intpt,1)
    eta = xi_int(intpt,2)
    zeta = xi_int(intpt,3)

    ! The shape functions
    sh(1) = EIGHTH*(ONE - xi)*(ONE - eta)*(ONE - zeta)
    sh(2) = EIGHTH*(ONE + xi)*(ONE - eta)*(ONE - zeta)
    sh(3) = EIGHTH*(ONE + xi)*(ONE + eta)*(ONE - zeta)
    sh(4) = EIGHTH*(ONE - xi)*(ONE + eta)*(ONE - zeta)
    sh(5) = EIGHTH*(ONE - xi)*(ONE - eta)*(ONE + zeta)
    sh(6) = EIGHTH*(ONE + xi)*(ONE - eta)*(ONE + zeta)
    sh(7) = EIGHTH*(ONE + xi)*(ONE + eta)*(ONE + zeta)
    sh(8) = EIGHTH*(ONE - xi)*(ONE + eta)*(ONE + zeta)

    ! The first derivatives
    dshxi(1,1) = -EIGHTH*(ONE - eta)*(ONE - zeta)
    dshxi(1,2) = -EIGHTH*(ONE - xi)*(ONE - zeta)
    dshxi(1,3) = -EIGHTH*(ONE - xi)*(ONE - eta)
    dshxi(2,1) = EIGHTH*(ONE - eta)*(ONE - zeta)
    dshxi(2,2) = -EIGHTH*(ONE + xi)*(ONE - zeta)
    dshxi(2,3) = -EIGHTH*(ONE + xi)*(ONE - eta)
    dshxi(3,1) = EIGHTH*(ONE + eta)*(ONE - zeta)
    dshxi(3,2) = EIGHTH*(ONE + xi)*(ONE - zeta)
    dshxi(3,3) = -EIGHTH*(ONE + xi)*(ONE + eta)
    dshxi(4,1) = -EIGHTH*(ONE + eta)*(ONE - zeta)
    dshxi(4,2) = EIGHTH*(ONE - xi)*(ONE - zeta)
    dshxi(4,3) = -EIGHTH*(ONE - xi)*(ONE + eta)
    dshxi(5,1) = -EIGHTH*(ONE - eta)*(ONE + zeta)
    dshxi(5,2) = -EIGHTH*(ONE - xi)*(ONE + zeta)
    dshxi(5,3) = EIGHTH*(ONE - xi)*(ONE - eta)
    dshxi(6,1) = EIGHTH*(ONE - eta)*(ONE + zeta)
    dshxi(6,2) = -EIGHTH*(ONE + xi)*(ONE + zeta)
    dshxi(6,3) = EIGHTH*(ONE + xi)*(ONE - eta)
    dshxi(7,1) = EIGHTH*(ONE + eta)*(ONE + zeta)
    dshxi(7,2) = EIGHTH*(ONE + xi)*(ONE + zeta)
    dshxi(7,3) = EIGHTH*(ONE + xi)*(ONE + eta)
    dshxi(8,1) = -EIGHTH*(ONE + eta)*(ONE + zeta)
    dshxi(8,2) = EIGHTH*(ONE - xi)*(ONE + zeta)
    dshxi(8,3) = EIGHTH*(ONE - xi)*(ONE + eta)

  end subroutine calcShape3DLinear

  !> Map derivatives of shape functions from the xi-eta-zeta domain to the
  !> x-y-z domain. Returns stat=0 (without modifying dsh) if the mapping
  !> Jacobian determinant is non-positive, so the caller can cut back.
  subroutine mapShape3D(nNode, dshxi, coords, dsh, detMapJ, stat)
    integer, intent(in)   :: nNode
    real(dp), intent(in)  :: dshxi(nNode,3), coords(3,nNode)
    real(dp), intent(out) :: dsh(nNode,3), detMapJ
    integer, intent(out)  :: stat

    integer  :: i, j, k
    real(dp) :: mapJ(3,3), mapJ_inv(3,3)

    stat = 1

    ! Calculate the mapping Jacobian matrix
    mapJ = ZERO
    do i = 1, 3
      do j = 1, 3
        do k = 1, nNode
          mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
        end do
      end do
    end do

    ! Determinant check: a non-positive mapping Jacobian flags a bad
    ! (inverted) element configuration
    call mdet(mapJ, detMapJ)
    if (detMapJ <= ZERO) then
      write(*,*) 'WARNING: mapShape3D: det of mapping Jacobian =', detMapJ
      stat = 0
      return
    end if

    ! Calculate the inverse of the Jacobian (reused from mod_tensor)
    call matinv3d(mapJ, mapJ_inv)

    ! Calculate first derivatives wrt x, y, z
    dsh = transpose(matmul(mapJ_inv, transpose(dshxi)))

  end subroutine mapShape3D

  !> Determinant of a 3x3 matrix.
  subroutine mdet(A, det)
    real(dp), intent(in)  :: A(3,3)
    real(dp), intent(out) :: det

    det = A(1,1)*A(2,2)*A(3,3) &
        + A(1,2)*A(2,3)*A(3,1) &
        + A(1,3)*A(2,1)*A(3,2) &
        - A(3,1)*A(2,2)*A(1,3) &
        - A(3,2)*A(2,3)*A(1,1) &
        - A(3,3)*A(2,1)*A(1,2)

  end subroutine mdet

  !> 3x3 identity matrix.
  subroutine identity3(iden)
    real(dp), intent(out) :: iden(3,3)
    integer :: i, j

    do i = 1, 3
      do j = 1, 3
        if (i == j) then
          iden(i,j) = ONE
        else
          iden(i,j) = ZERO
        end if
      end do
    end do

  end subroutine identity3

end module uel_shape
