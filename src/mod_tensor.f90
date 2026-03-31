!> @brief Tensor algebra operations for 3D continuum mechanics.
!>
!> Contains identity tensors, contractions, tensor products,
!> matrix inversion, eigenvalue decomposition, and push-forward/pull-back.
module mod_tensor
  use mod_constants, only: dp, ZERO, ONE, TWO, THREE, HALF
  implicit none
  private
  public :: onem, contraction22, contraction42, contraction24, contraction44
  public :: tensor_product_2, matinv3d
  public :: push2, push4, pull2, pull4
  public :: spectral, jacobi_eigen, eigsrt

contains

  !> Build 2nd-order identity, 4th-order identity, and 4th-order symmetric identity.
  subroutine onem(unit2, unit4, unit4s)
    real(dp), intent(out) :: unit2(3,3), unit4(3,3,3,3), unit4s(3,3,3,3)
    integer :: i, j, k, l

    unit2  = ZERO
    unit4  = ZERO
    unit4s = ZERO

    do i = 1, 3
      unit2(i,i) = ONE
    end do

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            unit4(i,j,k,l)  = unit2(i,k) * unit2(j,l)
            unit4s(i,j,k,l) = HALF * (unit2(i,k)*unit2(j,l) + unit2(i,l)*unit2(j,k))
          end do
        end do
      end do
    end do
  end subroutine onem

  !> Double contraction of two 2nd-order tensors: S = LT : RT
  subroutine contraction22(s, lt, rt)
    real(dp), intent(out) :: s
    real(dp), intent(in)  :: lt(3,3), rt(3,3)
    integer :: i, j

    s = ZERO
    do j = 1, 3
      do i = 1, 3
        s = s + lt(i,j) * rt(i,j)
      end do
    end do
  end subroutine contraction22

  !> Double contraction: 4th-order : 2nd-order -> 2nd-order. S = LT : RT
  subroutine contraction42(s, lt, rt)
    real(dp), intent(out) :: s(3,3)
    real(dp), intent(in)  :: lt(3,3,3,3), rt(3,3)
    integer :: i, j, k, l
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        aux = ZERO
        do k = 1, 3
          do l = 1, 3
            aux = aux + lt(i,j,k,l) * rt(k,l)
          end do
        end do
        s(i,j) = aux
      end do
    end do
  end subroutine contraction42

  !> Double contraction: 2nd-order : 4th-order -> 2nd-order. S = LT : RT
  subroutine contraction24(s, lt, rt)
    real(dp), intent(out) :: s(3,3)
    real(dp), intent(in)  :: lt(3,3), rt(3,3,3,3)
    integer :: i, j, k, l
    real(dp) :: aux

    do k = 1, 3
      do l = 1, 3
        aux = ZERO
        do i = 1, 3
          do j = 1, 3
            aux = aux + lt(i,j) * rt(i,j,k,l)
          end do
        end do
        s(k,l) = aux
      end do
    end do
  end subroutine contraction24

  !> Double contraction of two 4th-order tensors: S = LT : RT
  subroutine contraction44(s, lt, rt)
    real(dp), intent(out) :: s(3,3,3,3)
    real(dp), intent(in)  :: lt(3,3,3,3), rt(3,3,3,3)
    integer :: i, j, k, l, m, n
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            aux = ZERO
            do m = 1, 3
              do n = 1, 3
                aux = aux + lt(i,j,m,n) * rt(m,n,k,l)
              end do
            end do
            s(i,j,k,l) = aux
          end do
        end do
      end do
    end do
  end subroutine contraction44

  !> Dyadic (tensor) product of two 2nd-order tensors: C(i,j,k,l) = A(i,j)*B(k,l)
  subroutine tensor_product_2(a, b, c)
    real(dp), intent(in)  :: a(3,3), b(3,3)
    real(dp), intent(out) :: c(3,3,3,3)
    integer :: i, j, k, l

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            c(i,j,k,l) = a(i,j) * b(k,l)
          end do
        end do
      end do
    end do
  end subroutine tensor_product_2

  !> 3x3 matrix inversion using cofactor formula.
  subroutine matinv3d(a, ainv)
    real(dp), intent(in)  :: a(3,3)
    real(dp), intent(out) :: ainv(3,3)
    real(dp) :: det

    det = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))

    ainv(1,1) =  (a(2,2)*a(3,3) - a(2,3)*a(3,2)) / det
    ainv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2)) / det
    ainv(1,3) =  (a(1,2)*a(2,3) - a(1,3)*a(2,2)) / det
    ainv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1)) / det
    ainv(2,2) =  (a(1,1)*a(3,3) - a(1,3)*a(3,1)) / det
    ainv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1)) / det
    ainv(3,1) =  (a(2,1)*a(3,2) - a(2,2)*a(3,1)) / det
    ainv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1)) / det
    ainv(3,3) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1)) / det
  end subroutine matinv3d

  !> Push-forward of a 2nd-order tensor: sig = (1/det) * F * S * F^T
  subroutine push2(sig, pk, f, det)
    real(dp), intent(out) :: sig(3,3)
    real(dp), intent(in)  :: pk(3,3), f(3,3)
    real(dp), intent(in)  :: det
    integer :: i, j, ii, jj
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        aux = ZERO
        do ii = 1, 3
          do jj = 1, 3
            aux = aux + f(i,ii) * f(j,jj) * pk(ii,jj)
          end do
        end do
        sig(i,j) = aux / det
      end do
    end do
  end subroutine push2

  !> Push-forward of a 4th-order tensor: spatial = (1/det) * F_iI F_jJ F_kK F_lL * MAT_IJKL
  subroutine push4(spatial, mat, f, det)
    real(dp), intent(out) :: spatial(3,3,3,3)
    real(dp), intent(in)  :: mat(3,3,3,3), f(3,3)
    real(dp), intent(in)  :: det
    integer :: i, j, k, l, ii, jj, kk, ll
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            aux = ZERO
            do ii = 1, 3
              do jj = 1, 3
                do kk = 1, 3
                  do ll = 1, 3
                    aux = aux + f(i,ii)*f(j,jj)*f(k,kk)*f(l,ll)*mat(ii,jj,kk,ll)
                  end do
                end do
              end do
            end do
            spatial(i,j,k,l) = aux / det
          end do
        end do
      end do
    end do
  end subroutine push4

  !> Pull-back of a 2nd-order tensor (times det): PK = det * Finv * sig * Finv^T
  subroutine pull2(pk, sig, finv, det)
    real(dp), intent(out) :: pk(3,3)
    real(dp), intent(in)  :: sig(3,3), finv(3,3)
    real(dp), intent(in)  :: det
    integer :: i, j, ii, jj
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        aux = ZERO
        do ii = 1, 3
          do jj = 1, 3
            aux = aux + det * finv(i,ii) * finv(j,jj) * sig(ii,jj)
          end do
        end do
        pk(i,j) = aux
      end do
    end do
  end subroutine pull2

  !> Pull-back of a 4th-order tensor (times det)
  subroutine pull4(mat, spatial, finv, det)
    real(dp), intent(out) :: mat(3,3,3,3)
    real(dp), intent(in)  :: spatial(3,3,3,3), finv(3,3)
    real(dp), intent(in)  :: det
    integer :: i, j, k, l, ii, jj, kk, ll
    real(dp) :: aux

    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          do l = 1, 3
            aux = ZERO
            do ii = 1, 3
              do jj = 1, 3
                do kk = 1, 3
                  do ll = 1, 3
                    aux = aux + det * finv(i,ii)*finv(j,jj) &
                              * finv(k,kk)*finv(l,ll)*spatial(ii,jj,kk,ll)
                  end do
                end do
              end do
            end do
            mat(i,j,k,l) = aux
          end do
        end do
      end do
    end do
  end subroutine pull4

  !> Spectral decomposition of a 3x3 symmetric matrix.
  !> Returns eigenvalues (d) sorted ascending and eigenvectors as columns of v.
  subroutine spectral(a, d, v)
    real(dp), intent(in)  :: a(3,3)
    real(dp), intent(out) :: d(3), v(3,3)
    real(dp) :: acopy(3,3)
    integer  :: nrot

    acopy = a
    call jacobi_eigen(acopy, d, v, nrot)
    call eigsrt(d, v)
  end subroutine spectral

  !> Jacobi eigenvalue algorithm for 3x3 symmetric matrices.
  subroutine jacobi_eigen(a, d, v, nrot)
    real(dp), intent(inout) :: a(3,3)
    real(dp), intent(out)   :: d(3), v(3,3)
    integer,  intent(out)   :: nrot

    integer, parameter :: nmax = 3, max_sweeps = 50
    integer  :: i, ip, iq, p
    real(dp) :: c, g, h, s, sm, t, tau, theta, tresh
    real(dp) :: b(nmax), z(nmax)

    ! Initialize eigenvectors to identity
    v = ZERO
    do ip = 1, 3
      v(ip,ip) = ONE
      b(ip) = a(ip,ip)
      d(ip) = b(ip)
      z(ip) = ZERO
    end do

    nrot = 0
    do i = 1, max_sweeps
      sm = ZERO
      do ip = 1, 2
        do iq = ip + 1, 3
          sm = sm + abs(a(ip,iq))
        end do
      end do
      if (sm == ZERO) return

      if (i < 4) then
        tresh = 0.2_dp * sm / 9.0_dp
      else
        tresh = ZERO
      end if

      do ip = 1, 2
        do iq = ip + 1, 3
          g = 100.0_dp * abs(a(ip,iq))
          if (i > 4 .and. abs(d(ip)) + g == abs(d(ip)) &
                    .and. abs(d(iq)) + g == abs(d(iq))) then
            a(ip,iq) = ZERO
          else if (abs(a(ip,iq)) > tresh) then
            h = d(iq) - d(ip)
            if (abs(h) + g == abs(h)) then
              t = a(ip,iq) / h
            else
              theta = HALF * h / a(ip,iq)
              t = ONE / (abs(theta) + sqrt(ONE + theta*theta))
              if (theta < ZERO) t = -t
            end if
            c = ONE / sqrt(ONE + t*t)
            s = t * c
            tau = s / (ONE + c)
            h = t * a(ip,iq)
            z(ip) = z(ip) - h
            z(iq) = z(iq) + h
            d(ip) = d(ip) - h
            d(iq) = d(iq) + h
            a(ip,iq) = ZERO

            do p = 1, ip - 1
              g = a(p,ip); h = a(p,iq)
              a(p,ip) = g - s*(h + g*tau)
              a(p,iq) = h + s*(g - h*tau)
            end do
            do p = ip + 1, iq - 1
              g = a(ip,p); h = a(p,iq)
              a(ip,p) = g - s*(h + g*tau)
              a(p,iq) = h + s*(g - h*tau)
            end do
            do p = iq + 1, 3
              g = a(ip,p); h = a(iq,p)
              a(ip,p) = g - s*(h + g*tau)
              a(iq,p) = h + s*(g - h*tau)
            end do

            do p = 1, 3
              g = v(p,ip); h = v(p,iq)
              v(p,ip) = g - s*(h + g*tau)
              v(p,iq) = h + s*(g - h*tau)
            end do
            nrot = nrot + 1
          end if
        end do
      end do

      do ip = 1, 3
        b(ip) = b(ip) + z(ip)
        d(ip) = b(ip)
        z(ip) = ZERO
      end do
    end do
  end subroutine jacobi_eigen

  !> Sort eigenvalues ascending and reorder eigenvectors.
  subroutine eigsrt(d, v)
    real(dp), intent(inout) :: d(3), v(3,3)
    integer  :: i, j, k
    real(dp) :: p

    do i = 1, 2
      k = i
      p = d(i)
      do j = i + 1, 3
        if (d(j) < p) then
          k = j
          p = d(j)
        end if
      end do
      if (k /= i) then
        d(k) = d(i)
        d(i) = p
        do j = 1, 3
          p = v(j,i)
          v(j,i) = v(j,k)
          v(j,k) = p
        end do
      end if
    end do
  end subroutine eigsrt

end module mod_tensor
