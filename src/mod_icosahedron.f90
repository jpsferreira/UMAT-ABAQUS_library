!> @brief Icosahedron-based spherical integration for angular integration (AI)
!>        network models.
!>
!> Provides geometry routines for subdividing an icosahedron into spherical
!> triangles, computing their areas (Girard's formula), and projecting
!> barycentric coordinates onto the unit sphere.
!>
!> Based on John Burkardt's sphere triangulation code (GNU LGPL license).
module mod_icosahedron
  use mod_constants, only: dp, ZERO, ONE, TWO, HALF
  implicit none
  private

  ! Icosahedron constants
  integer, parameter, public :: ICOS_POINT_NUM = 12
  integer, parameter, public :: ICOS_EDGE_NUM  = 30
  integer, parameter, public :: ICOS_FACE_NUM  = 20
  integer, parameter, public :: ICOS_FACE_ORDER_MAX = 3

  public :: icos_shape
  public :: sphere01_triangle_project
  public :: sphere01_triangle_vertices_to_area

contains

  !> Set the icosahedron geometry: vertex coordinates, edges, and faces.
  !> Vertices lie on the unit sphere.
  subroutine icos_shape(point_coord, edge_point, face_order, face_point)
    real(dp), intent(out) :: point_coord(3, ICOS_POINT_NUM)
    integer,  intent(out) :: edge_point(2, ICOS_EDGE_NUM)
    integer,  intent(out) :: face_order(ICOS_FACE_NUM)
    integer,  intent(out) :: face_point(ICOS_FACE_ORDER_MAX, ICOS_FACE_NUM)

    real(dp) :: phi, a, b, z

    phi = HALF * (sqrt(5.0_dp) + ONE)
    a = phi / sqrt(ONE + phi*phi)
    b = ONE / sqrt(ONE + phi*phi)
    z = ZERO

    point_coord(1:3, 1:ICOS_POINT_NUM) = reshape( (/ &
         a,  b,  z, &
         a, -b,  z, &
         b,  z,  a, &
         b,  z, -a, &
         z,  a,  b, &
         z,  a, -b, &
         z, -a,  b, &
         z, -a, -b, &
        -b,  z,  a, &
        -b,  z, -a, &
        -a,  b,  z, &
        -a, -b,  z /), (/ 3, ICOS_POINT_NUM /) )

    edge_point(1:2, 1:ICOS_EDGE_NUM) = reshape( (/ &
         1,  2, &
         1,  3, &
         1,  4, &
         1,  5, &
         1,  6, &
         2,  3, &
         2,  4, &
         2,  7, &
         2,  8, &
         3,  5, &
         3,  7, &
         3,  9, &
         4,  6, &
         4,  8, &
         4, 10, &
         5,  6, &
         5,  9, &
         5, 11, &
         6, 10, &
         6, 11, &
         7,  8, &
         7,  9, &
         7, 12, &
         8, 10, &
         8, 12, &
         9, 11, &
         9, 12, &
        10, 11, &
        10, 12, &
        11, 12 /), (/ 2, ICOS_EDGE_NUM /) )

    face_order(1:ICOS_FACE_NUM) = 3

    face_point(1:ICOS_FACE_ORDER_MAX, 1:ICOS_FACE_NUM) = reshape( (/ &
         1,  2,  4, &
         1,  3,  2, &
         1,  4,  6, &
         1,  5,  3, &
         1,  6,  5, &
         2,  3,  7, &
         2,  7,  8, &
         2,  8,  4, &
         3,  5,  9, &
         3,  9,  7, &
         4,  8, 10, &
         4, 10,  6, &
         5,  6, 11, &
         5, 11,  9, &
         6, 10, 11, &
         7,  9, 12, &
         7, 12,  8, &
         8, 12, 10, &
         9, 11, 12, &
        10, 12, 11 /), (/ ICOS_FACE_ORDER_MAX, ICOS_FACE_NUM /) )

  end subroutine icos_shape

  !> Project barycentric coordinates (f1,f2,f3) on a planar triangle (A,B,C)
  !> onto the unit sphere.
  subroutine sphere01_triangle_project(a_xyz, b_xyz, c_xyz, f1, f2, f3, node_xyz)
    real(dp), intent(in)  :: a_xyz(3), b_xyz(3), c_xyz(3)
    integer,  intent(in)  :: f1, f2, f3
    real(dp), intent(out) :: node_xyz(3)
    real(dp) :: node_norm

    node_xyz(1:3) = ( real(f1, dp) * a_xyz(1:3)   &
                    + real(f2, dp) * b_xyz(1:3)   &
                    + real(f3, dp) * c_xyz(1:3) ) &
                    / real(f1 + f2 + f3, dp)

    node_norm = sqrt(sum(node_xyz(1:3)**2))
    node_xyz(1:3) = node_xyz(1:3) / node_norm

  end subroutine sphere01_triangle_project

  !> Compute the area of a spherical triangle on the unit sphere using
  !> Girard's formula: area = A + B + C - pi.
  subroutine sphere01_triangle_vertices_to_area(v1, v2, v3, area)
    real(dp), intent(in)  :: v1(3), v2(3), v3(3)
    real(dp), intent(out) :: area
    real(dp) :: as_side, bs_side, cs_side
    real(dp) :: a_angle, b_angle, c_angle
    real(dp), parameter :: pi = 3.141592653589793_dp

    ! Side lengths (geodesic distances)
    as_side = acos(max(-ONE, min(ONE, dot_product(v2, v3))))
    bs_side = acos(max(-ONE, min(ONE, dot_product(v3, v1))))
    cs_side = acos(max(-ONE, min(ONE, dot_product(v1, v2))))

    ! Spherical angles via half-angle tangent formula
    call sides_to_angles(as_side, bs_side, cs_side, a_angle, b_angle, c_angle)

    ! Girard's formula
    area = a_angle + b_angle + c_angle - pi

  end subroutine sphere01_triangle_vertices_to_area

  ! --------------------------------------------------------------------------
  ! Private helpers
  ! --------------------------------------------------------------------------

  !> Convert spherical triangle side lengths to angles.
  subroutine sides_to_angles(as_s, bs_s, cs_s, a, b, c)
    real(dp), intent(in)  :: as_s, bs_s, cs_s
    real(dp), intent(out) :: a, b, c
    real(dp) :: ssu, tan_a2, tan_b2, tan_c2

    ssu = (as_s + bs_s + cs_s) * HALF

    tan_a2 = sqrt( (sin(ssu - bs_s) * sin(ssu - cs_s)) / &
                   (sin(ssu) * sin(ssu - as_s)) )
    a = TWO * atan(tan_a2)

    tan_b2 = sqrt( (sin(ssu - as_s) * sin(ssu - cs_s)) / &
                   (sin(ssu) * sin(ssu - bs_s)) )
    b = TWO * atan(tan_b2)

    tan_c2 = sqrt( (sin(ssu - as_s) * sin(ssu - bs_s)) / &
                   (sin(ssu) * sin(ssu - cs_s)) )
    c = TWO * atan(tan_c2)

  end subroutine sides_to_angles

end module mod_icosahedron
