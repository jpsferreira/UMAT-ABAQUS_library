
subroutine icos_shape ( point_num, edge_num, face_num, face_order_max, &
  point_coord, edge_point, face_order, face_point )

!*****************************************************************************80
!
!! ICOS_SHAPE describes an icosahedron.
!
!  Discussion:
!
!    The input data required for this routine can be retrieved from ICOS_SIZE.
!
!    The vertices lie on the unit sphere.
!
!    The dual of an icosahedron is a dodecahedron.
!
!    The data has been rearranged from a previous assignment.  
!    The STRIPACK program refuses to triangulate data if the first
!    three nodes are "collinear" on the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points (12).
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges (30).
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces (20).
!
!    Input, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum number of 
!    vertices per face (3).
!
!    Output, real ( kind = 8 ) POINT_COORD(3,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) EDGE_POINT(2,EDGE_NUM), the points that 
!    make up each edge, listed in ascending order of their indexes.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer ( kind = 4 ) FACE_POINT(FACE_ORDER_MAX,FACE_NUM); 
!    FACE_POINT(I,J) is the index of the I-th point in the J-th face.  The
!    points are listed in the counter clockwise direction defined
!    by the outward normal at the face.  The nodes of each face are ordered 
!    so that the lowest index occurs first.  The faces are then sorted by
!    nodes.
!
  use global

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), parameter :: edge_order = 2
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) edge_point(edge_order,edge_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) face_point(face_order_max,face_num)
  real ( kind = 8 ) phi
  real ( kind = 8 ) point_coord(3,point_num)
  real ( kind = 8 ) z
!
!  Set the point coordinates.
!
  phi = 0.5D+00 * ( sqrt ( 5.0D+00 ) + 1.0D+00 )

  a = phi / sqrt ( 1.0D+00 + phi * phi )
  b = 1.0D+00 / sqrt ( 1.0D+00 + phi * phi )
  z = 0.0D+00
!
!  A*A + B*B + Z*Z = 1.
!
  point_coord(1:3,1:point_num) = reshape ( (/ &
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
     -a, -b,  z /), (/ 3, point_num /) )
!
!  Set the edges.
!
  edge_point(1:edge_order,1:edge_num) = reshape ( (/ &
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
    11, 12 /), (/ edge_order, edge_num /) )
!
!  Set the face orders.
!
  face_order(1:face_num) = (/ &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
!
!  Set the faces.
!
  face_point(1:face_order_max,1:face_num) = reshape ( (/ &
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
    10, 12, 11 /), (/ face_order_max, face_num /) )

  return
end
subroutine icos_size ( point_num, edge_num, face_num, face_order_max )
!*****************************************************************************80
!
!! ICOS_SIZE gives "sizes" for an icosahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) FACE_ORDER_MAX, the maximum order of any face.
!
  use global

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) point_num

  point_num = 12
  edge_num = 30
  face_num = 20
  face_order_max = 3

  return
end

subroutine sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
  node_xyz )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_PROJECT projects from plane to spherical triangle.
!
!  Discussion:
!
!    We assume that points A, B and C lie on the unit sphere, and they
!    thus define a spherical triangle.
!
!    They also, of course, define a planar triangle.
!
!    Let (F1,F2,F3) be the barycentric coordinates of a point in this 
!    planar triangle.
!
!    This function determines the coordinates of the point in the planar
!    triangle identified by the barycentric coordinates, and returns the
!    coordinates of the projection of that point onto the unit sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A_XYZ(3), B_XYZ(3), C_XYZ(3), the coordinates
!    of the points A, B, and C.
!
!    Input, integer ( kind = 4 ) F1, F2, F3, the barycentric coordinates
!    of a point in the triangle ABC.  Normally, these coordinates would
!    be real numbers, and would sum to 1.  For convenience, we allow these
!    to be integers which must be divided by F1+F2+F3.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3), the coordinates of the 
!    point on the unit sphere which is the projection of the point on the plane
!    whose barycentric coordinates with respect to A, B, and C is
!    (F1,F2,F3)/(F1+F2+F3).
!
  use global

  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) c_xyz(3)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  real ( kind = 8 ) node_norm
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ) r8vec_norm

  node_xyz(1:3) = &
    ( real ( f1,           kind = 8 ) * a_xyz(1:3)   &
    + real (      f2,      kind = 8 ) * b_xyz(1:3)   &
    + real (           f3, kind = 8 ) * c_xyz(1:3) ) &
    / real ( f1 + f2 + f3, kind = 8 )

  node_norm = r8vec_norm ( 3, node_xyz(1:3) )

  node_xyz(1:3) = node_xyz(1:3) / node_norm

  return
end

subroutine sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_SIDES computes spherical triangle sides.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the spherical
!    triangle.
!
!    Output, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
  use global

  real ( kind = 8 ) as
  real ( kind = 8 ) bs
  real ( kind = 8 ) cs
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  as = acos ( dot_product ( v2(1:3), v3(1:3) ) )
  bs = acos ( dot_product ( v3(1:3), v1(1:3) ) )
  cs = acos ( dot_product ( v1(1:3), v2(1:3) ) )

  return
end




subroutine sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_SIDES_TO_ANGLES computes spherical triangle angles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AS, BS, CS, the (geodesic) length of the 
!    sides of the triangle.
!
!    Output, real ( kind = 8 ) A, B, C, the spherical angles of the triangle.
!    Angle A is opposite the side of length AS, and so on.
!
  use global

  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) asu
  real ( kind = 8 ) b
  real ( kind = 8 ) bs
  real ( kind = 8 ) bsu
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) csu
  real ( kind = 8 ) ssu
  real ( kind = 8 ) tan_a2
  real ( kind = 8 ) tan_b2
  real ( kind = 8 ) tan_c2

  asu = as
  bsu = bs
  csu = cs
  ssu = ( asu + bsu + csu ) / 2.0D+00

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - asu )     ) )

  a = 2.0D+00 * atan ( tan_a2 )

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) )

  b = 2.0D+00 * atan ( tan_b2 )

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) / &
                  ( sin ( ssu ) * sin ( ssu - csu )     ) )

  c = 2.0D+00 * atan ( tan_c2 )

  return
end
subroutine sphere01_triangle_vertices_to_area ( v1, v2, v3, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_VERTICES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vertices of the triangle.
!
!    Output, real ( kind = 8 ) AREA, the area of the sphere.
!
  use global

  real ( kind = 8 ) area
  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) b
  real ( kind = 8 ) bs
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)
!
!  Compute the lengths of the sides of the spherical triangle.
!
  call sphere01_triangle_vertices_to_sides ( v1, v2, v3, as, bs, cs )
!
!  Get the spherical angles.
!
  call sphere01_triangle_sides_to_angles ( as, bs, cs, a, b, c )
!
!  Get the area.
!
  call sphere01_triangle_angles_to_area ( a, b, c, area )

  return
end


subroutine sphere01_triangle_angles_to_area ( a, b, c, area )

!*****************************************************************************80
!
!! SPHERE01_TRIANGLE_ANGLES_TO_AREA computes the area of a spherical triangle.
!
!  Discussion:
!
!    A unit sphere in 3D satisfies the equation:
!
!      X^2 + Y^2 + Z^2 = 1
!
!    A spherical triangle is specified by three points on the surface
!    of the sphere.
!
!    The area formula is known as Girard's formula.
!
!    The area of a spherical triangle is:
!
!      AREA = ( A + B + C - PI )
!
!    where A, B and C are the (surface) angles of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the angles of the triangle.
!
!    Output, real ( kind = 8 ) AREA, the area of the sphere.
!
  use global

  real ( kind = 8 ) area
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
!
!  Apply Girard's formula.
!
  area = a + b + c - pi

  return
end

subroutine polyterm_exponent ( action, e )

!*****************************************************************************80
!
!! POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 3 ) ACTION.
!    'GET' asks the routine to return the current values in E.
!    'SET' asks the routine to set the current values to E.
!
!    Input/output, integer ( kind = 4 ) E(3), storage used to set or get values.
!
  use global

  character ( len = * )  action
  integer   ( kind = 4 ) e(3)
  integer   ( kind = 4 ), save, dimension ( 3 ) :: e_save = (/ 0, 0, 0 /)
  character ( len = 80 ) text
  character ( len = 80 ) text2

  if ( action(1:1) == 'G' ) then

    e(1:3) = e_save(1:3)

  else if ( action(1:1) == 'P' ) then

    write ( *, '(a)' ) ' '

    if ( all ( e_save(1:3) == 0 ) ) then

      text = 'P(X,Y,Z) = 1'

    else

      text = 'P(X,Y,Z) = '

      if ( e_save(1) == 0 ) then

      else if ( e_save(1) == 1 ) then

        call s_cat ( text, ' X', text )

      else

        call s_cat ( text, ' X^', text )

        write ( text2, '(i2)' ) e_save(1)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if

      if ( e_save(2) == 0 ) then

      else if ( e_save(2) == 1 ) then

        call s_cat ( text, ' Y', text )

      else

        call s_cat ( text, ' Y^', text )

        write ( text2, '(i2)' ) e_save(2)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
       
      if ( e_save(3) == 0 ) then

      else if ( e_save(3) == 1 ) then

        call s_cat ( text, ' Z', text )

      else

        call s_cat ( text, ' Z^', text )

        write ( text2, '(i2)' ) e_save(3)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
 
    end if

    write ( *, '(a)' ) trim ( text )
    
  else if ( action(1:1) == 'S' ) then

    e_save(1:3) = e(1:3)

  end if

  return
end
subroutine polyterm_value_3d ( n, x, f )

!*****************************************************************************80
!
!! POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
!
!  Discussion:
!
!    The polynomial term has the form:
!
!      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
!
!    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(3,N), the points where the polynomial term 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) F(N), the value of the polynomial term.
!
  use global

  integer ( kind = 4 ) n

  integer ( kind = 4 ) e(3)
  real ( kind = 8 ) f(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(3,n)

  call polyterm_exponent ( 'GET', e )

  f(1:n) = 1.0D+00

  do i = 1, 3

    if ( e(i) /= 0 ) then
      f(1:n) = f(1:n) * x(i,1:n)**e(i)
    end if

  end do
  
  return
end

function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  use global

  real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real ( kind = 8 ), parameter :: eps = 2.22D-16
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) sum
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 171.624D+00
  real ( kind = 8 ) xden
  real ( kind = 8 ), parameter :: xinf = 1.0D+30
  real ( kind = 8 ), parameter :: xminin = 2.23D-308
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) ysq
  real ( kind = 8 ) z

  parity = .false.
  fact = 1.0D+00
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= 0.0D+00 ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= 0.0D+00 ) then

      if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + 1.0D+00

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = 1.0D+00 / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < 12.0D+00 ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < 1.0D+00 ) then

      z = y
      y = y + 1.0D+00
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - 1.0D+00

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = 0.0D+00
    xden = 1.0D+00
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + 1.0D+00
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + 1.0D+00
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - 0.5D+00 ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= 1.0D+00 ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  use global

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8vec_norm ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM, the L2 norm of A.
!
  use global

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) r8vec_norm

  r8vec_norm = sqrt ( sum ( a(1:n)**2 ) )

  return
end
subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the vector to be polarized.
!
!    Input, real ( kind = 8 ) P(N), the polarizing direction.
!
!    Output, real ( kind = 8 ) A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  use global

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_dot_p
  real ( kind = 8 ) a_normal(n)
  real ( kind = 8 ) a_parallel(n)
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
    return
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  use global

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if

  return
end
subroutine sphere01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SPHERE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit sphere.
!
!  Discussion:
!
!    The integration region is 
!
!      X^2 + Y^2 + Z^2 = 1.
!
!    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Academic Press, 1984, page 263.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) E(3), the exponents of X, Y and Z in the 
!    monomial.  Each exponent must be nonnegative.
!
!    Output, real ( kind = 8 ) INTEGRAL, the integral.
!
  use global

  integer ( kind = 4 ) e(3)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma

  if ( any ( e(1:3) < 0 ) ) then
    integral = - huge ( integral )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    write ( *, '(a,i8)' ) '  E(1) = ', e(1)
    write ( *, '(a,i8)' ) '  E(2) = ', e(2)
    write ( *, '(a,i8)' ) '  E(3) = ', e(3)
    stop
  end if

  if ( all ( e(1:3) == 0 ) ) then

    integral = 2.0D+00 * sqrt ( pi**3 ) / r8_gamma ( 1.5D+00 )

  else if ( any ( mod ( e(1:3), 2 ) == 1 ) ) then

    integral = 0.0D+00

  else

    integral = 2.0D+00

    do i = 1, 3
      integral = integral * r8_gamma ( 0.5D+00 * real ( e(i) + 1, kind = 8 ) )
    end do

    integral = integral &
      / r8_gamma ( 0.5D+00 * ( real ( sum ( e(1:3) + 1 ), kind = 8 ) ) )

  end if

  return
end
