SUBROUTINE naffnetfic_discrete(sfic,cfic,f,filprops,naffprops,  &
        det,factor,ndi)
!
!>    NON-AFFINE NETWORK:'FICTICIUOUS' CAUCHY STRESS AND ELASTICITY TNSR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: naffprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: h(ndi,ndi),hh(ndi,ndi,ndi,ndi),  &
    hi(ndi,ndi),hhi(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dw,ddw,aux1,aux2,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda,lambda0,n,pp
DOUBLE PRECISION :: r0f,r0c,etac



! INTEGRATION SCHEME
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) a
  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) ai !area o triangle i
  real ( kind = 8 ) area_total
  integer ( kind = 4 ) b
  real ( kind = 8 ) b_xyz(3)
  real ( kind = 8 ) b2_xyz(3)
  integer ( kind = 4 ) c
  real ( kind = 8 ) c_xyz(3)
  real ( kind = 8 ) c2_xyz(3)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_point
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) f3
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: face_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: face_point
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ) factor
  !external             fun
  real ( kind = 8 ) node_xyz(3)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) rr, aa
  real ( kind = 8 ) v



!  Size the icosahedron.
!
  call icos_size ( point_num, edge_num, face_num, face_order_max )
!
!  Set the icosahedron.
!
  allocate ( point_coord(1:3,1:point_num) )
  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )
!
!  Initialize the integral data.
!
  rr = 0.0D+00
  area_total = 0.0D+00
  node_num = 0
  
  
!     FILAMENT
l       = filprops(1)
r0f     = filprops(2)
r0c     = filprops(3)
etac    = filprops(4)
mu0     = filprops(5)
beta    = filprops(6)
b0      = filprops(7)
lambda0 = filprops(8)
!     NETWORK
n       = naffprops(1)
pp      = naffprops(2)


lambda=zero
h=zero
hh=zero
hi=zero
hhi=zero


r0=r0f+r0c



 !  Pick a face of the icosahedron, and identify its vertices as A, B, C.
!
  do face = 1, face_num

    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)

    a_xyz(1:3) = point_coord(1:3,a)
    b_xyz(1:3) = point_coord(1:3,b)
    c_xyz(1:3) = point_coord(1:3,c)
!
!  Some subtriangles will have the same direction as the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 1, 3 * factor - 2, 3
      do f2 = 1, 3 * factor - f3 - 1, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 2, f2 - 1, f3 - 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 + 2, f3 - 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 1, f2 - 1, f3 + 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, ai )

        !direction of the sphere triangle barycenter - direction i
        mf0i=node_xyz
        CALL deffil(lambdai,mfi,mf0i,f,ndi)
        lambda=lambda+(lambdai**pp)*ai
  
        CALL hfilfic(hi,hhi,pp,lambdai,mfi,ai,ndi)
          
        DO j1=1,ndi
          DO k1=1,ndi
             h(j1,k1)=h(j1,k1)+hi(j1,k1)
             DO l1=1,ndi
               DO m1=1,ndi
                 hh(j1,k1,l1,m1)=hh(j1,k1,l1,m1)+hhi(j1,k1,l1,m1)
               END DO
             END DO
          END DO
        END DO
        
        node_num = node_num + 1
        area_total = area_total + ai

      end do
    end do
!
!  The other subtriangles have the opposite direction from the face.
!  Generate each in turn, by determining the barycentric coordinates
!  of the centroid (F1,F2,F3), from which we can also work out the barycentric
!  coordinates of the vertices of the subtriangle.
!
    do f3 = 2, 3 * factor - 4, 3
      do f2 = 2, 3 * factor - f3 - 2, 3

        f1 = 3 * factor - f3 - f2

        call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
          node_xyz )

        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
        call sphere01_triangle_project ( &
          a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

        call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, ai )

              !direction of the sphere triangle barycenter - direction i
        mf0i=node_xyz
        CALL deffil(lambdai,mfi,mf0i,f,ndi)
        lambda=lambda+(lambdai**pp)*ai
  
        CALL hfilfic(hi,hhi,pp,lambdai,mfi,ai,ndi)
          
        DO j1=1,ndi
          DO k1=1,ndi
             h(j1,k1)=h(j1,k1)+hi(j1,k1)
             DO l1=1,ndi
               DO m1=1,ndi
                 hh(j1,k1,l1,m1)=hh(j1,k1,l1,m1)+hhi(j1,k1,l1,m1)
               END DO
             END DO
          END DO
        END DO
        
        node_num = node_num + 1
        area_total = area_total + ai

      end do
    end do

  end do
!

lambda = lambda/area_total
lambda=lambda**(pp**(-one))

CALL fil(fi,ffi,dw,ddw,lambda,lambda0,l,r0,mu0,beta,b0)

cfic=zero
sfic=zero
aux=n*(det**(-one))/area_total
aux1=aux*dw*lambda**(one-pp)
aux2=aux*(ddw*(lambda**(two*(one-pp)))-(pp-one)*dw*(lambda**(one-two*pp)))

DO j1=1,ndi
  DO k1=1,ndi
    sfic(j1,k1)=aux1*h(j1,k1)
    DO l1=1,ndi
      DO m1=1,ndi
        cfic(j1,k1,l1,m1)=aux1*hh(j1,k1,l1,m1)+aux2*h(j1,k1)*h(l1,m1)
      END DO
    END DO
  END DO
END DO


!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )


RETURN


END SUBROUTINE naffnetfic_discrete
