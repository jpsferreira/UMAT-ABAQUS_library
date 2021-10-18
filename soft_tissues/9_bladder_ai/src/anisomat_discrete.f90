SUBROUTINE anisomat_discrete(w,sfic,cfic,f,props,  &
        efi,noel,npt,kinc,det,factor,prefdir,ndi)

!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> DISCRETE ANGULAR INTEGRATION SCHEME (icosahedron)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: props(8)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel,npt,kinc
DOUBLE PRECISION, INTENT(IN OUT)         :: det

INTEGER :: j1,k1,l1,m1
DOUBLE PRECISION :: sfibfic(ndi,ndi), cfibfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dwi,ddwi,ddwi1,ddwi2,ddwi3
DOUBLE PRECISION :: bdisp,ang,w,wi,rho,aux2,lambdain,lf
DOUBLE PRECISION :: aic,ais
DOUBLE PRECISION :: avga,maxa,suma,dirmax(ndi),kk1,kk2,ei
DOUBLE PRECISION :: prefdir(nelem,4)

! INTEGRATION SCHEME
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) a
  real ( kind = 8 ) a_xyz(3)
  real ( kind = 8 ) a2_xyz(3)
  real ( kind = 8 ) ai !area of triangle i
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
  

!! initialize the model data
  !     fibers
kk1      = props(4)
kk2      = props(5)
bdisp    = props(6)
lf       = props(7)

  aux=two*(det**(-one))
  aux2=four*(det**(-four/three))
  cfic=zero
  sfic=zero

  aa = zero
  rr = zero
  avga=zero
  maxa=zero
  suma=zero
  dirmax=zero
  
!  Pick a face of the icosahedron, and identify its vertices as A, B, C.
!
  do face = 1, face_num
!
    a = face_point(1,face)
    b = face_point(2,face)
    c = face_point(3,face)
!
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
        CALL deffib(lambdai,mfi,mf0i,f,ndi)
  
        CALL bangle(ang,f,mfi,noel,prefdir,ndi)
        CALL density(rho,ang,bdisp,efi)
        !scaled weight
        ai = ai/(two*pi)
        !rho=rho/(four*pi)
        !strain-like of fiber i
        !lambdai=lambdai*lambdai
        lambdai=lambdai/lf
        ei = lambdai-one
         !calculate fiber sef and sef derivatives values
        if (ei .ge. zero) then
          !fiber sef
          wi   = (kk1/(two*kk2))*(dexp(kk2*ei*ei)-one)
          ! fiber derivatives
          !gho
          !dwi  = kk1*ei*dexp(kk2*ei*ei)
          !ddwi = kk1*dexp(kk2*ei*ei)*(two*kk2*ei*ei+one)
          !bladder model
          dwi=(kk1/(two*sqrt(lambdai)))*((sqrt(lambdai)-one)*(dexp(kk2*(sqrt(lambdai)-one)*(sqrt(lambdai)-one))))
          ddwi1=(kk1*kk2*(sqrt(lambdai)-one)*(sqrt(lambdai)-one)*(dexp(kk2*(sqrt(lambdai)-one)*(sqrt(lambdai)-one))))/(two*lambdai)
          ddwi2=(kk1*(dexp(kk2*(sqrt(lambdai)-one)*(sqrt(lambdai)-one))))/(four*lambdai)
          ddwi3=(-kk1*(sqrt(lambdai)-one)*(dexp(kk2*(sqrt(lambdai)-one)*(sqrt(lambdai)-one))))/(four*(lambdai**(three/two)))
          ddwi=ddwi1+ddwi2+ddwi3
          !update weight
          ais=ai*lambdai**(-one)*lambdai**(-two)
          !stress and material  tangent
        CALL sigfibfic(sfibfic,rho,dwi,mfi,ais,ndi)
        ! 
          aic=ais*lambdai**(-one)*lambdai**(-two)
        CALL csfibfic(cfibfic,rho,dwi,ddwi,mfi,aic,ndi)
!
        DO j1=1,ndi
           DO k1=1,ndi
              sfic(j1,k1)=sfic(j1,k1)+aux*sfibfic(j1,k1)
              DO l1=1,ndi
                DO m1=1,ndi
                  cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux2*cfibfic(j1,k1,l1,m1)
                END DO
              END DO
           END DO
         END DO
!
        w=w+rho*ai*wi
        aa=aa+1
       endif
        rr = rr +  rho*ai
        node_num = node_num + 1  
        area_total = area_total + ai
        !write(*,*) node_num,ang, rho
      end do
    end do
! !
! !  The other subtriangles have the opposite direction from the face.
! !  Generate each in turn, by determining the barycentric coordinates
! !  of the centroid (F1,F2,F3), from which we can also work out the barycentric
! !  coordinates of the vertices of the subtriangle.
! !
!     do f3 = 2, 3 * factor - 4, 3
!       do f2 = 2, 3 * factor - f3 - 2, 3

!         f1 = 3 * factor - f3 - f2

!         call sphere01_triangle_project ( a_xyz, b_xyz, c_xyz, f1, f2, f3, &
!           node_xyz )

!         call sphere01_triangle_project ( &
!           a_xyz, b_xyz, c_xyz, f1 - 2, f2 + 1, f3 + 1, a2_xyz )
!         call sphere01_triangle_project ( &
!           a_xyz, b_xyz, c_xyz, f1 + 1, f2 - 2, f3 + 1, b2_xyz )
!         call sphere01_triangle_project ( &
!           a_xyz, b_xyz, c_xyz, f1 + 1, f2 + 1, f3 - 2, c2_xyz )

!         call sphere01_triangle_vertices_to_area ( a2_xyz, b2_xyz, c2_xyz, ai )

!         !direction of the sphere triangle barycenter - direction i
!         mf0i=node_xyz
!         CALL deffib(lambdai,mfi,mf0i,f,ndi)
  
!         CALL bangle(ang,f,mfi,noel,ndi)
!         CALL density(rho,ang,bdisp,efi)
!         !scaled weight
!         ai = ai/(four*pi)
!         !rho=rho/(four*pi)
!         !strain-like of fiber i
!         lambdai=lambdai*lambdai
!         ei = lambdai-one
!          !calculate fiber sef and sef derivatives values
!         if (ei .ge. zero) then
!           !fiber sef
!           wi   = (kk1/(two*kk2))*(dexp(kk2*ei*ei)-one)
!           ! fiber derivatives
!           dwi  = kk1*ei*dexp(kk2*ei*ei)
!           ddwi = kk1*dexp(kk2*ei*ei)*(two*kk2*ei*ei+one)
!           !stress and material  tangent 
!         CALL sigfibfic(sfibfic,rho,dwi,mfi,ai,ndi)
!         ! 
!         CALL csfibfic(cfibfic,rho,dwi,ddwi,mfi,ai,ndi)
!         !
!         DO j1=1,ndi
!            DO k1=1,ndi
!               sfic(j1,k1)=sfic(j1,k1)+aux*sfibfic(j1,k1)
!               DO l1=1,ndi
!                 DO m1=1,ndi
!                   cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux2*cfibfic(j1,k1,l1,m1)
!                 END DO
!               END DO
!            END DO
!          END DO
! !
!         w=w+rho*ai*wi
!         aa=aa+1
!        endif
!         rr = rr +  rho*ai
!         node_num = node_num + 1  
!         area_total = area_total + ai
!         write(*,*) node_num,ang, rho
!       end do
!     end do

   end do
! !
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

!write(*,*) w,rr,area_total
RETURN
END SUBROUTINE anisomat_discrete
