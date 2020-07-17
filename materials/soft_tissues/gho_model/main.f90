!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'ABA_PARAM.INC'
INCLUDE 'PARAM_UMAT.INC'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 23, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),STATEVP(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


CHARACTER*8 filename
integer un


i=1.0d0
j=1.0d0
DO i=1,NTENS
    DO j=1,NTENS
        DDSDDE(i,j)=0.0D0
    ENDDO
    STRESS(i)=0.0D0
ENDDO


!open(unit=20,file='output.txt',status='replace',action="write", &
!            position="append",iostat=ierr,access='sequential')

!
! DEFORMATION GRADIENT
 !PURE SHEAR
! DFGRD1(1,1)=   1.2D0
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.D0
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.D0/DFGRD1(1,1)
!
!
!  !ARBITRARY
! DFGRD1(1,1)=   1.10000000000000
! DFGRD1(1,2)=  0.200000000000000E+000
! DFGRD1(1,3)=  0.300000000000000E+000
! DFGRD1(2,1)=  0.400000000000000E+000
! DFGRD1(2,2)=   2.50000000000000
! DFGRD1(2,3)=  0.600000000000000E+000
! DFGRD1(3,1)=  0.700000000000000E+000
! DFGRD1(3,2)=  0.800000000000000E+000
! DFGRD1(3,3)=   3.90000000000000

!
!      !DILATATION
! DFGRD1(1,1)=   1.20000000000000
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.1d0
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.1d0
!
!   !EXTENSION
! DFGRD1(1,1)=   1.1d0
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.d0/SQRT(DFGRD1(1,1))
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.d0/SQRT(DFGRD1(1,1))

  DFGRD1(1,1)=   1.0D0
 !DFGRD1(1,2)=  0.3d0
 DFGRD1(1,3)=  0.0E+000
 DFGRD1(2,1)=  0.0D0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0E+000
 DFGRD1(3,1)=  0.0E+000
 DFGRD1(3,2)=  0.0E+000
 DFGRD1(3,3)=   1.0D0
!!!
time(1)=0.d0
time(2)=0.d0
call UEXTERNALDB(0,0,time,0.D0,0,0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=1000.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ISOTROPIC MATRIX PARAMS
! C10=
PROPS(2)=1.00d0
! C01
PROPS(3)=0.00d0

PROPS(4)=1.00d0
PROPS(5)=1.0d0
PROPS(6)=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(7)=0
!tau1 %
PROPS(8)=2.0d0
!teta1
PROPS(9)=0.835d0
!tau2 %
PROPS(10)=1.2d0
!teta2
PROPS(11)=7.0d0
!tau3 %
PROPS(12)=12.d0
!teta3
PROPS(13)=2.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STATEV=0.D0
STATEVP=0.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
erf=0.d0
RHO=0.D0
!
OPEN (UNIT=150, FILE='data.out', STATUS='UNKNOWN')
        rewind(150)
!
!
 DFGRD1(1,1)=  1.0D0
 DFGRD1(1,2)=  0.0D0
 DFGRD1(1,3)=  0.0d0
 DFGRD1(2,1)=  0.0d0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0d0
 DFGRD1(3,1)=  0.0d0
 DFGRD1(3,2)=  0.0d0
 DFGRD1(3,3)=  1.0D0

!################################################################################################!
!!     TENSILE MONOTONIC LOAD TEST

 DFGRD1(1,1)=  1.6D0
 DFGRD1(1,2)=  0.0D0
 DFGRD1(1,3)=  0.0d0
 DFGRD1(2,1)=  0.0d0
 DFGRD1(2,2)=  1/sqrt(DFGRD1(1,1))
 DFGRD1(2,3)=  0.0d0
 DFGRD1(3,1)=  0.0d0
 DFGRD1(3,2)=  0.0d0
 DFGRD1(3,3)=  1/sqrt(DFGRD1(1,1))


 CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

 write(*,*) stress
  write(*,*)
   write(*,*) DDSDDE

!################################################################################################!
close(150)
!
END PROGRAM
