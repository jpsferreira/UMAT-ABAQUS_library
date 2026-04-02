PROGRAM TEST_GENERAL_UMAT
use,intrinsic :: ISO_Fortran_env
use global
INCLUDE 'aba_param.inc'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK
      COMMON /KFIL/MF0
      COMMON /KFILR/RW
      COMMON /KFILP/PREFDIR
      COMMON /KFILC/CHEM
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),PREFDIR(NELEM,4),CHEM(2,7)

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 18, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)
!
CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
i=1.0d0
j=1.0d0
DO i=1,NTENS
    DO j=1,NTENS
        DDSDDE(i,j)=0.0D0
    ENDDO
    STRESS(i)=0.0D0
ENDDO
!
! DEFORMATION GRADIENT
DFGRD1(1,1)= 1.1D0
DFGRD1(1,2)= 0.0D0
DFGRD1(1,3)= 0.0D0
DFGRD1(2,1)= 0.0D0
DFGRD1(2,2)= 1.0D0/DFGRD1(1,1)
DFGRD1(2,3)= 0.0D0
DFGRD1(3,1)= 0.0D0
DFGRD1(3,2)= 0.0D0
DFGRD1(3,3)= 1.0D0/DFGRD1(1,1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ISOTROPIC MATRIX PARAMS
! C10=
PROPS(2)=0.00d0
! C01
PROPS(3)=0.00d0
!PHI....
PROPS(4)=1.0000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SINGLE FILAMENT PARAMS
!L
PROPS(5)=0.988d0
!R0F
PROPS(6)=0.818d0-0.014d0
!R0C
PROPS(7)=0.014d0
!ETAC
PROPS(8)=0.667d0
!mu0
PROPS(9)=38600.0d0
!PROPS(9)=111111111111111111111111111138600.0d0
!beta
PROPS(10)=0.438d0
!PROPS(10)=0.5d0
!B0 = tk*lp*k0
PROPS(11)=294.d0*16.d0*1.38d-5
!lambda0.
PROPS(12)=1.007d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NON-AFFINE NETWORK PARAMS
!n - isotropic filaments per unit volume
PROPS(13)=15.2d0
!P....
PROPS(14)=3.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!AFFINE CONTRACTILE NETWORK PARAMS
!NA - fil in affine net per unit volume
PROPS(15)=0.2d0
!B - dispersion
PROPS(16)=0.001d0
!FRIC - activation
PROPS(17)=1.1e1
!PROPS(15)=6.0e4
!FFMAX - activation
PROPS(18)=1.1e1
! !
STATEV=0.D0
!
erf=0.d0
RHO=0.D0
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
!
!################################################################################################!
!!     TENSILE MONOTONIC LOAD TEST
 DFGRD1(1,1)=  1.0D0
 DFGRD1(1,2)=  0.1D0
 DFGRD1(1,3)=  0.0d0
 DFGRD1(2,1)=  0.0d0
 DFGRD1(2,2)=  1/sqrt(DFGRD1(1,1))
 DFGRD1(2,3)=  0.0d0
 DFGRD1(3,1)=  0.0d0
 DFGRD1(3,2)=  0.0d0
 DFGRD1(3,3)=  1/sqrt(DFGRD1(1,1))
!
time(1)=0.d0
time(2)=0.d0
dtime = 0.1d0
kstep = 1
dtime=0.1d0
call UEXTERNALDB(0,0,time,0.D0,0,0)
!
 CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

 write(*,*) STRESS
 write(*,*)
 write(*,*) DDSDDE
close(150)
!################################################################################################!
!
END PROGRAM
