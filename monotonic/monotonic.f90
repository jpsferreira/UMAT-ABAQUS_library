!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'ABA_PARAM.INC'
INCLUDE 'PARAM_UMAT.INC'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 13, NDI=3, NSHR=3)
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
!!!
time(1)=ZERO
time(2)=ZERO
call UEXTERNALDB(0,0,time,ZERO,0,0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=2.d0/1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELASTIC PARAM
PROPS(2)=1.00d0  ! C10
PROPS(3)=0.00d0 ! C01
PROPS(4)=0.00d0 !K1
PROPS(5)=0.1d0 !K2
PROPS(6)=0.1d0 !kdisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(7)=1
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
STATEV=ZERO
STATEVP=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
erf=ZERO
RHO=ZERO
!
!################################################################################################!
NSTEPS=400 !NUMBER OF POINTS PER CYCLE
DTIME       = 0.1d0
!
STRETCH_MAX = 2.D0
STRETCH_INITIAL = 0.5d0
DSTRETCH    = (STRETCH_MAX-STRETCH_INITIAL)/NSTEPS
!
GAMMA_MAX = 0.6D0
GAMMA_INITIAL = -.6D0
DGAMMA   = (GAMMA_MAX-GAMMA_INITIAL)/NSTEPS
!################################### |||   NO DEFORMATION  ||| ##################################!
!  DFGRD1(1,1)=  ONE
!  DFGRD1(1,2)=  ZERO
!  DFGRD1(1,3)=  ZERO
!  DFGRD1(2,1)=  ZERO
!  DFGRD1(2,2)=  ONE
!  DFGRD1(2,3)=  ZERO
!  DFGRD1(3,1)=  ZERO
!  DFGRD1(3,2)=  ZERO
!  DFGRD1(3,3)=  ONE
!################################## |||   UNIAXIAL   ||| ################################!
 CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=150, FILE='stress_curves/uniaxial.out', STATUS='UNKNOWN')
        rewind(150)
TIME(1)     = 0.d0
!
STRETCH = STRETCH_INITIAL
DFGRD1(1,1)=STRETCH
DFGRD1(2,2)=ONE/SQRT(STRETCH)
DFGRD1(3,3)=ONE/SQRT(STRETCH)
!
DO KSTEP=1,NSTEPS
!
    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
    write(150,*) TIME(1),STRETCH,STRESS(1)
    TIME(1)     = TIME(1)+DTIME
    STRETCH = STRETCH + DSTRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=ONE/SQRT(STRETCH)
    DFGRD1(3,3)=ONE/SQRT(STRETCH)
!   
ENDDO
! CALL SYSTEM('gnuplot -p data_uniaxial.plt')
close(150)
!################################## |||   BIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=151, FILE='stress_curves/biaxial.out', STATUS='UNKNOWN')
        rewind(151)
TIME(1)     = 0.d0
!
STRETCH = STRETCH_INITIAL
DFGRD1(1,1)=STRETCH
DFGRD1(2,2)=STRETCH
DFGRD1(3,3)=ONE/(STRETCH*STRETCH)
!
DO KSTEP=1,NSTEPS
!
    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
    write(151,*) TIME(1),STRETCH,STRESS(1)
    TIME(1)     = TIME(1)+DTIME
    STRETCH = STRETCH + DSTRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=STRETCH
    DFGRD1(3,3)=ONE/(STRETCH*STRETCH)
!   
ENDDO
close(151)
! 
!################################## |||   SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=152, FILE='stress_curves/shear.out', STATUS='UNKNOWN')
        rewind(152)
TIME(1)     = 0.d0
!
GAMMA = GAMMA_INITIAL
DFGRD1(1,2)=GAMMA
DFGRD1(2,1)=GAMMA
!
DO KSTEP=1,NSTEPS
!
    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
    write(152,*) TIME(1),GAMMA,STRESS(4)
    TIME(1)     = TIME(1)+DTIME
    GAMMA = GAMMA + DGAMMA
    DFGRD1(1,2)=GAMMA
    DFGRD1(2,1)=GAMMA
!   
ENDDO
close(152)
!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=153, FILE='stress_curves/s_shear.out', STATUS='UNKNOWN')
        rewind(153)
TIME(1)     = 0.d0
!
GAMMA = GAMMA_INITIAL
DFGRD1(1,2)=GAMMA
!
DO KSTEP=1,NSTEPS
!
    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
    write(153,*) TIME(1),GAMMA,STRESS(4)
    TIME(1)     = TIME(1)+DTIME
    GAMMA = GAMMA + DGAMMA
    DFGRD1(1,2)=GAMMA
!   
ENDDO
close(153)
! !################################################################################################!
CALL SYSTEM('gnuplot -p data_monotonic.plt')
!!################################################################################################!
END PROGRAM
