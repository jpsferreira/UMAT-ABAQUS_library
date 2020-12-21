!>********************************************************************
!> Record of revisions:                                              |
!>        Date        Programmer        Description of change        |
!>        ====        ==========        =====================        |
!>     05/11/2016    Joao Ferreira      full network model           |
!>--------------------------------------------------------------------
!>     Description:
!>     UMAT: USER MATERIAL FOR THE FULL NETWORK MODEL.
!>                AFFINE DEFORMATIONS
!>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!>--------------------------------------------------------------------
!>---------------------------------------------------------------------

SUBROUTINE umat(stress,statev,ddsdde,sse,spd,scd, rpl,ddsddt,drplde,drpldt,  &
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,  &
    ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt,  &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!
use global  
IMPLICIT NONE
!     FILAMENTS DIRECTION
COMMON /kfil/mf0
!     FILAMENTS WEIGHT
COMMON /kfilr/rw
!     PREFERED DIRETION
COMMON /kfilp/prefdir
!     CHEMICAL DYNAMICS MATRIX
COMMON /kfilf/frac0
COMMON /kfilk/kch
!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
INTEGER, INTENT(IN OUT)                  :: noel
INTEGER, INTENT(IN OUT)                  :: npt
INTEGER, INTENT(IN OUT)                  :: layer
INTEGER, INTENT(IN OUT)                  :: kspt
INTEGER, INTENT(IN OUT)                  :: kstep
INTEGER, INTENT(IN OUT)                  :: kinc
INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN OUT)                  :: nshr
INTEGER, INTENT(IN OUT)                  :: ntens
INTEGER, INTENT(IN OUT)                  :: nstatev
INTEGER, INTENT(IN OUT)                  :: nprops
DOUBLE PRECISION, INTENT(IN OUT)         :: sse
DOUBLE PRECISION, INTENT(IN OUT)         :: spd
DOUBLE PRECISION, INTENT(IN OUT)         :: scd
DOUBLE PRECISION, INTENT(IN OUT)         :: rpl
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: drpldt
DOUBLE PRECISION, INTENT(IN OUT)         :: temp
DOUBLE PRECISION, INTENT(IN OUT)         :: dtemp
CHARACTER (LEN=8), INTENT(IN OUT)        :: cmname
DOUBLE PRECISION, INTENT(IN OUT)         :: pnewdt
DOUBLE PRECISION, INTENT(IN OUT)         :: celent

DOUBLE PRECISION, INTENT(IN OUT)         :: stress(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: statev(nstatev)
DOUBLE PRECISION, INTENT(IN OUT)         :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: ddsddt(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: drplde(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: stran(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: dstran(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: time(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: predef(1)
DOUBLE PRECISION, INTENT(IN OUT)         :: dpred(1)
DOUBLE PRECISION, INTENT(IN)             :: props(nprops)
DOUBLE PRECISION, INTENT(IN OUT)         :: coords(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: drot(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd0(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd1(3,3)

!
!     FLAGS
!      INTEGER FLAG1
!     UTILITY TENSORS
DOUBLE PRECISION :: unit2(ndi,ndi),unit4(ndi,ndi,ndi,ndi),  &
    unit4s(ndi,ndi,ndi,ndi), proje(ndi,ndi,ndi,ndi),projl(ndi,ndi,ndi,ndi)
!     KINEMATICS
DOUBLE PRECISION :: distgr(ndi,ndi),c(ndi,ndi),b(ndi,ndi),  &
    cbar(ndi,ndi),bbar(ndi,ndi),distgrinv(ndi,ndi),  &
    ubar(ndi,ndi),vbar(ndi,ndi),rot(ndi,ndi), dfgrd1inv(ndi,ndi)
DOUBLE PRECISION :: det,cbari1,cbari2
!     VOLUMETRIC CONTRIBUTION
DOUBLE PRECISION :: pkvol(ndi,ndi),svol(ndi,ndi),  &
    cvol(ndi,ndi,ndi,ndi),cmvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: k,pv,ppv,ssev
!     ISOCHORIC CONTRIBUTION
DOUBLE PRECISION :: siso(ndi,ndi),pkiso(ndi,ndi),pk2(ndi,ndi),  &
    ciso(ndi,ndi,ndi,ndi),cmiso(ndi,ndi,ndi,ndi),  &
    sfic(ndi,ndi),cfic(ndi,ndi,ndi,ndi), pkfic(ndi,ndi),cmfic(ndi,ndi,ndi,ndi)
!     ISOCHORIC ISOTROPIC CONTRIBUTION
DOUBLE PRECISION :: c10,c01,sseiso,diso(5),pkmatfic(ndi,ndi),  &
    smatfic(ndi,ndi),sisomatfic(ndi,ndi), cmisomatfic(ndi,ndi,ndi,ndi),  &
    cisomatfic(ndi,ndi,ndi,ndi)
!     FILAMENTS NETWORK CONTRIBUTION
DOUBLE PRECISION :: mf0(nwp,3),rw(nwp),filprops(6), affprops(2)
DOUBLE PRECISION :: ll,lambda0,mu0,beta,nn,mm,b0,bb
DOUBLE PRECISION :: phi,p,r0
DOUBLE PRECISION :: pknetfic(ndi,ndi),cmnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: snetfic(ndi,ndi),cnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: pknetficaf(ndi,ndi),pknetficnaf(ndi,ndi)
DOUBLE PRECISION :: snetficaf(ndi,ndi),snetficnaf(ndi,ndi)
DOUBLE PRECISION :: cmnetficaf(ndi,ndi,ndi,ndi), cmnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: cnetficaf(ndi,ndi,ndi,ndi), cnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: efi
INTEGER :: nterm,factor
!     CONTRACTILE FILAMENT
DOUBLE PRECISION :: fric,ffmax,frac0(4),frac(4),kch(7),ru0(720),  &
    prefdir(nelem,4),varact,dirmax(ndi)

!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
DOUBLE PRECISION :: cjr(ndi,ndi,ndi,ndi)
!     CAUCHY STRESS AND ELASTICITY TENSOR
DOUBLE PRECISION :: sigma(ndi,ndi),ddsigdde(ndi,ndi,ndi,ndi),  &
    ddpkdde(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: stest(ndi,ndi), ctest(ndi,ndi,ndi,ndi)
!----------------------------------------------------------------------
!-------------------------- INITIALIZATIONS ---------------------------
!----------------------------------------------------------------------
!     IDENTITY AND PROJECTION TENSORS
unit2=zero
unit4=zero
unit4s=zero
proje=zero
projl=zero
!     KINEMATICS
distgr=zero
c=zero
b=zero
cbar=zero
bbar=zero
ubar=zero
vbar=zero
rot=zero
det=zero
cbari1=zero
cbari2=zero
!     VOLUMETRIC
pkvol=zero
svol=zero
cvol=zero
k=zero
pv=zero
ppv=zero
ssev=zero
!     ISOCHORIC
siso=zero
pkiso=zero
pk2=zero
ciso=zero
cfic=zero
sfic=zero
pkfic=zero
!     ISOTROPIC
c10=zero
c01=zero
sseiso=zero
diso=zero
pkmatfic=zero
smatfic=zero
sisomatfic=zero
cmisomatfic=zero
cisomatfic=zero
!     FILAMENTS NETWORK
snetfic=zero
cnetfic=zero
pknetfic=zero
pknetficaf=zero
pknetficnaf=zero
snetficaf=zero
snetficnaf=zero
cmnetfic=zero
cmnetficaf=zero
cmnetficnaf=zero
cnetficaf=zero
cnetficnaf=zero

ru0=zero
!     JAUMANN RATE
cjr=zero
!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
sigma=zero
ddsigdde=zero
!----------------------------------------------------------------------
!------------------------ IDENTITY TENSORS ----------------------------
!----------------------------------------------------------------------
CALL onem(unit2,unit4,unit4s,ndi)
!----------------------------------------------------------------------
!------------------- MATERIAL CONSTANTS AND DATA ----------------------
!----------------------------------------------------------------------
!     VOLUMETRIC
k        = props(1)
!     ISOCHORIC ISOTROPIC
c10      = props(2)
c01      = props(3)
phi      = props(4)
!     FILAMENT
ll       = props(5)
r0      = props(6)
mu0      = props(7)
beta     = props(8)
b0       = props(9)
lambda0  = props(10)
filprops = props(5:10)
!     NONAFFINE NETWORK
nn       = props(11)
bb        = props(12)
affprops= props(11:12)


!     NUMERICAL COMPUTATIONS
nterm    = 60

!        STATE VARIABLES AND CHEMICAL PARAMETERS

IF ((time(1) == zero).AND.(kstep == 1)) THEN
  CALL initialize(statev)
END IF
!        READ STATEV
CALL sdvread(frac0,ru0,statev)
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
CALL fslip(dfgrd1,distgr,det,ndi)
!     INVERSE OF DISTORTION GRADIENT
CALL matinv3d(dfgrd1,dfgrd1inv,ndi)
!     INVERSE OF DISTORTION GRADIENT
CALL matinv3d(distgr,distgrinv,ndi)
!     CAUCHY-GREEN DEFORMATION TENSORS
CALL deformation(dfgrd1,c,b,ndi)
CALL deformation(distgr,cbar,bbar,ndi)
!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
CALL invariants(cbar,cbari1,cbari2,ndi)
!     STRETCH TENSORS
CALL stretch(cbar,bbar,ubar,vbar,ndi)
!     ROTATION TENSORS
CALL rotation(distgr,rot,ubar,ndi)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!     DEVIATORIC PROJECTION TENSORS
CALL projeul(unit2,unit4s,proje,ndi)

CALL projlag(c,unit4,projl,ndi)

!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY
CALL vol(ssev,pv,ppv,k,det)

!---- ISOCHORIC ISOTROPIC ---------------------------------------------
IF (phi < one) THEN
!     STRAIN-ENERGY
  CALL isomat(sseiso,diso,c10,c01,cbari1,cbari2)
!     PK2 'FICTICIOUS' STRESS TENSOR
  CALL pk2isomatfic(pkmatfic,diso,cbar,cbari1,unit2,ndi)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
  CALL sigisomatfic(sisomatfic,pkmatfic,distgr,det,ndi)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
  CALL cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2, diso,unit2,unit4,det,ndi)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
  CALL csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)
  
END IF
factor = 6
!---- FILAMENTS NETWORK -----------------------------------------------
!     IMAGINARY ERROR FUNCTION BASED ON DISPERSION PARAMETER
CALL erfi(efi,bb,nterm)
!     'FICTICIOUS' PK2 STRESS AND MATERIAL ELASTICITY TENSORS
!------------ AFFINE NETWORK --------------
IF (nn > zero) THEN
  CALL affnetfic_discrete(snetficaf,cnetficaf,distgr,filprops,  &
      affprops,efi,noel,det,factor,ndi)
END IF

!      PKNETFIC=PKNETFICNAF+PKNETFICAF
snetfic=snetficnaf+snetficaf
!      CMNETFIC=CMNETFICNAF+CMNETFICAF
cnetfic=cnetficnaf+cnetficaf
!----------------------------------------------------------------------
!     STRAIN-ENERGY
!      SSE=SSEV+SSEISO
!     PK2 'FICTICIOUS' STRESS
pkfic=(one-phi)*pkmatfic+pknetfic
!     CAUCHY 'FICTICIOUS' STRESS
sfic=(one-phi)*sisomatfic+snetfic
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
cmfic=(one-phi)*cmisomatfic+cmnetfic
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
cfic=(one-phi)*cisomatfic+cnetfic

!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
CALL pk2vol(pkvol,pv,c,ndi)
!      CAUCHY STRESS
CALL sigvol(svol,pv,unit2,ndi)

!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
CALL pk2iso(pkiso,pkfic,projl,det,ndi)
!      CAUCHY STRESS
CALL sigiso(siso,sfic,proje,ndi)
!      ACTIVE CAUCHY STRESS
!      CALL SIGISO(SACTISO,SNETFICAF,PROJE,NDI)

!      CALL SPECTRAL(SACTISO,SACTVL,SACTVC)

!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
pk2 = pkvol + pkiso
!      CAUCHY STRESS
sigma = svol + siso

!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------

!      CALL METVOL(CMVOL,C,PV,PPV,DET,NDI)

!---- ISOCHORIC -------------------------------------------------------

!      CALL METISO(CMISO,CMFIC,PROJL,PKISO,PKFIC,C,UNIT2,DET,NDI)

!----------------------------------------------------------------------

!      DDPKDDE=CMVOL+CMISO

!----------------------------------------------------------------------
!--------------------- SPATIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------

CALL setvol(cvol,pv,ppv,unit2,unit4s,ndi)

!---- ISOCHORIC -------------------------------------------------------

CALL setiso(ciso,cfic,proje,siso,sfic,unit2,ndi)

!-----JAUMMAN RATE ----------------------------------------------------

CALL setjr(cjr,sigma,unit2,ndi)

!----------------------------------------------------------------------

!     ELASTICITY TENSOR
ddsigdde=cvol+ciso+cjr


!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SIMMETRY IMPOSED
CALL indexx(stress,ddsdde,sigma,ddsigdde,ntens,ndi)

!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
CALL sdvwrite(frac,ru0,det,varact,dirmax,statev)
!     END DO
!----------------------------------------------------------------------
RETURN
END SUBROUTINE umat
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------- AUXILIAR SUBROUTINES -------------------------
!----------------------------------------------------------------------
!                         INPUT FILES
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!                         KINEMATIC QUANTITIES
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                         STRESS TENSORS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                   LINEARISED ELASTICITY TENSORS
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------- UTILITY SUBROUTINES --------------------------
!----------------------------------------------------------------------

