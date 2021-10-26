!>********************************************************************
!> Record of revisions:                                              |
!>        Date        Programmer        Description of change        |
!>        ====        ==========        =====================        |
!>     15/11/2017    Joao Ferreira      cont mech general eqs        |
!>     01/11/2018    Joao Ferreira      comments added               |
!>--------------------------------------------------------------------
!>     Description:
!>     UMAT: IMPLEMENTATION OF THE CONSTITUTIVE EQUATIONS BASED UPON
!>     A STRAIN-ENERGY FUNCTION (SEF).
!>     THIS CODE, AS IS, EXPECTS A SEF BASED ON THE INVARIANTS OF THE
!>     CAUCHY-GREEN TENSORS. A VISCOELASTIC COMPONENT IS ALSO
!>     INCLUDED IF NEEDED.
!>     YOU CAN CHOOSE TO COMPUTE AT THE MATERIAL FRAME AND THEN
!>     PUSHFORWARD OR  COPUTE AND THE SPATIAL FRAME DIRECTLY.
!>--------------------------------------------------------------------
!>     IF YOU WANT TO ADAPT THE CODE ACCORDING TO YOUR SEF:
!>    ISOMAT - DERIVATIVES OF THE SEF IN ORDER TO THE INVARIANTS
!>    ADD OTHER CONTRIBUTIONS: STRESS AND TANGENT MATRIX
!>--------------------------------------------------------------------
!      STATE VARIABLES: CHECK ROUTINES - INITIALIZE, WRITESDV, READSDV.
!>--------------------------------------------------------------------
!>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!>--------------------------------------------------------------------
!>---------------------------------------------------------------------

SUBROUTINE umat(stress,statev,ddsdde,sse,spd,scd, rpl,ddsddt,drplde,drpldt,  &
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,  &
    ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt,  &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

use global 
implicit none
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

DOUBLE PRECISION, INTENT(IN OUT)         :: stress(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: statev(nstatev)
DOUBLE PRECISION, INTENT(IN OUT)         :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(OUT)            :: sse
DOUBLE PRECISION, INTENT(IN OUT)         :: spd
DOUBLE PRECISION, INTENT(IN OUT)         :: scd
DOUBLE PRECISION, INTENT(IN OUT)         :: rpl
DOUBLE PRECISION, INTENT(IN OUT)         :: ddsddt(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: drplde(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: drpldt
DOUBLE PRECISION, INTENT(IN OUT)         :: stran(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: dstran(ntens)
DOUBLE PRECISION, INTENT(IN OUT)         :: time(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: temp
DOUBLE PRECISION, INTENT(IN OUT)         :: dtemp
DOUBLE PRECISION, INTENT(IN OUT)         :: predef(1)
DOUBLE PRECISION, INTENT(IN OUT)         :: dpred(1)
CHARACTER (LEN=8), INTENT(IN OUT)        :: cmname
DOUBLE PRECISION, INTENT(IN)             :: props(nprops)
DOUBLE PRECISION, INTENT(IN OUT)         :: coords(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: drot(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: pnewdt
DOUBLE PRECISION, INTENT(IN OUT)         :: celent
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd0(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd1(3,3)
!
COMMON /kfilp/prefdir
!
DOUBLE PRECISION :: prefdir(nelem,4)
!
INTEGER :: nterm

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
DOUBLE PRECISION :: kbulk,pv,ppv,ssev
!     ISOCHORIC CONTRIBUTION
DOUBLE PRECISION :: siso(ndi,ndi),pkiso(ndi,ndi),pk2(ndi,ndi),  &
    ciso(ndi,ndi,ndi,ndi),cmiso(ndi,ndi,ndi,ndi),  &
    sfic(ndi,ndi),cfic(ndi,ndi,ndi,ndi), pkfic(ndi,ndi),cmfic(ndi,ndi,ndi,ndi)
!     ISOCHORIC ISOTROPIC CONTRIBUTION
DOUBLE PRECISION :: c10,c01,sseiso,diso(5),pkmatfic(ndi,ndi),  &
    smatfic(ndi,ndi),sisomatfic(ndi,ndi), cmisomatfic(ndi,ndi,ndi,ndi),  &
    cisomatfic(ndi,ndi,ndi,ndi)
!     ISOCHORIC ANISOTROPIC CONTRIBUTION
DOUBLE PRECISION :: k1,k2,bdisp,sseaniso, daniso(4),  &
    pkmatficaniso(ndi,ndi), sanisomatfic(ndi,ndi),  &
    cmanisomatfic(ndi,ndi,ndi,ndi), canisomatfic(ndi,ndi,ndi,ndi),  &
    lambda,barlambda, cbari4,efi,lr
DOUBLE PRECISION :: vorif(3),vd(3),m0(3,3),mm(3,3),  &
    vorif2(3),vd2(3),n0(3,3),nn(3,3)
!     LIST VARS OF OTHER CONTRIBUTIONS HERE
INTEGER factor
!     VISCOUS PROPERTIES (GENERALIZED MAXWEL DASHPOTS)
DOUBLE PRECISION :: vscprops(6)
INTEGER :: vv
!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
DOUBLE PRECISION :: cjr(ndi,ndi,ndi,ndi)
!     CAUCHY STRESS AND ELASTICITY TENSOR
DOUBLE PRECISION :: sigma(ndi,ndi),ddsigdde(ndi,ndi,ndi,ndi),  &
    ddpkdde(ndi,ndi,ndi,ndi)
!     MICROSTRUCTURE STATE 
DOUBLE PRECISION :: frac_pld,prefang
!     TESTING/DEBUG VARS
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
kbulk=zero
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
!     INITIALIZE OTHER CONT HERE

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
kbulk    = props(1)
!     ISOCHORIC ISOTROPIC NEO HOOKE
c10      = props(2)
c01      = props(3)
!     ISOCHORIC ANISOTROPIC GHO
k1      = props(4)
k2      = props(5)
bdisp   = props(6)
lr      = props(7)
!
factor  = props(8)
!     VISCOUS EFFECTS: MAXWELL ELEMENTS (MAX:3)
!      VV       = INT(PROPS(7))
!      VSCPROPS = PROPS(8:13)
!     NUMERICAL COMPUTATIONS
nterm    = 60

!     STATE VARIABLES

IF ((time(1) == zero).AND.(kstep == 1)) THEN
  CALL initialize(statev)
END IF
!        READ STATEV
CALL sdvread(statev)

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
!     FIBER UNIT VECTOR AND STRUCTURAL TENSOR
CALL fibdir(prefdir,m0,mm,nelem,noel,ndi,vorif,vd,distgr,dfgrd1)

!     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
CALL invariants(cbar,cbari1,cbari2,ndi)
!
CALL pinvariants(cbar,cbari4,ndi,m0,lambda,barlambda,det)

!     STRETCH TENSORS
CALL stretch(cbar,bbar,ubar,vbar,ndi)
!     ROTATION TENSORS
CALL rotation(distgr,rot,ubar,ndi)
!     DEVIATORIC PROJECTION TENSORS
CALL projeul(unit2,unit4s,proje,ndi)

CALL projlag(c,unit4,projl,ndi)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY AND DERIVATIVES (CHANGE HERE ACCORDING TO YOUR MODEL)
CALL vol(ssev,pv,ppv,kbulk,det)
!---- ISOCHORIC -------------------------------------------------------
CALL isomat(sseiso,diso,c10,c01,cbari1,cbari2)
!CALL anisomat(sseaniso,daniso,diso,k1,k2,kdisp,cbari4,cbari1)
!
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
!     PK2 'FICTICIOUS' STRESS TENSOR
CALL pk2isomatfic(pkmatfic,diso,cbar,cbari1,unit2,ndi)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
CALL sigisomatfic(sisomatfic,pkmatfic,distgr,det,ndi)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
CALL cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2, diso,unit2,unit4,det,ndi)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
CALL csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)
!
!---- FIBERS (ONE FAMILY)   -------------------------------------------
!
!CALL pk2anisomatfic(pkmatficaniso,daniso,cbar,cbari4,m0,ndi)
!CALL push2(sanisomatfic,pkmatficaniso,distgr,det,ndi)
!     IMAGINARY ERROR FUNCTION BASED ON DISPERSION PARAMETER
CALL erfi(efi,bdisp,nterm)
!
! material configuration
CALL manisomat_discrete(sseaniso,pkmatficaniso,cmanisomatfic,distgr,props, &
    efi,noel, npt, kinc, det, factor, prefdir, prefang, frac_pld, ndi )
!
! spatial configuration
CALL anisomat_discrete(sseaniso,sanisomatfic,canisomatfic,distgr,props, &
    efi,noel, npt, kinc, det, factor, prefdir, prefang, frac_pld, ndi )
!CALL cmatanisomatfic(cmanisomatfic,m0,daniso,unit2,det,ndi)
!CALL push4(canisomatfic,cmanisomatfic,distgr,det,ndi)!
!----------------------------------------------------------------------
!     SUM OF ALL ELASTIC CONTRIBUTIONS
!----------------------------------------------------------------------
!     STRAIN-ENERGY
sse=ssev+sseiso+sseaniso
!     PK2   'FICTICIOUS' STRESS
pkfic=pkmatfic+pkmatficaniso
!     CAUCHY 'FICTICIOUS' STRESS
sfic=sisomatfic+sanisomatfic
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
cmfic=cmisomatfic+cmanisomatfic
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
cfic=cisomatfic+canisomatfic
!
!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
CALL pk2vol(pkvol,pv,c,ndi)
!      CAUCHY STRESS
CALL sigvol(svol,pv,unit2,ndi)
!
!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
CALL pk2iso(pkiso,pkfic,projl,det,ndi)
!      CAUCHY STRESS
CALL sigiso(siso,sfic,proje,ndi)
!
!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
pk2   = pkvol + pkiso
!      CAUCHY STRESS
sigma = svol  + siso
!
!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------
!
!---- VOLUMETRIC ------------------------------------------------------
!
CALL metvol(cmvol,c,pv,ppv,det,ndi)
!
!---- ISOCHORIC -------------------------------------------------------

CALL metiso(cmiso,cmfic,projl,pkiso,pkfic,c,unit2,det,ndi)

!----------------------------------------------------------------------

ddpkdde=  cmvol + cmiso

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


!lets test pk2...
! call push4(ctest,ddpkdde,distgr,det,ndi)
! write(*,*) ctest-ciso-cvol
! write(*,*) '*********************************************'
! call push2(stest,pkiso,distgr,det,ndi)
! write(*,*) stest-siso
! write(*,*) '*********************************************'
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
CALL sdvwrite(det,lambda,frac_pld,prefang,statev)
!     END DO
!----------------------------------------------------------------------
RETURN
END SUBROUTINE umat
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------

