!>********************************************************************
!> Record of revisions:                                              |
!>        Date        Programmer        Description of change        |
!>        ====        ==========        =====================        |
!>     05/11/2016    Joao Ferreira      full network model           |
!>     05/11/2016    Joao Ferreira      nonaffine network            |
!>     06/11/2016    Joao Ferreira      affine network               |
!>     16/02/2018    Joao Ferreira      compliant cross-linkers      |
!>--------------------------------------------------------------------
!>     Description:
!>     UMAT: USER MATERIAL FOR THE FULL NETWORK MODEL.
!>                AFFINE AND NON-AFFINE DEFORMATIONS
!>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
!>--------------------------------------------------------------------
!>---------------------------------------------------------------------

SUBROUTINE umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, drplde,drpldt,stran,     &
dstran,time,dtime,temp,dtemp,predef,dpred,cmname,ndi,nshr,ntens,nstatev,props,  &
nprops,coords,drot,pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

!----------------------------------------------------------------------
!--------------------------- DECLARATIONS -----------------------------
!----------------------------------------------------------------------
use global
IMPLICIT NONE


INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN OUT)                  :: nshr
INTEGER, INTENT(IN OUT)                  :: ntens
INTEGER, INTENT(IN OUT)                  :: nstatev
INTEGER, INTENT(IN OUT)                  :: nprops
DOUBLE PRECISION, INTENT(IN OUT)         :: pnewdt
INTEGER, INTENT(IN OUT)                  :: noel
INTEGER, INTENT(IN OUT)                  :: npt
INTEGER, INTENT(IN OUT)                  :: kstep
INTEGER, INTENT(IN OUT)                  :: kinc
DOUBLE PRECISION, INTENT(IN)             :: props(nprops)
DOUBLE PRECISION, INTENT(IN OUT)         :: coords(3)
DOUBLE PRECISION, INTENT(OUT)            :: sigma(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: statev(nstatev)
DOUBLE PRECISION, INTENT(OUT)            :: ddsigdde(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd0(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: dfgrd1(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: time(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: predef(1)
!     FILAMENTS DIRECTION
COMMON /kfil/mf0
!     FILAMENTS WEIGHT
COMMON /kfilr/rw
!     PREFERED DIRETION
COMMON /kfilp/prefdir
!
CHARACTER (LEN=8) :: cmname
INTEGER :: layer, kspt
DOUBLE PRECISION :: stress(ntens),  &
    ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),  &
    stran(ntens),dstran(ntens), dpred(1),  drot(3,3)
!
DOUBLE PRECISION :: sse, spd, scd, rpl, drpldt, temp, dtemp, celent
!     FLAGS
!      INTEGER FLAG1
!     UTILITY TENSORS
DOUBLE PRECISION :: unit2(ndi,ndi),unit4(ndi,ndi,ndi,ndi),  &
    unit4s(ndi,ndi,ndi,ndi), proje(ndi,ndi,ndi,ndi),projl(ndi,ndi,ndi,ndi)
!     KINEMATICS
DOUBLE PRECISION :: distgr(ndi,ndi),distgr0(ndi,ndi),  &
    c(ndi,ndi),b(ndi,ndi), cbar(ndi,ndi),bbar(ndi,ndi),distgrinv(ndi,ndi),  &
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
DOUBLE PRECISION :: mf0(nwp,3),rw(nwp),filprops(8), naffprops(2),affprops(2)
DOUBLE PRECISION :: ll,r0,lambda0,mu0,beta,nn,mm,b0,bb,phi,p
DOUBLE PRECISION :: nna,r0f,r0c,etac
DOUBLE PRECISION :: pknetfic(ndi,ndi),cmnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: snetfic(ndi,ndi),cnetfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: pknetficaf(ndi,ndi),pknetficnaf(ndi,ndi)
DOUBLE PRECISION :: snetficaf(ndi,ndi),snetficnaf(ndi,ndi),  &
    snetficafp(ndi,ndi)
DOUBLE PRECISION :: cmnetficaf(ndi,ndi,ndi,ndi), cmnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: cnetficaf(ndi,ndi,ndi,ndi),  &
    cnetficafp(ndi,ndi,ndi,ndi), cnetficnaf(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: efi
INTEGER :: nterm,stat
!     CONTRACTILE FILAMENT
DOUBLE PRECISION :: fric,ffmax,frac0(nch),frac(nch),kch(7),ru0(nwp),  &
    prefdir(nelem,4),varact,sactiso(ndi,ndi),  &
    sactvl(ndi),sactvc(ndi,ndi),dirmax(ndi)
!     VISCOUS PROPERTIES (GENERALIZED MAXWEL DASHPOTS)
DOUBLE PRECISION :: vscprops(6)
INTEGER :: vv
!     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
!      DOUBLE PRECISION CJR(NDI,NDI,NDI,NDI)
!     CAUCHY STRESS AND ELASTICITY TENSOR
DOUBLE PRECISION :: ddpkdde(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: stest(ndi,ndi), ctest(ndi,ndi,ndi,ndi)
INTEGER :: i1,j1, factor
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
c=zero
b=zero
cbar=zero
bbar=zero
ubar=zero
vbar=zero
rot=zero
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
snetficafp=zero
snetficnaf=zero
cmnetfic=zero
cmnetficaf=zero
cmnetficnaf=zero
cnetficaf=zero
cnetficnaf=zero
cnetficafp=zero
!     DIFFUSION
dpdt=zero
dphidmu=zero
dphidotdmu=zero
mob=zero
dmdmu=zero
dmdj=zero
dduddc=one
ddcddu=one
dete=one
dets=one
!     JAUMANN RATE
!      CJR=ZERO
!     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
sigma=zero
ddsigdde=zero

ru0=zero
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
r0f      = props(6)
r0c      = props(7)
etac     = props(8)
mu0      = props(9)
beta     = props(10)
b0       = props(11)
lambda0  = props(12)
filprops = props(5:12)
!     NONAFFINE NETWORK
nn       = props(13)
p        = props(14)
naffprops= props(13:14)
!     AFFINE NETWORK
nna      = props(15)
bb       = props(16)
affprops = props(15:16)

!     NUMERICAL COMPUTATIONS
nterm    = 60

!        STATE VARIABLES: INITIALIZATION/FROM LAST INCREMENT
IF ((kinc <= 1).AND.(kstep == 1)) THEN
  CALL initialize(statev,phi0)  
END IF
!        READ STATEV
CALL sdvread(statev,phit,cr)
!if (npt==1) then
!write(*,*) kinc,statev
!endif
!----------------------------------------------------------------------
!---------------------------- KINEMATICS ------------------------------
!----------------------------------------------------------------------
!     DISTORTION GRADIENT
CALL fslip(dfgrd1,distgr,det,ndi)
!     INVERSE OF DEFORMATION GRADIENT
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
!     ROTATION TENSOR
CALL rotation(distgr,rot,ubar,ndi)
!----------------------------------------------------------------------
!--------------------- CONSTITUTIVE RELATIONS  ------------------------
!----------------------------------------------------------------------
!     DEVIATORIC PROJECTION TENSORS
CALL projeul(unit2,unit4s,proje,ndi)
CALL projlag(c,unit4,projl,ndi)

!---- CHEMICAL POTENTIAL ---------------------------------------------
!---- CHEMICAL INTERACTION ---------------------------------------------
      PHIT=ONE
     call chempot(phit,phitau,dpdt,dphidmu,dphidotdmu,mob,dmdmu,dmdj, &
                   dduddc,ddcddu,vmol,mutau,theta0,   &
                    dete,dets,det,cr,k,dffprops,unit2,dtime,ndi)
      PHITAU=PHIT
!     TIME DERIVATIVE OF THE FRACTION CONCENTRATION (FINITE DIF)
      DPDT= (PHITAU - PHIT)/DTIME
!     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
      DPHIDMU=ZERO
!     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
      DPHIDMU=ZERO
      DDUDDC=zero
      DDCDDU=zero
!      VMOL = ONE
!---- VOLUMETRIC ------------------------------------------------------
!     STRAIN-ENERGY
CALL vol(ssev,pv,ppv,k,det,dete,dets)
!----------------------------------------------------------------------
!---- ISOCHORIC ISOTROPIC ---------------------------------------------
!----------------------------------------------------------------------
IF (phi < one) THEN
!     ISOTROPIC STRAIN-ENERGY
  CALL isomat(sseiso,diso,c10,c01,cbari1,cbari2)
!     PK2 'FICTICIOUS' STRESS TENSOR
  CALL pk2isomatfic(pkmatfic,diso,cbar,cbari1,unit2,ndi)
!     CAUCHY 'FICTICIOUS' STRESS TENSOR
  CALL sigisomatfic(sisomatfic,pkmatfic,distgr,det,ndi)
!     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
  CALL cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2,  &
      diso,unit2,unit4,det,ndi)
!     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
  CALL csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)
  
END IF
!----------------------------------------------------------------------
!---- FILAMENTOUS NETWORK ---------------------------------------------
!----------------------------------------------------------------------
!     IMAGINARY ERROR FUNCTION BASED ON DISPERSION PARAMETER
CALL erfi(efi,bb,nterm)
!----------------------------------------------------------------------
!       ONLY SPATIAL TENSORS ARE CORRECT ...
!------------ NONAFFINE NETWORK --------------
factor = 6 !720 points
IF (nn > zero) THEN
!     'FICTICIOUS' PK2 STRESS AND MATERIAL ELASTICITY TENSORS
!      CALL NAFFMATNETFIC(PKNETFICNAF,CMNETFICNAF,DISTGR,MF0,RW,FILPROPS,
!     1                               NAFFPROPS,DET,NDI)
!     'FICTICIOUS' CAUCHY STRESS AND SPATIAL ELASTICITY TENSORS
!  CALL naffnetfic_bazant(snetficnaf,cnetficnaf,distgr,mf0,rw,filprops,  &
!      naffprops,det,ndi)
   CALL naffnetfic_discrete(snetficnaf,cnetficnaf,distgr,mf0,rw,filprops,  &
      naffprops,det,factor,ndi)
END IF
!------------ AFFINE NETWORK --------------
IF (nna > zero) THEN
!      CALL AFFMATCLNETFIC(PKNETFICAF,CMNETFICAF,DISTGR,MF0,RW,FILPROPS,
!     1                        AFFPROPS,EFI,NOEL,DET,NDI)
 ! CALL affclnetfic_bazant(snetficafp,cnetficafp,distgr,mf0,rw,filprops,  &
 !     affprops,efi,noel,det,ndi)
!      
  CALL affclnetfic_discrete(snetficafp,cnetficafp,distgr,mf0,rw,filprops,  &
      affprops,efi,noel,det,factor,ndi)
END IF
!----------------------------------------------------------------------
!      PKNETFIC=PKNETFICAF+PKNETFICNAF
snetfic=snetficnaf+snetficafp
!      CMNETFIC=CMNETFICAF+CMNETFICNAF
cnetfic=cnetficnaf+cnetficafp
!----------------------------------------------------------------------
!     STRAIN-ENERGY
!      SSE=SSEV+SSEISO
!     PK2 'FICTICIOUS' STRESS
!      PKFIC=(ONE-PHI)*PKMATFIC+PKNETFIC
!     CAUCHY 'FICTICIOUS' STRESS
sfic=(one-phi)*sisomatfic+snetfic
!     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
!      CMFIC=(ONE-PHI)*CMISOMATFIC+CMNETFIC
!     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
cfic=(one-phi)*cisomatfic+cnetfic

!----------------------------------------------------------------------
!-------------------------- STRESS MEASURES ---------------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------
!      PK2 STRESS
!      CALL PK2VOL(PKVOL,PV,C,NDI)
!      CAUCHY STRESS
CALL sigvol(svol,pv,unit2,ndi)

!---- ISOCHORIC -------------------------------------------------------
!      PK2 STRESS
!      CALL PK2ISO(PKISO,PKFIC,PROJL,DETE,NDI)
!      CAUCHY STRESS
CALL sigiso(siso,sfic,proje,ndi)

!---- VOLUMETRIC + ISOCHORIC ------------------------------------------
!      PK2 STRESS
!      PK2   = PKVOL + PKISO
!      CAUCHY STRESS
sigma = svol + siso
!
!----------------------------------------------------------------------
!-------------------- MATERIAL ELASTICITY TENSOR ----------------------
!----------------------------------------------------------------------

!---- VOLUMETRIC ------------------------------------------------------

!      CALL METVOL(CMVOL,C,PV,PPV,DETE,NDI)

!---- ISOCHORIC -------------------------------------------------------

!      CALL METISO(CMISO,CMFIC,PROJL,PKISO,PKFIC,C,UNIT2,DETE,NDI)

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

!      CALL SETJR(CJR,SIGMA,UNIT2,NDI)

!----------------------------------------------------------------------

!     ELASTICITY TENSOR
ddsigdde=cvol+ciso!+CJR

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!------------------------- INDEX ALLOCATION ---------------------------
!----------------------------------------------------------------------
!     VOIGT NOTATION  - FULLY SIMMETRY IMPOSED
!      CALL INDEXX(STRESS,DDSDDE,SIGMA,DDSIGDDE,NTENS,NDI)
!----------------------------------------------------------------------
!--------------------------- STATE VARIABLES --------------------------
!----------------------------------------------------------------------
!     DO K1 = 1, NTENS
!      STATEV(1:27) = VISCOUS TENSORS
CALL sdvwrite(statev,phitau,cr,sigma(1,2))
!     END DO
!----------------------------------------------------------------------
!if(npt == 1) then 
!write(*,*) cbari1, statev
!endif
      
RETURN
END SUBROUTINE material
!----------------------------------------------------------------------
!--------------------------- END OF UMAT ------------------------------
!----------------------------------------------------------------------

