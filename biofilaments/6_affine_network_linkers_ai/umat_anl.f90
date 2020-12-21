SUBROUTINE getoutdir(outdir, lenoutdir)



!>     GET CURRENT WORKING DIRECTORY
INCLUDE 'aba_param.inc'


CHARACTER (LEN=256), INTENT(IN OUT)      :: outdir
INTEGER, INTENT(OUT)                     :: lenoutdir



CALL getcwd(outdir)
!        OUTDIR=OUTDIR(1:SCAN(OUTDIR,'\',BACK=.TRUE.)-1)
lenoutdir=len_trim(outdir)

RETURN
END SUBROUTINE getoutdir
module global

! This module is used to transfer SDV's from the UEL
!  to the UVARM so that SDV's can be visualized on a
!  dummy mesh
!
!  globalSdv(X,Y,Z)
!   X - element pointer
!   Y - integration point pointer
!   Z - SDV pointer
!
!  numElem
!   Total number of elements in the real mesh, the dummy
!   mesh needs to have the same number of elements, and
!   the dummy mesh needs to have the same number of integ
!   points.  You must set that parameter value here.
!
!  ElemOffset
!   Offset between element numbers on the real mesh and
!    dummy mesh.  That is set in the input file, and
!    that value must be set here the same.

INTEGER :: numelem,elemoffset,ERR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the number of UEL elements used here
PARAMETER(numelem=1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the offset here for UVARM plotting, must match input file!
PARAMETER(elemoffset=1000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8), allocatable :: globalsdv(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the control parameters to run the material-related routines

INTEGER NWP,NELEM,NCH,NCNT,NSDV,NVSC,NMECH,NQUEM
PARAMETER (NELEM=1)
PARAMETER (NWP=21,NVSC=27,NMECH=1,NCH=4,NCNT=4,NQUEM=5)
!      PARAMETER (NSDV=NWP+NCH+NVSC+NMECH+NCNT+NQUEM)
PARAMETER (NSDV=3)
DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
PARAMETER (ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
PARAMETER (THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
DOUBLE PRECISION HALF,THIRD
PARAMETER (HALF=0.5d0,THIRD=1.d0/3.d0)
CHARACTER(256) DIR1,DIR2,DIR3
PARAMETER (DIR1='sphere_int21c.inp',DIR2='prefdir.inp')
PARAMETER (DIR3='initial_cond.inp')


END module global
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

SUBROUTINE affclnetfic_bazant(sfic,cfic,f,mf0,rw,filprops,affprops,  &
        efi,noel,det,ndi)

!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> BAZANT SPHERICAL INTEGRATION SCHEME
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,k1,l1,m1, im1
DOUBLE PRECISION :: sfilfic(ndi,ndi), cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,r0f,r0,mu0,b0,beta,lambda0,rho,n,fi,ffi,dtime
DOUBLE PRECISION :: r0c,etac,lambdaif
DOUBLE PRECISION :: b,fric,ffmax,ang, frac(4),ru0(nwp),ru
DOUBLE PRECISION :: vara,avga,maxa,aux0,ffic,suma,rho0,dirmax(ndi)
DOUBLE PRECISION :: aa,rr
INTEGER          :: nn,factor
real ( kind = 8 ), external :: polyterm_value_3d
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
n       = affprops(1)
b       = affprops(2)

pi=four*ATAN(one)
aux=n*(det**(-one))*four*pi
cfic=zero
sfic=zero

rho=one
r0=r0f+r0c

aa = zero
avga=zero
maxa=zero
suma=zero
dirmax=zero
!       CALL DENSITY(RHO0,ZERO,B,EFI)

!             OPEN (UNIT=20,FILE="projfil.out",action="write",
!     1 status="replace")

!        LOOP OVER THE INTEGRATION DIRECTIONS (BAZANT INTEGRATION SCHEME)
!bazant integration scheme
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  rho=one
  
  IF((etac > zero).AND.(etac < one))THEN
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaic=zero
  END IF
  
  CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,rwi,ndi)
  
  
  DO j1=1,ndi
    DO k1=1,ndi
      sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
  !bazant integration scheme
 ! aa = aa +rwi*dwi
END DO !end of discretization scheme

!!discrete angular integration test
!factor = 4
!nn = 0
!rr = zero
!call sphere01_quad_icos1c ( factor, polyterm_value_3d, nn, rr )
!write(*,*) 'bazant', four*pi*aa

RETURN
END SUBROUTINE affclnetfic_bazant
SUBROUTINE affclnetfic_discrete(sfic,cfic,f,mf0,rw,filprops,affprops,  &
        efi,noel,det,factor,ndi)

!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
!> DISCRETE ANGULAR INTEGRATION SCHEME (icosahedron)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det

INTEGER :: i1,j1,k1,l1,m1, im1
DOUBLE PRECISION :: sfilfic(ndi,ndi), cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,r0f,r0,mu0,b0,beta,lambda0,rho,n,fi,ffi,dtime
DOUBLE PRECISION :: r0c,etac,lambdaif
DOUBLE PRECISION :: bdisp,fric,ffmax,ang, frac(4),ru0(nwp),ru
DOUBLE PRECISION :: vara,avga,maxa,aux0,ffic,suma,rho0,dirmax(ndi)

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
  

!! initialize the model data
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
n       = affprops(1)
bdisp   = affprops(2)

  aux=n*(det**(-one))
  cfic=zero
  sfic=zero

  rho=one
  r0=r0f+r0c

  aa = zero
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
        CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
        CALL bangle(ang,f,mfi,noel,ndi)
  
        CALL density(rho,ang,bdisp,efi)
        rho=one
  
        IF((etac > zero).AND.(etac < one))THEN
            lambdaif=etac*(r0/r0f)*(lambdai-one)+one
            lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
        ELSE
            lambdaif=lambdai
            lambdaic=zero
        END IF
  
        CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
        
        CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,ai,ndi)
  
        CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,ai,ndi)
  
  
        DO j1=1,ndi
           DO k1=1,ndi
              sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
              DO l1=1,ndi
                DO m1=1,ndi
                  cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
                END DO
              END DO
           END DO
         END DO
        
        !v=dwi
        node_num = node_num + 1
        !rr = rr + ai * v
        !area_total = area_total + ai

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
  
        CALL bangle(ang,f,mfi,noel,ndi)
  
        CALL density(rho,ang,bdisp,efi)
        rho=one
  
        IF((etac > zero).AND.(etac < one))THEN
            lambdaif=etac*(r0/r0f)*(lambdai-one)+one
            lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
        ELSE
            lambdaif=lambdai
            lambdaic=zero
        END IF
  
        CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
        
        CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,ai,ndi)
  
        CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,ai,ndi)
  
  
        DO j1=1,ndi
           DO k1=1,ndi
              sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
              DO l1=1,ndi
                DO m1=1,ndi
                  cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
                END DO
              END DO
           END DO
         END DO
        
        
        !v=dwi
        node_num = node_num + 1  
        !rr = rr + ai * v
        !area_total = area_total + ai

      end do
    end do

  end do
!
!  Discard allocated memory.
!
  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )


RETURN
END SUBROUTINE affclnetfic_discrete
SUBROUTINE affmatclnetfic(pkfic,cfic,f,mf0,rw,filprops,affprops,  &
        efi,noel,det,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' PK STRESS AND ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: sfilfic(ndi,ndi),   &
     cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,fi,ffi,r0
DOUBLE PRECISION :: r0c,etac,lambdaif,lambdaic
DOUBLE PRECISION :: l,r0f,mu0,b0,beta,lambda0,m,rho
DOUBLE PRECISION :: b,fric,ffmax,ang

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
m       = affprops(1)
b       = affprops(2)

pi=four*ATAN(one)
aux=m*(det**(-one))*four*pi
cfic=zero
pkfic=zero
r0=r0f+r0c

rho=one

!        LOOP OVER THE INTEGRATION DIRECTIONS
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
!       CALL FIL(FI,FFI,DWI,DDWI,LAMBDAI,LAMBDA0,L,R0,MU0,BETA,B0)
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  rho = one
  
  IF((etac > zero).AND.(etac < one))THEN
    
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaic=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaic=zero
  END IF
  
  CALL fil(fi,ffi,dwi,ddwi,lambdaif,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdai,dwi,mf0i,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdai,dwi,ddwi,mf0i,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      pkfic(j1,k1)=pkfic(j1,k1)+aux*sfilfic(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
END DO

RETURN
END SUBROUTINE affmatclnetfic
SUBROUTINE affmatnetfic(pkfic,cmfic,f,mf0,rw,filprops,affprops,  &
        ru0,dtime,frac,efi,noel,det,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' PK STRESS AND ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(6)
DOUBLE PRECISION, INTENT(IN)             :: affprops(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: ru0
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: frac(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: sfilfic(ndi,ndi),   &
     cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda0,m,rho
DOUBLE PRECISION :: b,fric,ffmax,ang

!     FILAMENT
l       = filprops(1)
r0      = filprops(2)
mu0     = filprops(3)
beta    = filprops(4)
b0      = filprops(5)
lambda0 = filprops(6)
!     NETWORK
m       = affprops(1)
b       = affprops(2)
fric    = affprops(3)
ffmax   = affprops(4)

pi=four*ATAN(one)
aux=m*(det**(-one))*four*pi
cmfic=zero
pkfic=zero

rho=one

!        LOOP OVER THE INTEGRATION DIRECTIONS
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  CALL fil(fi,ffi,dwi,ddwi,lambdai,lambda0,l,r0,mu0,beta,b0)
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  
  CALL sigfilfic(sfilfic,rho,lambdai,dwi,mf0i,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdai,dwi,ddwi,mf0i,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      pkfic(j1,k1)=pkfic(j1,k1)+aux*sfilfic(j1,k1)
      DO l1=1,ndi
        DO m1=1,ndi
          cmfic(j1,k1,l1,m1)=cmfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
        END DO
      END DO
    END DO
  END DO
  
END DO

RETURN
END SUBROUTINE affmatnetfic
SUBROUTINE affnetfic(sfic,cfic,f,mf0,rw,filprops,affprops,ru0,  &
        dtime,frac,efi,noel,vara,dirmax,det,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: affprops(5)
DOUBLE PRECISION, INTENT(IN OUT)         :: ru0(nwp)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: frac(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: efi
INTEGER, INTENT(IN OUT)                  :: noel
DOUBLE PRECISION, INTENT(OUT)            :: vara
DOUBLE PRECISION, INTENT(OUT)            :: dirmax(ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1, im1
DOUBLE PRECISION :: sfilfic(ndi,ndi),   &
     cfilfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,lambdaic
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda0,rho,m,fi,ffi
DOUBLE PRECISION :: r0f,r0c,etac,lambdaif,lambdaicl
DOUBLE PRECISION :: b,fric,ffmax,ang, ru
DOUBLE PRECISION :: avga,maxa,aux0,ffic,suma,rho0

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
m       = affprops(2)
b       = affprops(3)
fric    = affprops(4)
ffmax   = affprops(5)

pi=four*ATAN(one)
aux=m*(det**(-one))*four*pi
cfic=zero
sfic=zero

rho=one
r0=r0f+r0c

avga=zero
maxa=zero
suma=zero
dirmax=zero
!       CALL DENSITY(RHO0,ZERO,B,EFI)

!             OPEN (UNIT=20,FILE="projfil.out",action="write",
!     1 status="replace")

!        LOOP OVER THE INTEGRATION DIRECTIONS
DO i1=1,nwp
  
  mfi=zero
  mf0i=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  IF((etac > zero).AND.(etac < one))THEN
    
    lambdaif=etac*(r0/r0f)*(lambdai-one)+one
    lambdaicl=(lambdai*r0-lambdaif*r0f)/r0c
  ELSE
    lambdaif=lambdai
    lambdaicl=zero
  END IF
  
  ru=ru0(i1)
  CALL contractile(fi,ffi,dwi,ddwi,ffic,ru,ru,lambdaic,lambdaif,  &
      lambda0,l,r0,mu0,beta,b0,ffmax,fric,frac,dtime)
  ru0(i1)=ru
  
  CALL bangle(ang,f,mfi,noel,ndi)
  
  CALL density(rho,ang,b,efi)
  
!        AUX0=(FFIC/FFMAX)*(RHO)
  aux0=ru*rho
  
  IF (ru > zero) THEN
!        AVERAGE CONTRACTION LEVEL
    avga=avga+aux0
    suma=suma+one
    IF (aux0 > maxa) THEN
!        MAXIMUM CONTRACTION LEVEL
      maxa = aux0
      dirmax=mfi
      im1=i1
    END IF
  END IF
  
  CALL sigfilfic(sfilfic,rho,lambdaif,dwi,mfi,rwi,ndi)
  
  CALL csfilfic(cfilfic,rho,lambdaif,dwi,ddwi,mfi,rwi,ndi)
  
  
  IF(ru > zero)THEN
    
    DO j1=1,ndi
      DO k1=1,ndi
        sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
        DO l1=1,ndi
          DO m1=1,ndi
            cfic(j1,k1,l1,m1)=cfic(j1,k1,l1,m1)+aux*cfilfic(j1,k1,l1,m1)
          END DO
        END DO
      END DO
    END DO
    
  END IF
  
END DO
!      close(20)

IF (suma > zero) THEN
  avga=avga/nwp
END IF
vara=(maxa-avga)*((maxa)**(-one))
!        WRITE(*,*) VARA,MAXA,IM1,SUMA

RETURN
END SUBROUTINE affnetfic
SUBROUTINE assembleelement(ndim,nnode,ndofel,  &
        ru,rc,kuu,kuc,kcu,kcc,  &
        rhs,amatrx)

!
! Subroutine to assemble the local elements residual and tangent
!
use global

INTEGER, INTENT(IN)                      :: ndim
INTEGER, INTENT(IN)                      :: nnode
INTEGER, INTENT(IN)                      :: ndofel
real(8), INTENT(IN)                       :: ru(ndim*nnode,1)
real(8), INTENT(IN)                       :: rc(nnode,1)
real(8), INTENT(IN)                       :: kuu(ndim*nnode,ndim*nnode)
real(8), INTENT(IN)                       :: kuc(ndim*nnode,nnode)
real(8), INTENT(IN)                       :: kcu(nnode,ndim*nnode)
real(8), INTENT(IN)                       :: kcc(nnode,nnode)
real(8), INTENT(OUT)                      :: rhs(ndofel,1)
real(8), INTENT(OUT)                      :: amatrx(ndofel,ndofel)


INTEGER :: i,j,k,l,m,n,a11,a12,b11,b12, ndofn

! Total number of degrees of freedom per node
!
ndofn = ndofel/nnode
! init
!
rhs = 0.d0
amatrx = 0.d0
IF(ndim == 2) THEN
!
! Assemble the element level residual
!
  DO i=1,nnode
    a11 = ndofn*(i-1)+1
    a12 = ndim*(i-1)+1
!
! displacement
!
    rhs(a11,1) = ru(a12,1)
    rhs(a11+1,1) = ru(a12+1,1)
!
! chemical potential
!
    rhs(a11+2,1) = rc(i,1)
!            ! temperature
!            !
!            rhs(A11+3,1) = Rt(i,1)
  END DO
!
! Assemble the element level tangent matrix
!
  DO i=1,nnode
    DO j=1,nnode
      a11 = ndofn*(i-1)+1
      a12 = ndim*(i-1)+1
      b11 = ndofn*(j-1)+1
      b12 = ndim*(j-1)+1
!
! displacement
!
      amatrx(a11,b11) = kuu(a12,b12)
      amatrx(a11,b11+1) = kuu(a12,b12+1)
      amatrx(a11+1,b11) = kuu(a12+1,b12)
      amatrx(a11+1,b11+1) = kuu(a12+1,b12+1)
!
! chemical potential
!
      amatrx(a11+2,b11+2) = kcc(i,j)
!               ! temperature
!               !
!               amatrx(A11+3,B11+3) = Ktt(i,j)
!
! displacement - chemical potential
!
      amatrx(a11,b11+2) = kuc(a12,j)
      amatrx(a11+1,b11+2) = kuc(a12+1,j)
!               ! displacement - temperature
!               !
!               amatrx(A11,B11+3) = Kut(A12,j)
!               amatrx(A11+1,B11+3) = Kut(A12+1,j)
!
! chemical potential - displacement
!
      amatrx(a11+2,b11) = kcu(i,b12)
      amatrx(a11+2,b11+1) = kcu(i,b12+1)
!
!               ! chemical potential - temperature
!               !
!               amatrx(A11+2,B11+3) = Kct(i,j)
!               !
!               ! temperature - displacement
!               !
!               amatrx(A11+3,B11) = Ktu(i,B12)
!               amatrx(A11+3,B11+1) = Ktu(i,B12+1)
!               !
!               ! temperature - chemical potential
!               !
!               amatrx(A11+3,B11+2) = Ktc(i,j)
    END DO
  END DO
!
ELSE IF(ndim == 3) THEN
!
! Assemble the element level residual
!
  DO i=1,nnode
    a11 = ndofn*(i-1)+1
    a12 = ndim*(i-1)+1
!
! displacement
!
    rhs(a11,1)   = ru(a12,1)
    rhs(a11+1,1) = ru(a12+1,1)
    rhs(a11+2,1) = ru(a12+2,1)
!
! chemical potential
!
    rhs(a11+3,1) = rc(i,1)
!
!             ! temperature
!            !
!            rhs(A11+4,1) = Rt(i,1)
  END DO
!
! Assembly the element level tangent matrix
!
  DO i=1,nnode
    DO j=1,nnode
      a11 = ndofn*(i-1)+1
      a12 = ndim*(i-1)+1
      b11 = ndofn*(j-1)+1
      b12 = ndim*(j-1)+1
!
! displacement
!
      amatrx(a11,b11)     = kuu(a12,b12)
      amatrx(a11,b11+1)   = kuu(a12,b12+1)
      amatrx(a11,b11+2)   = kuu(a12,b12+2)
      amatrx(a11+1,b11)   = kuu(a12+1,b12)
      amatrx(a11+1,b11+1) = kuu(a12+1,b12+1)
      amatrx(a11+1,b11+2) = kuu(a12+1,b12+2)
      amatrx(a11+2,b11)   = kuu(a12+2,b12)
      amatrx(a11+2,b11+1) = kuu(a12+2,b12+1)
      amatrx(a11+2,b11+2) = kuu(a12+2,b12+2)
!
! chemical potential
!
      amatrx(a11+3,b11+3) = kcc(i,j)
!                ! temperature
!               !
!               amatrx(A11+4,B11+4) = Ktt(i,j)
!
! displacement - chemical potential
!
      amatrx(a11,b11+3) = kuc(a12,j)
      amatrx(a11+1,b11+3) = kuc(a12+1,j)
      amatrx(a11+2,b11+3) = kuc(a12+2,j)
      
!               ! displacement - temperature
!               !
!               amatrx(A11,B11+4) = Kut(A12,j)
!               amatrx(A11+1,B11+4) = Kut(A12+1,j)
!               amatrx(A11+2,B11+4) = Kut(A12+2,j)
!
! chemical potential - displacement
!
      amatrx(a11+3,b11) = kcu(i,b12)
      amatrx(a11+3,b11+1) = kcu(i,b12+1)
      amatrx(a11+3,b11+2) = kcu(i,b12+2)
!
!
!               ! chemical potential - temperature
!               !
!               amatrx(A11+3,B11+4) = Kct(i,j)
!               !
!               ! temperature - displacement
!               !
!               amatrx(A11+4,B11) = Ktu(i,B12)
!               amatrx(A11+4,B11+1) = Ktu(i,B12+1)
!               amatrx(A11+4,B11+2) = Ktu(i,B12+2)
!               !
!               ! temperature - chemical potential
!               !
!               amatrx(A11+4,B11+3) = Ktc(i,j)
    END DO
  END DO
!
ELSE
  WRITE(*,*) 'How did you get nDim=',ndim
  CALL xit()
END IF

RETURN
END SUBROUTINE assembleelement
SUBROUTINE bangle(ang,f,mf,noel,ndi)



!>    ANGLE BETWEEN FILAMENT AND PREFERED DIRECTION

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ang
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf(ndi)
INTEGER, INTENT(IN OUT)                  :: noel



COMMON /kfilp/prefdir
DOUBLE PRECISION :: prefdir(nelem,4)

INTEGER :: inoel,i,j
DOUBLE PRECISION :: dnorm,pdir(ndi), mfa(ndi),aux
DOUBLE PRECISION :: c(ndi,ndi),egvc(ndi,ndi),egvl(ndi)

inoel=0
i=0
DO i=1,nelem
!               ELEMENT IDENTIFICATION
  IF(noel == INT(prefdir(i,1))) THEN
    inoel=i
  END IF
END DO

DO i=1,ndi
  j=i+1
!       PREFERED ORIENTATION  ORIENTATION NORMALIZED
  pdir(i)=prefdir(inoel,j)
END DO
!        ALTERNATIVE APPROACH: BUNDLES FOLLOW PRINCIPAL DIRECTIONS
c=matmul(transpose(f),f)
CALL spectral(c,egvl,egvc)
!       WRITE(*,*) EGVC
pdir(1)=egvc(1,1)
pdir(2)=egvc(2,1)
pdir(3)=egvc(3,1)
!        END OF ALTERNATIVE

!     PREFERED ORIENTATION
dnorm=dot_product(pdir,pdir)
dnorm=DSQRT(dnorm)

!       PREFERED ORIENTATION  NORMALIZED
pdir=pdir/dnorm

!       FILAMENT ORIENTATION
mfa=mf
dnorm=dot_product(mfa,mfa)
dnorm=DSQRT(dnorm)

!       FILAMENT ORIENTATION  NORMALIZED
mfa=mfa/dnorm
!        ANGLE BETWEEN PREFERED ORIENTATION AND FILAMENT - BANGLE
aux=dot_product(mfa,pdir)

!        if AUX.GT.ONE
!        endif
!        write(*,*) aux
ang=ACOS(aux)

RETURN
END SUBROUTINE bangle

SUBROUTINE chemicalstat(frac,frac0,k,dtime)


use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: frac(4)
DOUBLE PRECISION, INTENT(IN)             :: frac0(4)
DOUBLE PRECISION, INTENT(IN)             :: k(7)
DOUBLE PRECISION, INTENT(IN)             :: dtime


DOUBLE PRECISION :: stiff(4,4)

DOUBLE PRECISION :: aux
INTEGER :: i1,j1

frac=zero
stiff=zero
stiff(1,1)=-k(1)
stiff(1,2)=k(2)
stiff(1,4)=k(7)
stiff(2,1)=k(1)
stiff(2,2)=-k(2)-k(3)
stiff(2,3)=k(4)
stiff(3,2)=k(3)
stiff(3,3)=-k(4)-k(5)
stiff(3,4)=k(6)
stiff(4,3)=k(5)
stiff(4,4)=-k(6)-k(7)

DO i1=1,4
  aux=zero
  DO j1=1,4
    aux=aux+stiff(i1,j1)*frac0(j1)
  END DO
  frac(i1)=aux*dtime+frac0(i1)
END DO

!      FRAC0=FRAC

RETURN
END SUBROUTINE chemicalstat
SUBROUTINE chempot(phi_t,phi_tau,dpdt,dphidmu,dphidotdmu,mfluid,  &
        dmdmu,dmdj,spucmod,spcumodfac,vmol,mu_tau,theta,  &
        detfe,detfs,detf,cr,kbulk,dffprops,iden,dtime,ndi)



!>    ASDASDASD
use global
IMPLICIT NONE

INTEGER :: nargs,ndi
PARAMETER(nargs=8)
DOUBLE PRECISION :: phi_tau,phi_t,phi_per,phi_m,dtime
DOUBLE PRECISION :: d,mu_tau,mu0,rgas,theta,chi,vmol,kbulk,detf,cr
DOUBLE PRECISION :: dpdt,dpdt_per,dpdt_m,dphidmu,dphidj,dphidotdmu
DOUBLE PRECISION :: detfe,detfs,deltamu,mfluid,dmdmu,dmdj,dmdphi
DOUBLE PRECISION :: spcumodfac(ndi,ndi),spucmod(ndi,ndi),  &
                    iden(ndi,ndi),args(nargs)
DOUBLE PRECISION :: dffprops(5),cbari1

!      ! COMPUTE THE POLYMER VOLUME FRACTION
!      !
chi      = dffprops(1)
d        = dffprops(2)
mu0      = dffprops(3)
vmol     = dffprops(4)
rgas     = dffprops(5)
!
args(1)  = mu_tau
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
!
CALL solvephi(phi_tau,args,nargs,phi_t)
! COMPUTE THE ELASTIC VOLUME RATIO, DETFE
!
cr=(one/vmol)*(one*(phi_tau**(-one))-one)
detfe = detf*phi_tau
!
! COMPUTE THE SWELLING VOLUME RATIO, DETFS
!
detfs = one/phi_tau
! COMPUTE THE TIME RATE OF THE POLYMER VOLUME FRACTION USING
!  A FINITE DIFFERENCE IN TIME
!
dpdt = (phi_tau - phi_t)/dtime


! COMPUTE THE DERIVATIVE OF THE POLYMER VOLUME FRACTION WITH
!  RESPECT TO THE CHEMICAL POTENTIAL.  COMPUTED VIA IMPLICIT
!  DIFFERENTIATION ON THE CHEMICAL POTENTIAL EQUATION.
!
dphidmu =  (one/(rgas*theta))/ (  &
    (one/(phi_tau - one)) + one + two*chi*phi_tau  &
    - ((vmol*kbulk)/(rgas*theta*phi_tau))  &
    + ((vmol*kbulk)/(rgas*theta*phi_tau))*DLOG(detf*phi_tau) )
! COMPUTE THE DERIVATIVE OF THE POLYMER VOLUME FRACTION WITH
!  RESPECT TO THE CHEMICAL POTENTIAL.  COMPUTED VIA IMPLICIT
!  DIFFERENTIATION ON THE CHEMICAL POTENTIAL EQUATION.
!
dphidj = ( (vmol*kbulk)/(rgas*theta*detf)  &
    - ((vmol*kbulk)/(rgas*theta*detf))*DLOG(detf*phi_tau) )/  &
    ( (one/(phi_tau - one)) + one + two*chi*phi_tau  &
    - ((vmol*kbulk)/(rgas*theta*phi_tau))  &
    + ((vmol*kbulk)/(rgas*theta*phi_tau))*DLOG(detf*phi_tau) )
! COMPUTE THE PERTURBATION ON THE CHEMICAL POTENTIAL
!
IF(DABS(mu_tau) > one) THEN
  deltamu = DABS(mu_tau)*1.d-8
ELSE
  deltamu = 1.d-8
END IF

! COMPUTE A PERTURBED POLYMER VOLUME FRACTION
!
args(1)  = mu_tau + deltamu
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
CALL solvephi(phi_per,args,nargs,phi_t)

! COMPUTE A PERTURBED POLYMER VOLUME FRACTION
!
args(1)  = mu_tau - deltamu
args(2)  = mu0
args(3)  = rgas
args(4)  = theta
args(5)  = chi
args(6)  = vmol
args(7)  = kbulk
args(8)  = detf
CALL solvephi(phi_m,args,nargs,phi_t)


! COMPUTE THE DERIVATIVE OF THE TIME RATE OF CHANGE OF THE
!  POLYMER VOLUME FRACTION WITH RESPECT TO THE CHEMICAL POTENTIAL
!
dpdt_per = (phi_per - phi_t)/dtime
dpdt_m   = (phi_m - phi_t)/dtime
dphidotdmu = (dpdt_per - dpdt_m)/(two*deltamu)

! COMPUTE THE FLUID PERMEABILITY AT THIS INTEG. POINT
!
!MFLUID = D/(VMOL*RGAS*THETA)
!
! TO DO M = (D*CR)/(R*T), USE THE FOLLOWING LINE
!MFLUID = (D*(ONE/PHI_TAU - ONE))/(VMOL*RGAS*THETA)
!
! TO DO M = (D*C)/(R*T), USE THE FOLLOWING LINE
!   MFLUID = (cbari1**(cbari1)/3.d0)*
! 1                (D*(ONE/PHI_TAU - ONE))/(DETF*VMOL*RGAS*THETA)
mfluid = (d*(one/phi_tau - one))/(detf*vmol*rgas*theta)
!write(*,*) cbari1
! COMPUTE THE TANGENTS OF THE FLUID MOBILITY
!
!DMDPHI = ZERO
!
! TO DO M = (D*CR)/(R*T), USE THE FOLLOWING LINE
!DMDPHI = -(D/(VMOL*PHI_TAU*PHI_TAU*RGAS*THETA))
!
! TO DO M = (D*C)/(R*T), USE THE FOLLOWING LINE
dmdphi = -(d/(detf*vmol*phi_tau*phi_tau*rgas*theta))
!
dmdmu = dmdphi*dphidmu
dmdj  = dmdphi*dphidj
! COMPUTE THE DISPLACEMENT - CHEMICAL POTENTIAL MODULUS
!
spucmod = (kbulk/(detfe*phi_tau))*iden*dphidmu
! COMPUTE THE CHEMICAL POTENTIAL - DISPLACEMENT MODULUS
!
spcumodfac = mfluid*iden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PHIT=ONE
!      PHITAU=PHIT
!C     TIME DERIVATIVE OF THE FRACTION CONCENTRATION (FINITE DIF)
!      DPDT= (PHITAU - PHIT)/DTIME
!C     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
!      DPHIDMU=ZERO
!C     DERIVATIVE OF CONCENTRATION WRT CHEMICAL POTENTIAL
!      DPHIDMU=ZERO
!      DDUDDC=ONE
!      DDCDDU=ONE

!      CHI    = PROPS(3)
!      D      = PROPS(4)
!      MU0    = PROPS(5)
!      VMOL   = PROPS(6)
!      RGAS   = PROPS(7)
!
!
!        DDUDD
!        DDCDDU
!         MUTAU
!        THETA0
!        PHIT
!        PHITAU
!        DPDT
!        DPHIDMU
!        DPHITDMU,
!         MOB
!        DMDMU
!        DMDJ
!        VMOL = ONE



END SUBROUTINE
SUBROUTINE csisomatfic(cisomatfic,cmisomatfic,distgr,det,ndi)



!>    ISOTROPIC MATRIX: SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: distgr(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



call push4(cisomatfic,cmisomatfic,distgr,det,ndi)

RETURN
END SUBROUTINE csisomatfic
SUBROUTINE cmatisomatfic(cmisomatfic,cbar,cbari1,cbari2,  &
        diso,unit2,unit4,det,ndi)



!>    ISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: cmisomatfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
    

DOUBLE PRECISION :: dudi1,dudi2,d2ud2i1,d2ud2i2,d2udi1i2
DOUBLE PRECISION :: aux,aux1,aux2,aux3,aux4
DOUBLE PRECISION :: uij,ukl,cij,ckl

dudi1=diso(1)
dudi2=diso(2)
d2ud2i1=diso(3)
d2ud2i2=diso(4)
d2udi1i2=diso(5)

aux1=four*(d2ud2i1+two*cbari1*d2udi1i2+ dudi2+cbari1*cbari1*d2ud2i2)
aux2=-four*(d2udi1i2+cbari1*d2ud2i2)
aux3=four*d2ud2i2
aux4=-four*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        uij=unit2(i1,j1)
        ukl=unit2(k1,l1)
        cij=cbar(i1,j1)
        ckl=cbar(k1,l1)
        aux=aux1*uij*ukl+ aux2*(uij*ckl+cij*ukl)+aux3*cij*ckl+  &
            aux4*unit4(i1,j1,k1,l1)
        cmisomatfic(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE cmatisomatfic
SUBROUTINE contractile(fc,ffc,dw,ddw,ffcc,ru,ru0,lambdac,lambdap,  &
        lambda0,l,r0,mu0,beta,b0,ffcmax,fric,frac,dtime)


use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(IN OUT)         :: fc
DOUBLE PRECISION, INTENT(IN OUT)         :: ffc
DOUBLE PRECISION, INTENT(IN OUT)         :: dw
DOUBLE PRECISION, INTENT(IN OUT)         :: ddw
DOUBLE PRECISION, INTENT(IN OUT)         :: ffcc
DOUBLE PRECISION, INTENT(IN OUT)         :: ru
DOUBLE PRECISION, INTENT(IN OUT)             :: ru0
DOUBLE PRECISION, INTENT(OUT)            :: lambdac
DOUBLE PRECISION, INTENT(IN OUT)             :: lambdap
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: l
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION, INTENT(IN OUT)         :: ffcmax
DOUBLE PRECISION, INTENT(IN OUT)         :: fric
DOUBLE PRECISION, INTENT(IN OUT)         :: frac(nch)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime





!      CHECK STRETCH
IF (ru0 > zero) THEN
  lambdac=lambdap+ru0
ELSE
  lambdac=lambdap
END IF
!        IF(LAMBDAC.LT.LAMBDAP)THEN
!          LAMBDAC=LAMBDAP
!           RU=ZERO
!        ENDIF
!        FORCE WITH CONTRACTION
CALL fil(fc,ffc,dw,ddw,lambdac,lambda0,l,r0,mu0,beta,b0)
!        PASSIVE FORCE
!          CALL FIL(FC,FFC,DW,DDW,LAMBDAP,LAMBDA0,L,R0,MU0,BETA,B0)
!        RELATIVE SLIDING
CALL sliding(ffcc,ru,ffc,ru0,ffcmax,fric,frac,dtime)

RETURN

END SUBROUTINE contractile

SUBROUTINE contraction24(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - RIGHT 2ND ORDER TENSOR
!>       RT - LEFT  4TH ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: lt(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO k1=1,ndi
  DO l1=1,ndi
    aux=zero
    DO i1=1,ndi
      DO j1=1,ndi
        aux=aux+lt(k1,l1)*rt(i1,j1,k1,l1)
      END DO
    END DO
    s(k1,l1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction24
SUBROUTINE contraction42(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
!>      INPUT:
!>       LT - left 4TH ORDER TENSOR
!>       RT - right  2ND ODER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi)


INTEGER :: i1,j1,k1,l1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO k1=1,ndi
      DO l1=1,ndi
        aux=aux+LT(i1,j1,k1,l1)*rt(k1,l1)
      END DO
    END DO
    s(i1,j1)=aux
  END DO
END DO
RETURN
END SUBROUTINE contraction42
SUBROUTINE contraction44(s,LT,rt,ndi)



!>       DOUBLE CONTRACTION BETWEEN 4TH ORDER TENSORS
!>      INPUT:
!>       LT - RIGHT 4TH ORDER TENSOR
!>       RT - LEFT  4TH ORDER TENSOR
!>      OUTPUT:
!>       S - DOUBLE CONTRACTED TENSOR (4TH ORDER)
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: s(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: LT(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rt(ndi,ndi,ndi,ndi)



INTEGER :: i1,j1,k1,l1,m1,n1


DOUBLE PRECISION :: aux



DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO m1=1,ndi
          DO n1=1,ndi
            aux=aux+LT(i1,j1,m1,n1)*rt(m1,n1,k1,l1)
          END DO
        END DO
        s(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE contraction44
SUBROUTINE csfilfic(cfic,rho,lambda,dw,ddw,m,rw,ndi)



!>    AFFINE NETWORK: 'FICTICIOUS' ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rho
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN)             :: dw
DOUBLE PRECISION, INTENT(IN)             :: ddw
DOUBLE PRECISION, INTENT(IN)             :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw



INTEGER :: i1,j1,k1,l1

DOUBLE PRECISION :: aux, aux0

aux0=ddw-(lambda**(-one))*dw
aux=rho*aux0*rw*(lambda**(-two))
DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        cfic(i1,j1,k1,l1)=aux*m(i1)*m(j1)*m(k1)*m(l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE csfilfic
SUBROUTINE deffil(lambda,m,m0,f,ndi)



!>      SINGLE FILAMENT: STRETCH AND DEFORMED DIRECTION
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: lambda
DOUBLE PRECISION, INTENT(OUT)            :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: m0(ndi)
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)


INTEGER :: i1,j1

DOUBLE PRECISION :: aux

lambda=zero
DO i1=1,ndi
  aux=zero
  DO j1=1,ndi
    aux=aux+f(i1,j1)*m0(j1)
  END DO
  m(i1)=aux
END DO
lambda=dot_product(m,m)
lambda=SQRT(lambda)

RETURN
END SUBROUTINE deffil
SUBROUTINE deformation(f,c,b,ndi)



!>     RIGHT AND LEFT CAUCHY-GREEN DEFORMATION TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: b(ndi,ndi)


!     RIGHT CAUCHY-GREEN DEFORMATION TENSOR
c=matmul(transpose(f),f)
!     LEFT CAUCHY-GREEN DEFORMATION TENSOR
b=matmul(f,transpose(f))
RETURN
END SUBROUTINE deformation
SUBROUTINE density(rho,ang,bb,erfi)



!>    SINGLE FILAMENT: DENSITY FUNCTION VALUE
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: rho
DOUBLE PRECISION, INTENT(IN OUT)         :: ang
DOUBLE PRECISION, INTENT(IN OUT)         :: bb
DOUBLE PRECISION, INTENT(IN OUT)         :: erfi



DOUBLE PRECISION :: pi,aux1,aux2

pi=four*ATAN(one)
aux1=SQRT(bb/(two*pi))
aux2=DEXP(bb*(COS(two*ang)+one))
rho=four*aux1*aux2*(erfi**(-one))
!      RHO=RHO*((FOUR*PI)**(-ONE))

RETURN
END SUBROUTINE density
!****************************************************************************



!     Element subroutines
!****************************************************************************

SUBROUTINE xint2d1pt(xi,w,ninttpt)
!
! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 1 gauss point for integration
!
!  xi(nInttPt,2): xi,eta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
!
use global

real(8), INTENT(OUT)                      :: xi(1,2)
real(8), INTENT(OUT)                      :: w(1)
INTEGER, INTENT(OUT)                     :: ninttpt

!
INTEGER :: ndim
!


! Initialize
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 1


! Gauss weights
!
w = 4.d0


! Gauss pt location in master element
!
xi(1,1) = 0.d0
xi(1,2) = 0.d0


RETURN
END SUBROUTINE xint2d1pt

!************************************************************************

SUBROUTINE xint2d4pt(xi,w,ninttpt)
!
! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 4 gauss points for integration
!
!  xi(nInttPt,2): xi,eta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
!
use global

real(8), INTENT(OUT)                      :: xi(4,2)
real(8), INTENT(OUT)                      :: w(4)
INTEGER, INTENT(OUT)                     :: ninttpt

!
INTEGER :: ndim
!



! Initialize
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 4


! Gauss weights
!
w(1) = 1.d0
w(2) = 1.d0
w(3) = 1.d0
w(4) = 1.d0


! Gauss pt locations in master element
!
xi(1,1) = -DSQRT(1.d0/3.d0)
xi(1,2) = -DSQRT(1.d0/3.d0)
xi(2,1) = DSQRT(1.d0/3.d0)
xi(2,2) = -DSQRT(1.d0/3.d0)
xi(3,1) = -DSQRT(1.d0/3.d0)
xi(3,2) = DSQRT(1.d0/3.d0)
xi(4,1) = DSQRT(1.d0/3.d0)
xi(4,2) = DSQRT(1.d0/3.d0)


RETURN
END SUBROUTINE xint2d4pt

!***********************************************************************

SUBROUTINE xint3d1pt(xi,w,ninttpt)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using a 2 gauss points for integration
!
!  xi(nInttPt,3): xi,eta,zeta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights

use global

real(8), INTENT(OUT)                      :: xi(1,3)
real(8), INTENT(OUT)                      :: w(1)
INTEGER, INTENT(OUT)                     :: ninttpt


INTEGER :: ndim




! Init
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 1


! Gauss weights
!
w(1) = 8.d0


! Gauss pt locations in master element
!
xi(1,1) = 0.d0
xi(1,2) = 0.d0
xi(1,3) = 0.d0

RETURN
END SUBROUTINE xint3d1pt

!***********************************************************************

SUBROUTINE xint3d8pt(xi,w,ninttpt)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 8 gauss points for integration
!
!  xi(nInttPt,3): xi,eta,zeta coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

real(8), INTENT(OUT)                      :: xi(8,3)
real(8), INTENT(OUT)                      :: w(8)
INTEGER, INTENT(OUT)                     :: ninttpt


INTEGER :: ndim




! Init
!
w = 0.d0
xi = 0.d0


! Number of Gauss points
!
ninttpt = 8


! Gauss weights
!
w(1) = 1.d0
w(2) = 1.d0
w(3) = 1.d0
w(4) = 1.d0
w(5) = 1.d0
w(6) = 1.d0
w(7) = 1.d0
w(8) = 1.d0


! Gauss pt locations in master element
!
xi(1,1) = -DSQRT(1.d0/3.d0)
xi(1,2) = -DSQRT(1.d0/3.d0)
xi(1,3) = -DSQRT(1.d0/3.d0)
xi(2,1) = DSQRT(1.d0/3.d0)
xi(2,2) = -DSQRT(1.d0/3.d0)
xi(2,3) = -DSQRT(1.d0/3.d0)
xi(3,1) = -DSQRT(1.d0/3.d0)
xi(3,2) = DSQRT(1.d0/3.d0)
xi(3,3) = -DSQRT(1.d0/3.d0)
xi(4,1) = DSQRT(1.d0/3.d0)
xi(4,2) = DSQRT(1.d0/3.d0)
xi(4,3) = -DSQRT(1.d0/3.d0)
xi(5,1) = -DSQRT(1.d0/3.d0)
xi(5,2) = -DSQRT(1.d0/3.d0)
xi(5,3) = DSQRT(1.d0/3.d0)
xi(6,1) = DSQRT(1.d0/3.d0)
xi(6,2) = -DSQRT(1.d0/3.d0)
xi(6,3) = DSQRT(1.d0/3.d0)
xi(7,1) = -DSQRT(1.d0/3.d0)
xi(7,2) = DSQRT(1.d0/3.d0)
xi(7,3) = DSQRT(1.d0/3.d0)
xi(8,1) = DSQRT(1.d0/3.d0)
xi(8,2) = DSQRT(1.d0/3.d0)
xi(8,3) = DSQRT(1.d0/3.d0)


RETURN
END SUBROUTINE xint3d8pt

!***********************************************************************

SUBROUTINE xintsurf2d1pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(1)
real(8), INTENT(OUT)                      :: ylocal(1)
real(8), INTENT(OUT)                      :: w(1)



! Gauss weights
!
w(1) = two


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = zero
  ylocal(1) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = zero
ELSE IF(face == 3) THEN
  xlocal(1) = zero
  ylocal(1) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = zero
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d1pt

!***********************************************************************

SUBROUTINE xintsurf2d2pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(2)
real(8), INTENT(OUT)                      :: ylocal(2)
real(8), INTENT(OUT)                      :: w(2)



! Gauss weights
!
w(1) = one
w(2) = one


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(one/three)
  xlocal(2) = one
  ylocal(2) = DSQRT(one/three)
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = DSQRT(one/three)
  xlocal(2) = -one
  ylocal(2) = -DSQRT(one/three)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d2pt

!***********************************************************************

SUBROUTINE xintsurf2d3pt(face,xlocal,ylocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 2D elements
!  using 2 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(3)
real(8), INTENT(OUT)                      :: ylocal(3)
real(8), INTENT(OUT)                      :: w(3)



real(8), PARAMETER :: five=5.d0
real(8), PARAMETER :: eight=8.d0
real(8), PARAMETER :: nine=9.d0


! Gauss weights
!
w(1) = five/nine
w(2) = eight/nine
w(3) = five/nine


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(three/five)
  ylocal(1) = -one
  xlocal(2) = zero
  ylocal(2) = -one
  xlocal(2) = DSQRT(three/five)
  ylocal(2) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(three/five)
  xlocal(2) = one
  ylocal(2) = zero
  xlocal(3) = one
  ylocal(3) = DSQRT(three/five)
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(three/five)
  ylocal(1) = one
  xlocal(2) = zero
  ylocal(2) = one
  xlocal(3) = DSQRT(three/five)
  ylocal(3) = one
ELSE IF(face == 4) THEN
  xlocal(1) = -one
  ylocal(1) = DSQRT(three/five)
  xlocal(2) = -one
  ylocal(2) = zero
  xlocal(3) = -one
  ylocal(3) = -DSQRT(three/five)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4'
! !write(80,*) 'face.ne.1,2,3,4'
  CALL xit()
END IF

END SUBROUTINE xintsurf2d3pt

!***********************************************************************

SUBROUTINE xintsurf3d1pt(face,xlocal,ylocal,zlocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 1 gauss point for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  zLocal(nInttPt): z coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(1)
real(8), INTENT(OUT)                      :: ylocal(1)
real(8), INTENT(OUT)                      :: zlocal(1)
real(8), INTENT(OUT)                      :: w(1)




! Gauss weights
!
w(1) = four


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = zero
  ylocal(1) = zero
  zlocal(1) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = zero
  ylocal(1) = zero
  zlocal(1) = one
ELSE IF(face == 3) THEN
  xlocal(1) = zero
  ylocal(1) = -one
  zlocal(1) = zero
ELSE IF(face == 4) THEN
  xlocal(1) = one
  ylocal(1) = zero
  zlocal(1) = zero
ELSE IF(face == 5) THEN
  xlocal(1) = zero
  ylocal(1) = one
  zlocal(1) = zero
ELSE IF(face == 6) THEN
  xlocal(1) = -one
  ylocal(1) = zero
  zlocal(1) = zero
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4,5,6'
! !write(80,*) 'face.ne.1,2,3,4,5,6'
  CALL xit()
END IF

END SUBROUTINE xintsurf3d1pt

!***********************************************************************

SUBROUTINE xintsurf3d4pt(face,xlocal,ylocal,zlocal,w)

! This subroutine will get the integration point locations
!  and corresponding gauss quadrature weights for 3D elements
!  using 4 gauss points for surface integration
!
!  xLocal(nInttPt): x coordinates for the integration pts
!  yLocal(nInttPt): y coordinates for the integration pts
!  yLocal(nInttPt): z coordinates for the integration pts
!  w(nInttPt):    corresponding integration weights
use global

INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(OUT)                      :: xlocal(4)
real(8), INTENT(OUT)                      :: ylocal(4)
real(8), INTENT(OUT)                      :: zlocal(4)
real(8), INTENT(OUT)                      :: w(4)


! Gauss weights
!
w(1) = one
w(2) = one
w(3) = one
w(4) = one


! Gauss pt locations in master element
!
IF(face == 1) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -DSQRT(one/three)
  zlocal(2) = -one
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = -one
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = DSQRT(one/three)
  zlocal(4) = -one
ELSE IF(face == 2) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = one
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -DSQRT(one/three)
  zlocal(2) = one
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = one
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = DSQRT(one/three)
  zlocal(4) = one
ELSE IF(face == 3) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = -one
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = -one
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = -one
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = -one
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 4) THEN
  xlocal(1) = one
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = one
  ylocal(2) = DSQRT(one/three)
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = one
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = one
  ylocal(4) = -DSQRT(one/three)
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 5) THEN
  xlocal(1) = -DSQRT(one/three)
  ylocal(1) = one
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = DSQRT(one/three)
  ylocal(2) = one
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = DSQRT(one/three)
  ylocal(3) = one
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -DSQRT(one/three)
  ylocal(4) = one
  zlocal(4) = DSQRT(one/three)
ELSE IF(face == 6) THEN
  xlocal(1) = -one
  ylocal(1) = -DSQRT(one/three)
  zlocal(1) = -DSQRT(one/three)
  xlocal(2) = -one
  ylocal(2) = DSQRT(one/three)
  zlocal(2) = -DSQRT(one/three)
  xlocal(3) = -one
  ylocal(3) = DSQRT(one/three)
  zlocal(3) = DSQRT(one/three)
  xlocal(4) = -one
  ylocal(4) = -DSQRT(one/three)
  zlocal(4) = DSQRT(one/three)
ELSE
  WRITE(*,*) 'face.ne.1,2,3,4,5,6'
! !write(80,*) 'face.ne.1,2,3,4,5,6'
  CALL xit()
END IF

END SUBROUTINE xintsurf3d4pt

!************************************************************************

SUBROUTINE calcshape2dlinear(ninttpt,xi_int,intpt,sh,dshxi)
!
! Calculate the shape functions and their derivatives at the
! given integration point in the master element


! Calculate the shape functions and their derivatives at the
! given integration point in the master element
!
!                          eta
!   4-----------3          |
!   |           |          |
!   |           |          |
!   |           |          |
!   |           |          |
!   |           |          O--------- xi
!   1-----------2        origin at center
!
!
! sh(i) = shape function of node i at the intpt.
! dshxi(i,j) = derivative wrt j direction of shape fn of node i
!
use global
INTEGER, INTENT(IN OUT)                  :: ninttpt
real(8), INTENT(IN)                       :: xi_int(ninttpt,2)
INTEGER, INTENT(IN OUT)                  :: intpt
real(8), INTENT(OUT)                      :: sh(4)
real(8), INTENT(OUT)                      :: dshxi(4,2)

!
INTEGER :: ndim
!
real(8)  xi,eta
!
real(8), PARAMETER :: fourth=1.d0/4.d0


! Location in the master element
!
xi = xi_int(intpt,1)
eta = xi_int(intpt,2)


! The shape functions
!
sh(1) = fourth*(one - xi)*(one - eta)
sh(2) = fourth*(one + xi)*(one - eta)
sh(3) = fourth*(one + xi)*(one + eta)
sh(4) = fourth*(one - xi)*(one + eta)


! The first derivatives
!
dshxi(1,1) = -fourth*(one - eta)
dshxi(1,2) = -fourth*(one - xi)
dshxi(2,1) = fourth*(one - eta)
dshxi(2,2) = -fourth*(one + xi)
dshxi(3,1) = fourth*(one + eta)
dshxi(3,2) = fourth*(one + xi)
dshxi(4,1) = -fourth*(one + eta)
dshxi(4,2) = fourth*(one - xi)


RETURN
END SUBROUTINE calcshape2dlinear

!***********************************************************************

SUBROUTINE calcshape3dlinear(ninttpt,xi_int,intpt,sh,dshxi)
!
!
! Calculate the shape functions and their derivatives at the
! given integration point in the master element
!
! This subroutine uses a 8-node linear 3D element as shown
!
!      8-----------7
!     /|          /|       zeta
!    / |         / |
!   5-----------6  |       |     eta
!   |  |        |  |       |   /
!   |  |        |  |       |  /
!   |  4--------|--3       | /
!   | /         | /        |/
!   |/          |/         O--------- xi
!   1-----------2        origin at cube center
!
!
! sh(i) = shape function of node i at the intpt.
! dshxi(i,j) = derivative wrt j direction of shape fn of node i
! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i
use global

INTEGER, INTENT(IN OUT)                  :: ninttpt
real(8), INTENT(IN)                       :: xi_int(ninttpt,3)
INTEGER, INTENT(IN OUT)                  :: intpt
real(8), INTENT(OUT)                      :: sh(8)
real(8), INTENT(OUT)                      :: dshxi(8,3)


INTEGER :: ndim, i,j


real(8) d2shxi(8,3,3),xi,eta,zeta

real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Location in the master element
!
xi = xi_int(intpt,1)
eta = xi_int(intpt,2)
zeta = xi_int(intpt,3)
!
! The shape functions
!
sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
!
! The first derivatives
!
dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
dshxi(1,3) = -eighth*(one - xi)*(one - eta)
dshxi(2,1) = eighth*(one - eta)*(one - zeta)
dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
dshxi(2,3) = -eighth*(one + xi)*(one - eta)
dshxi(3,1) = eighth*(one + eta)*(one - zeta)
dshxi(3,2) = eighth*(one + xi)*(one - zeta)
dshxi(3,3) = -eighth*(one + xi)*(one + eta)
dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
dshxi(4,2) = eighth*(one - xi)*(one - zeta)
dshxi(4,3) = -eighth*(one - xi)*(one + eta)
dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
dshxi(5,3) = eighth*(one - xi)*(one - eta)
dshxi(6,1) = eighth*(one - eta)*(one + zeta)
dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
dshxi(6,3) = eighth*(one + xi)*(one - eta)
dshxi(7,1) = eighth*(one + eta)*(one + zeta)
dshxi(7,2) = eighth*(one + xi)*(one + zeta)
dshxi(7,3) = eighth*(one + xi)*(one + eta)
dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
dshxi(8,2) = eighth*(one - xi)*(one + zeta)
dshxi(8,3) = eighth*(one - xi)*(one + eta)
!
! The second derivatives
!
d2shxi = zero
d2shxi(1,1,2) = eighth*(one - zeta)
d2shxi(1,2,1) = d2shxi(1,1,2)
d2shxi(1,1,3) = eighth*(one - eta)
d2shxi(1,3,1) = d2shxi(1,1,3)
d2shxi(1,2,3) = eighth*(one - xi)
d2shxi(1,3,2) = d2shxi(1,2,3)
d2shxi(2,1,2) = -eighth*(one - zeta)
d2shxi(2,2,1) = d2shxi(2,1,2)
d2shxi(2,1,3) = -eighth*(one - eta)
d2shxi(2,3,1) = d2shxi(2,1,3)
d2shxi(2,2,3) = eighth*(one + xi)
d2shxi(2,3,2) = d2shxi(2,2,3)
d2shxi(3,1,2) = eighth*(one - zeta)
d2shxi(3,2,1) = d2shxi(2,1,2)
d2shxi(3,1,3) = -eighth*(one + eta)
d2shxi(3,3,1) = d2shxi(2,1,3)
d2shxi(3,2,3) = -eighth*(one + xi)
d2shxi(3,3,2) = d2shxi(2,2,3)
d2shxi(4,1,2) = -eighth*(one - zeta)
d2shxi(4,2,1) = d2shxi(2,1,2)
d2shxi(4,1,3) = eighth*(one + eta)
d2shxi(4,3,1) = d2shxi(2,1,3)
d2shxi(4,2,3) = -eighth*(one - xi)
d2shxi(4,3,2) = d2shxi(2,2,3)
d2shxi(5,1,2) = eighth*(one + zeta)
d2shxi(5,2,1) = d2shxi(2,1,2)
d2shxi(5,1,3) = -eighth*(one - eta)
d2shxi(5,3,1) = d2shxi(2,1,3)
d2shxi(5,2,3) = -eighth*(one - xi)
d2shxi(5,3,2) = d2shxi(2,2,3)
d2shxi(6,1,2) = eighth*(one + zeta)
d2shxi(6,2,1) = d2shxi(2,1,2)
d2shxi(6,1,3) = eighth*(one - eta)
d2shxi(6,3,1) = d2shxi(2,1,3)
d2shxi(6,2,3) = -eighth*(one + xi)
d2shxi(6,3,2) = d2shxi(2,2,3)
d2shxi(7,1,2) = eighth*(one + zeta)
d2shxi(7,2,1) = d2shxi(2,1,2)
d2shxi(7,1,3) = eighth*(one + eta)
d2shxi(7,3,1) = d2shxi(2,1,3)
d2shxi(7,2,3) = eighth*(one + xi)
d2shxi(7,3,2) = d2shxi(2,2,3)
d2shxi(8,1,2) = -eighth*(one + zeta)
d2shxi(8,2,1) = d2shxi(2,1,2)
d2shxi(8,1,3) = -eighth*(one + eta)
d2shxi(8,3,1) = d2shxi(2,1,3)
d2shxi(8,2,3) = eighth*(one - xi)
d2shxi(8,3,2) = d2shxi(2,2,3)

RETURN
END SUBROUTINE calcshape3dlinear

!************************************************************************


SUBROUTINE computesurf(xlocal,ylocal,face,coords,sh,ds)

! This subroutine computes the shape functions, derivatives
!  of shape functions, and the length ds, so that one can
!  do the numerical integration on the boundary for fluxes
!  on the 4-node quadrilateral elements

use global

real(8), INTENT(IN OUT)                   :: xlocal
real(8), INTENT(IN OUT)                   :: ylocal
INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(IN)                       :: coords(2,4)
real(8), INTENT(OUT)                      :: sh(4)
real(8), INTENT(OUT)                      :: ds



real(8)  dshxi(4,2), dxdxi,dxdeta,dydxi
real(8) dydeta, shape,normal(2,1)
real(8), PARAMETER :: fourth=1.d0/4.d0

sh(1) = fourth*(one - xlocal)*(one - ylocal)
sh(2) = fourth*(one + xlocal)*(one - ylocal)
sh(3) = fourth*(one + xlocal)*(one + ylocal)
sh(4) = fourth*(one - xlocal)*(one + ylocal)

dshxi(1,1) = -fourth*(one - ylocal)
dshxi(1,2) = -fourth*(one - xlocal)
dshxi(2,1) = fourth*(one - ylocal)
dshxi(2,2) = -fourth*(one + xlocal)
dshxi(3,1) = fourth*(one + ylocal)
dshxi(3,2) = fourth*(one + xlocal)
dshxi(4,1) = -fourth*(one + ylocal)
dshxi(4,2) = fourth*(one - xlocal)

dxdxi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)  &
    + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
dxdeta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)  &
    + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
dydxi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)  &
    + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
dydeta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)  &
    + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)


! Jacobian of the mapping
!
IF((face == 2).OR.(face == 4)) THEN
  ds = DSQRT(dxdeta*dxdeta + dydeta*dydeta)
ELSE IF((face == 1).OR.(face == 3)) THEN
  ds = DSQRT(dxdxi*dxdxi + dydxi*dydxi)
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF


! Surface normal, outward pointing in this case. Useful for
!  ``follower'' type loads. The normal is referential or spatial
!  depending on which coords were supplied to this subroutine
!  (NOT fully tested)
!
IF((face == 2).OR.(face == 4)) THEN
  normal(1,1) = dydeta/DSQRT(dxdeta*dxdeta + dydeta*dydeta)
  normal(2,1) = -dxdeta/DSQRT(dxdeta*dxdeta + dydeta*dydeta)
  IF(face == 4) normal = -normal
ELSE IF((face == 1).OR.(face == 3)) THEN
  normal(1,1) = dydxi/DSQRT(dxdxi*dxdxi + dydxi*dydxi)
  normal(2,1) = -dxdxi/DSQRT(dxdxi*dxdxi + dydxi*dydxi)
  IF(face == 3) normal = -normal
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF

RETURN
END SUBROUTINE computesurf

!***********************************************************************

SUBROUTINE computesurf3d(xlocal,ylocal,zlocal,face,coords,sh,da)

! This subroutine computes the shape functions, derivatives
!  of shape functions, and the area dA, so that one can
!  do the numerical integration on the boundary for fluxes
!  on the 8-node brick elements

use global

real(8), INTENT(IN OUT)                   :: xlocal
real(8), INTENT(IN OUT)                   :: ylocal
real(8), INTENT(IN OUT)                   :: zlocal
INTEGER, INTENT(IN OUT)                  :: face
real(8), INTENT(IN)                       :: coords(3,8)
real(8), INTENT(OUT)                      :: sh(8)
real(8), INTENT(OUT)                      :: da

INTEGER :: stat,i,j,k

real(8)  dshxi(8,3), dsh(8,3)
real(8)  mapj(3,3),mag,normal(3,1)

real(8) dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta
real(8) dzdzeta


real(8), PARAMETER :: eighth=1.d0/8.d0


! The shape functions
!
sh(1) = eighth*(one - xlocal)*(one - ylocal)*(one - zlocal)
sh(2) = eighth*(one + xlocal)*(one - ylocal)*(one - zlocal)
sh(3) = eighth*(one + xlocal)*(one + ylocal)*(one - zlocal)
sh(4) = eighth*(one - xlocal)*(one + ylocal)*(one - zlocal)
sh(5) = eighth*(one - xlocal)*(one - ylocal)*(one + zlocal)
sh(6) = eighth*(one + xlocal)*(one - ylocal)*(one + zlocal)
sh(7) = eighth*(one + xlocal)*(one + ylocal)*(one + zlocal)
sh(8) = eighth*(one - xlocal)*(one + ylocal)*(one + zlocal)


! Shape function derivatives
!
dshxi(1,1) = -eighth*(one - ylocal)*(one - zlocal)
dshxi(1,2) = -eighth*(one - xlocal)*(one - zlocal)
dshxi(1,3) = -eighth*(one - xlocal)*(one - ylocal)
dshxi(2,1) = eighth*(one - ylocal)*(one - zlocal)
dshxi(2,2) = -eighth*(one + xlocal)*(one - zlocal)
dshxi(2,3) = -eighth*(one + xlocal)*(one - ylocal)
dshxi(3,1) = eighth*(one + ylocal)*(one - zlocal)
dshxi(3,2) = eighth*(one + xlocal)*(one - zlocal)
dshxi(3,3) = -eighth*(one + xlocal)*(one + ylocal)
dshxi(4,1) = -eighth*(one + ylocal)*(one - zlocal)
dshxi(4,2) = eighth*(one - xlocal)*(one - zlocal)
dshxi(4,3) = -eighth*(one - xlocal)*(one + ylocal)
dshxi(5,1) = -eighth*(one - ylocal)*(one + zlocal)
dshxi(5,2) = -eighth*(one - xlocal)*(one + zlocal)
dshxi(5,3) = eighth*(one - xlocal)*(one - ylocal)
dshxi(6,1) = eighth*(one - ylocal)*(one + zlocal)
dshxi(6,2) = -eighth*(one + xlocal)*(one + zlocal)
dshxi(6,3) = eighth*(one + xlocal)*(one - ylocal)
dshxi(7,1) = eighth*(one + ylocal)*(one + zlocal)
dshxi(7,2) = eighth*(one + xlocal)*(one + zlocal)
dshxi(7,3) = eighth*(one + xlocal)*(one + ylocal)
dshxi(8,1) = -eighth*(one + ylocal)*(one + zlocal)
dshxi(8,2) = eighth*(one - xlocal)*(one + zlocal)
dshxi(8,3) = eighth*(one - xlocal)*(one + ylocal)


dxdxi = zero
dxdeta = zero
dxdzeta = zero
dydxi = zero
dydeta = zero
dydzeta = zero
dzdxi = zero
dzdeta = zero
dzdzeta = zero
DO k=1,8
  dxdxi = dxdxi + dshxi(k,1)*coords(1,k)
  dxdeta = dxdeta + dshxi(k,2)*coords(1,k)
  dxdzeta = dxdzeta + dshxi(k,3)*coords(1,k)
  dydxi = dydxi + dshxi(k,1)*coords(2,k)
  dydeta = dydeta + dshxi(k,2)*coords(2,k)
  dydzeta = dydzeta + dshxi(k,3)*coords(2,k)
  dzdxi = dzdxi + dshxi(k,1)*coords(3,k)
  dzdeta = dzdeta + dshxi(k,2)*coords(3,k)
  dzdzeta = dzdzeta + dshxi(k,3)*coords(3,k)
END DO


! Jacobian of the mapping
!
IF((face == 1).OR.(face == 2)) THEN
! zeta = constant on this face
  da = DSQRT( (dydxi*dzdeta - dydeta*dzdxi)**two  &
      + (dxdxi*dzdeta - dxdeta*dzdxi)**two + (dxdxi*dydeta - dxdeta*dydxi)**two  &
      )
ELSE IF((face == 3).OR.(face == 5)) THEN
! eta = constant on this face
  da = DSQRT( (dydxi*dzdzeta - dydzeta*dzdxi)**two  &
      + (dxdxi*dzdzeta - dxdzeta*dzdxi)**two  &
      + (dxdxi*dydzeta - dxdzeta*dydxi)**two )
ELSE IF((face == 4).OR.(face == 6)) THEN
! xi = constant on this face
  da = DSQRT( (dydeta*dzdzeta - dydzeta*dzdeta)**two  &
      + (dxdeta*dzdzeta - dxdzeta*dzdeta)**two  &
      + (dxdeta*dydzeta - dxdzeta*dydeta)**two )
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF


! Surface normal, outward pointing in this case. Useful for
!  ``follower'' type loads. The normal is referential or spatial
!  depending on which coords were supplied to this subroutine
!  (NOT fully tested)
!
IF((face == 1).OR.(face == 2)) THEN
! zeta = constant on this face
  normal(1,1) = dydxi*dzdeta - dydeta*dzdxi
  normal(2,1) = dxdxi*dzdeta - dxdeta*dzdxi
  normal(3,1) = dxdxi*dydeta - dxdeta*dydxi
  IF(face == 1) normal = -normal
ELSE IF((face == 3).OR.(face == 5)) THEN
! eta = constant on this face
  normal(1,1) = dydxi*dzdzeta - dydzeta*dzdxi
  normal(2,1) = dxdxi*dzdzeta - dxdzeta*dzdxi
  normal(3,1) = dxdxi*dydzeta - dxdzeta*dydxi
  IF(face == 5) normal = -normal
ELSE IF((face == 4).OR.(face == 6)) THEN
! xi = constant on this face
  normal(1,1) = dydeta*dzdzeta - dydzeta*dzdeta
  normal(2,1) = dxdeta*dzdzeta - dxdzeta*dzdeta
  normal(3,1) = dxdeta*dydzeta - dxdzeta*dydeta
  IF(face == 6) normal = -normal
ELSE
  WRITE(*,*) 'never should get here'
  CALL xit()
END IF
mag = DSQRT(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
normal(1,1) = normal(1,1)/mag
normal(2,1) = normal(2,1)/mag
normal(3,1) = normal(3,1)/mag

END SUBROUTINE computesurf3d

!***********************************************************************

SUBROUTINE mapshape2d(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.
!
use global
INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,2)
real(8), INTENT(IN)                       :: coords(3,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,2)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat

!
INTEGER :: i,j,k, ieror
!
real(8)  mapj(2,2), mapj_inv(2,2)
!
real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,2
  DO j=1,2
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv2d(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape2D'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


RETURN
END SUBROUTINE mapshape2d

!*************************************************************************

SUBROUTINE mapshape2da(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.
!
! This subroutine is exactly the same as the regular mapShape2D
!with the exception that coords(2,nNode) here and coords(3,nNode)
!in the regular.  I have noticed that a "heat transfer" and
!"static" step uses MCRD=2,
!but for "coupled-temperature-displacement"
!you will get MCRD=3, even for a plane analysis.
use global

INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,2)
real(8), INTENT(IN)                       :: coords(2,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,2)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat

!
INTEGER :: i,j,k, ieror
!
real(8)  mapj(2,2), mapj_inv(2,2)
!
real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,2
  DO j=1,2
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv2d(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape2Da'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


RETURN
END SUBROUTINE mapshape2da

!***********************************************************************

SUBROUTINE mapshape3d(nnode,dshxi,coords,dsh,detmapj,stat)
!
! Map derivatives of shape fns from xi-eta-zeta domain
!  to x-y-z domain.  This subroutine works for both 8-node
!  linear and 20-node quadratic 3D elements.
!
use global
INTEGER, INTENT(IN)                      :: nnode
real(8), INTENT(IN)                       :: dshxi(nnode,3)
real(8), INTENT(IN)                       :: coords(3,nnode)
real(8), INTENT(OUT)                      :: dsh(nnode,3)
real(8), INTENT(IN OUT)                   :: detmapj
INTEGER, INTENT(IN OUT)                  :: stat


INTEGER :: i,j,k, ieror


real(8) mapj(3,3),mapj_inv(3,3)

real(8), PARAMETER :: fourth=0.25D0
real(8), PARAMETER :: eighth=1.d0/8.d0


! Calculate the mapping Jacobian matrix:
!
mapj = zero
DO i=1,3
  DO j=1,3
    DO k=1,nnode
      mapj(i,j) = mapj(i,j) + dshxi(k,i)*coords(j,k)
    END DO
  END DO
END DO


! Calculate the inverse and the determinant of Jacobian
!
CALL matinv3dd(mapj,mapj_inv,detmapj,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero in mapShape3D'
  CALL xit()
END IF


! Calculate first derivatives wrt x, y, z
!
dsh = transpose(matmul(mapj_inv,transpose(dshxi)))


! The second derivatives may be calculated.
!

RETURN
END SUBROUTINE mapshape3d
SUBROUTINE erfi(erf,b,nterm)



!>    IMAGINARY ERROR FUNCTION OF SQRT(B); B IS THE DISPERSION PARAM
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: erf
DOUBLE PRECISION, INTENT(IN OUT)         :: b
INTEGER, INTENT(IN)                      :: nterm


DOUBLE PRECISION :: pi
DOUBLE PRECISION :: aux,aux1,aux2,aux3,aux4,fact
INTEGER :: i1,j1

pi=four*ATAN(one)
aux=SQRT(two*b)
aux1=two*aux
aux2=(two/three)*(aux**three)
aux4=zero
DO j1=3,nterm
  i1=j1-1
  CALL factorial (fact,i1)
  aux3=two*j1-one
  aux4=aux4+(aux**aux3)/(half*aux3*fact)
END DO

erf=pi**(-one/two)*(aux1+aux2+aux4)
RETURN
END SUBROUTINE erfi
SUBROUTINE evalg(g,f,lambda,lambda0,l,r0,mu0,beta,b0)



!>     ESTABLISHMENT OF G(F)=LHS-RHS(F) THAT RELATES
!>       STRETCHFORCE RELATIONSHIP OF A SINGLE EXNTESIBLE FILAMENT
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: g
DOUBLE PRECISION, INTENT(IN)             :: f
DOUBLE PRECISION, INTENT(IN)             :: lambda
DOUBLE PRECISION, INTENT(IN)             :: lambda0
DOUBLE PRECISION, INTENT(IN)             :: l
DOUBLE PRECISION, INTENT(IN)             :: r0
DOUBLE PRECISION, INTENT(IN)             :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0


DOUBLE PRECISION :: lhs,rhs


DOUBLE PRECISION :: aux0,aux1,aux,aux2,aux3,aux4,pi

pi=four*ATAN(one)
aux0=one-r0/l
aux1=l*l*((pi*pi*b0)**(-one))
aux=f/mu0
aux2=one+aux
aux3=one+two*aux
aux4=one+f*aux1+f*aux*aux1

rhs=one+aux-aux0*(aux2**beta)*aux3*(aux4**(-beta))
lhs=lambda*lambda0*r0*(l**(-one))

g=lhs-rhs

RETURN
END SUBROUTINE evalg
SUBROUTINE factorial(fact,term)



!>    FACTORIAL
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: fact
INTEGER, INTENT(IN)                      :: term



INTEGER :: m

fact = 1

DO  m = 1, term
  fact = fact * m
END DO

RETURN
END SUBROUTINE factorial
SUBROUTINE fil(f,ff,dw,ddw,lambda,lambda0,ll,r0,mu0,beta,b0)



!>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: f
DOUBLE PRECISION, INTENT(OUT)            :: ff
DOUBLE PRECISION, INTENT(OUT)            :: dw
DOUBLE PRECISION, INTENT(OUT)            :: ddw
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: ll
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0




DOUBLE PRECISION :: a,b,machep,t
DOUBLE PRECISION :: aux, pi,alpha
DOUBLE PRECISION :: aux0,aux1,aux2,aux3,aux4,aux5,aux6,y

a=zero
b=1.0E09
machep=2.2204E-16
t=1.0E-6
f=zero

CALL pullforce(f, a, b, machep, t, lambda,lambda0,ll,r0,mu0,beta,b0)

pi=four*ATAN(one)
ff=f*ll*(pi*pi*b0)**(-one)
alpha=pi*pi*b0*(ll*ll*mu0)**(-one)

aux0=beta/alpha
aux=alpha*ff
aux1=one+ff+aux*ff
aux2=one+two*aux
aux3=one+aux
aux4=lambda0*r0*r0*mu0*(ll**(-one))
aux5=((one+aux)*(aux1**(-one)))**beta
aux6=one-r0*((ll)**(-one))

y=aux0*(aux2*aux2*(aux1**(-one)))-beta*(aux2*(aux3**(-one)))-two

dw=lambda0*r0*f
ddw=aux4*((one+y*aux5*aux6)**(-one))

RETURN
END SUBROUTINE fil
SUBROUTINE fslip(f,fbar,det,ndi)



!>     DISTORTION GRADIENT
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: fbar(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: det


INTEGER :: i1,j1

DOUBLE PRECISION :: scale1

!     JACOBIAN DETERMINANT
det = f(1,1) * f(2,2) * f(3,3) - f(1,2) * f(2,1) * f(3,3)

IF (ndi == 3) THEN
  det = det + f(1,2) * f(2,3) * f(3,1) + f(1,3) * f(3,2) * f(2,1)  &
      - f(1,3) * f(3,1) * f(2,2) - f(2,3) * f(3,2) * f(1,1)
END IF

scale1=det**(-one /three)

DO i1=1,ndi
  DO j1=1,ndi
    fbar(i1,j1)=scale1*f(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE fslip
SUBROUTINE hfilfic(h,hh,pp,lambda,m,rw,ndi)



!>      NON-AFFINE NETWORK: STRUCTURE TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: h(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: hh(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pp
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN)             :: m(ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw



INTEGER :: i1,j1,k1,l1

DOUBLE PRECISION :: aux0,aux, pi,aux1

pi=four*ATAN(one)
aux0=four*pi
aux=(lambda**(pp-two))*rw
aux1=(pp-two)*(lambda**(pp-four))*rw

DO i1=1,ndi
  DO j1=1,ndi
    h(i1,j1)=aux*m(i1)*m(j1)
    DO k1=1,ndi
      DO l1=1,ndi
        hh(i1,j1,k1,l1)=aux1*m(i1)*m(j1)*m(k1)*m(l1)
      END DO
    END DO
  END DO
END DO

RETURN

END SUBROUTINE hfilfic
SUBROUTINE hvread(hv,statev,v1,ndi)



!>    VISCOUS DISSIPATION: READ STATE VARS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                  :: ndi

DOUBLE PRECISION, INTENT(OUT)            :: hv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
INTEGER, INTENT(IN)                      :: v1



INTEGER :: pos


pos=9*v1-9
hv(1,1)=statev(1+pos)
hv(1,2)=statev(2+pos)
hv(1,3)=statev(3+pos)
hv(2,1)=statev(4+pos)
hv(2,2)=statev(5+pos)
hv(2,3)=statev(6+pos)
hv(3,1)=statev(7+pos)
hv(3,2)=statev(8+pos)
hv(3,3)=statev(9+pos)

RETURN

END SUBROUTINE hvread
SUBROUTINE hvwrite(statev,hv,v1,ndi)



!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                  :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
DOUBLE PRECISION, INTENT(IN)             :: hv(ndi,ndi)
INTEGER, INTENT(IN)                      :: v1



INTEGER :: pos


pos=9*v1-9
statev(1+pos)=hv(1,1)
statev(2+pos)=hv(1,2)
statev(3+pos)=hv(1,3)
statev(4+pos)=hv(2,1)
statev(5+pos)=hv(2,2)
statev(6+pos)=hv(2,3)
statev(7+pos)=hv(3,1)
statev(8+pos)=hv(3,2)
statev(9+pos)=hv(3,3)

RETURN

END SUBROUTINE hvwrite
SUBROUTINE onem(a,aa,aas,ndi)



!>      THIS SUBROUTINE GIVES:
!>          2ND ORDER IDENTITY TENSORS - A
!>          4TH ORDER IDENTITY TENSOR - AA
!>          4TH ORDER SYMMETRIC IDENTITY TENSOR - AAS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aas(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    IF (i == j) THEN
      a(i,j) = one
    ELSE
      a(i,j) = zero
    END IF
  END DO
END DO

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        aa(i,j,k,l)=a(i,k)*a(j,l)
        aas(i,j,k,l)=(one/two)*(a(i,k)*a(j,l)+a(i,l)*a(j,k))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE onem
SUBROUTINE indexx(stress,ddsdde,sig,tng,ntens,ndi)



!>    INDEXATION: FULL SIMMETRY  IN STRESSES AND ELASTICITY TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
INTEGER, INTENT(IN)                      :: ntens
DOUBLE PRECISION, INTENT(OUT)            :: stress(ntens)
DOUBLE PRECISION, INTENT(OUT)            :: ddsdde(ntens,ntens)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: tng(ndi,ndi,ndi,ndi)



INTEGER :: ii1(6),ii2(6), i1,j1


DOUBLE PRECISION :: pp1,pp2

ii1(1)=1
ii1(2)=2
ii1(3)=3
ii1(4)=1
ii1(5)=1
ii1(6)=2

ii2(1)=1
ii2(2)=2
ii2(3)=3
ii2(4)=2
ii2(5)=3
ii2(6)=3

DO i1=1,ntens
!       STRESS VECTOR
  stress(i1)=sig(ii1(i1),ii2(i1))
  DO j1=1,ntens
!       DDSDDE - FULLY SIMMETRY IMPOSED
    pp1=tng(ii1(i1),ii2(i1),ii1(j1),ii2(j1))
    pp2=tng(ii1(i1),ii2(i1),ii2(j1),ii1(j1))
    ddsdde(i1,j1)=(one/two)*(pp1+pp2)
  END DO
END DO

RETURN

END SUBROUTINE indexx
SUBROUTINE initialize(statev,phi0)


use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
DOUBLE PRECISION, INTENT(IN)            :: phi0

statev(1)=phi0

RETURN

END SUBROUTINE initialize
SUBROUTINE integ(props,nprops,dtime,  &
        f_tau,mu_tau,phi_t,theta,  &
        t_tau,sptanmod,  &
        phi_tau)



! This subroutine computes everything required for the time integration
! of the problem.
!
! Inputs:
!  1) material parameters, props(nprops)
!  2) time increment, dtime
!  3) deformation gradient, F_tau(3,3)
!  4) chemical potential, mu_tau
!  5) old polymer volume fraction, phi_t
!  6) temperature, theta
!
! Outputs:
!  1) Cauchy stress, T_tau(3,3)
!  2) spatial tangent modulus, SpTanMod(3,3,3,3)
!  3) polymer volume fraction, phi_tau
!  4) time rate of polymer volume fraction, dPdt
!  5) derivative of the phi with mu, DphiDmu
!  6) derivative of the time rate of phi with mu, DphidotDmu
!  7) scalar fluid permeability, Mfluid
!  8) derivative of permeability with chemical potential, DmDmu
!  9) volume of a mole of fluid, Vmol
! 10) displacement - chemical potential modulus terms
! 11) chemical potential - displacement modulus terms
use global


real(8), INTENT(IN OUT)                   :: props(nprops)
INTEGER, INTENT(IN OUT)                  :: nprops
real(8), INTENT(IN OUT)                   :: dtime
real(8), INTENT(IN OUT)                   :: f_tau(3,3)
real(8), INTENT(IN OUT)                   :: mu_tau
real(8), INTENT(IN OUT)                   :: phi_t
real(8), INTENT(IN OUT)                   :: theta
real(8), INTENT(OUT)                      :: t_tau(3,3)
real(8), INTENT(OUT)                      :: sptanmod(3,3,3,3)
real(8), INTENT(OUT)                      :: phi_tau

INTEGER :: i,j,k,l,m,n, stat
INTEGER, PARAMETER :: nargs=8

real(8) iden(3,3)
real(8)  tr_tau(3,3), dtrdf(3,3,3,3),gshear,kbulk
real(8)  chi,d,mu0,vmol,rgas,detf,finvt(3,3),dpdt
real(8) b_tau(3,3),trb_tau,c_tau(3,3),trc_tau,args(nargs),detfe
real(8) deltamu,dphidmu,dpdt_per,dpdt_m,dphidotdmu,mfluid,finv(3,3)
real(8) phi_per,phi_m, dmdmu,dphidj,spucmod(3,3),dmdphi,dmdj
real(8) spcumodfac(3,3),detfs




! Identity tensor
!
CALL onem0(iden)



! Compute the inverse of F, its determinant, and its transpose
!
CALL matinv3dd(f_tau,finv,detf,stat)
IF(stat == 0) THEN
  WRITE(*,*) 'Problem: detF.lt.zero'
  CALL xit()
END IF
finvt = transpose(finv)


! Compute the left Cauchy-Green tensor and its trace
!
b_tau = matmul(f_tau,transpose(f_tau))
trb_tau = b_tau(1,1) + b_tau(2,2) + b_tau(3,3)


! Compute the right Cauchy-Green tensor and its trace
!
c_tau = matmul(transpose(f_tau),f_tau)
trc_tau = c_tau(1,1) + c_tau(2,2) + c_tau(3,3)


! Compute the Cauchy stress
!
t_tau = (gshear*(b_tau-iden) + kbulk*DLOG(detf)*iden)/detf

! Compute the 1st Piola stress
!
tr_tau = gshear*(f_tau - finvt) + kbulk*DLOG(detf)*finvt


! Compute dTRdF, the so-called material tangent modulus
!
dtrdf = zero
DO i=1,3
  DO j = 1,3
    DO k = 1,3
      DO l = 1,3
        dtrdf(i,j,k,l) = dtrdf(i,j,k,l) + gshear*iden(i,k)*iden(j,l)  &
            + gshear*finv(l,i)*finv(j,k) + kbulk*finv(j,i)*finv(l,k)  &
            - kbulk*DLOG(detf)*finv(l,i)*finv(j,k)
      END DO
    END DO
  END DO
END DO
!
! Calculate the so-called spatial tangent modulus, based
!  on the push forward of the material tangent modulus
!
sptanmod = zero
DO i=1,3
  DO j=1,3
    DO k=1,3
      DO l=1,3
        DO m=1,3
          DO n=1,3
            sptanmod(i,j,k,l) = sptanmod(i,j,k,l) +  &
                (dtrdf(i,m,k,n)*f_tau(j,m)*f_tau(l,n))/detf
          END DO
        END DO
      END DO
    END DO
  END DO
END DO
phi_tau=detf
RETURN
END SUBROUTINE integ
SUBROUTINE invariants(a,inv1,inv2,ndi)



!>    1ST AND 2ND INVARIANTS OF A TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: inv1
DOUBLE PRECISION, INTENT(OUT)            :: inv2



INTEGER :: i1
DOUBLE PRECISION :: aa(ndi,ndi)
DOUBLE PRECISION :: inv1aa

inv1=zero
inv1aa=zero
aa=matmul(a,a)
DO i1=1,ndi
  inv1=inv1+a(i1,i1)
  inv1aa=inv1aa+aa(i1,i1)
END DO
inv2=(one/two)*(inv1*inv1-inv1aa)

RETURN
END SUBROUTINE invariants
SUBROUTINE isomat(sseiso,diso,c10,c01,cbari1,cbari2)



!>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: sseiso
DOUBLE PRECISION, INTENT(OUT)            :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: c10
DOUBLE PRECISION, INTENT(IN)             :: c01
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari1
DOUBLE PRECISION, INTENT(IN OUT)         :: cbari2


sseiso=c10*(cbari1-three)+c01*(cbari2-three)

diso(1)=c10
diso(2)=c01
diso(3)=zero
diso(4)=zero
diso(5)=zero

RETURN
END SUBROUTINE isomat
SUBROUTINE metiso(cmiso,cmfic,pl,pkiso,pkfic,c,unit2,det,ndi)



!>    ISOCHORIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cmiso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi),  &
    plt(ndi,ndi,ndi,ndi),cinv(ndi,ndi), pll(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz, aux,aux1

CALL matinv3d(c,cinv,ndi)
cisoaux1=zero
cisoaux=zero
CALL contraction44(cisoaux1,pl,cmfic,ndi)
DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        plt(i1,j1,k1,l1)=pl(k1,l1,i1,j1)
      END DO
    END DO
  END DO
END DO

CALL contraction44(cisoaux,cisoaux1,plt,ndi)

trfic=zero
aux=det**(-two/three)
aux1=aux**two
DO i1=1,ndi
  trfic=trfic+aux*pkfic(i1,i1)
END DO

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=aux1*cisoaux(i1,j1,k1,l1)
        pll(i1,j1,k1,l1)=(one/two)*(cinv(i1,k1)*cinv(j1,l1)+  &
            cinv(i1,l1)*cinv(j1,k1))- (one/three)*cinv(i1,j1)*cinv(k1,l1)
        yy=trfic*pll(i1,j1,k1,l1)
        zz=pkiso(i1,j1)*cinv(k1,l1)+cinv(i1,j1)*pkiso(k1,l1)
        
        cmiso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metiso
SUBROUTINE metvol(cvol,c,pv,ppv,det,ndi)



!>    VOLUMETRIC MATERIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= det*ppv*cinv(i1,j1)*cinv(k1,l1)  &
            -det*pv*(cinv(i1,k1)*cinv(j1,l1) +cinv(i1,l1)*cinv(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE metvol
SUBROUTINE matinv3d(a,a_inv,ndi)



!>    INVERSE OF A 3X3 MATRIX
!     RETURN THE INVERSE OF A(3,3) - A_INV
use global

INTEGER, INTENT(IN)                  :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: a_inv(ndi,ndi)




DOUBLE PRECISION :: det_a,det_a_inv

det_a = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -  &
    a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) +  &
    a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3))

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
  WRITE(*,*) 'WARNING: DET OF MAT=',det_a
  RETURN
END IF

det_a_inv = 1.d0/det_a

a_inv(1,1) = det_a_inv*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
a_inv(1,2) = det_a_inv*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
a_inv(1,3) = det_a_inv*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
a_inv(2,1) = det_a_inv*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
a_inv(2,2) = det_a_inv*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
a_inv(2,3) = det_a_inv*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
a_inv(3,1) = det_a_inv*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
a_inv(3,2) = det_a_inv*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
a_inv(3,3) = det_a_inv*(a(1,1)*a(2,2)-a(2,1)*a(1,2))

RETURN
END SUBROUTINE matinv3d
SUBROUTINE naffmatnetfic(pkfic,cmfic,f,mf0,rw,filprops,naffprops,  &
        det,ndi)

!>    NON-AFFINE NETWORK:'FICTICIOUS' PK2 STRESS AND ELASTICITY TNSR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cmfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(6)
DOUBLE PRECISION, INTENT(IN)             :: naffprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: h(ndi,ndi),hh(ndi,ndi,ndi,ndi),  &
    hi(ndi,ndi),hhi(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dw,ddw,rwi,aux1,aux2,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda,lambda0,n,pp

!     FILAMENT
l      = filprops(1)
r0      = filprops(2)
mu0     = filprops(3)
beta    = filprops(4)
b0      = filprops(5)
lambda0 = filprops(6)
!     NETWORK
n       = naffprops(1)
pp      = naffprops(2)

lambda=zero
h=zero
hh=zero
hi=zero
hhi=zero
aux=n*(det**(-one))

DO i1=1,nwp
  
  mfi=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  lambda=lambda+(lambdai**pp)*rwi
  
  CALL hfilfic(hi,hhi,pp,lambdai,mf0i,rwi,ndi)
  
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
  
END DO

lambda=lambda**(pp**(-one))

CALL fil(fi,ffi,dw,ddw,lambda,lambda0,l,r0,mu0,beta,b0)

pkfic=zero
cmfic=zero
aux1=aux*dw*lambda**(one-pp)
aux2=aux*(ddw*(lambda**(two*(one-pp)))- (pp-one)*dw*(lambda**(one-two*pp)))

DO j1=1,ndi
  DO k1=1,ndi
    pkfic(j1,k1)=aux1*h(j1,k1)
    DO l1=1,ndi
      DO m1=1,ndi
        cmfic(j1,k1,l1,m1)=aux1*hh(j1,k1,l1,m1)+aux2*h(j1,k1)*h(l1,m1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE naffmatnetfic
SUBROUTINE naffnetfic_bazant(sfic,cfic,f,mf0,rw,filprops,naffprops,  &
        det,ndi)
!
!>    NON-AFFINE NETWORK:'FICTICIUOUS' CAUCHY STRESS AND ELASTICITY TNSR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: naffprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: h(ndi,ndi),hh(ndi,ndi,ndi,ndi),  &
    hi(ndi,ndi),hhi(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dw,ddw,rwi,aux1,aux2,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta,lambda,lambda0,n,pp
DOUBLE PRECISION :: r0f,r0c,etac

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
aux=n*(det**(-one))

r0=r0f+r0c

DO i1=1,nwp
  
  mfi=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  lambda=lambda+(lambdai**pp)*rwi
  
  CALL hfilfic(hi,hhi,pp,lambdai,mfi,rwi,ndi)
  
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
  
END DO

lambda=lambda**(pp**(-one))

CALL fil(fi,ffi,dw,ddw,lambda,lambda0,l,r0,mu0,beta,b0)

cfic=zero
sfic=zero
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

RETURN
END SUBROUTINE naffnetfic_bazant
SUBROUTINE naffnetfic_discrete(sfic,cfic,f,mf0,rw,filprops,naffprops,  &
        det,factor,ndi)
!
!>    NON-AFFINE NETWORK:'FICTICIUOUS' CAUCHY STRESS AND ELASTICITY TNSR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(8)
DOUBLE PRECISION, INTENT(IN)             :: naffprops(2)
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1,l1,m1
DOUBLE PRECISION :: h(ndi,ndi),hh(ndi,ndi,ndi,ndi),  &
    hi(ndi,ndi),hhi(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,lambdai,dw,ddw,rwi,aux1,aux2,fi,ffi
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
SUBROUTINE pk2iso(pkiso,pkfic,pl,det,ndi)



!>    ISOCHORIC PK2 STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pl(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1

DOUBLE PRECISION :: scale2

CALL contraction42(pkiso,pl,pkfic,ndi)

scale2=det**(-two/three)
DO i1=1,ndi
  DO j1=1,ndi
    pkiso(i1,j1)=scale2*pkiso(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2iso
SUBROUTINE pk2isomatfic(fic,diso,cbar,cbari1,unit2,ndi)



!>     ISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
!      INPUT:
!       DISO - STRAIN-ENERGY DERIVATIVES
!       CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
!       CBARI1,CBARI2 - CBAR INVARIANTS
!       UNIT2 - 2ND ORDER IDENTITY TENSOR
!      OUTPUT:
!       FIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: fic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: diso(5)
DOUBLE PRECISION, INTENT(IN)             :: cbar(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cbari1
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1

DOUBLE PRECISION :: dudi1,dudi2
DOUBLE PRECISION :: aux1,aux2

dudi1=diso(1)
dudi2=diso(2)

aux1=two*(dudi1+cbari1*dudi2)
aux2=-two*dudi2

DO i1=1,ndi
  DO j1=1,ndi
    fic(i1,j1)=aux1*unit2(i1,j1)+aux2*cbar(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2isomatfic
SUBROUTINE pk2vol(pkvol,pv,c,ndi)



!>    VOLUMETRIC PK2 STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pkvol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)


INTEGER :: i1,j1
DOUBLE PRECISION :: cinv(ndi,ndi)


CALL matinv3d(c,cinv,ndi)

DO i1=1,ndi
  DO j1=1,ndi
    pkvol(i1,j1)=pv*cinv(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE pk2vol
SUBROUTINE projeul(a,aa,pe,ndi)



!>    EULERIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!      OUTPUTS:
!          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: a(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pe(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l



DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pe(i,j,k,l)=aa(i,j,k,l)-(one/three)*(a(i,j)*a(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projeul
SUBROUTINE projlag(c,aa,pl,ndi)



!>    LAGRANGIAN PROJECTION TENSOR
!      INPUTS:
!          IDENTITY TENSORS - A, AA
!          ISOCHORIC LEFT CAUCHY GREEN TENSOR - C
!          INVERSE OF C - CINV
!      OUTPUTS:
!          4TH ORDER SYMMETRIC LAGRANGIAN PROJECTION TENSOR - PL
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(IN)             :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: aa(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: pl(ndi,ndi,ndi,ndi)



INTEGER :: i,j,k,l

DOUBLE PRECISION :: cinv(ndi,ndi)

CALL matinv3d(c,cinv,ndi)

DO i=1,ndi
  DO j=1,ndi
    DO k=1,ndi
      DO l=1,ndi
        pl(i,j,k,l)=aa(i,j,k,l)-(one/three)*(cinv(i,j)*c(k,l))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE projlag
SUBROUTINE pull2(pk,sig,finv,det,ndi)



!>       PULL-BACK TIMES DET OF A 2ND ORDER TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+det*finv(i1,ii1)*finv(j1,jj1)*sig(ii1,jj1)
      END DO
    END DO
    pk(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE pull2
SUBROUTINE pull4(mat,spatial,finv,det,ndi)



!>        PULL-BACK TIMES DET OF 4TH ORDER TENSOR

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: finv(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: det



INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+det* finv(i1,ii1)*finv(j1,jj1)*  &
                    finv(k1,kk1)*finv(l1,ll1)*spatial(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        mat(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE pull4
SUBROUTINE pullforce(zero0, a, b, machep, t,  &
        lambda,lambda0,l,r0,mu0,beta,b0)



!>    SINGLE FILAMENT: COMPUTES PULLING FORCE FOR A GIVEN STRETCH
!*********************************************************************72

!     ZERO SEEKS THE ROOT OF A FUNCTION F(X) IN AN INTERVAL [A,B].

!     DISCUSSION:

!     THE INTERVAL [A,B] MUST BE A CHANGE OF SIGN INTERVAL FOR F.
!     THAT IS, F(A) AND F(B) MUST BE OF OPPOSITE SIGNS.  THEN
!     ASSUMING THAT F IS CONTINUOUS IMPLIES THE EXISTENCE OF AT LEAST
!     ONE VALUE C BETWEEN A AND B FOR WHICH F(C) = 0.

!     THE LOCATION OF THE ZERO IS DETERMINED TO WITHIN AN ACCURACY
!     OF 6 * MACHEPS * ABS ( C ) + 2 * T.


!     LICENSING:

!     THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.

!     MODIFIED:

!     11 FEBRUARY 2013

!     AUTHOR:

!     RICHARD BRENT
!     MODIFICATIONS BY JOHN BURKARDT

!     REFERENCE:

!     RICHARD BRENT,
!     ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES,
!     DOVER, 2002,
!     ISBN: 0-486-41998-3,
!     LC: QA402.5.B74.

!     PARAMETERS:

!     INPUT, DOUBLE PRECISION A, B, THE ENDPOINTS OF THE CHANGE OF SIGN
!     INTERVAL.
!     INPUT, DOUBLE PRECISION MACHEP, AN ESTIMATE FOR THE RELATIVE
!     MACHINE PRECISION.

!     INPUT, DOUBLE PRECISION T, A POSITIVE ERROR TOLERANCE.

!     INPUT, EXTERNAL DOUBLE PRECISION F, THE NAME OF A USER-SUPPLIED
!     FUNCTION, OF THE FORM "FUNCTION G ( F )", WHICH EVALUATES THE
!     FUNCTION WHOSE ZERO IS BEING SOUGHT.

!     OUTPUT, DOUBLE PRECISION ZERO, THE ESTIMATED VALUE OF A ZERO OF
!     THE FUNCTION G.
use global

DOUBLE PRECISION, INTENT(OUT)            :: zero0
DOUBLE PRECISION, INTENT(IN)             :: a
DOUBLE PRECISION, INTENT(IN)             :: b
DOUBLE PRECISION, INTENT(IN)             :: machep
DOUBLE PRECISION, INTENT(IN)             :: t
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda
DOUBLE PRECISION, INTENT(IN OUT)         :: lambda0
DOUBLE PRECISION, INTENT(IN OUT)         :: l
DOUBLE PRECISION, INTENT(IN OUT)         :: r0
DOUBLE PRECISION, INTENT(IN OUT)         :: mu0
DOUBLE PRECISION, INTENT(IN OUT)         :: beta
DOUBLE PRECISION, INTENT(IN OUT)         :: b0
DOUBLE PRECISION :: c
DOUBLE PRECISION :: d
DOUBLE PRECISION :: e
DOUBLE PRECISION :: fa
DOUBLE PRECISION :: fb
DOUBLE PRECISION :: fc
DOUBLE PRECISION :: m

DOUBLE PRECISION :: p
DOUBLE PRECISION :: q
DOUBLE PRECISION :: r
DOUBLE PRECISION :: s
DOUBLE PRECISION :: sa
DOUBLE PRECISION :: sb

DOUBLE PRECISION :: tol




!     MAKE LOCAL COPIES OF A AND B.

sa = a
sb = b
CALL evalg(fa,sa,lambda,lambda0,l,r0,mu0,beta,b0)
CALL evalg(fb,sb,lambda,lambda0,l,r0,mu0,beta,b0)
!      FA = F ( SA )
!      FB = F ( SB )

10    CONTINUE

c = sa
fc = fa
e = sb - sa
d = e

20    CONTINUE

IF ( ABS ( fc ) < ABS ( fb ) ) THEN
  sa = sb
  sb = c
  c = sa
  fa = fb
  fb = fc
  fc = fa
END IF

30    CONTINUE

tol = 2.0D+00 * machep * ABS ( sb ) + t
m = 0.5D+00 * ( c - sb )
IF ( ABS ( m ) <= tol .OR. fb == 0.0D+00 ) GO TO 140
IF ( ABS ( e ) >= tol .AND. ABS ( fa ) > ABS ( fb ) ) GO TO 40

e = m
d = e
GO TO 100

40    CONTINUE

s = fb / fa
IF ( sa /= c ) GO TO 50

p = 2.0D+00 * m * s
q = 1.0D+00 - s
GO TO 60

50    CONTINUE

q = fa / fc
r = fb / fc
p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

60    CONTINUE

IF ( p <= 0.0D+00 ) GO TO 70

q = - q
GO TO 80

70    CONTINUE

p = - p

80    CONTINUE

s = e
e = d
IF ( 2.0D+00 * p >= 3.0D+00 * m * q - ABS ( tol * q ) .OR.  &
    p >= ABS ( 0.5D+00 * s * q ) ) GO TO 90

d = p / q
GO TO 100

90    CONTINUE

e = m
d = e

100   CONTINUE

sa = sb
fa = fb
IF ( ABS ( d ) <= tol ) GO TO 110
sb = sb + d
GO TO 130

110   CONTINUE

IF ( m <= 0.0D+00 ) GO TO 120
sb = sb + tol
GO TO 130

120   CONTINUE

sb = sb - tol

130   CONTINUE

!      FB = F ( SB )
CALL evalg(fb,sb,lambda,lambda0,l,r0,mu0,beta,b0)
IF ( fb > 0.0D+00 .AND. fc > 0.0D+00 ) GO TO 10
IF ( fb <= 0.0D+00 .AND. fc <= 0.0D+00 ) GO TO 10
GO TO 20

140   CONTINUE

zero0 = sb

RETURN
END SUBROUTINE pullforce

!*********************************************************************72
SUBROUTINE push2(sig,pk,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       PK - 2ND PIOLA KIRCHOOF STRESS TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SIG - CAUCHY STRESS TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sig(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,ii1,jj1


DOUBLE PRECISION :: aux

DO i1=1,ndi
  DO j1=1,ndi
    aux=zero
    DO ii1=1,ndi
      DO jj1=1,ndi
        aux=aux+(det**(-one))*f(i1,ii1)*f(j1,jj1)*pk(ii1,jj1)
      END DO
    END DO
    sig(i1,j1)=aux
  END DO
END DO

RETURN
END SUBROUTINE push2
SUBROUTINE push4(spatial,mat,f,det,ndi)



!>        PIOLA TRANSFORMATION
!>      INPUT:
!>       MAT - MATERIAL ELASTICITY TENSOR
!>       F - DEFORMATION GRADIENT
!>       DET - DEFORMATION DETERMINANT
!>      OUTPUT:
!>       SPATIAL - SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: spatial(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: mat(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det


INTEGER :: i1,j1,k1,l1,ii1,jj1,kk1,ll1


DOUBLE PRECISION :: aux


DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        aux=zero
        DO ii1=1,ndi
          DO jj1=1,ndi
            DO kk1=1,ndi
              DO ll1=1,ndi
                aux=aux+(det**(-one))* f(i1,ii1)*f(j1,jj1)*  &
                    f(k1,kk1)*f(l1,ll1)*mat(ii1,jj1,kk1,ll1)
              END DO
            END DO
          END DO
        END DO
        spatial(i1,j1,k1,l1)=aux
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE push4
SUBROUTINE relax(qv,hv,aux1,hv0,pkiso,dtime,tau,teta,ndi)



!>    VISCOUS DISSIPATION: STRESS RELAXATION TENSORS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: qv(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: hv(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: aux1
DOUBLE PRECISION, INTENT(IN)             :: hv0(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: dtime
DOUBLE PRECISION, INTENT(IN OUT)         :: tau
DOUBLE PRECISION, INTENT(IN)             :: teta



INTEGER :: i1,j1

DOUBLE PRECISION :: aux

qv=zero
hv=zero

aux=DEXP(-dtime*((two*tau)**(-one)))
aux1=teta*aux
DO i1=1,ndi
  DO j1=1,ndi
    qv(i1,j1)=hv0(i1,j1)+aux1*pkiso(i1,j1)
    hv(i1,j1)=aux*(aux*qv(i1,j1)-teta*pkiso(i1,j1))
  END DO
END DO

RETURN
END SUBROUTINE relax
SUBROUTINE rotation(f,r,u,ndi)



!>    COMPUTES ROTATION TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: r(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: u(ndi,ndi)




DOUBLE PRECISION :: uinv(ndi,ndi)

CALL matinv3d(u,uinv,ndi)

r = matmul(f,uinv)
RETURN
END SUBROUTINE rotation
SUBROUTINE sdvread(statev,phi_t,cr)



!>     READ STATE VARS
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)             :: statev(nsdv)
DOUBLE PRECISION, INTENT(OUT)            :: cr,phi_t
!DOUBLE PRECISION, INTENT(IN OUT)         :: dmudx(3,1)


phi_t=statev(1)
cr=statev(2)

RETURN

END SUBROUTINE sdvread
SUBROUTINE sdvwrite(statev,phi_tau,cr,stress)



!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
IMPLICIT NONE

!DOUBLE PRECISION, INTENT(IN OUT)         :: det
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!DOUBLE PRECISION, INTENT(IN OUT)         :: dirmax(3)
DOUBLE PRECISION, INTENT(IN)             :: phi_tau,cr,stress
!DOUBLE PRECISION, INTENT(IN OUT)         :: dmudx(3,1)

statev(1)=phi_tau
statev(2)=cr
statev(3)=stress

RETURN

END SUBROUTINE sdvwrite
SUBROUTINE setiso(ciso,cfic,pe,siso,sfic,unit2,ndi)


use global
IMPLICIT NONE

!>    ISOCHORIC SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: ciso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cfic(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pe(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1
DOUBLE PRECISION :: cisoaux(ndi,ndi,ndi,ndi), cisoaux1(ndi,ndi,ndi,ndi)
DOUBLE PRECISION :: trfic,xx,yy,zz

cisoaux1=zero
cisoaux=zero

CALL contraction44(cisoaux1,pe,cfic,ndi)
CALL contraction44(cisoaux,cisoaux1,pe,ndi)

trfic=zero
DO i1=1,ndi
  trfic=trfic+sfic(i1,i1)
END DO

DO i1=1,ndi
  DO j1=1,ndi
    DO k1=1,ndi
      DO l1=1,ndi
        xx=cisoaux(i1,j1,k1,l1)
        yy=trfic*pe(i1,j1,k1,l1)
        zz=siso(i1,j1)*unit2(k1,l1)+unit2(i1,j1)*siso(k1,l1)
        
        ciso(i1,j1,k1,l1)=xx+(two/three)*yy-(two/three)*zz
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setiso
SUBROUTINE setjr(cjr,sigma,unit2,ndi)


use global
IMPLICIT NONE
!>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cjr(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sigma(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)



INTEGER :: i1,j1,k1,l1


DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        
        cjr(i1,j1,k1,l1)= (one/two)*(unit2(i1,k1)*sigma(j1,l1)  &
            +sigma(i1,k1)*unit2(j1,l1)+unit2(i1,l1)*sigma(j1,k1)  &
            +sigma(i1,l1)*unit2(j1,k1))
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setjr
SUBROUTINE setvol(cvol,pv,ppv,unit2,unit4s,ndi)



!>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: cvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pv
DOUBLE PRECISION, INTENT(IN OUT)         :: ppv
DOUBLE PRECISION, INTENT(IN OUT)         :: unit2(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: unit4s(ndi,ndi,ndi,ndi)


INTEGER :: i1,j1,k1,l1



DO i1 = 1, ndi
  DO j1 = 1, ndi
    DO k1 = 1, ndi
      DO l1 = 1, ndi
        cvol(i1,j1,k1,l1)= ppv*unit2(i1,j1)*unit2(k1,l1)  &
            -two*pv*unit4s(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE setvol
SUBROUTINE sigfilfic(sfic,rho,lambda,dw,m,rw,ndi)



!>    SINGLE FILAMENT:  'FICTICIOUS' CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi !number of dimensions
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi) !ficticious cauchy stress 
DOUBLE PRECISION, INTENT(IN)             :: rho  !angular density at m
DOUBLE PRECISION, INTENT(IN)             :: lambda !filament stretch
DOUBLE PRECISION, INTENT(IN)             :: dw !derivative of filament strain energy
DOUBLE PRECISION, INTENT(IN)             :: m(ndi) !direction vector
DOUBLE PRECISION, INTENT(IN)             :: rw ! integration weights


INTEGER :: i1,j1

DOUBLE PRECISION :: aux

aux=rho*lambda**(-one)*rw*dw
DO i1=1,ndi
  DO j1=1,ndi
    sfic(i1,j1)=aux*m(i1)*m(j1)
  END DO
END DO

RETURN
END SUBROUTINE sigfilfic
SUBROUTINE sigiso(siso,sfic,pe,ndi)



!>    ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: siso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pe(ndi,ndi,ndi,ndi)


CALL contraction42(siso,pe,sfic,ndi)

RETURN
END SUBROUTINE sigiso
SUBROUTINE sigisomatfic(sfic,pkfic,f,det,ndi)



!>    ISOTROPIC MATRIX:  ISOCHORIC CAUCHY STRESS
use global
IMPLICIT NONE


INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(IN OUT)         :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: pkfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: det





CALL push2(sfic,pkfic,f,det,ndi)

RETURN
END SUBROUTINE sigisomatfic
SUBROUTINE signetfic(sfic,f,mf0,rw,filprops,netprops,  &
        rho,lambda0,n,det,ndi)



!>    AFFINE NETWORK:  'FICTICIUOUS' CAUCHY STRESS
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: sfic(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: f(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: mf0(nwp,ndi)
DOUBLE PRECISION, INTENT(IN)             :: rw(nwp)
DOUBLE PRECISION, INTENT(IN)             :: filprops(5)
DOUBLE PRECISION, INTENT(IN)             :: netprops(3)
DOUBLE PRECISION, INTENT(OUT)            :: rho
DOUBLE PRECISION, INTENT(OUT)            :: lambda0
DOUBLE PRECISION, INTENT(OUT)            :: n
DOUBLE PRECISION, INTENT(IN OUT)         :: det



INTEGER :: i1,j1,k1
DOUBLE PRECISION :: sfilfic(ndi,ndi)
DOUBLE PRECISION :: mfi(ndi),mf0i(ndi)
DOUBLE PRECISION :: aux,pi,lambdai,dwi,ddwi,rwi,fi,ffi
DOUBLE PRECISION :: l,r0,mu0,b0,beta

!     FILAMENT
l      = filprops(1)
r0      = filprops(2)
mu0     = filprops(3)
beta    = filprops(4)
b0      = filprops(5)
!     NETWORK
n      = netprops(1)
lambda0 = netprops(2)
rho     = netprops(3)

pi=four*ATAN(one)
aux=n*(det**(-one))*four*pi
sfic=zero

DO i1=1,nwp
  
  mfi=zero
  DO j1=1,ndi
    mf0i(j1)=mf0(i1,j1)
  END DO
  rwi=rw(i1)
  
  
  CALL deffil(lambdai,mfi,mf0i,f,ndi)
  
  CALL fil(fi,ffi,dwi,ddwi,lambdai,lambda0,l,r0,mu0,beta,b0)
  
  CALL sigfilfic(sfilfic,rho,lambdai,dwi,mfi,rwi,ndi)
  
  DO j1=1,ndi
    DO k1=1,ndi
      sfic(j1,k1)=sfic(j1,k1)+aux*sfilfic(j1,k1)
    END DO
  END DO
  
END DO

RETURN
END SUBROUTINE signetfic
SUBROUTINE sigvol(svol,pv,unit2,ndi)



!>    VOLUMETRIC CAUCHY STRESS

use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: svol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pv
DOUBLE PRECISION, INTENT(IN)             :: unit2(ndi,ndi)



INTEGER :: i1,j1



DO i1=1,ndi
  DO j1=1,ndi
    svol(i1,j1)=pv*unit2(i1,j1)
  END DO
END DO

RETURN
END SUBROUTINE sigvol
SUBROUTINE sliding(ffc,ru,ffc0,ru0,ffcmax,fric,frac,dtime)



use global
IMPLICIT NONE


DOUBLE PRECISION, INTENT(OUT)            :: ffc
DOUBLE PRECISION, INTENT(OUT)            :: ru
DOUBLE PRECISION, INTENT(IN)             :: ffc0
DOUBLE PRECISION, INTENT(IN OUT)         :: ru0
DOUBLE PRECISION, INTENT(IN)             :: ffcmax
DOUBLE PRECISION, INTENT(IN OUT)         :: fric
DOUBLE PRECISION, INTENT(IN)             :: frac(4)
DOUBLE PRECISION, INTENT(IN)             :: dtime




DOUBLE PRECISION :: aux0,aux1,aux2, arg
!      INTEGER STATE

aux1=frac(3)
aux0=aux1+frac(4)
aux2=aux0
arg=ffc0/ffcmax

IF(arg < aux1) THEN
  ffc=aux1*ffcmax
ELSE IF (arg > aux2)THEN
  ffc=aux2*ffcmax
ELSE
  ffc=ffc0
END IF

ru=ru0+dtime*(fric**(-one))*(ffc-ffc0)
ru0=ru

RETURN

END SUBROUTINE sliding
SUBROUTINE solvephi(root,args,nargs,rootold)



! This subroutine will numerically solve for the polymer
!  volume fraction based on the current osmotic pressure
!  and the previous state.  See numerical recipies RTSAFE.

use global

real(8), INTENT(OUT)                      :: root
real(8), INTENT(IN)                       :: args(nargs)
INTEGER, INTENT(IN OUT)                  :: nargs
real(8), INTENT(IN OUT)                   :: rootold

INTEGER :: j

real(8) f,df,fl,fh,xl,xh,x1,x2,swap, dxold
real(8) dx,  temp,rootmax,rootmin

INTEGER, PARAMETER :: maxit=50
real(8), PARAMETER :: xacc=1.d-6

rootmax = 0.9999D0 ! corresponds to nearly 100% dry polymer
rootmin = 0.05D0   ! corresponds to nearly 100% fluid

x1 = rootmin
x2 = rootmax
CALL phifunc(x1,fl,df,args,nargs)
CALL phifunc(x2,fh,df,args,nargs)

IF(fl*fh >= zero) THEN
  root = rootold
  WRITE(*,*) 'FYI, root not bracketed on phi'
  WRITE(*,*) 'fl=',fl
  WRITE(*,*) 'fh=',fh
  WRITE(*,*) 'rootOld=',rootold
!write(80,*) 'FYI, the root is not bracketed on phi'
!write(80,*) 'fl=',fl
!write(80,*) 'fh=',fh
!write(80,*) 'rootOld=',rootOld
  
  WRITE(*,*) 'mu =',args(1)
  WRITE(*,*) 'mu0=',args(2)
  WRITE(*,*) 'Rgas=',args(3)
  WRITE(*,*) 'theta=',args(4)
  WRITE(*,*) 'chi=',args(5)
  WRITE(*,*) 'Vmol=',args(6)
  WRITE(*,*) 'Kbulk=',args(7)
  WRITE(*,*) 'detF=',args(8)
  
  CALL xit()
  RETURN
END IF


!  ORIENT THE SEARCH SO THAT F(XL) < 0.

IF( fl < 0.d0 ) THEN
  xl = x1
  xh = x2
ELSE
  xh = x1
  xl = x2
  swap = fl
  fl = fh
  fh = swap
END IF

!  INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
!  BEFORE LAST'', AND THE LAST STEP

IF(rootold < rootmin) rootold = rootmin
IF(rootold > rootmax) rootold = rootmax
root = rootold !0.5D0 *( X1 + X2)
dxold = DABS(x2 - x1)
dx = dxold

CALL phifunc(root,f,df,args,nargs)


!   LOOP OVER ALLOWED ITERATIONS

DO  j = 1,maxit
  
!   BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
!   FAST ENOUGH.
  
  IF( ((root-xh)*df - f)*((root - xl)*df -f) >= 0.d0  &
        .OR. DABS(2.d0*f) > DABS(dxold*df) ) THEN
    
    dxold = dx
    dx = 0.5D0*(xh-xl)
    root = xl + dx
    IF( xl == root ) THEN
      
!   CHANGE IN ROOT IS NEGLIGIBLE
      
      RETURN
    END IF
    
  ELSE
    
!   NEWTON STEP IS ACCEPTABLE. TAKE IT.
    
    dxold = dx
    dx = f/df
    temp = root
    root = root - dx
    IF( temp == root) THEN
      
!    CHANGE IN ROOT IS NEGLIGIBLE
      
      RETURN
    END IF
    
  END IF
  
!  CONVERVEGENCE CRITERION
  
  IF( DABS(dx) < xacc) RETURN
  
  
!   THE ONE NEW FUNCTION EVALUATION PER ITERATION
  
  CALL phifunc(root,f,df,args,nargs)
  
  
!  MAINTAIN THE BRACKET ON THE ROOT
  
  IF( f < 0.d0) THEN
    xl = root
    fl = f
  ELSE
    xh = root
    fh = f
  END IF
  
END DO

WRITE(*,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'
WRITE(80,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'

RETURN
END SUBROUTINE solvephi

!***************************************************************************
SUBROUTINE spectral(a,d,v)



!>    EIGENVALUES AND EIGENVECTOR OF A 3X3 MATRIX
!     THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
!     A SYMMETRIC 3X3 MATRIX A.

!     THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
!     EIGENVALUES IN ASCENDING ORDER, AND A MATRIX V WHOSE
!     COLUMNS CONTAIN THE CORRESPONDING EIGENVECTORS.

use global

DOUBLE PRECISION, INTENT(IN OUT)         :: a(3,3)
DOUBLE PRECISION                         :: e(3,3)
DOUBLE PRECISION, INTENT(IN OUT)         :: d(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(3,3)

INTEGER :: nrot
INTEGER :: np=3



e = a

CALL jacobi(e,3,np,d,v,nrot)
CALL eigsrt(d,v,3,np)

RETURN
END SUBROUTINE spectral

!***********************************************************************

SUBROUTINE jacobi(a,n,np,d,v,nrot)

! COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
!  MATRIX A, WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL
!  NP BY NP ARRAY.  ON OUTPUT, ELEMENTS OF A ABOVE THE DIAGONAL
!  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
!  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
!  VECTOR D RETURNS THE EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
!  V IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
!  A WHOSE COLUMNS CONTAIN, UPON OUTPUT, THE NORMALIZED
!  EIGENVECTORS OF A.  NROT RETURNS THE NUMBER OF JACOBI ROTATION
!  WHICH WERE REQUIRED.

! THIS SUBROUTINE IS TAKEN FROM 'NUMERICAL RECIPES.'
use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: np
DOUBLE PRECISION, INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(OUT)            :: d(np)
DOUBLE PRECISION, INTENT(OUT)            :: v(np,np)
INTEGER, INTENT(OUT)                     :: nrot

INTEGER :: ip,iq,  i,j
INTEGER, PARAMETER :: nmax=100

DOUBLE PRECISION :: b(nmax),z(nmax), sm,tresh,g,t,h,theta,s,c,tau


! INITIALIZE V TO THE IDENTITY MATRIX
DO i=1,3
  v(i,i)=one
  DO j=1,3
    IF (i /= j)THEN
      v(i,j)=zero
    END IF
  END DO
END DO
! INITIALIZE B AND D TO THE DIAGONAL OF A, AND Z TO ZERO.
!  THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
!  IN EQUATION (11.1.14)

DO ip = 1,n
  b(ip) = a(ip,ip)
  d(ip) = b(ip)
  z(ip) = 0.d0
END DO


! BEGIN ITERATION

nrot = 0
DO i=1,50
  
!         SUM OFF-DIAGONAL ELEMENTS
  
  sm = 0.d0
  DO ip=1,n-1
    DO iq=ip+1,n
      sm = sm + DABS(a(ip,iq))
    END DO
  END DO
  
!          IF SM = 0., THEN RETURN.  THIS IS THE NORMAL RETURN,
!          WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE
!          UNDERFLOW.
  
  IF (sm == 0.d0) RETURN
  
!          IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
!           |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE,
!           SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
  
  IF (i < 4) THEN
    tresh = 0.2D0*sm/n**2
  ELSE
    tresh = 0.d0
  END IF
  
  DO ip=1,n-1
    DO iq=ip+1,n
      g = 100.d0*DABS(a(ip,iq))
      
!              AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
!               OFF-DIAGONAL ELEMENT IS SMALL.
      
      IF ((i > 4).AND.(DABS(d(ip))+g == DABS(d(ip)))  &
            .AND.(DABS(d(iq))+g == DABS(d(iq)))) THEN
        a(ip,iq) = 0.d0
      ELSE IF (DABS(a(ip,iq)) > tresh) THEN
        h = d(iq) - d(ip)
        IF (DABS(h)+g == DABS(h)) THEN
          
!                  T = 1./(2.*THETA), EQUATION (11.1.10)
          
          t =a(ip,iq)/h
        ELSE
          theta = 0.5D0*h/a(ip,iq)
          t =1.d0/(DABS(theta)+DSQRT(1.d0+theta**2.d0))
          IF (theta < 0.d0) t = -t
        END IF
        c = 1.d0/DSQRT(1.d0 + t**2.d0)
        s = t*c
        tau = s/(1.d0 + c)
        h = t*a(ip,iq)
        z(ip) = z(ip) - h
        z(iq) = z(iq) + h
        d(ip) = d(ip) - h
        d(iq) = d(iq) + h
        a(ip,iq) = 0.d0
        
!               CASE OF ROTATIONS 1 <= J < P
        
        DO j=1,ip-1
          g = a(j,ip)
          h = a(j,iq)
          a(j,ip) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                CASE OF ROTATIONS P < J < Q
        
        DO j=ip+1,iq-1
          g = a(ip,j)
          h = a(j,iq)
          a(ip,j) = g - s*(h + g*tau)
          a(j,iq) = h + s*(g - h*tau)
        END DO
        
!                 CASE OF ROTATIONS Q < J <= N
        
        DO j=iq+1,n
          g = a(ip,j)
          h = a(iq,j)
          a(ip,j) = g - s*(h + g*tau)
          a(iq,j) = h + s*(g - h*tau)
        END DO
        DO j = 1,n
          g = v(j,ip)
          h = v(j,iq)
          v(j,ip) = g - s*(h + g*tau)
          v(j,iq) = h + s*(g - h*tau)
        END DO
        nrot = nrot + 1
      END IF
    END DO
  END DO
  
!          UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
  
  DO ip=1,n
    b(ip) = b(ip) + z(ip)
    d(ip) = b(ip)
    z(ip) = 0.d0
  END DO
END DO

! IF THE ALGORITHM HAS REACHED THIS STAGE, THEN THERE
!  ARE TOO MANY SWEEPS.  PRINT A DIAGNOSTIC AND CUT THE
!  TIME INCREMENT.

WRITE (*,'(/1X,A/)') '50 ITERATIONS IN JACOBI SHOULD NEVER HAPPEN'

RETURN
END SUBROUTINE jacobi

!**********************************************************************

SUBROUTINE eigsrt(d,v,n,np)

!     GIVEN THE EIGENVALUES D AND EIGENVECTORS V AS OUTPUT FROM
!     JACOBI, THIS SUBROUTINE SORTS THE EIGENVALUES INTO ASCENDING
!     ORDER AND REARRANGES THE COLMNS OF V ACCORDINGLY.

!     THE SUBROUTINE WAS TAKEN FROM 'NUMERICAL RECIPES.'
use global

DOUBLE PRECISION, INTENT(IN OUT)         :: d(np)
DOUBLE PRECISION, INTENT(IN OUT)         :: v(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: np


INTEGER :: i,j,k

DOUBLE PRECISION :: p

DO i=1,n-1
  k = i
  p = d(i)
  DO j=i+1,n
    IF (d(j) >= p) THEN
      k = j
      p = d(j)
    END IF
  END DO
  IF (k /= i) THEN
    d(k) = d(i)
    d(i) = p
    DO j=1,n
      p = v(j,i)
      v(j,i) = v(j,k)
      v(j,k) = p
    END DO
  END IF
END DO

RETURN
END SUBROUTINE eigsrt

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
SUBROUTINE stretch(c,b,u,v,ndi)



!>    STRETCH TENSORS

use global
IMPLICIT NONE

INTEGER, INTENT(IN OUT)                  :: ndi

DOUBLE PRECISION, INTENT(IN OUT)         :: c(ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: b(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: u(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: v(ndi,ndi)




DOUBLE PRECISION :: eigval(ndi),omega(ndi),eigvec(ndi,ndi)

CALL spectral(c,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))

u(1,1) = eigval(1)
u(2,2) = eigval(2)
u(3,3) = eigval(3)

u = matmul(matmul(eigvec,u),transpose(eigvec))

CALL spectral(b,omega,eigvec)

eigval(1) = DSQRT(omega(1))
eigval(2) = DSQRT(omega(2))
eigval(3) = DSQRT(omega(3))
!      write(*,*) eigvec(1,1),eigvec(2,1),eigvec(3,1)

v(1,1) = eigval(1)
v(2,2) = eigval(2)
v(3,3) = eigval(3)

v = matmul(matmul(eigvec,v),transpose(eigvec))
RETURN
END SUBROUTINE stretch
SUBROUTINE uexternaldb(lop,lrestart,time,dtime,kstep,kinc)



!>    READ FILAMENTS ORIENTATION AND PREFERED DIRECTIONS
use global
INCLUDE 'aba_param.inc'
!       this subroutine get the directions and weights for
!      the numerical integration

!     UEXTERNAL just called once; work in parallel computing

INTEGER, INTENT(IN OUT)                  :: lop
INTEGER, INTENT(IN OUT)                  :: lrestart
REAL, INTENT(IN OUT)                     :: time(2)
real(8), INTENT(IN OUT)                   :: dtime
INTEGER, INTENT(IN OUT)                  :: kstep
INTEGER, INTENT(IN OUT)                  :: kinc

COMMON /kfil/mf0
COMMON /kfilr/rw
COMMON /kfilp/prefdir
COMMON /kinit/init


DOUBLE PRECISION :: sphere(nwp,5),mf0(nwp,3),rw(nwp),  &
                    prefdir(nelem,4), init(2)
CHARACTER (LEN=256) ::  filename, jobdir
INTEGER :: lenjobdir,i,j,k

!     LOP=0 --> START OF THE ANALYSIS
IF(lop == 0.OR.lop == 4) THEN
  
  CALL getoutdir(jobdir,lenjobdir)
  
  !direction and weights of bazant integration 
  filename=jobdir(:lenjobdir)//'/'//dir1
  OPEN(15,FILE=filename)
  DO i=1,nwp
    READ(15,*) (sphere(i,j),j=1,5)
  END DO
  CLOSE(15)
  
  DO i=1,nwp
    DO j=1,3
      k=j+1
      mf0(i,j)=sphere(i,k)
    END DO
!         ORTHOGONAL SYMMETRIES
    rw(i)=two*sphere(i,5)
!        NO ORTHOGONAL SYMMETRIES
!          RW(I)=SPHERE(I,5)
  END DO
  
  !preferential direction
  filename=jobdir(:lenjobdir)//'/'//dir2
  OPEN(16,FILE=filename)
  DO i=1,nelem
    READ(16,*) (prefdir(i,j),j=1,4)
  END DO
  CLOSE(16)
  !initial conditions
  filename=jobdir(:lenjobdir)//'/'//dir3
  OPEN(18,FILE=filename)
  READ(18,*) (init(j),j=1,2)
  CLOSE(18)
  
END IF

RETURN

END SUBROUTINE uexternaldb
!****************************************************************************



!     Utility subroutines
!****************************************************************************
!***************************************************************************

SUBROUTINE matinv3dd(a,a_inv,det_a,istat)
!
! Returns A_inv, the inverse and det_A, the determinant
! Note that the det is of the original matrix, not the
! inverse
!
use global

real(8), INTENT(IN)                       :: a(3,3)
real(8), INTENT(OUT)                      :: a_inv(3,3)
real(8), INTENT(OUT)                      :: det_a
INTEGER, INTENT(OUT)                     :: istat

!

!
real(8)  det_a_inv
!
istat = 1

det_a = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -  &
    a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) +  &
    a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3))

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: subroutine matInv3Dd:'
  WRITE(*,*) 'WARNING: det of mat=',det_a
  istat = 0
  RETURN
END IF
!
det_a_inv = 1.d0/det_a
!
a_inv(1,1) = det_a_inv*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
a_inv(1,2) = det_a_inv*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
a_inv(1,3) = det_a_inv*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
a_inv(2,1) = det_a_inv*(a(3,1)*a(2,3)-a(2,1)*a(3,3))
a_inv(2,2) = det_a_inv*(a(1,1)*a(3,3)-a(3,1)*a(1,3))
a_inv(2,3) = det_a_inv*(a(2,1)*a(1,3)-a(1,1)*a(2,3))
a_inv(3,1) = det_a_inv*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
a_inv(3,2) = det_a_inv*(a(3,1)*a(1,2)-a(1,1)*a(3,2))
a_inv(3,3) = det_a_inv*(a(1,1)*a(2,2)-a(2,1)*a(1,2))
!
RETURN
END SUBROUTINE matinv3dd
!!!

SUBROUTINE matinv2d(a,a_inv,det_a,istat)
!
! Returns A_inv, the inverse, and det_A, the determinant
! Note that the det is of the original matrix, not the
! inverse
!
use global

real(8), INTENT(IN)                       :: a(2,2)
real(8), INTENT(OUT)                      :: a_inv(2,2)
real(8), INTENT(OUT)                      :: det_a
INTEGER, INTENT(OUT)                     :: istat

!

!
real(8)  det_a_inv


istat = 1

det_a = a(1,1)*a(2,2) - a(1,2)*a(2,1)

IF (det_a <= 0.d0) THEN
  WRITE(*,*) 'WARNING: subroutine matInv2D:'
  WRITE(*,*) 'WARNING: det of mat=',det_a
  istat = 0
  RETURN
END IF

det_a_inv = 1.d0/det_a

a_inv(1,1) =  det_a_inv*a(2,2)
a_inv(1,2) = -det_a_inv*a(1,2)
a_inv(2,1) = -det_a_inv*a(2,1)
a_inv(2,2) =  det_a_inv*a(1,1)


RETURN
END SUBROUTINE matinv2d

!****************************************************************************

SUBROUTINE mdet(a,det)
!
! This subroutine calculates the determinant
! of a 3 by 3 matrix [A]
!
use global

real(8), INTENT(IN)                       :: a(3,3)
real(8), INTENT(OUT)                      :: det

!



det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)  &
    + a(1,3)*a(2,1)*a(3,2) - a(3,1)*a(2,2)*a(1,3)  &
    - a(3,2)*a(2,3)*a(1,1) - a(3,3)*a(2,1)*a(1,2)


RETURN
END SUBROUTINE mdet

!****************************************************************************

SUBROUTINE onem0(a)
!
! This subroutine stores the identity matrix in the
! 3 by 3 matrix [A]
!
use global

real(8), INTENT(OUT)                      :: a(3,3)

!
INTEGER :: i,j
!



DO i=1,3
  DO j=1,3
    IF (i == j) THEN
      a(i,j) = 1.0
    ELSE
      a(i,j) = 0.0
    END IF
  END DO
END DO


RETURN
END SUBROUTINE onem0
!***************************************************************************
subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname, &
        nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord, &
        jmac,jmatyp,matlayo,laccfla)

! this subroutine is used to transfer sdv's from the uel
!  onto the dummy mesh for viewing.  note that an offset of
!  elemoffset is used between the real mesh and the dummy mesh.
!  if your model has more than elemoffset uel elements, then
!  this will need to be modified.

use global

include 'aba_param.inc'

character*80 cmname,orname
character*3 flgray(15)
dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)
integer i1

!     the dimensions of the variables flgray, array and jarray
!     must be set equal to or greater than 15.

!      uvar(1) = globalsdv(noel-elemoffset,npt,1)
!      for example
!      uvar(2) = globalsdv(noel-elemoffset,npt,2)
do i1=1,nsdv
uvar(i1) = globalsdv(noel-elemoffset,npt,i1)
enddo
!      uvar(3) = globalsdv(noel-elemoffset,npt,3)
!      uvar(4) = globalsdv(noel-elemoffset,npt,4)

return
end subroutine uvarm
 
SUBROUTINE visco(pk,cmat,vv,pkvol,pkiso,cmatvol,cmatiso,dtime,  &
        vscprops,statev,ndi)



!>    VISCOUS DISSIPATION: MAXWELL SPRINGS AND DASHPOTS SCHEME
use global
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: pk(ndi,ndi)
DOUBLE PRECISION, INTENT(OUT)            :: cmat(ndi,ndi,ndi,ndi)
INTEGER, INTENT(IN)                      :: vv
DOUBLE PRECISION, INTENT(IN)             :: pkvol(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: pkiso(ndi,ndi)
DOUBLE PRECISION, INTENT(IN)             :: cmatvol(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: cmatiso(ndi,ndi,ndi,ndi)
DOUBLE PRECISION, INTENT(IN OUT)         :: dtime
DOUBLE PRECISION, INTENT(IN)             :: vscprops(6)
DOUBLE PRECISION, INTENT(IN OUT)         :: statev(nsdv)



INTEGER :: i1,j1,k1,l1, v1

DOUBLE PRECISION :: q(ndi,ndi),qv(ndi,ndi),hv(ndi,ndi), hv0(ndi,ndi)
DOUBLE PRECISION :: teta,tau,aux,auxc

q=zero
qv=zero
hv=zero
auxc=zero

!     ( GENERAL MAXWELL DASHPOTS)
DO v1=1,vv
  
  tau=vscprops(2*v1-1)
  teta=vscprops(2*v1)
  
!      READ STATE VARIABLES
  CALL hvread(hv,statev,v1,ndi)
  hv0=hv
!        RALAXATION TENSORS
  CALL relax(qv,hv,aux,hv0,pkiso,dtime,tau,teta,ndi)
  auxc=auxc+aux
!        WRITE STATE VARIABLES
  CALL hvwrite(statev,hv,v1,ndi)
  
  q=q+qv
  
END DO

auxc=one+auxc
pk=pkvol+pkiso


DO i1=1,ndi
  DO j1=1,ndi
    pk(i1,j1)=pk(i1,j1)+q(i1,j1)
    DO k1=1,ndi
      DO l1=1,ndi
        cmat(i1,j1,k1,l1)= cmatvol(i1,j1,k1,l1)+ auxc*cmatiso(i1,j1,k1,l1)
      END DO
    END DO
  END DO
END DO



RETURN
END SUBROUTINE visco
SUBROUTINE vol(ssev,pv,ppv,k,det,dete,dets)



!>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
use global
IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: ssev
DOUBLE PRECISION, INTENT(OUT)            :: pv
DOUBLE PRECISION, INTENT(OUT)            :: ppv
DOUBLE PRECISION, INTENT(IN)             :: k
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN OUT)         :: dete
DOUBLE PRECISION, INTENT(IN OUT)         :: dets



DOUBLE PRECISION :: g, aux

g=(one/four)*(det*det-one-two*log(det))
ssev=k*g

!no sweelling
pv=k*(one/two)*(det-one/det)
aux=k*(one/two)*(one+one/(det*det))
ppv=pv+det*aux
dete=det
dets=det

!sweelling
!pv=dets*k*(dlog(dete))/det
!aux=dets*k*(one/two)*(one+one/(dete*dete))
!ppv=pv+dete*aux

RETURN
END SUBROUTINE vol
SUBROUTINE xit()



CALL EXIT()

END SUBROUTINE
