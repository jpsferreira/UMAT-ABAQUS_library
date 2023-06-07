C>********************************************************************
C> Record of revisions:                                              |
C>        Date        Programmer        Description of change        |
C>        ====        ==========        =====================        |
C>     05/11/2016    Joao Ferreira      full network model           |
C>--------------------------------------------------------------------
C>     Description:
C>     UMAT: USER MATERIAL FOR THE FULL NETWORK MODEL.
C>                AFFINE DEFORMATIONS
C>     UEXTERNALDB: READ FILAMENTS ORIENTATION AND PREFERED DIRECTION
C>--------------------------------------------------------------------
C>---------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
C----------------------------------------------------------------------
C--------------------------- DECLARATIONS -----------------------------
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C     FILAMENTS DIRECTION
      COMMON /KFIL/MF0
C     FILAMENTS WEIGHT      
      COMMON /KFILR/RW
C     PREFERED DIRETION
      COMMON /KFILP/PREFDIR
C
      DOUBLE PRECISION PREFDIR(NELEM,4)
      CHARACTER*8 CMNAME
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1        LAYER, KSPT, KSTEP, KINC
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1                 DTEMP,PNEWDT,CELENT
C
C     FLAGS
C      INTEGER FLAG1
C     UTILITY TENSORS
      DOUBLE PRECISION UNIT2(NDI,NDI),UNIT4(NDI,NDI,NDI,NDI),
     1                 UNIT4S(NDI,NDI,NDI,NDI),
     2                 PROJE(NDI,NDI,NDI,NDI),PROJL(NDI,NDI,NDI,NDI)
C     KINEMATICS
      DOUBLE PRECISION DISTGR(NDI,NDI),C(NDI,NDI),B(NDI,NDI),
     1                 CBAR(NDI,NDI),BBAR(NDI,NDI),DISTGRINV(NDI,NDI),
     2                 UBAR(NDI,NDI),VBAR(NDI,NDI),ROT(NDI,NDI),
     3                 DFGRD1INV(NDI,NDI)
      DOUBLE PRECISION DET,CBARI1,CBARI2
C     VOLUMETRIC CONTRIBUTION
      DOUBLE PRECISION PKVOL(NDI,NDI),SVOL(NDI,NDI),
     1                 CVOL(NDI,NDI,NDI,NDI),CMVOL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION K,PV,PPV,SSEV
C     ISOCHORIC CONTRIBUTION
      DOUBLE PRECISION SISO(NDI,NDI),PKISO(NDI,NDI),PK2(NDI,NDI),
     1                 CISO(NDI,NDI,NDI,NDI),CMISO(NDI,NDI,NDI,NDI),
     2                 SFIC(NDI,NDI),CFIC(NDI,NDI,NDI,NDI),
     3                 PKFIC(NDI,NDI),CMFIC(NDI,NDI,NDI,NDI)    
C     ISOCHORIC ISOTROPIC CONTRIBUTION
      DOUBLE PRECISION C10,C01,SSEISO,DISO(5),PKMATFIC(NDI,NDI),
     1                 SMATFIC(NDI,NDI),SISOMATFIC(NDI,NDI),
     2                 CMISOMATFIC(NDI,NDI,NDI,NDI),
     3                 CISOMATFIC(NDI,NDI,NDI,NDI)
C     FILAMENTS NETWORK CONTRIBUTION
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),FILPROPS(6),
     1                 NAFFPROPS(2),AFFPROPS(2)
      DOUBLE PRECISION LL,R0,LAMBDA0,MU0,BETA,NN,MM,B0,BB,PHI,P
      DOUBLE PRECISION PKNETFIC(NDI,NDI),CMNETFIC(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION SNETFIC(NDI,NDI),CNETFIC(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION PKNETFICAF(NDI,NDI),PKNETFICNAF(NDI,NDI)
      DOUBLE PRECISION SNETFICAF(NDI,NDI),SNETFICNAF(NDI,NDI)
      DOUBLE PRECISION CMNETFICAF(NDI,NDI,NDI,NDI),
     1                 CMNETFICNAF(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION CNETFICAF(NDI,NDI,NDI,NDI),
     1                 CNETFICNAF(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION EFI
      INTEGER NTERM
C
C     VISCOUS PROPERTIES (GENERALIZED MAXWEL DASHPOTS)
      DOUBLE PRECISION VSCPROPS(6)
      INTEGER VV   
C     JAUMMAN RATE CONTRIBUTION (REQUIRED FOR ABAQUS UMAT)
      DOUBLE PRECISION CJR(NDI,NDI,NDI,NDI)
C     CAUCHY STRESS AND ELASTICITY TENSOR
      DOUBLE PRECISION SIGMA(NDI,NDI),DDSIGDDE(NDI,NDI,NDI,NDI),
     1                                 DDPKDDE(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION STEST(NDI,NDI), CTEST(NDI,NDI,NDI,NDI)
C----------------------------------------------------------------------
C-------------------------- INITIALIZATIONS ---------------------------
C----------------------------------------------------------------------
C     IDENTITY AND PROJECTION TENSORS
      UNIT2=ZERO
      UNIT4=ZERO
      UNIT4S=ZERO
      PROJE=ZERO
      PROJL=ZERO
C     KINEMATICS
      DISTGR=ZERO
      C=ZERO
      B=ZERO
      CBAR=ZERO
      BBAR=ZERO
      UBAR=ZERO
      VBAR=ZERO
      ROT=ZERO
      DET=ZERO
      CBARI1=ZERO
      CBARI2=ZERO
C     VOLUMETRIC
      PKVOL=ZERO
      SVOL=ZERO
      CVOL=ZERO
      K=ZERO
      PV=ZERO
      PPV=ZERO
      SSEV=ZERO
C     ISOCHORIC
      SISO=ZERO
      PKISO=ZERO
      PK2=ZERO
      CISO=ZERO
      CFIC=ZERO
      SFIC=ZERO
      PKFIC=ZERO
C     ISOTROPIC
      C10=ZERO
      C01=ZERO
      SSEISO=ZERO
      DISO=ZERO
      PKMATFIC=ZERO
      SMATFIC=ZERO
      SISOMATFIC=ZERO
      CMISOMATFIC=ZERO
      CISOMATFIC=ZERO
C     FILAMENTS NETWORK
      SNETFIC=ZERO
      CNETFIC=ZERO
      PKNETFIC=ZERO
      PKNETFICAF=ZERO
      PKNETFICNAF=ZERO
      SNETFICAF=ZERO
      SNETFICNAF=ZERO
      CMNETFIC=ZERO
      CMNETFICAF=ZERO
      CMNETFICNAF=ZERO
      CNETFICAF=ZERO
      CNETFICNAF=ZERO
C     JAUMANN RATE
      CJR=ZERO
C     TOTAL CAUCHY STRESS AND ELASTICITY TENSORS
      SIGMA=ZERO
      DDSIGDDE=ZERO   
C----------------------------------------------------------------------
C------------------------ IDENTITY TENSORS ----------------------------
C----------------------------------------------------------------------
            CALL ONEM(UNIT2,UNIT4,UNIT4S,NDI)
C----------------------------------------------------------------------
C------------------- MATERIAL CONSTANTS AND DATA ----------------------
C----------------------------------------------------------------------
C     VOLUMETRIC
      K        = PROPS(1)
C     ISOCHORIC ISOTROPIC
      C10      = PROPS(2)
      C01      = PROPS(3)
      PHI      = PROPS(4)
C     FILAMENT
      LL       = PROPS(5)
      R0       = PROPS(6)
      MU0      = PROPS(7)
      BETA     = PROPS(8)
      B0       = PROPS(9)
      LAMBDA0  = PROPS(10)
      FILPROPS = PROPS(5:10)
C     NETWORK
      NN       = PROPS(11)
      BB        = PROPS(12)   
      AFFPROPS = PROPS(11:12)     
C
C     VISCOUS EFFECTS
      VV       = INT(PROPS(13))
      VSCPROPS = PROPS(14:19)
C     NUMERICAL COMPUTATIONS
      NTERM    = 60      
C 
C        STATE VARIABLES AND CHEMICAL PARAMETERS
C
      IF ((TIME(1).EQ.ZERO).AND.(KSTEP.EQ.1)) THEN
      CALL INITIALIZE(STATEV,VV)   
      ENDIF
C        READ STATEV
      CALL SDVREAD(STATEV,VV)       
C----------------------------------------------------------------------
C---------------------------- KINEMATICS ------------------------------
C----------------------------------------------------------------------
C     DISTORTION GRADIENT
      CALL FSLIP(DFGRD1,DISTGR,DET,NDI)
C     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DFGRD1,DFGRD1INV,NDI)      
C     INVERSE OF DISTORTION GRADIENT
      CALL MATINV3D(DISTGR,DISTGRINV,NDI)
C     CAUCHY-GREEN DEFORMATION TENSORS
      CALL DEFORMATION(DFGRD1,C,B,NDI)
      CALL DEFORMATION(DISTGR,CBAR,BBAR,NDI)
C     INVARIANTS OF DEVIATORIC DEFORMATION TENSORS
      CALL INVARIANTS(CBAR,CBARI1,CBARI2,NDI)
C     STRETCH TENSORS
      CALL STRETCH(CBAR,BBAR,UBAR,VBAR,NDI)
C     ROTATION TENSORS
      CALL ROTATION(DISTGR,ROT,UBAR,NDI)
C----------------------------------------------------------------------
C--------------------- CONSTITUTIVE RELATIONS  ------------------------
C----------------------------------------------------------------------
C     DEVIATORIC PROJECTION TENSORS
      CALL PROJEUL(UNIT2,UNIT4S,PROJE,NDI)
C
      CALL PROJLAG(C,UNIT4,PROJL,NDI)
C
C---- VOLUMETRIC ------------------------------------------------------
C     STRAIN-ENERGY
      CALL VOL(SSEV,PV,PPV,K,DET)
C
C---- ISOCHORIC ISOTROPIC ---------------------------------------------
C     STRAIN-ENERGY
      CALL ISOMAT(SSEISO,DISO,C10,C01,CBARI1,CBARI2)
C     PK2 'FICTICIOUS' STRESS TENSOR
      CALL PK2ISOMATFIC(PKMATFIC,DISO,CBAR,CBARI1,UNIT2,NDI)
C     CAUCHY 'FICTICIOUS' STRESS TENSOR
      CALL SIGISOMATFIC(SISOMATFIC,PKMATFIC,DISTGR,DET,NDI)
C     'FICTICIOUS' MATERIAL ELASTICITY TENSOR
      CALL CMATISOMATFIC(CMISOMATFIC,CBAR,CBARI1,CBARI2,
     1                          DISO,UNIT2,UNIT4,DET,NDI)
C     'FICTICIOUS' SPATIAL ELASTICITY TENSOR
      CALL CSISOMATFIC(CISOMATFIC,CMISOMATFIC,DISTGR,DET,NDI)
C
C---- FILAMENTS NETWORK -----------------------------------------------
C     IMAGINARY ERROR FUNCTION BASED ON DISPERSION PARAMETER
      CALL ERFI(EFI,BB,NTERM)
C
C     'FICTICIOUS' PK2 STRESS AND MATERIAL ELASTICITY TENSORS
C      CALL AFFMATNETFIC(PKNETFICAF,CMNETFICAF,DISTGR,MF0,RW,FILPROPS,
C     1                         AFFPROPS,EFI,NOEL,DET,NDI)
C     
      CALL AFFNETFIC(SNETFICAF,CNETFICAF,DISTGR,MF0,RW,FILPROPS,
     1          AFFPROPS,EFI,NOEL,DET,NDI)
C
C      PKNETFIC=PKNETFICAF
      SNETFIC=SNETFICAF
C      CMNETFIC=CMNETFICAF
      CNETFIC=CNETFICAF
C----------------------------------------------------------------------
C     STRAIN-ENERGY
C      SSE=SSEV+SSEISO
C     PK2 'FICTICIOUS' STRESS
      PKFIC=(ONE-PHI)*PKMATFIC+PKNETFIC
C     CAUCHY 'FICTICIOUS' STRESS
      SFIC=(ONE-PHI)*SISOMATFIC+SNETFIC
C     MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      CMFIC=(ONE-PHI)*CMISOMATFIC+CMNETFIC
C     SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      CFIC=(ONE-PHI)*CISOMATFIC+CNETFIC
C
C----------------------------------------------------------------------
C-------------------------- STRESS MEASURES ---------------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C      PK2 STRESS
      CALL PK2VOL(PKVOL,PV,C,NDI)
C      CAUCHY STRESS
      CALL SIGVOL(SVOL,PV,UNIT2,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C      PK2 STRESS
      CALL PK2ISO(PKISO,PKFIC,PROJL,DET,NDI)
C      CAUCHY STRESS
      CALL SIGISO(SISO,SFIC,PROJE,NDI)
C      ACTIVE CAUCHY STRESS
C      CALL SIGISO(SACTISO,SNETFICAF,PROJE,NDI)
C     
C      CALL SPECTRAL(SACTISO,SACTVL,SACTVC)
C
C---- VOLUMETRIC + ISOCHORIC ------------------------------------------
C      PK2 STRESS
      PK2 = PKVOL + PKISO
C      CAUCHY STRESS
      SIGMA = SVOL + SISO
C
C----------------------------------------------------------------------
C-------------------- MATERIAL ELASTICITY TENSOR ----------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C
C      CALL METVOL(CMVOL,C,PV,PPV,DET,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C
C      CALL METISO(CMISO,CMFIC,PROJL,PKISO,PKFIC,C,UNIT2,DET,NDI)
C
C----------------------------------------------------------------------
C
C      DDPKDDE=CMVOL+CMISO
C                                                 
C----------------------------------------------------------------------
C--------------------- SPATIAL ELASTICITY TENSOR ----------------------
C----------------------------------------------------------------------
C
C---- VOLUMETRIC ------------------------------------------------------
C
      CALL SETVOL(CVOL,PV,PPV,UNIT2,UNIT4S,NDI)
C
C---- ISOCHORIC -------------------------------------------------------
C
      CALL SETISO(CISO,CFIC,PROJE,SISO,SFIC,UNIT2,NDI)
C
C-----JAUMMAN RATE ----------------------------------------------------
C
      CALL SETJR(CJR,SIGMA,UNIT2,NDI)
C
C----------------------------------------------------------------------
C
C     ELASTICITY TENSOR
      DDSIGDDE=CVOL+CISO+CJR
C
C----------------------------------------------------------------------
C-------------------------- VISCOUS PART ------------------------------
C----------------------------------------------------------------------
C      PULLBACK OF STRESS AND ELASTICITY TENSORS
      CALL PULL2(PKVOL,SVOL,DFGRD1INV,DET,NDI)
      CALL PULL2(PKISO,SISO,DFGRD1INV,DET,NDI)
      CALL PULL4(CMVOL,CVOL,DFGRD1INV,DET,NDI)
      CALL PULL4(CMISO,CISO,DFGRD1INV,DET,NDI)
C      VISCOUS DAMPING 
      CALL VISCO(PK2,DDPKDDE,VV,PKVOL,PKISO,CMVOL,CMISO,DTIME,
     1            VSCPROPS,STATEV,NDI) 
C      PUSH FORWARD OF STRESS AND ELASTICITY TENSOR
      CALL PUSH2(SIGMA,PK2,DFGRD1,DET,NDI)
C
      CALL PUSH4(DDSIGDDE,DDPKDDE,DFGRD1,DET,NDI)
C----------------------------------------------------------------------
C------------------------- INDEX ALLOCATION ---------------------------
C----------------------------------------------------------------------
C     VOIGT NOTATION  - FULLY SIMMETRY IMPOSED
      CALL INDEXX(STRESS,DDSDDE,SIGMA,DDSIGDDE,NTENS,NDI)
C
C----------------------------------------------------------------------
C--------------------------- STATE VARIABLES --------------------------
C----------------------------------------------------------------------
C     DO K1 = 1, NTENS
C      STATEV(1:27) = VISCOUS TENSORS
       CALL SDVWRITE(STATEV,DET,VV)
C     END DO
C----------------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------------
C--------------------------- END OF UMAT ------------------------------
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C----------------------- AUXILIAR SUBROUTINES -------------------------
C----------------------------------------------------------------------
C                         INPUT FILES
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C                         KINEMATIC QUANTITIES
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                         STRESS TENSORS
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C                   LINEARISED ELASTICITY TENSORS
C----------------------------------------------------------------------
C
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------- UTILITY SUBROUTINES --------------------------
C----------------------------------------------------------------------
C
