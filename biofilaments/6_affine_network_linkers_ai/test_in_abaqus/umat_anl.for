      SUBROUTINE FACTORIAL(FACT,TERM)
C>    FACTORIAL
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION FACT
      INTEGER M,TERM
C
      FACT = 1
C
      DO  M = 1, TERM
           FACT = FACT * M
      ENDDO
C
      RETURN
      END SUBROUTINE FACTORIAL
      SUBROUTINE PULL2(PK,SIG,FINV,DET,NDI)
C>       PULL-BACK TIMES DET OF A 2ND ORDER TENSOR
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,II1,JJ1,NDI
       DOUBLE PRECISION PK(NDI,NDI),FINV(NDI,NDI)
       DOUBLE PRECISION SIG(NDI,NDI)
       DOUBLE PRECISION AUX,DET
C
       DO I1=1,NDI
        DO J1=1,NDI
          AUX=ZERO
         DO II1=1,NDI
          DO JJ1=1,NDI
            AUX=AUX+DET*FINV(I1,II1)*FINV(J1,JJ1)*SIG(II1,JJ1)
         END DO
        END DO
           PK(I1,J1)=AUX
        END DO
       END DO
C
       RETURN
      END SUBROUTINE PULL2
      SUBROUTINE METVOL(CVOL,C,PV,PPV,DET,NDI)
C>    VOLUMETRIC MATERIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION C(NDI,NDI),CINV(NDI,NDI),
     1                 CVOL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION PV,PPV,DET

C
      CALL MATINV3D(C,CINV,NDI)
C
      DO I1 = 1, NDI
        DO J1 = 1, NDI
         DO K1 = 1, NDI
           DO L1 = 1, NDI
             CVOL(I1,J1,K1,L1)=
     1                 DET*PPV*CINV(I1,J1)*CINV(K1,L1)
     2           -DET*PV*(CINV(I1,K1)*CINV(J1,L1)
     3                      +CINV(I1,L1)*CINV(J1,K1))
           END DO
         END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE METVOL
      SUBROUTINE CONTRACTION42(S,LT,RT,NDI)
C>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
C>      INPUT:
C>       LT - RIGHT 4TH ORDER TENSOR
C>       RT - LEFT  2ND ODER TENSOR
C>      OUTPUT:
C>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,K1,L1,NDI
C
       DOUBLE PRECISION RT(NDI,NDI),LT(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION S(NDI,NDI)
       DOUBLE PRECISION AUX
C
       DO I1=1,NDI
        DO J1=1,NDI
          AUX=ZERO
         DO K1=1,NDI
          DO L1=1,NDI
            AUX=AUX+LT(I1,J1,K1,L1)*RT(K1,L1)
         END DO
        END DO
           S(I1,J1)=AUX
       END DO
      END DO
       RETURN
      END SUBROUTINE CONTRACTION42
       SUBROUTINE EVALG(G,F,LAMBDA,LAMBDA0,L,R0,MU0,BETA,B0)
C>     ESTABLISHMENT OF G(F)=LHS-RHS(F) THAT RELATES
C>       STRETCHFORCE RELATIONSHIP OF A SINGLE EXNTESIBLE FILAMENT
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION G,F,LHS,RHS
C
      DOUBLE PRECISION L,R0,MU0,B0,BETA,LAMBDA,LAMBDA0
      DOUBLE PRECISION AUX0,AUX1,AUX,AUX2,AUX3,AUX4,PI
C
      PI=FOUR*ATAN(ONE)
      AUX0=ONE-R0/L
      AUX1=L*L*((PI*PI*B0)**(-ONE))
      AUX=F/MU0
      AUX2=ONE+AUX
      AUX3=ONE+TWO*AUX
      AUX4=ONE+F*AUX1+F*AUX*AUX1
C
      RHS=ONE+AUX-AUX0*(AUX2**BETA)*AUX3*(AUX4**(-BETA))
      LHS=LAMBDA*LAMBDA0*R0*(L**(-ONE))
C
      G=LHS-RHS
C
      RETURN
      END SUBROUTINE EVALG
      SUBROUTINE PROJEUL(A,AA,PE,NDI)
C>    EULERIAN PROJECTION TENSOR
C      INPUTS:
C          IDENTITY TENSORS - A, AA
C      OUTPUTS:
C          4TH ORDER SYMMETRIC EULERIAN PROJECTION TENSOR - PE
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER I,J,K,L,NDI
C
      DOUBLE PRECISION A(NDI,NDI),AA(NDI,NDI,NDI,NDI),
     1                 PE(NDI,NDI,NDI,NDI)
C
      DO I=1,NDI
         DO J=1,NDI
          DO K=1,NDI
             DO L=1,NDI
              PE(I,J,K,L)=AA(I,J,K,L)-(ONE/THREE)*(A(I,J)*A(K,L))
           END DO
          END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PROJEUL
      SUBROUTINE SIGVOL(SVOL,PV,UNIT2,NDI)
C>    VOLUMETRIC CAUCHY STRESS 
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION UNIT2(NDI,NDI),SVOL(NDI,NDI)
      DOUBLE PRECISION PV
C
      DO I1=1,NDI
        DO J1=1,NDI
          SVOL(I1,J1)=PV*UNIT2(I1,J1)
        END DO
      END DO
C
      RETURN
      END SUBROUTINE SIGVOL
      SUBROUTINE ISOMAT(SSEISO,DISO,C10,C01,CBARI1,CBARI2)
C>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION SSEISO,DISO(5)
      DOUBLE PRECISION C10,C01,CBARI1,CBARI2
C
      SSEISO=C10*(CBARI1-THREE)+C01*(CBARI2-THREE)
C
      DISO(1)=C10
      DISO(2)=C01
      DISO(3)=ZERO
      DISO(4)=ZERO
      DISO(5)=ZERO
C
      RETURN
      END SUBROUTINE ISOMAT
      SUBROUTINE SDVREAD(STATEV)
C>    VISCOUS DISSIPATION: READ STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER VV,POS1,POS2,POS3,I1
      DOUBLE PRECISION STATEV(NSDV)
C

C
      RETURN
C
      END SUBROUTINE SDVREAD
      SUBROUTINE PK2VOL(PKVOL,PV,C,NDI)
C>    VOLUMETRIC PK2 STRESS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION PKVOL(NDI,NDI),C(NDI,NDI),CINV(NDI,NDI)
      DOUBLE PRECISION PV
C
      CALL MATINV3D(C,CINV,NDI)
C
      DO I1=1,NDI
        DO J1=1,NDI
          PKVOL(I1,J1)=PV*CINV(I1,J1)
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PK2VOL
      SUBROUTINE PULL4(MAT,SPATIAL,FINV,DET,NDI)
C>        PULL-BACK TIMES DET OF 4TH ORDER TENSOR
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,K1,L1,II1,JJ1,KK1,LL1,NDI
       DOUBLE PRECISION MAT(NDI,NDI,NDI,NDI),FINV(NDI,NDI)
       DOUBLE PRECISION SPATIAL(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION AUX,DET
C
       DO I1=1,NDI
        DO J1=1,NDI
         DO K1=1,NDI
          DO L1=1,NDI
           AUX=ZERO
           DO II1=1,NDI
            DO JJ1=1,NDI
             DO KK1=1,NDI
              DO LL1=1,NDI
              AUX=AUX+DET*
     +        FINV(I1,II1)*FINV(J1,JJ1)*
     +        FINV(K1,KK1)*FINV(L1,LL1)*SPATIAL(II1,JJ1,KK1,LL1)
              END DO
             END DO
            END DO
           END DO
           MAT(I1,J1,K1,L1)=AUX
          END DO
         END DO
        END DO
       END DO
C
       RETURN
      END SUBROUTINE PULL4
      SUBROUTINE DEFFIL(LAMBDA,M,M0,F,NDI)
C>      SINGLE FILAMENT: STRETCH AND DEFORMED DIRECTION
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION F(NDI,NDI),M(NDI),M0(NDI)
      DOUBLE PRECISION AUX,LAMBDA
C
      LAMBDA=ZERO
      DO I1=1,NDI
        AUX=ZERO
       DO J1=1,NDI
         AUX=AUX+F(I1,J1)*M0(J1)
       END DO
         M(I1)=AUX
      END DO
         LAMBDA=DOT_PRODUCT(M,M)
         LAMBDA=SQRT(LAMBDA)
C
      RETURN
      END SUBROUTINE DEFFIL
      SUBROUTINE SETJR(CJR,SIGMA,UNIT2,NDI)
C>    JAUMAN RATE CONTRIBUTION FOR THE SPATIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),
     1                 CJR(NDI,NDI,NDI,NDI),SIGMA(NDI,NDI)
C
      DO I1 = 1, NDI
        DO J1 = 1, NDI
         DO K1 = 1, NDI
           DO L1 = 1, NDI
           
              CJR(I1,J1,K1,L1)=
     1             (ONE/TWO)*(UNIT2(I1,K1)*SIGMA(J1,L1)
     2             +SIGMA(I1,K1)*UNIT2(J1,L1)+UNIT2(I1,L1)*SIGMA(J1,K1)
     3             +SIGMA(I1,L1)*UNIT2(J1,K1))
           END DO
         END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE SETJR
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
      INCLUDE 'PARAM_UMAT.INC'
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
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),FILPROPS(8),
     1                 NAFFPROPS(2),AFFPROPS(2)
      DOUBLE PRECISION LL,R0,LAMBDA0,MU0,BETA,NN,MM,B0,BB
      DOUBLE PRECISION PHI,P,R0C,R0F,ETAC
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
      R0F      = PROPS(6)
      R0C      = PROPS(7)
      ETAC     = PROPS(8)
      MU0      = PROPS(9)
      BETA     = PROPS(10)
      B0       = PROPS(11)
      LAMBDA0  = PROPS(12)
      FILPROPS = PROPS(5:12)
C     NETWORK
      NN       = PROPS(13)
      BB        = PROPS(14)   
      AFFPROPS = PROPS(13:14)    
C
C     NUMERICAL COMPUTATIONS
      NTERM    = 60      
C 
C        STATE VARIABLES AND CHEMICAL PARAMETERS
C
      IF ((TIME(1).EQ.ZERO).AND.(KSTEP.EQ.1)) THEN
      CALL INITIALIZE(STATEV)   
      ENDIF
C        READ STATEV
      CALL SDVREAD(STATEV)       
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
C      CALL AFFCLMATNETFIC(PKNETFICAF,CMNETFICAF,DISTGR,MF0,RW,FILPROPS,
C     1                         AFFPROPS,EFI,NOEL,DET,NDI)
C     
      CALL AFFCLNETFIC(SNETFICAF,CNETFICAF,DISTGR,MF0,RW,FILPROPS,
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
C
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
       CALL SDVWRITE(DET,STATEV)
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
      SUBROUTINE CONTRACTION44(S,LT,RT,NDI)
C>       DOUBLE CONTRACTION BETWEEN 4TH ORDER TENSORS
C>      INPUT:
C>       LT - RIGHT 4TH ORDER TENSOR
C>       RT - LEFT  4TH ORDER TENSOR
C>      OUTPUT:
C>       S - DOUBLE CONTRACTED TENSOR (4TH ORDER)
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,K1,L1,M1,N1,NDI
C
       DOUBLE PRECISION LT(NDI,NDI,NDI,NDI),RT(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION S(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION AUX

C                     
       DO I1=1,NDI
        DO J1=1,NDI
         DO K1=1,NDI
          DO L1=1,NDI
           AUX=ZERO
           DO M1=1,NDI
            DO N1=1,NDI
                AUX=AUX+LT(I1,J1,M1,N1)*RT(M1,N1,K1,L1)
            END DO
           END DO
           S(I1,J1,K1,L1)=AUX
          END DO
         END DO
        END DO
       END DO
C
       RETURN
      END SUBROUTINE CONTRACTION44
       SUBROUTINE INITIALIZE(STATEV)
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
C      DOUBLE PRECISION TIME(2),KSTEP
      INTEGER I1,POS,POS1,POS2,POS3,VV
      DOUBLE PRECISION STATEV(NSDV)
C
        POS1=0
C       DETERMINANT        
          STATEV(POS1+1)=ONE
C        CONTRACTION VARIANCE
          STATEV(POS1+2)=ZERO          
C     
      RETURN
C
      END SUBROUTINE INITIALIZE
      SUBROUTINE SETVOL(CVOL,PV,PPV,UNIT2,UNIT4S,NDI)
C>    VOLUMETRIC SPATIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),UNIT4S(NDI,NDI,NDI,NDI),
     1                 CVOL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION PV,PPV
C
      DO I1 = 1, NDI
        DO J1 = 1, NDI
         DO K1 = 1, NDI
           DO L1 = 1, NDI
             CVOL(I1,J1,K1,L1)=
     1                 PPV*UNIT2(I1,J1)*UNIT2(K1,L1)
     2                 -TWO*PV*UNIT4S(I1,J1,K1,L1)
           END DO
         END DO
        END DO
      END DO
C      
      RETURN
      END SUBROUTINE SETVOL
      SUBROUTINE PK2ISOMATFIC(FIC,DISO,CBAR,CBARI1,UNIT2,NDI)
C>     ISOTROPIC MATRIX: 2PK 'FICTICIOUS' STRESS TENSOR
C      INPUT:
C       DISO - STRAIN-ENERGY DERIVATIVES
C       CBAR - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
C       CBARI1,CBARI2 - CBAR INVARIANTS
C       UNIT2 - 2ND ORDER IDENTITY TENSOR
C      OUTPUT:
C       FIC - 2ND PIOLA KIRCHOOF 'FICTICIOUS' STRESS TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER I1,J1,NDI
      DOUBLE PRECISION FIC(NDI,NDI),DISO(5),CBAR(NDI,NDI),UNIT2(NDI,NDI)
      DOUBLE PRECISION DUDI1,DUDI2,CBARI1
      DOUBLE PRECISION AUX1,AUX2

C
      DUDI1=DISO(1)
      DUDI2=DISO(2)
C
      AUX1=TWO*(DUDI1+CBARI1*DUDI2)
      AUX2=-TWO*DUDI2
C
      DO I1=1,NDI
       DO J1=1,NDI
        FIC(I1,J1)=AUX1*UNIT2(I1,J1)+AUX2*CBAR(I1,J1)
       END DO
      END DO
C
      RETURN
      END SUBROUTINE PK2ISOMATFIC
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C>    READ FILAMENTS ORIENTATION AND PREFERED DIRECTIONS
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'PARAM_UMAT.INC'
C       this subroutine get the directions and weights for
C      the numerical integration
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KFIL/MF0
      COMMON /KFILR/RW
      COMMON /KFILP/PREFDIR
C      
      DIMENSION TIME(2)
      DOUBLE PRECISION  SPHERE(NWP,5),MF0(NWP,3),RW(NWP),
     1                  PREFDIR(NELEM,4)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
       CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
         OPEN(15,FILE=FILENAME)
         DO I=1,NWP
            READ(15,*) (SPHERE(I,J),J=1,5)
         END DO
          CLOSE(15)
C
         DO I=1,NWP
           DO J=1,3
            K=J+1
           MF0(I,J)=SPHERE(I,K)
          END DO
          RW(I)=TWO*SPHERE(I,5)
         END DO
C
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR2
         OPEN(16,FILE=FILENAME)
         DO I=1,NELEM
            READ(16,*) (PREFDIR(I,J),J=1,4)
         END DO
         CLOSE(16)
C
C         
      END IF    
C
      RETURN
C
      END SUBROUTINE UEXTERNALDB
      SUBROUTINE SIGISOMATFIC(SFIC,PKFIC,F,DET,NDI)
C>    ISOTROPIC MATRIX:  ISOCHORIC CAUCHY STRESS 
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION SFIC(NDI,NDI),F(NDI,NDI),
     1                 PKFIC(NDI,NDI)
      DOUBLE PRECISION DET
C
      CALL PUSH2(SFIC,PKFIC,F,DET,NDI)
C
      RETURN
      END SUBROUTINE SIGISOMATFIC
      SUBROUTINE FIL(F,FF,DW,DDW,LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0)
C>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION DW,DDW
      DOUBLE PRECISION LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0
      DOUBLE PRECISION A,B,MACHEP,T
      DOUBLE PRECISION AUX,F,FF,PI,ALPHA
      DOUBLE PRECISION AUX0,AUX1,AUX2,AUX3,AUX4,AUX5,AUX6,Y
C
      A=ZERO
      B=1.0E09
      MACHEP=2.2204e-16
      T=1.0E-6
      F=ZERO
C
      CALL PULLFORCE(F, A, B, MACHEP, T,
     1                   LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0)
C
      PI=FOUR*ATAN(ONE)
      FF=F*LL*(PI*PI*B0)**(-ONE)
      ALPHA=PI*PI*B0*(LL*LL*MU0)**(-ONE)
C
      AUX0=BETA/ALPHA
      AUX=ALPHA*FF
      AUX1=ONE+FF+AUX*FF
      AUX2=ONE+TWO*AUX
      AUX3=ONE+AUX
      AUX4=LAMBDA0*R0*R0*MU0*(LL**(-ONE))
      AUX5=((ONE+AUX)*(AUX1**(-ONE)))**BETA
      AUX6=ONE-R0*((LL)**(-ONE))
C      
      Y=AUX0*(AUX2*AUX2*(AUX1**(-ONE)))-BETA*(AUX2*(AUX3**(-ONE)))-TWO
C      
      DW=LAMBDA0*R0*F
      DDW=AUX4*((ONE+Y*AUX5*AUX6)**(-ONE))
C
      RETURN
      END SUBROUTINE FIL
      SUBROUTINE PULLFORCE(ZERO, A, B, MACHEP, T,
     1                   LAMBDA,LAMBDA0,L,R0,MU0,BETA,B0)
C>    SINGLE FILAMENT: COMPUTES PULLING FORCE FOR A GIVEN STRETCH     
C*********************************************************************72
C
C     ZERO SEEKS THE ROOT OF A FUNCTION F(X) IN AN INTERVAL [A,B].
C
C     DISCUSSION:
C
C     THE INTERVAL [A,B] MUST BE A CHANGE OF SIGN INTERVAL FOR F.
C     THAT IS, F(A) AND F(B) MUST BE OF OPPOSITE SIGNS.  THEN
C     ASSUMING THAT F IS CONTINUOUS IMPLIES THE EXISTENCE OF AT LEAST
C     ONE VALUE C BETWEEN A AND B FOR WHICH F(C) = 0.
C
C     THE LOCATION OF THE ZERO IS DETERMINED TO WITHIN AN ACCURACY
C     OF 6 * MACHEPS * ABS ( C ) + 2 * T.
C
C
C     LICENSING:
C
C     THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
C
C     MODIFIED:
C
C     11 FEBRUARY 2013
C
C     AUTHOR:
C
C     RICHARD BRENT
C     MODIFICATIONS BY JOHN BURKARDT
C
C     REFERENCE:
C
C     RICHARD BRENT,
C     ALGORITHMS FOR MINIMIZATION WITHOUT DERIVATIVES,
C     DOVER, 2002,
C     ISBN: 0-486-41998-3,
C     LC: QA402.5.B74.
C
C     PARAMETERS:
C
C     INPUT, DOUBLE PRECISION A, B, THE ENDPOINTS OF THE CHANGE OF SIGN
C     INTERVAL.
C     INPUT, DOUBLE PRECISION MACHEP, AN ESTIMATE FOR THE RELATIVE
C     MACHINE PRECISION.
C
C     INPUT, DOUBLE PRECISION T, A POSITIVE ERROR TOLERANCE.
C
C     INPUT, EXTERNAL DOUBLE PRECISION F, THE NAME OF A USER-SUPPLIED
C     FUNCTION, OF THE FORM "FUNCTION G ( F )", WHICH EVALUATES THE
C     FUNCTION WHOSE ZERO IS BEING SOUGHT.
C
C     OUTPUT, DOUBLE PRECISION ZERO, THE ESTIMATED VALUE OF A ZERO OF
C     THE FUNCTION G.
C
      IMPLICIT NONE
C
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION C
      DOUBLE PRECISION D
      DOUBLE PRECISION E
      DOUBLE PRECISION FA
      DOUBLE PRECISION FB
      DOUBLE PRECISION FC
      DOUBLE PRECISION M
      DOUBLE PRECISION MACHEP
      DOUBLE PRECISION P
      DOUBLE PRECISION Q
      DOUBLE PRECISION R
      DOUBLE PRECISION S
      DOUBLE PRECISION SA
      DOUBLE PRECISION SB
      DOUBLE PRECISION T
      DOUBLE PRECISION TOL
      DOUBLE PRECISION ZERO
C
      DOUBLE PRECISION L,R0,MU0,B0,BETA,LAMBDA,LAMBDA0
C      
C     MAKE LOCAL COPIES OF A AND B.
C
      SA = A
      SB = B
      CALL EVALG(FA,SA,LAMBDA,LAMBDA0,L,R0,MU0,BETA,B0)
      CALL EVALG(FB,SB,LAMBDA,LAMBDA0,L,R0,MU0,BETA,B0)
C      FA = F ( SA )
C      FB = F ( SB )
C
10    CONTINUE
C
      C = SA
      FC = FA
      E = SB - SA
      D = E

20    CONTINUE
C
      IF ( ABS ( FC ) .LT. ABS ( FB ) ) THEN
        SA = SB
        SB = C
        C = SA
        FA = FB
        FB = FC
        FC = FA
      END IF
C
30    CONTINUE
C
      TOL = 2.0D+00 * MACHEP * ABS ( SB ) + T
      M = 0.5D+00 * ( C - SB )
      IF ( ABS ( M ) .LE. TOL .OR. FB .EQ. 0.0D+00 ) GO TO 140
      IF ( ABS ( E ) .GE. TOL .AND. ABS ( FA ) .GT. ABS ( FB ) )
     &  GO TO 40
C
      E = M
      D = E
      GO TO 100
C
40    CONTINUE
C
      S = FB / FA
      IF ( SA .NE. C ) GO TO 50
C
      P = 2.0D+00 * M * S
      Q = 1.0D+00 - S
      GO TO 60
C
50    CONTINUE
C
      Q = FA / FC
      R = FB / FC
      P = S *
     &  ( 2.0D+00 * M * Q * ( Q - R ) - ( SB - SA ) * ( R - 1.0D+00 ) )
      Q = ( Q - 1.0D+00 ) * ( R - 1.0D+00 ) * ( S - 1.0D+00 )
C
60    CONTINUE
C
      IF ( P .LE. 0.0D+00 ) GO TO 70
C
      Q = - Q
      GO TO 80
C
70    CONTINUE
C
      P = - P
C
80    CONTINUE
C
      S = E
      E = D
      IF ( 2.0D+00 * P .GE. 3.0D+00 * M * Q - ABS ( TOL * Q ) .OR.
     &  P .GE. ABS ( 0.5D+00 * S * Q ) ) GO TO 90
C
      D = P / Q
      GO TO 100
C
90    CONTINUE
C
      E = M
      D = E
C
100   CONTINUE
C
      SA = SB
      FA = FB
      IF ( ABS ( D ) .LE. TOL ) GO TO 110
      SB = SB + D
      GO TO 130
C
110   CONTINUE
C
      IF ( M .LE. 0.0D+00 ) GO TO 120
      SB = SB + TOL
      GO TO 130
C
120   CONTINUE
C
      SB = SB - TOL
C
130   CONTINUE
C
C      FB = F ( SB )
      CALL EVALG(FB,SB,LAMBDA,LAMBDA0,L,R0,MU0,BETA,B0)
      IF ( FB .GT. 0.0D+00 .AND. FC .GT. 0.0D+00 ) GO TO 10
      IF ( FB .LE. 0.0D+00 .AND. FC .LE. 0.0D+00 ) GO TO 10
      GO TO 20
C
140   CONTINUE
C
      ZERO = SB
C
      RETURN
      END SUBROUTINE PULLFORCE
C
C*********************************************************************72
      SUBROUTINE STRETCH(C,B,U,V,NDI)
C>    STRETCH TENSORS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION C(NDI,NDI),B(NDI,NDI),U(NDI,NDI),V(NDI,NDI)
      DOUBLE PRECISION EIGVAL(NDI),OMEGA(NDI),EIGVEC(NDI,NDI)
C
      CALL SPECTRAL(C,OMEGA,EIGVEC)
C
      EIGVAL(1) = DSQRT(OMEGA(1))
      EIGVAL(2) = DSQRT(OMEGA(2))
      EIGVAL(3) = DSQRT(OMEGA(3))
C
      U(1,1) = EIGVAL(1)
      U(2,2) = EIGVAL(2)
      U(3,3) = EIGVAL(3)
C
      U = MATMUL(MATMUL(EIGVEC,U),TRANSPOSE(EIGVEC))
C
      CALL SPECTRAL(B,OMEGA,EIGVEC)
C
      EIGVAL(1) = DSQRT(OMEGA(1))
      EIGVAL(2) = DSQRT(OMEGA(2))
      EIGVAL(3) = DSQRT(OMEGA(3))      
C      write(*,*) eigvec(1,1),eigvec(2,1),eigvec(3,1)
C
      V(1,1) = EIGVAL(1)
      V(2,2) = EIGVAL(2)
      V(3,3) = EIGVAL(3)
C
      V = MATMUL(MATMUL(EIGVEC,V),TRANSPOSE(EIGVEC))
      RETURN
      END SUBROUTINE STRETCH
      SUBROUTINE AFFCLNETFIC(SFIC,CFIC,F,MF0,RW,FILPROPS,AFFPROPS,
     1                        EFI,NOEL,DET,NDI)
C>    AFFINE NETWORK: 'FICTICIOUS' CAUCHY STRESS AND ELASTICITY TENSOR
C>    MOBILE LINKERS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1,M1,NOEL,IM1
      DOUBLE PRECISION SFIC(NDI,NDI),SFILFIC(NDI,NDI),
     1                  CFIC(NDI,NDI,NDI,NDI),F(NDI,NDI),MF0(NWP,NDI),
     2                 RW(NWP),CFILFIC(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION FILPROPS(8),AFFPROPS(2),MFI(NDI),MF0I(NDI)
      DOUBLE PRECISION DET,AUX,PI,LAMBDAI,DWI,DDWI,RWI,LAMBDAIC
      DOUBLE PRECISION L,R0F,R0,MU0,B0,BETA,LAMBDA0,RHO,N,FI,FFI,DTIME
      DOUBLE PRECISION R0C,ETAC,LAMBDAIF
      DOUBLE PRECISION B,FRIC,FFMAX,ANG,EFI,FRAC(4),RU0(NWP),RU
      DOUBLE PRECISION VARA,AVGA,MAXA,AUX0,FFIC,SUMA,RHO0,DIRMAX(NDI)
C
C     FILAMENT
      L       = FILPROPS(1)
      R0F     = FILPROPS(2)
      R0C     = FILPROPS(3)
      ETAC    = FILPROPS(4)      
      MU0     = FILPROPS(5)
      BETA    = FILPROPS(6)
      B0      = FILPROPS(7)
      LAMBDA0 = FILPROPS(8)
C     NETWORK
      N       = AFFPROPS(1)
      B       = AFFPROPS(2)      
C
      PI=FOUR*ATAN(ONE)
      AUX=N*(DET**(-ONE))*FOUR*PI
      CFIC=ZERO
      SFIC=ZERO
C      
       RHO=ONE
       R0=R0F+R0C
C       
C       CALL DENSITY(RHO0,ZERO,B,EFI)
C        
C             OPEN (UNIT=20,FILE="projfil.out",action="write",
C     1 status="replace")
     
C        LOOP OVER THE INTEGRATION DIRECTIONS
      DO I1=1,NWP
C
       MFI=ZERO
       MF0I=ZERO
       DO J1=1,NDI
        MF0I(J1)=MF0(I1,J1)
       END DO
       RWI=RW(I1)
C
       CALL DEFFIL(LAMBDAI,MFI,MF0I,F,NDI)
C        
       CALL BANGLE(ANG,F,MFI,NOEL,NDI)
C      
       CALL DENSITY(RHO,ANG,B,EFI)
C       
       IF((ETAC.GT.ZERO).AND.(ETAC.LT.ONE))THEN
C
        LAMBDAIF=ETAC*(R0/R0F)*(LAMBDAI-ONE)+ONE
        LAMBDAIC=(LAMBDAI*R0-LAMBDAIF*R0F)/R0C
       ELSE
        LAMBDAIF=LAMBDAI
        LAMBDAIC=ZERO
       ENDIF   
        

C     
       CALL FIL(FI,FFI,DWI,DDWI,LAMBDAIF,LAMBDA0,L,R0,MU0,BETA,B0)          
C       
       CALL SIGFILFIC(SFILFIC,RHO,LAMBDAIF,DWI,MFI,RWI,NDI)
C
       CALL CSFILFIC(CFILFIC,RHO,LAMBDAIF,DWI,DDWI,MFI,RWI,NDI)
C     
C
       DO J1=1,NDI
        DO K1=1,NDI
         SFIC(J1,K1)=SFIC(J1,K1)+AUX*SFILFIC(J1,K1)
         DO L1=1,NDI
          DO M1=1,NDI
           CFIC(J1,K1,L1,M1)=CFIC(J1,K1,L1,M1)+AUX*CFILFIC(J1,K1,L1,M1)
          END DO
         END DO
        END DO
       END DO  
C             
      END DO  
C        
      RETURN
      END SUBROUTINE AFFCLNETFIC
      SUBROUTINE MATINV3D(A,A_INV,NDI)
C>    INVERSE OF A 3X3 MATRIX
C     RETURN THE INVERSE OF A(3,3) - A_INV
      IMPLICIT NONE
C
      INTEGER NDI
C
      DOUBLE PRECISION A(NDI,NDI),A_INV(NDI,NDI),DET_A,DET_A_INV
C
      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))

      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) = DET_A_INV*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_INV(1,2) = DET_A_INV*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_INV(1,3) = DET_A_INV*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_INV(2,1) = DET_A_INV*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_INV(2,2) = DET_A_INV*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_INV(2,3) = DET_A_INV*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_INV(3,1) = DET_A_INV*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_INV(3,2) = DET_A_INV*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_INV(3,3) = DET_A_INV*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
C
      RETURN
      END SUBROUTINE MATINV3D
      SUBROUTINE SPECTRAL(A,D,V)
C>    EIGENVALUES AND EIGENVECTOR OF A 3X3 MATRIX
C     THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C     A SYMMETRIC 3X3 MATRIX A.
C
C     THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C     EIGENVALUES IN ASCENDING ORDER, AND A MATRIX V WHOSE
C     COLUMNS CONTAIN THE CORRESPONDING EIGENVECTORS.
C
      IMPLICIT NONE
C
      INTEGER NP,NROT
      PARAMETER(NP=3)
C
      DOUBLE PRECISION D(3),V(3,3),A(3,3),E(3,3)
C
      E = A
C
      CALL JACOBI(E,3,NP,D,V,NROT)
      CALL EIGSRT(D,V,3,NP)
C
      RETURN
      END SUBROUTINE SPECTRAL

C***********************************************************************

      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
C
C COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C  MATRIX A, WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL
C  NP BY NP ARRAY.  ON OUTPUT, ELEMENTS OF A ABOVE THE DIAGONAL
C  ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C  AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C  VECTOR D RETURNS THE EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
C  V IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C  A WHOSE COLUMNS CONTAIN, UPON OUTPUT, THE NORMALIZED
C  EIGENVECTORS OF A.  NROT RETURNS THE NUMBER OF JACOBI ROTATION
C  WHICH WERE REQUIRED.
C
C THIS SUBROUTINE IS TAKEN FROM 'NUMERICAL RECIPES.'
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER IP,IQ,N,NMAX,NP,NROT,I,J
      PARAMETER (NMAX=100)
C
      DOUBLE PRECISION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX),
     +  SM,TRESH,G,T,H,THETA,S,C,TAU


C INITIALIZE V TO THE IDENTITY MATRIX
      DO I=1,3
          V(I,I)=ONE
        DO J=1,3
          IF (I.NE.J)THEN
           V(I,J)=ZERO
         ENDIF
       END DO
      END DO
C INITIALIZE B AND D TO THE DIAGONAL OF A, AND Z TO ZERO.
C  THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C  IN EQUATION (11.1.14)
C
      DO IP = 1,N
        B(IP) = A(IP,IP)
        D(IP) = B(IP)
        Z(IP) = 0.D0
      END DO


C BEGIN ITERATION
C
      NROT = 0
      DO I=1,50
C
C         SUM OFF-DIAGONAL ELEMENTS
C
          SM = 0.D0
          DO IP=1,N-1
            DO IQ=IP+1,N
              SM = SM + DABS(A(IP,IQ))
            END DO
          END DO
C
C          IF SM = 0., THEN RETURN.  THIS IS THE NORMAL RETURN,
C          WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE
C          UNDERFLOW.
C
          IF (SM.EQ.0.D0) RETURN
C
C          IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C           |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE,
C           SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
C
          IF (I.LT.4) THEN
            TRESH = 0.2D0*SM/N**2
          ELSE
            TRESH = 0.D0
          END IF
C
          DO IP=1,N-1
            DO IQ=IP+1,N
              G = 100.D0*DABS(A(IP,IQ))
C
C              AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE
C               OFF-DIAGONAL ELEMENT IS SMALL.
C
              IF ((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP)))
     +            .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF (DABS(A(IP,IQ)).GT.TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G.EQ.DABS(H)) THEN
C
C                  T = 1./(2.*THETA), EQUATION (11.1.10)
C
                  T =A(IP,IQ)/H
                ELSE
                  THETA = 0.5D0*H/A(IP,IQ)
                  T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2.D0))
                  IF (THETA.LT.0.D0) T = -T
                END IF
                C = 1.D0/DSQRT(1.D0 + T**2.D0)
                S = T*C
                TAU = S/(1.D0 + C)
                H = T*A(IP,IQ)
                Z(IP) = Z(IP) - H
                Z(IQ) = Z(IQ) + H
                D(IP) = D(IP) - H
                D(IQ) = D(IQ) + H
                A(IP,IQ) = 0.D0
C
C               CASE OF ROTATIONS 1 <= J < P
C
                DO J=1,IP-1
                  G = A(J,IP)
                  H = A(J,IQ)
                  A(J,IP) = G - S*(H + G*TAU)
                  A(J,IQ) = H + S*(G - H*TAU)
                END DO
C
C                CASE OF ROTATIONS P < J < Q
C
                DO J=IP+1,IQ-1
                  G = A(IP,J)
                  H = A(J,IQ)
                  A(IP,J) = G - S*(H + G*TAU)
                  A(J,IQ) = H + S*(G - H*TAU)
                END DO
C
C                 CASE OF ROTATIONS Q < J <= N
C
                DO J=IQ+1,N
                  G = A(IP,J)
                  H = A(IQ,J)
                  A(IP,J) = G - S*(H + G*TAU)
                  A(IQ,J) = H + S*(G - H*TAU)
                END DO
                DO J = 1,N
                  G = V(J,IP)
                  H = V(J,IQ)
                  V(J,IP) = G - S*(H + G*TAU)
                  V(J,IQ) = H + S*(G - H*TAU)
                END DO
                NROT = NROT + 1
             END IF
               END DO
             END DO
C
C          UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
C
       DO IP=1,N
         B(IP) = B(IP) + Z(IP)
         D(IP) = B(IP)
         Z(IP) = 0.D0
       END DO
      END DO
C
C IF THE ALGORITHM HAS REACHED THIS STAGE, THEN THERE
C  ARE TOO MANY SWEEPS.  PRINT A DIAGNOSTIC AND CUT THE
C  TIME INCREMENT.
C
      WRITE (*,'(/1X,A/)') '50 ITERATIONS IN JACOBI SHOULD NEVER HAPPEN'
C
      RETURN
      END SUBROUTINE JACOBI

C**********************************************************************
      SUBROUTINE EIGSRT(D,V,N,NP)
C
C     GIVEN THE EIGENVALUES D AND EIGENVECTORS V AS OUTPUT FROM
C     JACOBI, THIS SUBROUTINE SORTS THE EIGENVALUES INTO ASCENDING
C     ORDER AND REARRANGES THE COLMNS OF V ACCORDINGLY.
C
C     THE SUBROUTINE WAS TAKEN FROM 'NUMERICAL RECIPES.'
C
      IMPLICIT NONE
C
      INTEGER N,NP,I,J,K
C
      DOUBLE PRECISION D(NP),V(NP,NP),P
C
      DO I=1,N-1
              K = I
              P = D(I)
              DO J=I+1,N
               IF (D(J).GE.P) THEN
                K = J
                P = D(J)
               END IF
              END DO
              IF (K.NE.I) THEN
               D(K) = D(I)
               D(I) = P
               DO J=1,N
                P = V(J,I)
                V(J,I) = V(J,K)
                V(J,K) = P
               END DO
              END IF
      END DO
C
      RETURN
      END SUBROUTINE EIGSRT
      SUBROUTINE CSFILFIC(CFIC,RHO,LAMBDA,DW,DDW,M,RW,NDI)
C>    AFFINE NETWORK: 'FICTICIOUS' ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION CFIC(NDI,NDI,NDI,NDI),M(NDI)
      DOUBLE PRECISION RHO,AUX,DW,DDW,RW,LAMBDA,AUX0
C
      AUX0=DDW-(LAMBDA**(-ONE))*DW
      AUX=RHO*AUX0*RW*(LAMBDA**(-TWO))
      DO I1=1,NDI
       DO J1=1,NDI
        DO K1=1,NDI
         DO L1=1,NDI
          CFIC(I1,J1,K1,L1)=AUX*M(I1)*M(J1)*M(K1)*M(L1)
         END DO
        END DO
       END DO
      END DO
C
      RETURN
      END SUBROUTINE CSFILFIC
      SUBROUTINE CONTRACTION24(S,LT,RT,NDI)
C>       DOUBLE CONTRACTION BETWEEN 4TH ORDER AND 2ND ORDER  TENSOR
C>      INPUT:
C>       LT - RIGHT 2ND ORDER TENSOR
C>       RT - LEFT  4TH ODER TENSOR
C>      OUTPUT:
C>       S - DOUBLE CONTRACTED TENSOR (2ND ORDER)
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,K1,L1,NDI
C
       DOUBLE PRECISION LT(NDI,NDI),RT(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION S(NDI,NDI)
       DOUBLE PRECISION AUX
C
      DO K1=1,NDI
       DO L1=1,NDI
         AUX=ZERO
        DO I1=1,NDI
         DO J1=1,NDI
           AUX=AUX+LT(K1,L1)*RT(I1,J1,K1,L1)
        END DO
       END DO
          S(K1,L1)=AUX
      END DO
      END DO
       RETURN
      END SUBROUTINE CONTRACTION24
      SUBROUTINE INVARIANTS(A,INV1,INV2,NDI)
C>    1ST AND 2ND INVARIANTS OF A TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1
      DOUBLE PRECISION A(NDI,NDI),AA(NDI,NDI)
      DOUBLE PRECISION INV1,INV1AA, INV2
C
      INV1=ZERO
      INV1AA=ZERO
      AA=MATMUL(A,A)
      DO I1=1,NDI
         INV1=INV1+A(I1,I1)
         INV1AA=INV1AA+AA(I1,I1)
      END DO
         INV2=(ONE/TWO)*(INV1*INV1-INV1AA)
C
      RETURN
      END SUBROUTINE INVARIANTS
      SUBROUTINE SIGISO(SISO,SFIC,PE,NDI)
C>    ISOCHORIC CAUCHY STRESS 
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION SISO(NDI,NDI),
     1                 PE(NDI,NDI,NDI,NDI),SFIC(NDI,NDI)
C
      CALL CONTRACTION42(SISO,PE,SFIC,NDI)
C
      RETURN
      END SUBROUTINE SIGISO
      SUBROUTINE VOL(SSEV,PV,PPV,K,DET)
C>     VOLUMETRIC CONTRIBUTION :STRAIN ENERGY FUNCTION AND DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION SSEV,PV,PPV
      DOUBLE PRECISION K,G,DET,AUX
C
      G=(ONE/FOUR)*(DET*DET-ONE-TWO*LOG(DET))
C
      SSEV=K*G
C
      PV=K*(ONE/TWO)*(DET-ONE/DET)
      AUX=K*(ONE/TWO)*(ONE+ONE/(DET*DET))
      PPV=PV+DET*AUX
C
      RETURN
      END SUBROUTINE VOL
      SUBROUTINE PUSH4(SPATIAL,MAT,F,DET,NDI)
C>        PIOLA TRANSFORMATION
C>      INPUT:
C>       MAT - MATERIAL ELASTICITY TENSOR
C>       F - DEFORMATION GRADIENT
C>       DET - DEFORMATION DETERMINANT
C>      OUTPUT:
C>       SPATIAL - SPATIAL ELASTICITY TENSOR
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,K1,L1,II1,JJ1,KK1,LL1,NDI
       DOUBLE PRECISION MAT(NDI,NDI,NDI,NDI),F(NDI,NDI)
       DOUBLE PRECISION SPATIAL(NDI,NDI,NDI,NDI)
       DOUBLE PRECISION AUX,DET
C
       DO I1=1,NDI
        DO J1=1,NDI
         DO K1=1,NDI
          DO L1=1,NDI
           AUX=ZERO
           DO II1=1,NDI
            DO JJ1=1,NDI
             DO KK1=1,NDI
              DO LL1=1,NDI
              AUX=AUX+(DET**(-ONE))*
     +        F(I1,II1)*F(J1,JJ1)*
     +        F(K1,KK1)*F(L1,LL1)*MAT(II1,JJ1,KK1,LL1)
              END DO
             END DO
            END DO
           END DO
           SPATIAL(I1,J1,K1,L1)=AUX
          END DO
         END DO
        END DO
       END DO
C
       RETURN
      END SUBROUTINE PUSH4
      SUBROUTINE SIGFILFIC(SFIC,RHO,LAMBDA,DW,M,RW,NDI)
C>    SINGLE FILAMENT:  'FICTICIUOUS' CAUCHY STRESS 
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION SFIC(NDI,NDI),M(NDI)
      DOUBLE PRECISION RHO,AUX,DW,RW,LAMBDA
C
      AUX=RHO*LAMBDA**(-ONE)*RW*DW
      DO I1=1,NDI
       DO J1=1,NDI
        SFIC(I1,J1)=AUX*M(I1)*M(J1)
       END DO
      END DO
C
      RETURN
      END SUBROUTINE SIGFILFIC
      SUBROUTINE PK2ISO(PKISO,PKFIC,PL,DET,NDI)
C>    ISOCHORIC PK2 STRESS TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION PKISO(NDI,NDI),
     1                 PL(NDI,NDI,NDI,NDI),PKFIC(NDI,NDI)
      DOUBLE PRECISION DET,SCALE2

C
      CALL CONTRACTION42(PKISO,PL,PKFIC,NDI)
C
      SCALE2=DET**(-TWO/THREE)
      DO I1=1,NDI
        DO J1=1,NDI
          PKISO(I1,J1)=SCALE2*PKISO(I1,J1)
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PK2ISO
      SUBROUTINE ERFI(ERF,B,NTERM)
C>    IMAGINARY ERROR FUNCTION OF SQRT(B); B IS THE DISPERSION PARAM    
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER I1,J1,NTERM
      DOUBLE PRECISION ERF,B,PI
      DOUBLE PRECISION AUX,AUX1,AUX2,AUX3,AUX4,HALF,FACT
C
      PI=FOUR*ATAN(ONE)
      AUX=SQRT(TWO*B)
      AUX1=TWO*AUX
      AUX2=(TWO/THREE)*(AUX**THREE)
      HALF=ONE/TWO
      AUX4=ZERO
      DO J1=3,NTERM
        I1=J1-1
       CALL FACTORIAL(FACT,I1)
       AUX3=TWO*J1-ONE
       AUX4=AUX4+(AUX**AUX3)/(HALF*AUX3*FACT)
      ENDDO
C
      ERF=PI**(-ONE/TWO)*(AUX1+AUX2+AUX4)
      RETURN
      END SUBROUTINE ERFI
      SUBROUTINE ONEM(A,AA,AAS,NDI)
C
C>      THIS SUBROUTINE GIVES:
C>          2ND ORDER IDENTITY TENSORS - A
C>          4TH ORDER IDENTITY TENSOR - AA
C>          4TH ORDER SYMMETRIC IDENTITY TENSOR - AAS
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER I,J,K,L,NDI
C
      DOUBLE PRECISION A(NDI,NDI),AA(NDI,NDI,NDI,NDI),
     1                 AAS(NDI,NDI,NDI,NDI)
C
      DO I=1,NDI
         DO J=1,NDI
            IF (I .EQ. J) THEN
              A(I,J) = ONE
            ELSE
              A(I,J) = ZERO
            END IF
         END DO
      END DO
C
      DO I=1,NDI
         DO J=1,NDI
          DO K=1,NDI
             DO L=1,NDI
              AA(I,J,K,L)=A(I,K)*A(J,L)
              AAS(I,J,K,L)=(ONE/TWO)*(A(I,K)*A(J,L)+A(I,L)*A(J,K))
           END DO
          END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE ONEM
      SUBROUTINE METISO(CMISO,CMFIC,PL,PKISO,PKFIC,C,UNIT2,DET,NDI)
C>    ISOCHORIC MATERIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),PL(NDI,NDI,NDI,NDI),
     1                 CMISO(NDI,NDI,NDI,NDI),PKISO(NDI,NDI),
     2                 CMFIC(NDI,NDI,NDI,NDI),PKFIC(NDI,NDI),
     3                 CISOAUX(NDI,NDI,NDI,NDI),
     4                 CISOAUX1(NDI,NDI,NDI,NDI),C(NDI,NDI),
     5                 PLT(NDI,NDI,NDI,NDI),CINV(NDI,NDI),
     6                 PLL(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION TRFIC,XX,YY,ZZ,DET,AUX,AUX1
C
      CALL MATINV3D(C,CINV,NDI)
      CISOAUX1=ZERO
      CISOAUX=ZERO
      CALL CONTRACTION44(CISOAUX1,PL,CMFIC,NDI)
      DO I1=1,NDI
        DO J1=1,NDI
           DO K1=1,NDI
              DO L1=1,NDI
                PLT(I1,J1,K1,L1)=PL(K1,L1,I1,J1)
              END DO
            END DO
         END DO
      END DO
C
      CALL CONTRACTION44(CISOAUX,CISOAUX1,PLT,NDI)
C
      TRFIC=ZERO
      AUX=DET**(-TWO/THREE)
      AUX1=AUX**TWO
      DO I1=1,NDI
         TRFIC=TRFIC+AUX*PKFIC(I1,I1)
      END DO
C
      DO I1=1,NDI
        DO J1=1,NDI
           DO K1=1,NDI
              DO L1=1,NDI
                XX=AUX1*CISOAUX(I1,J1,K1,L1)
                PLL(I1,J1,K1,L1)=(ONE/TWO)*(CINV(I1,K1)*CINV(J1,L1)+
     1                                      CINV(I1,L1)*CINV(J1,K1))-
     2                           (ONE/THREE)*CINV(I1,J1)*CINV(K1,L1)
                YY=TRFIC*PLL(I1,J1,K1,L1)
                ZZ=PKISO(I1,J1)*CINV(K1,L1)+CINV(I1,J1)*PKISO(K1,L1)
C
                CMISO(I1,J1,K1,L1)=XX+(TWO/THREE)*YY-(TWO/THREE)*ZZ
              END DO
           END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE METISO
      SUBROUTINE SETISO(CISO,CFIC,PE,SISO,SFIC,UNIT2,NDI)
C>    ISOCHORIC SPATIAL ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION UNIT2(NDI,NDI),PE(NDI,NDI,NDI,NDI),
     1                 CISO(NDI,NDI,NDI,NDI),SISO(NDI,NDI),
     2                 CFIC(NDI,NDI,NDI,NDI),SFIC(NDI,NDI),
     3                 CISOAUX(NDI,NDI,NDI,NDI),
     4                 CISOAUX1(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION TRFIC,XX,YY,ZZ
C
      CISOAUX1=ZERO
      CISOAUX=ZERO

      CALL CONTRACTION44(CISOAUX1,PE,CFIC,NDI)      
      CALL CONTRACTION44(CISOAUX,CISOAUX1,PE,NDI)
C
      TRFIC=ZERO
      DO I1=1,NDI
         TRFIC=TRFIC+SFIC(I1,I1)
      END DO
C
      DO I1=1,NDI
        DO J1=1,NDI
           DO K1=1,NDI
              DO L1=1,NDI
                XX=CISOAUX(I1,J1,K1,L1)
                YY=TRFIC*PE(I1,J1,K1,L1)
                ZZ=SISO(I1,J1)*UNIT2(K1,L1)+UNIT2(I1,J1)*SISO(K1,L1)
C                
                CISO(I1,J1,K1,L1)=XX+(TWO/THREE)*YY-(TWO/THREE)*ZZ
              END DO
           END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE SETISO
      SUBROUTINE PUSH2(SIG,PK,F,DET,NDI)
C>        PIOLA TRANSFORMATION
C>      INPUT:
C>       PK - 2ND PIOLA KIRCHOOF STRESS TENSOR
C>       F - DEFORMATION GRADIENT
C>       DET - DEFORMATION DETERMINANT
C>      OUTPUT:
C>       SIG - CAUCHY STRESS TENSOR
       IMPLICIT NONE
       INCLUDE 'PARAM_UMAT.INC'
C
       INTEGER I1,J1,II1,JJ1,NDI
       DOUBLE PRECISION PK(NDI,NDI),F(NDI,NDI)
       DOUBLE PRECISION SIG(NDI,NDI)
       DOUBLE PRECISION AUX,DET
C
       DO I1=1,NDI
        DO J1=1,NDI
          AUX=ZERO
         DO II1=1,NDI
          DO JJ1=1,NDI
            AUX=AUX+(DET**(-ONE))*F(I1,II1)*F(J1,JJ1)*PK(II1,JJ1)
         END DO
        END DO
           SIG(I1,J1)=AUX
        END DO
       END DO
C
       RETURN
      END SUBROUTINE PUSH2
      SUBROUTINE DEFORMATION(F,C,B,NDI)
C>     RIGHT AND LEFT CAUCHY-GREEN DEFORMATION TENSORS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION F(NDI,NDI),C(NDI,NDI),B(NDI,NDI)
C     RIGHT CAUCHY-GREEN DEFORMATION TENSOR
      C=MATMUL(TRANSPOSE(F),F)
C     LEFT CAUCHY-GREEN DEFORMATION TENSOR
      B=MATMUL(F,TRANSPOSE(F))
      RETURN
      END SUBROUTINE DEFORMATION
      SUBROUTINE AFFCLMATNETFIC(PKFIC,CMFIC,F,MF0,RW,FILPROPS,AFFPROPS,
     1                        EFI,NOEL,DET,NDI)
C>    AFFINE NETWORK: 'FICTICIOUS' PK STRESS AND ELASTICITY TENSOR 
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1,M1,NOEL
      DOUBLE PRECISION PKFIC(NDI,NDI),SFILFIC(NDI,NDI),
     1                  CMFIC(NDI,NDI,NDI,NDI),F(NDI,NDI),MF0(NWP,NDI),
     2                 RW(NWP),CFILFIC(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION FILPROPS(8),AFFPROPS(2),MFI(NDI),MF0I(NDI)
      DOUBLE PRECISION DET,AUX,PI,LAMBDAI,DWI,DDWI,RWI,LAMBDAIC
      DOUBLE PRECISION L,R0F,R0,MU0,B0,BETA,LAMBDA0,RHO,N,FI,FFI,DTIME
      DOUBLE PRECISION R0C,ETAC,LAMBDAIF
      DOUBLE PRECISION B,FRIC,FFMAX,ANG,EFI,FRAC(4),RU0(NWP),RU
      DOUBLE PRECISION VARA,AVGA,MAXA,AUX0,FFIC,SUMA,RHO0,DIRMAX(NDI)

C
C     FILAMENT
      L       = FILPROPS(1)
      R0F     = FILPROPS(2)
      R0C     = FILPROPS(3)
      ETAC    = FILPROPS(4)      
      MU0     = FILPROPS(5)
      BETA    = FILPROPS(6)
      B0      = FILPROPS(7)
      LAMBDA0 = FILPROPS(8)
C     NETWORK
      N       = AFFPROPS(1)
      B       = AFFPROPS(2)      
C
      PI=FOUR*ATAN(ONE)
      AUX=N*(DET**(-ONE))*FOUR*PI
      CMFIC=ZERO
      PKFIC=ZERO
C      
       RHO=ONE
       R0=R0F+R0C
C  
C        LOOP OVER THE INTEGRATION DIRECTIONS
      DO I1=1,NWP
C
       MFI=ZERO
       MF0I=ZERO
       DO J1=1,NDI
        MF0I(J1)=MF0(I1,J1)
       END DO
       RWI=RW(I1)
C
       CALL DEFFIL(LAMBDAI,MFI,MF0I,F,NDI)
C
       CALL BANGLE(ANG,F,MFI,NOEL,NDI)
C      
       CALL DENSITY(RHO,ANG,B,EFI)
C
       IF((ETAC.GT.ZERO).AND.(ETAC.LT.ONE))THEN
C
        LAMBDAIF=ETAC*(R0/R0F)*(LAMBDAI-ONE)+ONE
        LAMBDAIC=(LAMBDAI*R0-LAMBDAIF*R0F)/R0C
       ELSE
        LAMBDAIF=LAMBDAI
        LAMBDAIC=ZERO
       ENDIF  
C
       CALL FIL(FI,FFI,DWI,DDWI,LAMBDAIF,LAMBDA0,L,R0,MU0,BETA,B0)
C       
       CALL SIGFILFIC(SFILFIC,RHO,LAMBDAIF,DWI,MF0I,RWI,NDI)
C
       CALL CSFILFIC(CFILFIC,RHO,LAMBDAIF,DWI,DDWI,MF0I,RWI,NDI)
C
       DO J1=1,NDI
        DO K1=1,NDI
         PKFIC(J1,K1)=PKFIC(J1,K1)+AUX*SFILFIC(J1,K1)
         DO L1=1,NDI
          DO M1=1,NDI
          CMFIC(J1,K1,L1,M1)=CMFIC(J1,K1,L1,M1)+AUX*CFILFIC(J1,K1,L1,M1)
          END DO
         END DO
        END DO
       END DO
C
      END DO
C
      RETURN
      END SUBROUTINE AFFCLMATNETFIC
      SUBROUTINE FSLIP(F,FBAR,DET,NDI)
C>     DISTORTION GRADIENT
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1
      DOUBLE PRECISION F(NDI,NDI),FBAR(NDI,NDI)
      DOUBLE PRECISION DET,SCALE1
C     
C     JACOBIAN DETERMINANT
      DET = F(1,1) * F(2,2) * F(3,3)
     1    - F(1,2) * F(2,1) * F(3,3)
C
      IF (NDI .EQ. 3) THEN
          DET = DET + F(1,2) * F(2,3) * F(3,1)
     1              + F(1,3) * F(3,2) * F(2,1)
     2              - F(1,3) * F(3,1) * F(2,2)
     3              - F(2,3) * F(3,2) * F(1,1)
      END IF 
C
      SCALE1=DET**(-ONE /THREE)
C      
      DO I1=1,NDI
        DO J1=1,NDI
          FBAR(I1,J1)=SCALE1*F(I1,J1)
        END DO
      END DO
C
      RETURN      
      END SUBROUTINE FSLIP
      SUBROUTINE SDVWRITE(DET,STATEV)
C>    VISCOUS DISSIPATION: WRITE STATE VARS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER VV,POS1,POS2,POS3,I1
      DOUBLE PRECISION STATEV(NSDV),DET
C
        STATEV(1)=DET
      RETURN
C
      END SUBROUTINE SDVWRITE
      SUBROUTINE BANGLE(ANG,F,MF,NOEL,NDI)
C>    ANGLE BETWEEN FILAMENT AND PREFERED DIRECTION
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C      
      COMMON /KFILP/PREFDIR
      DOUBLE PRECISION PREFDIR(NELEM,4)
C
      INTEGER INOEL,I,J,NDI,NOEL
      DOUBLE PRECISION DNORM,PDIR(NDI),ANG,MF(NDI),MFA(NDI),AUX
      DOUBLE PRECISION F(NDI,NDI),C(NDI,NDI),EGVC(NDI,NDI),EGVL(NDI)
      
C
        INOEL=0
        I=0
        DO I=1,NELEM
C               ELEMENT IDENTIFICATION
            IF(NOEL.EQ.INT(PREFDIR(I,1))) THEN
                INOEL=I
            ENDIF
        ENDDO
C
       DO I=1,NDI
        J=I+1
C       PREFERED ORIENTATION  ORIENTATION NORMALIZED 
        PDIR(I)=PREFDIR(INOEL,J)
       END DO 
C        ALTERNATIVE APPROACH: BUNDLES FOLLOW PRINCIPAL DIRECTIONS       
       C=MATMUL(TRANSPOSE(F),F)
       CALL SPECTRAL(C,EGVL,EGVC)
C       WRITE(*,*) EGVC
       PDIR(1)=EGVC(1,1)
       PDIR(2)=EGVC(2,1)
       PDIR(3)=EGVC(3,1)
C        END OF ALTERNATIVE      
C       
C     PREFERED ORIENTATION 
        DNORM=DOT_PRODUCT(PDIR,PDIR)
        DNORM=DSQRT(DNORM)
C
C       PREFERED ORIENTATION  NORMALIZED 
        PDIR=PDIR/DNORM 
C
C       FILAMENT ORIENTATION 
        MFA=MF
        DNORM=DOT_PRODUCT(MFA,MFA)
        DNORM=DSQRT(DNORM)
C
C       FILAMENT ORIENTATION  NORMALIZED 
        MFA=MFA/DNORM
C        ANGLE BETWEEN PREFERED ORIENTATION AND FILAMENT - BANGLE        
        AUX=DOT_PRODUCT(MFA,PDIR)
        
!        if AUX.GT.ONE
!        endif
!        write(*,*) aux
        ANG=ACOS(AUX)
C
      RETURN
      END SUBROUTINE BANGLE

      SUBROUTINE INDEXX(STRESS,DDSDDE,SIG,TNG,NTENS,NDI)
C>    INDEXATION: FULL SIMMETRY  IN STRESSES AND ELASTICITY TENSORS
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER II1(6),II2(6),NTENS,NDI,I1,J1
      DOUBLE PRECISION STRESS(NTENS),DDSDDE(NTENS,NTENS)
      DOUBLE PRECISION SIG(NDI,NDI),TNG(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION PP1,PP2
C
      II1(1)=1
      II1(2)=2
      II1(3)=3
      II1(4)=1
      II1(5)=1
      II1(6)=2
C
      II2(1)=1
      II2(2)=2
      II2(3)=3
      II2(4)=2
      II2(5)=3
      II2(6)=3
C
      DO I1=1,NTENS
C       STRESS VECTOR
         STRESS(I1)=SIG(II1(I1),II2(I1))
         DO J1=1,NTENS
C       DDSDDE - FULLY SIMMETRY IMPOSED
            PP1=TNG(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=TNG(II1(I1),II2(I1),II2(J1),II1(J1))
            DDSDDE(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
C
      RETURN
C
      END SUBROUTINE INDEXX
      SUBROUTINE ROTATION(F,R,U,NDI)
C>    COMPUTES ROTATION TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION F(NDI,NDI),R(NDI,NDI),U(NDI,NDI),UINV(NDI,NDI)
C
      CALL MATINV3D(U,UINV,NDI)
C
      R = MATMUL(F,UINV)
      RETURN
      END SUBROUTINE ROTATION
      SUBROUTINE CMATISOMATFIC(CMISOMATFIC,CBAR,CBARI1,CBARI2,
     1                          DISO,UNIT2,UNIT4,DET,NDI)
C>    ISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI,I1,J1,K1,L1
      DOUBLE PRECISION CMISOMATFIC(NDI,NDI,NDI,NDI),UNIT2(NDI,NDI),
     1                 CBAR(NDI,NDI),DISO(5),
     2                 UNIT4(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION CBARI1,CBARI2
      DOUBLE PRECISION DUDI1,DUDI2,D2UD2I1,D2UD2I2,D2UDI1I2
      DOUBLE PRECISION AUX,AUX1,AUX2,AUX3,AUX4,DET
      DOUBLE PRECISION UIJ,UKL,CIJ,CKL
      
C
      DUDI1=DISO(1)
      DUDI2=DISO(2)
      D2UD2I1=DISO(3)
      D2UD2I2=DISO(4)
      D2UDI1I2=DISO(5)
C
      AUX1=FOUR*(D2UD2I1+TWO*CBARI1*D2UDI1I2+
     1           DUDI2+CBARI1*CBARI1*D2UD2I2)
      AUX2=-FOUR*(D2UDI1I2+CBARI1*D2UD2I2)
      AUX3=FOUR*D2UD2I2
      AUX4=-FOUR*DUDI2

      DO I1=1,NDI
        DO J1=1,NDI
           DO K1=1,NDI
              DO L1=1,NDI
                     UIJ=UNIT2(I1,J1)
                     UKL=UNIT2(K1,L1)
                     CIJ=CBAR(I1,J1)
                     CKL=CBAR(K1,L1)
                     AUX=AUX1*UIJ*UKL+
     1                   AUX2*(UIJ*CKL+CIJ*UKL)+AUX3*CIJ*CKL+
     3                   AUX4*UNIT4(I1,J1,K1,L1)
                     CMISOMATFIC(I1,J1,K1,L1)=AUX        
              END DO
           END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE CMATISOMATFIC
      SUBROUTINE PROJLAG(C,AA,PL,NDI)
C>    LAGRANGIAN PROJECTION TENSOR
C      INPUTS:
C          IDENTITY TENSORS - A, AA
C          ISOCHORIC LEFT CAUCHY GREEN TENSOR - C
C          INVERSE OF C - CINV
C      OUTPUTS:
C          4TH ORDER SYMMETRIC LAGRANGIAN PROJECTION TENSOR - PL
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER I,J,K,L,NDI
C
      DOUBLE PRECISION CINV(NDI,NDI),AA(NDI,NDI,NDI,NDI),
     1                 PL(NDI,NDI,NDI,NDI),C(NDI,NDI)
C
      CALL MATINV3D(C,CINV,NDI)
C
      DO I=1,NDI
         DO J=1,NDI
          DO K=1,NDI
             DO L=1,NDI
              PL(I,J,K,L)=AA(I,J,K,L)-(ONE/THREE)*(CINV(I,J)*C(K,L))
           END DO
          END DO
        END DO
      END DO
C
      RETURN
      END SUBROUTINE PROJLAG
      SUBROUTINE DENSITY(RHO,ANG,BB,ERFI)
C>    SINGLE FILAMENT: DENSITY FUNCTION VALUE
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION RHO,ANG,BB,ERFI
      DOUBLE PRECISION PI,AUX1,AUX2
C
      PI=FOUR*ATAN(ONE)
      AUX1=SQRT(BB/(TWO*PI))
      AUX2=DEXP(BB*(COS(TWO*ANG)+ONE))
      RHO=FOUR*AUX1*AUX2*(ERFI**(-ONE))
C      RHO=RHO*((FOUR*PI)**(-ONE))
C
      RETURN
      END SUBROUTINE DENSITY
      SUBROUTINE CSISOMATFIC(CISOMATFIC,CMISOMATFIC,DISTGR,DET,NDI)
C>    ISOTROPIC MATRIX: SPATIAL 'FICTICIOUS' ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      INTEGER NDI
      DOUBLE PRECISION CMISOMATFIC(NDI,NDI),DISTGR(NDI,NDI),
     1                 CISOMATFIC(NDI,NDI,NDI,NDI)
      DOUBLE PRECISION DET
      
C
      CALL PUSH4(CISOMATFIC,CMISOMATFIC,DISTGR,DET,NDI)
C
      RETURN
      END SUBROUTINE CSISOMATFIC
