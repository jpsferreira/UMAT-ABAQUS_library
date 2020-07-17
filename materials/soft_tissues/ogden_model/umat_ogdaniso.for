C********************************************************************
C Record of revisions:                                              |
C        Date        Programmer        Description of change        |
C        ====        ==========        =====================        |
C     15/02/2015    Joao Ferreira      adapted from Jabaren eqs     |
C     23/06/2015    Joao Ferreira      anisotropy added             |
C     04/12/2015    Joao Ferreira      ddsdde corrected             |
C     14/03/2016    Joao Ferreira      k dispersion added           |
C--------------------------------------------------------------------
C     Description:
C     UEXTERNALDB: mount database; fibers directions;.inc, .inp files
C     UMAT: ogden matrix + Holzapfel fibers + k dispersion
C--------------------------------------------------------------------
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'PARAM_UMAT.INC'
C       this subroutine get the fibers directions resorting the
C      .inc files and read the fibers direction
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
      DIMENSION TIME(2)
      REAL*8 FIBORI(NELEM,4),FIBORI2(NELEM,4)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
        CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI(I,J),J=1,4)
          END DO
           CLOSE(15)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR2
          OPEN(15,FILE=FILENAME)
          DO I=1,NELEM
             READ(15,*) (FIBORI2(I,J),J=1,4)
          END DO
           CLOSE(15)
      END IF
C
      RETURN
C
      END
C
C---------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      COMMON /KFIB/FIBORI
      COMMON /KFIB/FIBORI2
C
      CHARACTER*8 CMNAME
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     3 DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),DDSDDE(NTENS,NTENS),
     4 FIBORI(NELEM,4),FIBORI2(NELEM,4)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1   DTEMP,PNEWDT,CELENT
C
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1 LAYER, KSPT, KSTEP, KINC
C
C---------------------------------------------------------------------
C     LOCAL ARRAYS
C---------------------------------------------------------------------
C     STATEV - STATE VARIABLES ARRAY
C     UNIT   - IDENTITY TENSOR
C     II1    - POINTERS VECTOR
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C     C      -  RIGHT CAUCHY-GREEN TENSOR
C     B      -  LEFT CAUCHY-GREEN TENSOR
C     CBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     BBAR   - DEVIATORIC LEFT CAUCHY-GREEN TENSOR
C     B2BAR  - SQUARE OF THE DEVIATORIC LEFT CAUCHY-GREEN TENSOR
C     VORIF  - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 1)
C     VORIF2 - FIBER ORIENTATION IN UNDERFORMED CONFIGURATION (FAMILY 2)
C     M0,N0  - STRUCTURAL FIBERS TENSOR IN UNDEFORMED  CONFIGURATION
C     M,MM,N,NN - STRUCTURAL FIBERS TENSOR IN DEFORMED  CONFIGURATION
C     VD     - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 1)
C     VD2    - FIBER ORIENTATION IN DERFORMED CONFIGURATION (FAMILY 2)
C     SSX    - AUXILIAR TENSORS
C     SX     - CAUCHY STRESSES OF EACH CONTRIBUTION: SVOL, SM, SF, SDEV
C     SIGMA  - TOTAL CUCHY STRESS
C     DDSDDEXX - JACOBIAN MATERIAL CONTRIBUTIONS
C---------------------------------------------------------------------
C     LOCAL VARIABLES
C---------------------------------------------------------------------
C     DET    - DETERMINANT OF THE DEFORMATION GRADIENT
C     SCALEX - SCALES OF DET
C     CHECK  - STABILITY CHECK VARIABLE
C     CBARI1 - 1ST INVARIANT OF CBAR
C     CBARI2 - 2ND INVARIANT OF CBAR
C     DNORM  - NORM OF THE FIBER DIRECTION VECTOR
C     CBARI4 - SQUARE OF THE FIBER1 STRETCH
C     BARLAMBDA - DEVIATORIC STRETCH OF FIBER 1
C     LAMBDA - STRETCH OF FIBER 1
C     CBARI6 - SQUARE OF THE FIBER2 STRETCH
C     BARLAMBDA2 - DEVIATORIC STRETCH OF FIBER 2
C     LAMBDA2 - STRETCH OF FIBER 2
C     SEV    - VOLUMETRIC PART OF THE STRAIN-ENERGY
C     SEM    - ISOTROPIC PART OF THE STRAIN-ENERGY
C     SEF    - ANISOTROPIC PART OF THE STRAIN-ENERGY
C     SSE    - STRAIN-ENERGY
C     D2DUDXXX - 2ND DERIVATIVES OF SSE RESPECTIVE TO THE INVARIANTS
C     DUDXX  - 1ST DERIVATIVES OF SSE RESPECTIVE TO THE INVARIANTS
C              ***TO BE ADDED LATER ON....***
C     PJMMAX - EQUIVALENT STRAIN FOR MATRIX (SIMO TYPE)
C     PJFMAX - EQUIVALENT STRAIN FOR FIBERS (SIMO TYPE)
C     PJ     - CURRENT EQUIVALENT STRAIN
C     ZD     - DAMAGE CRITERIA
C     DSTAT  - DAMAGE FLAG
C     GBAR   - REDUCTION FACTOR
C     GDERM  - DERIVATIVE OF GBAR FOR MATRIX
C     GDERF  - DERIVATIVE OF GBAR FOR FIBERS
C     DMGM   - MATRIX DAMAGE
C     DMGF   - FIBERS DAMAGE
C---------------------------------------------------------------------
C
C                COUNTERS AND STATE INDICATORS
      INTEGER  II1(6),II2(6),ISTAT,INOEL,I,J,K,I1,J1,K1,K2,L,
     1         NEIG,NEGN
C
C                KINEMATIC ARRAY VARIABLES
      DOUBLE PRECISION  UNIT(3,3),C(3,3),B(3,3),CBAR(3,3),BBAR(3,3),
     1        C2BAR(3,3),B2BAR(3,3),DISTGR(3,3),VCBAR(6),PS(3),
     1        AN(3,3),DD(3),
C                FIBERS STRUCTURE VARIABLES
     2        VORIF(3),VD(3),M(3,3),M0(3,3),MM(3,3),
     3        VORIF2(3),VD2(3),N(3,3),N0(3,3),NN(3,3),
C                STRESS VARIABLES
     4        SS1(3,3),SS2(3,3),SS3(3,3),SS4(3,3),SS5(3,3),
     5        S1(3,3),S2(3,3),S3(3,3),S4(3,3),S5(3,3), SS6(3,3),
     6        S6(3,3),SS8(3,3),SS7(3,3),S7(3,3),S8(3,3),SIGMA(3,3),
C                MATERIAL TANGENT VARIABLES
     7        DDSIGDDE(3,3,3,3),DDSDDEJR(3,3,3,3),DDSDDEISO(3,3,3,3),
     8        DDSDDEVOL(3,3,3,3),DDSDDEFIB1(3,3,3,3),
     9        DDSDDEFIB2(3,3,3,3),DDSDDEFIB3(3,3,3,3)
C
C                MATERIAL PARAMETERS VARIABLES
      DOUBLE PRECISION  D1,C10,C01,C3,C4,C5,C6,KAPPA,
C                KINEMATIC SCALAR VARIABLES
     1        DET,SCALE,SCALE0,SCALE1,SCALE2,SCALE3,SUM,SUMB,
     2        CI1,CBARI1,C2BARI1,CBARI2,E1,EE1,EE2,EE3,PP1,PP2,A,
     3        LAMBDA,AUX0,AUX1,AUX2,AUX3,UNIT4,SUM1,LAMBDA2,
     4        DNORM,BARLAMBDA,CBARI4,DNORM2,SUM2,CBARI6,BARLAMBDA2,
C                STRAIN-ENERGY DERIVATIVES VARIABLES
     5        DUDI1,DUDI2,D2UD2I2,D2DUDI1DI2,D2DUDJDI1,DUDJ,D2UDJ2,
     6        D2DUDJDI2,DUDI4,DUDLAMBDA,D2UD2LAMBDA,D2DUDJDI4,
     7        D2DUDJDLAMBDA,D2UD2I1,D2DUDI1DLAMBDA,D2DUDI2DLAMBDA,
     8        D2DUDI4DLAMBDA,D2DUDI2DI4,D2DUDI1DI4,D2UD2I4,DUDI6,DUDI8,
     9        D2DUDJDI6,D2DUDJDI8,D2DUDJDLAMBDA2,D2DUDI1DI6,D2DUDI2DI6,
     9        D2DUDI1DLAMBDA2,D2DUDI2DLAMBDA2,D2UD2I6,DUDLAMBDA2,
     9        D2UD2LAMBDA2,D2DUDI6DLAMBDA2,D2UD2I8,D2UDJDI8,D2UDI1DI8,
     9	       D2UDI2DI8,D2DUDI1DI8,D2DUDI2DI8,
     9        CHECK12,CHECK13,CHECK23,CHECKD,CHECK12D,CHECK13D,CHECK23D,
     9        TOL,MINCHECK,AA,BB,COEF,DI1DA,DI1DB,DI2DA,DI2DB,D,
     9        DADI1,DADI2,D2DAD2I1,D2DAD2I2,D2DADI1DI2,DBDI1,DBDI2
     9        GAMMA,D2DI1DA2,D2DI1DADB,D2DI1DB2,D2DI2DA2,D2DI2DB2,
     9        D2DI2DADB,V1,V2,V3,V4,V5,V6,V7,D2DBDI1DI2,D2DBD2I1,
     9        D2DBD2I2,GAMMA,GAMMA1,GAMMA2,GAMMA3,GAMMA4,GAMMA01,
     9        GAMMA02,GAMMA03,GAMMA04,GAMMA05,DUDA,DUDB,D2DUDA,D2DUDADB,
     9        DBDI2,D2DUDB,C11,C20,C02,CI2
C----------------------------------------------------------------------
C     MATERIAL CONSTANTS
C----------------------------------------------------------------------
C      CXX = PROPS(X)
      D1  = PROPS(NPROPS)
      C3  = PROPS(NPROPS-5)
      C4  = PROPS(NPROPS-4)
      C5  = PROPS(NPROPS-3)
      C6  = PROPS(NPROPS-2)
      KAPPA = PROPS(NPROPS-1)
C----------------------------------------------------------------------
C     INITIALIZATIONS
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      ISTAT=1
      DET=ZERO
      SCALE=ZERO
      SCALE0=ZERO
      SCALE1=ZERO
      SCALE2=ZERO
      SCALE3=ZERO
      SUM=ZERO
      SUMB=ZERO
      CI1=ZERO
      CBARI1=ZERO
      C2BARI1=ZERO
      CBARI2=ZERO
      DUDI1=ZERO
      D2UD2I1=ZERO
      DUDI2=ZERO
      D2UD2I2=ZERO
      D2DUDI1DI2=ZERO
      D2DUDJDI1=ZERO
      D2DUDJDI2=ZERO
      DUDJ=ZERO
      D2UDJ2=ZERO
      AUX0=ZERO
      AUX1=ZERO
      AUX2=ZERO
      AUX3=ZERO
      UNIT4=ZERO
C
      LAMBDA=ONE
      LAMBDA2=ONE
      DNORM=ZERO
      DNORM2=ZERO
      BARLAMBDA=ZERO
      CBARI4=ZERO
      BARLAMBDA2=ZERO
      CBARI6=ZERO
      DUDLAMBDA=ZERO
      D2UD2LAMBDA=ZERO
      D2DUDJDI4=ZERO
      D2DUDJDLAMBDA=ZERO
      D2DUDI1DLAMBDA=ZERO
      D2DUDI2DLAMBDA=ZERO
      D2DUDI4DLAMBDA=ZERO
      D2DUDI2DI4=ZERO
      D2DUDI1DI4=ZERO
      D2UD2I4=ZERO
C
      PP1=ZERO
      PP2=ZERO
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     IDENTITY TENSOR DEFINITION                                       *
C----------------------------------------------------------------------
C
      CAll ONEM(UNIT)
C
C----------------------------------------------------------------------
C     POINTERS DEFINITION
C----------------------------------------------------------------------
C     POINTERS TO STORAGE STRESS AND DDSDDE ARRAYS
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
C----------------------------------------------------------------------
C     JACOBIAN AND DISTORTION TENSOR
C----------------------------------------------------------------------
C     JACOBIAN
      DET = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
C
      IF (NSHR .EQ. 3) THEN
          DET = DET + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     1              + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     2              - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     3              - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
      END IF
C     AUXILIAR VARIABLES
      SCALE0=DET**(-ONE)
      SCALE1=(ONE/DET)
      SCALE = DET**(-ONE /THREE)
      SCALE2 = SCALE**(TWO)
      SCALE3 = SCALE2**(TWO)
C
      DO K1 = 1, 3
        DO K2 = 1, 3
          DISTGR(K2,K1) = SCALE * DFGRD1(K2,K1)
        END DO
      END DO
C
C----------------------------------------------------------------------
C     SUBROUTINE STABILITY VERIFICATION
C----------------------------------------------------------------------
C
C	  CHECK=DABS(DFGRD1(1,1))+DABS(DFGRD1(1,2))+DABS(DFGRD1(1,3)) +
C     +      DABS(DFGRD1(2,1))+DABS(DFGRD1(2,2))+DABS(DFGRD1(2,3)) +
C     +      DABS(DFGRD1(3,1))+DABS(DFGRD1(3,2))+DABS(DFGRD1(3,3))
C
C	  IF(DET.GT.1.10D0.OR.DET.LT.0.90D0.OR.CHECK.GT.10.0D0) THEN
C	    WRITE(7,*) "NOEL=",NOEL
C	    WRITE(*,*) "DET=",DET,CHECK
C		WRITE(7,*) DDSDDE
C	    WRITE(7,*) "STRESS=",STRESS
C		DDSDDE=0.0D0
C		STRESS=0.0D0
C		PNEWDT=0.5D0
C		RETURN
C	  ELSE
C		CONTINUE
C	  ENDIF
C
C----------------------------------------------------------------------
C     RIGHT AND LEFT CAUCHY-GREEN TENSOR AND INVARIANTS
C----------------------------------------------------------------------
C
C     RIGHT AND LEFT CAUCHY-GREEN TENSORS
      DO I1=1,NDI
         DO J1=1,NDI
            SUM=ZERO
            SUMB=ZERO
            DO K1=1,NDI
               SUM=SUM+DFGRD1(K1,I1)*DFGRD1(K1,J1)
               SUMB=SUMB+DFGRD1(I1,K1)*DFGRD1(J1,K1)
            END DO
            C(I1,J1)=SUM
            B(I1,J1)=SUMB
            CBAR(I1,J1)=SCALE2*C(I1,J1)
            BBAR(I1,J1)=SCALE2*B(I1,J1)
         END DO
         CI1=CI1+C(I1,I1)
      END DO
C    FIRST INVARIANT OF CBAR
      CBARI1=SCALE2*CI1
C     SQUARE OF DEFORMATION TENSORS
      C2BAR=MATMUL(CBAR,CBAR)
      B2BAR=MATMUL(BBAR,BBAR)
      DO I1=1,NDI
        C2BARI1=C2BARI1+C2BAR(I1,I1)
      END DO
        CBARI2=(ONE/TWO)*(CBARI1*CBARI1-C2BARI1)
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE UNDEFORMED FIBERS
C----------------------------------------------------------------------
C
        INOEL=0
        I=0
        DO I=1,NELEM
C               ELEMENT IDENTIFICATION
            IF(NOEL.EQ.INT(FIBORI(I,1))) THEN
                INOEL=I
            ENDIF
        ENDDO
C
C     FIBORI - FIBER ORIENTATION - FAMILY 1
             DNORM=DSQRT(FIBORI(INOEL,2)*FIBORI(INOEL,2)+
     1                   FIBORI(INOEL,3)*FIBORI(INOEL,3)+
     2                   FIBORI(INOEL,4)*FIBORI(INOEL,4))
C
C      FIBORI2 - FIBER ORIENTATION - FAMILY 2
             DNORM2=DSQRT(FIBORI2(INOEL,2)*FIBORI2(INOEL,2)+
     1                   FIBORI2(INOEL,3)*FIBORI2(INOEL,3)+
     2                   FIBORI2(INOEL,4)*FIBORI2(INOEL,4))
C
C       UNDERFORMED FIBER ORIENTATION TENSOR
C
        DO I=1,NDI
        J=I+1
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 1
        VORIF(I)=FIBORI(INOEL,J)/DNORM
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 2
        VORIF2(I)=FIBORI2(INOEL,J)/DNORM2
        END DO
C
      DO I=1,NDI
       DO J=1,NDI
C       STRUCTURAL TENSOR - FAMILY 1
       M0(I,J)=VORIF(I)*VORIF(J)
C       STRUCTURAL TENSOR - FAMILY 2
       N0(I,J)=VORIF2(I)*VORIF2(J)
       END DO
      END DO
C
C----------------------------------------------------------------------
C     STRETCH RATIO IN THE FIBER DIRECTION
C----------------------------------------------------------------------
C
        CBARI4=ZERO
        CBARI6=ZERO
      DO I=1,NDI
        DO J=1, NDI
            CBARI4=CBARI4+CBAR(I,J)*M0(I,J)
            CBARI6=CBARI6+CBAR(I,J)*N0(I,J)
        ENDDO
      ENDDO
C     FAMILY 1
      BARLAMBDA=DSQRT(CBARI4)
      LAMBDA=BARLAMBDA/SCALE
C
C     FAMILY 2
      BARLAMBDA2=DSQRT(CBARI6)
      LAMBDA2=BARLAMBDA2/SCALE
C
C----------------------------------------------------------------------
C     UNIT VECTOR IN THE DIRECTION OF THE DEFORMED FIBERS
C----------------------------------------------------------------------
C
      DO I=1,NDI
         SUM1=ZERO
         SUM2=ZERO
         DO J=1,NDI
          SUM1=SUM1+DFGRD1(I,J)*VORIF(J)
          SUM2=SUM2+DFGRD1(I,J)*VORIF2(J)
         ENDDO
C     FIBER DIRECTIONS IN THE DEFORMED CONFIGURATION
C               -FAMILY 1
         VD(I)=SUM1
C
C               -FAMILY 2
         VD2(I)=SUM2
      ENDDO
      DNORM=DSQRT(VD(1)*VD(1)+
     1             VD(2)*VD(2)+
     2             VD(3)*VD(3))
      DNORM2=DSQRT(VD2(1)*VD2(1)+
     1             VD2(2)*VD2(2)+
     2             VD2(3)*VD2(3))
C           COSINE OF THE ANGLE BETWEEN FIBERS
      DO I=1,NDI
       VD(I)=VD(I)/DNORM
       VD2(I)=VD2(I)/DNORM2
       A=A+VD(I)*VD2(I)
      END DO
C
C
C       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 1
      M=MATMUL(M0,TRANSPOSE(DFGRD1))
      M=MATMUL(DFGRD1,M)
      MM=MATMUL(M0,TRANSPOSE(DISTGR))
      MM=MATMUL(DISTGR,MM)
C
C       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 2
      N=MATMUL(N0,TRANSPOSE(DFGRD1))
      N=MATMUL(DFGRD1,N)
      NN=MATMUL(N0,TRANSPOSE(DISTGR))
      NN=MATMUL(DISTGR,NN)
C
C----------------------------------------------------------------------
C     STRAIN-ENERGY DERIVATIVES WITH RESPECT TO INVARIANTS
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES UI - VOLUME-PRESERVING PART
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1 AND J
       D2DUDJDI1=ZERO
C  SECOND DERIVATIVE UI  WITH RESPECT TO I2 AND J
       D2DUDJDI2=ZERO
C  FIRST DERIVATIVE UJ  WITH RESPECT TO J
       DUDJ=(ONE/D1)*(DET-ONE/DET)
C  SECOND DERIVATIVE UJ  WITH RESPECT TO J
      D2UDJ2=(ONE/D1)*(ONE+ONE/(DET*DET))
      D2DUDJDI4=ZERO
      D2DUDJDLAMBDA=ZERO
      D2DUDJDI6=ZERO
      D2DUDJDLAMBDA2=ZERO
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES UI - ISOTROPIC CONTRIBUTION
C     EIGENVALUES AND EIGENVECTORS OF CBAR
       DO I1=1,6
         VCBAR(I1)=C(II1(I1),II2(I1))
       END DO
      CALL SPRIND(VCBAR, PS, AN, 2, NDI, NSHR)
C      CALL SPECTRAL(CBAR,PS,AN,1)
C
C
          CHECK12=ABS(PS(1)-PS(2))
          CHECK13=ABS(PS(1)-PS(3))
          CHECK23=ABS(PS(2)-PS(3))
C
          CHECKD=PS(1)*PS(1)+PS(2)*PS(2)+PS(3)*PS(3)
          CHECKD=DSQRT(CHECKD)
C
          CHECK12D=CHECK12/CHECKD
          CHECK13D=CHECK13/CHECKD
          CHECK23D=CHECK23/CHECKD
C
C           NEIG IS THE NUMBER OF EQUAL EIGENVALUES
          TOL=1.D-4
          MINCHECK=CHECK12
C
            IF((CHECK12D.LT.TOL).OR.
     &         (CHECK13D.LT.TOL).OR.
     &         (CHECK23D.LT.TOL)) THEN
C          CHECK FOR THREE EQUAL EIGENVALUES
            IF ((CHECK12D.LT.TOL).AND.
     &              (CHECK13D.LT.TOL).AND.
     &                  (CHECK23D.LT.TOL)) THEN
                NEIG=3
             ELSE
C          CHECK FOR TWO EQUAL EIGENVALUES
                NEIG=2
                 AA=(MINCHECK/TWO)**TWO
                 BB=(PS(1)+PS(2))/TWO
                IF (CHECK13D .LT. MINCHECK) THEN
                MINCHECK=CHECK13
                AA=(MINCHECK/TWO)**TWO
                BB=(PS(1)+PS(3))/TWO
                ELSEIF (CHECK23D .LT. MINCHECK) THEN
                MINCHECK=CHECK23
                AA=(MINCHECK/TWO)**TWO
                BB=(PS(2)+PS(3))/TWO
                END IF
            END IF
            ELSE
C          CHECK FOR THREE DIFFERENT EIGENVALUES
              NEIG=1
           ENDIF
C
C     FIRST AND SECOND DERIVATIVES UI
C
C     DO NOT FORGET TO ADAPT NEGN FOR EACH FIBER STRAIN-ENERGY FUNCTION!
      NEGN=(NPROPS-5-1)/2
C
       IF (NEIG==1) THEN
C---------------------------------
C      THREE DIFFERENT EIGENVALUES
C---------------------------------
      DD(1)=(PS(1)-PS(2))*(PS(1)-PS(3))
      DD(2)=(PS(2)-PS(1))*(PS(2)-PS(3))
      DD(3)=(PS(3)-PS(1))*(PS(3)-PS(2))
      DO I1=1,NEGN
        GAMMA=(ONE/TWO)*PROPS(2*I1)
        COEF=(TWO*PROPS(2*I1-1)/(FOUR*GAMMA**TWO))
        DO K1=1,NDI
C  FIRST DERIVATIVE UI WITH RESPECT TO I1
         DUDI1=DUDI1+COEF*GAMMA*((PS(K1)**(GAMMA+ONE))/DD(K1))
C  FIRST DERIVATIVE UI WITH RESPECT TO I2
       DUDI2=DUDI2-COEF*GAMMA*((PS(K1)**(GAMMA))/DD(K1))
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1
       D2UD2I1=D2UD2I1
     +      +COEF*GAMMA*(GAMMA+ONE)*
     +      ((PS(K1)**(GAMMA+TWO))/(DD(K1)**TWO))
     +      +COEF*TWO*GAMMA*((PS(K1)**(GAMMA+ONE))
     +      *(1-PS(K1)**THREE)/(DD(K1)**THREE))
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1 AND I2
       D2DUDI1DI2=D2DUDI1DI2
     +      -COEF*
     +       GAMMA*GAMMA*((PS(K1)**(GAMMA+ONE))/(DD(K1)**TWO))
     +      -COEF*TWO*GAMMA*((PS(K1)**(GAMMA))
     +      *(1-PS(K1)**THREE)/(DD(K1)**THREE)
     +              )
C  SECOND DERIVATIVE UI  WITH RESPECT TO I2
       D2UD2I2=D2UD2I2
     +  +COEF*
     +      GAMMA*GAMMA*((PS(K1)**(GAMMA))/(DD(K1)**TWO))
     +      +COEF*GAMMA*((PS(K1)**(GAMMA))
     +      *(CBARI2-THREE*PS(K1)**TWO)/(DD(K1)**THREE))
C
        END DO
      END DO
C---------------------------------
C
      ELSEIF (NEIG==2) THEN
C---------------------------------
C      TWO EQUAL EIGENVALUES
C---------------------------------
C       DERIVATIVES I1 AND I2 W.R.T A AND B
      DI1DA=BB**(-FOUR)+TWO*(BB**(-6.D0))*AA
      DI2DA=TWO*(BB**(-THREE))-ONE+FOUR*(BB**(-5.D0))*AA
      DI1DB=TWO-TWO*(BB**(-THREE))-FOUR*(BB**(-5.D0))*
     +          AA-6.D0*(BB**(-7.D0))*AA*AA
      DI2DB=TWO*BB-TWO*(BB**(-TWO))-6.D0*(BB**(-4.D0))*
     +          AA-10.D0*(BB**(-6.D0))*AA*AA
      D=DI1DA*DI2DB-DI1DB*DI2DA
C       DERIVATIVES A AND B W.R.T I1 AND I2
      DADI1=(ONE/D)*DI2DB
      DADI2=-(ONE/D)*DI1DB
      DBDI1=-(ONE/D)*DI2DA
      DBDI2=(ONE/D)*DI1DA
C      SECOND DERIVATIVES I1 AND I2 W.R.T A AND B
      D2DI1DA2=TWO*BB**(-6.D0)
      D2DI1DADB=-FOUR*BB**(-5.D0)-12.D0*BB**(-7.D0)*AA
      D2DI1DB2=6.D0*BB**(-FOUR)+20.D0*BB**(-6.D0)*AA
     +                      +42.D0*BB**(-8.D0)*AA*AA
      D2DI2DA2=FOUR*BB**(-5.D0)
      D2DI2DADB=-6.D0*BB**(-4.D0)-20.D0*BB**(-6.D0)*AA
      D2DI2DB2=TWO+4.D0*BB**(-THREE)+24.D0*BB**(-5.D0)*AA
     +                      +60.D0*BB**(-7.D0)*AA*AA
C      SECOND DERIVATIVES A AND B W.R.T I1 AND I2
      V1=DADI1*DADI2
      V2=DBDI1*DBDI2
      V3=DADI1*DBDI2+DBDI1*DADI2
      V4=D2DI1DA2*V1+D2DI1DADB*V3+D2DI1DB2*V2
      V5=D2DI2DA2*V1+D2DI2DADB*V3+D2DI2DB2*V2
      D2DADI1DI2=-DADI1*V4-DADI2*V5
      D2DBDI1DI2=-DBDI1*V4-DBDI2*V5
C
      V1=DADI1*DADI1
      V2=DBDI1*DBDI1
      V6=DADI1*DBDI1+DBDI1*DADI1
      V4=D2DI1DA2*V1+D2DI1DADB*V6+D2DI1DB2*V2
      V5=D2DI2DA2*V1+D2DI2DADB*V6+D2DI2DB2*V2
      D2DAD2I1=-DADI1*V4-DADI2*V5
      D2DBD2I1=-DBDI1*V4-DBDI2*V5
C
      V1=DADI2*DADI2
      V2=DBDI2*DBDI2
      V7=DADI2*DBDI2+DBDI2*DADI2
      V4=D2DI1DA2*V1+D2DI1DADB*V7+D2DI1DB2*V2
      V5=D2DI2DA2*V1+D2DI2DADB*V7+D2DI2DB2*V2
      D2DAD2I2=-DADI1*V4-DADI2*V5
      D2DBD2I2=-DBDI1*V4-DBDI2*V5
C
       DO I1=1,NEGN
C       AUXIALIAR VARS
        GAMMA=(ONE/TWO)*PROPS(2*I1)
        COEF=(TWO*PROPS(2*I1-1)/(FOUR*GAMMA**TWO))
        GAMMA1=GAMMA+ONE
        GAMMA2=GAMMA1+ONE
        GAMMA3=GAMMA2+ONE
        GAMMA4=GAMMA3+ONE
        GAMMA01=GAMMA-ONE
        GAMMA02=GAMMA01-ONE
        GAMMA03=GAMMA02-ONE
        GAMMA04=GAMMA03-ONE
        GAMMA05=GAMMA04-ONE
C
C       DERIVATIVES U W.R.T A AND B
        DUDA=GAMMA*BB**(-TWO*GAMMA1)+
     +      GAMMA*GAMMA01*BB**GAMMA02+
     +      (ONE/6.D0)*GAMMA*GAMMA01*GAMMA02*GAMMA03*
     +      BB**GAMMA04*AA+GAMMA*GAMMA1*
     +      BB**(-TWO*GAMMA2)*AA
        DUDB=TWO*GAMMA*BB**GAMMA01-TWO*GAMMA*BB**(-GAMMA1-GAMMA)-
     +       TWO*GAMMA*GAMMA1*BB**(-TWO*GAMMA-THREE)*AA+
     +       GAMMA*GAMMA01*GAMMA02*BB**GAMMA03*AA+
     +       (ONE/12.D0)*GAMMA*GAMMA01*GAMMA02*GAMMA03*
     +       GAMMA04*BB**GAMMA05*AA*AA-GAMMA*GAMMA1*
     +       GAMMA2*BB**(-TWO*GAMMA-5.D0)*AA*AA
C       SECOND DERIVATIVES U W.R.T A AND B
        D2DUDA=(ONE/6.D0)*GAMMA*GAMMA01*GAMMA02*GAMMA03*BB**GAMMA04+
     +          GAMMA*GAMMA1*BB**(-TWO*GAMMA2)
        D2DUDADB=-TWO*GAMMA*GAMMA1*BB**(-TWO*GAMMA-THREE)+
     +   GAMMA*GAMMA01*GAMMA02*BB**(GAMMA03)+
     +   (ONE/6.D0)*GAMMA*GAMMA01*GAMMA02*GAMMA03*GAMMA04*BB**GAMMA05*AA
     +   -TWO*GAMMA*GAMMA1*GAMMA2*BB**(-TWO*GAMMA2-ONE)*AA
        D2DUDB=TWO*GAMMA*GAMMA01*BB**GAMMA02+
     +   TWO*GAMMA*(TWO*GAMMA+ONE)*BB**(-TWO*GAMMA1)+
     +   TWO*GAMMA*GAMMA1*(TWO*GAMMA+THREE)*BB**(-TWO*GAMMA2)*AA+
     +   GAMMA*GAMMA01*GAMMA02*GAMMA03*BB**GAMMA04*AA+
     +   GAMMA*GAMMA1*GAMMA2*(TWO*GAMMA+5.D0)*BB**(-TWO*GAMMA3)*(AA*AA)+
     +   (ONE/12.D0)*GAMMA*GAMMA01*GAMMA02*GAMMA03*GAMMA04*GAMMA05*
     +    BB**(GAMMA-6.D0)*AA*AA
C       DERIVATIVES U W.R.T I1 AND I2
        DUDI1=DUDI1+COEF*(DUDA*DADI1+DUDB*DBDI1)
        DUDI2=DUDI2+COEF*(DUDA*DADI2+DUDB*DBDI2)
C       SECOND DERIVATIVES U W.R.T I1 AND I2
        D2DUDI1DI2=D2DUDI1DI2+COEF*
     +          (D2DUDA*DADI1*DADI2+D2DUDADB*V3+
     +          D2DUDB*DBDI1*DBDI2+DUDA*D2DADI1DI2+DUDB*D2DBDI1DI2)
        D2UD2I1=D2UD2I1+COEF*
     +          (D2DUDA*DADI1*DADI1+D2DUDADB*V6+
     +          D2DUDB*DBDI1*DBDI1+DUDA*D2DAD2I1+DUDB*D2DBD2I1)
        D2UD2I2=D2UD2I2+COEF*
     +          (D2DUDA*DADI2*DADI2+D2DUDADB*V7+
     +          D2DUDB*DBDI2*DBDI2+DUDA*D2DAD2I2+DUDB*D2DBD2I2)
C
      END DO
C---------------------------------
      ELSE
C---------------------------------
C      THREE EQUAL EIGENVALUES
C---------------------------------
      DO I1=1,NEGN
        GAMMA=(ONE/TWO)*PROPS(2*I1)
        COEF=(TWO*PROPS(2*I1-1)/(PROPS(2*I1)**TWO))
      C10=(ONE/TWO)*(ONE+GAMMA)
      C01=(ONE/TWO)*(ONE-GAMMA)
      C11=(-ONE/120.D0)*(GAMMA*GAMMA-ONE)*(GAMMA*GAMMA-FOUR)
      C20=(ONE/240.D0)*(GAMMA*GAMMA-ONE)*
     + (GAMMA*GAMMA+5.D0*GAMMA+6.D0)
      C02=(ONE/240.D0)*(GAMMA*GAMMA-ONE)*
     + (GAMMA*GAMMA-5.D0*GAMMA+6.D0)
C  FIRST DERIVATIVE UI WITH RESPECT TO I1
         DUDI1=DUDI1+COEF*
     + GAMMA**2*(C10+2*C20*(CI1-THREE)+C11*((CI2-THREE)-
     +             (ONE/8.D0)*(CI1+CI2-6.D0)**TWO))
C  FIRST DERIVATIVE UI WITH RESPECT TO I2
       DUDI2=DUDI2+COEF*
     + GAMMA**TWO*(C01+2*C02*(CI2-THREE)+C11*((CI1-THREE)-
     +             (ONE/8.D0)*(CI1+CI2-6.D0)**TWO))
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1
       D2UD2I1=D2UD2I1+COEF*
     + GAMMA**TWO*(TWO*C20-
     +             (ONE/FOUR)*C11*(CI1+CI2-6.D0))
C  SECOND DERIVATIVE UI  WITH RESPECT TO I2
       D2UD2I2=D2UD2I2+COEF*
     + GAMMA**TWO*(TWO*C02-
     +             (ONE/FOUR)*C11*(CI1+CI2-6.D0))
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1 AND I2
       D2DUDI1DI2=D2DUDI1DI2+COEF*
     + GAMMA**TWO*(C11*(ONE-
     +             (ONE/FOUR)*(CI1+CI2-6.D0)))
      END DO
C---------------------------------
      ENDIF
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES UI - FIBERS CONTRIBUTION
C
      D2DUDI1DI4=ZERO
      D2DUDI2DI4=ZERO
      D2DUDI1DLAMBDA=ZERO
      D2DUDI2DLAMBDA=ZERO
      DUDI4=ZERO
      D2UD2I4=ZERO
      DUDLAMBDA=ZERO
      D2UD2LAMBDA=ZERO
      D2DUDI4DLAMBDA=ZERO
C
      D2DUDI1DI6=ZERO
      D2DUDI2DI6=ZERO
      D2DUDI1DLAMBDA2=ZERO
      D2DUDI2DLAMBDA2=ZERO
      DUDI6=ZERO
      D2UD2I6=ZERO
      DUDLAMBDA2=ZERO
      D2UD2LAMBDA2=ZERO
      D2DUDI6DLAMBDA2=ZERO
C
C	FAMILY 1
C
      E1=CBARI4*(ONE-THREE*KAPPA)+CBARI1*KAPPA-ONE
      IF(E1.GT.ZERO) THEN
C
      EE2=DEXP(C4*E1*E1)
      EE3=(ONE+TWO*C4*E1*E1)
C
      DUDI4=C3*(ONE-THREE*KAPPA)*E1*EE2
      DUDI1=DUDI1+C3*KAPPA*E1*EE2
C
      D2UD2I4=C3*((ONE-THREE*KAPPA)**TWO)*EE3*EE2
      D2UD2I1=D2UD2I1+C3*KAPPA*KAPPA*EE3*EE2
      D2DUDI1DI4=C3*(ONE-THREE*KAPPA)*KAPPA*EE3*EE2
C
      ELSE
      DUDI4=ZERO
      D2UD2I4=ZERO
      DUDLAMBDA=ZERO
      D2UD2LAMBDA=ZERO
      D2DUDI4DLAMBDA=ZERO
C
      END IF
C
C	FAMILY 2
C
      EE1=CBARI6*(ONE-THREE*KAPPA)+CBARI1*KAPPA-ONE
      IF(EE1.GT.ZERO) THEN
C
      EE2=DEXP(C6*EE1*EE1)
      EE3=(ONE+TWO*C6*EE1*EE1)
C
      DUDI6=C5*(ONE-THREE*KAPPA)*EE1*EE2
      DUDI1=DUDI1+C5*KAPPA*EE1*EE2
C
      D2UD2I6=C5*((ONE-THREE*KAPPA)**TWO)*EE3*EE2
      D2UD2I1=D2UD2I1+C5*KAPPA*KAPPA*EE3*EE2
      D2DUDI1DI6=C5*(ONE-THREE*KAPPA)*KAPPA*EE3*EE2
      ELSE
C

      DUDI6=ZERO
      D2UD2I6=ZERO
      DUDLAMBDA2=ZERO
      D2UD2LAMBDA2=ZERO
      D2DUDI6DLAMBDA2=ZERO
C
      END IF
C
C	INTERACTION BETWEEN FAMILY 1 AND FAMILY 2 - ZERO
C
      DUDI8=ZERO
      D2UD2I8=ZERO
      D2UDJDI8=ZERO
      D2DUDI1DI8=ZERO
      D2DUDI2DI8=ZERO
      D2UDJDI8=ZERO
      D2UDI1DI8=ZERO
      D2UDI2DI8=ZERO
C
C----------------------------------------------------------------------
C     CAUCHY STRESS TENSOR
C----------------------------------------------------------------------
C
C     CAUCHY TENSOR SIGMA + AUXILIAR SS TENSORS
      DO I1 = 1, NDI
        DO J1=1,NDI
c	VOLUME PRESERVING CONTRIBUTION
               SS1(I1,J1)=UNIT(I1,J1)
C       ISOTROPIC CONTRIBUTION
               SS2(I1,J1)=TWO*SCALE1*(BBAR(I1,J1)-
     +                           (ONE/THREE)*CBARI1*UNIT(I1,J1))
               SS3(I1,J1)=TWO*SCALE1*(CBARI1*BBAR(I1,J1)-
     +                B2BAR(I1,J1)-(TWO/THREE)*CBARI2*UNIT(I1,J1))
               S1(I1,J1)=DUDJ*SS1(I1,J1)
               S2(I1,J1)=DUDI1*SS2(I1,J1)
               S3(I1,J1)=DUDI2*SS3(I1,J1)
C     FIBER CONTRIBUTION
C	FAMILY 1
               SS4(I1,J1)=TWO*SCALE1*(MM(I1,J1)
     +                       -(ONE/THREE)*CBARI4*UNIT(I1,J1)
     +                                        )
               SS5(I1,J1)=(ONE/LAMBDA)*SCALE1*M(I1,J1)
               S4(I1,J1)=DUDI4*SS4(I1,J1)
               S5(I1,J1)=DUDLAMBDA*SS5(I1,J1)
C	FAMILY 2
               SS6(I1,J1)=TWO*SCALE1*(NN(I1,J1)
     +                       -(ONE/THREE)*CBARI6*UNIT(I1,J1)
     +                                        )
               SS7(I1,J1)=(ONE/LAMBDA2)*SCALE1*N(I1,J1)
               S6(I1,J1)=DUDI6*SS6(I1,J1)
               S7(I1,J1)=DUDLAMBDA2*SS7(I1,J1)
C
C	FAMILY1 INTERACT WITH FAMILY 2
C
	       SS8(I1,J1)=TWO*SCALE1*(MM(I1,J1)+NN(I1,J1)
     +                       -(ONE/THREE)*CBARI6*UNIT(I1,J1)
     +					)
	       S8(I1,J1)=DUDI8*SS8(I1,J1)
C
C     CAUCHY STRESS TENSOR
               SIGMA(I1,J1) =
     +                           S1(I1,J1)
     +                          +S2(I1,J1)+S3(I1,J1)
     +       			+S4(I1,J1)+S5(I1,J1)
     +       			+S6(I1,J1)+S7(I1,J1)
     +       			+S8(I1,J1)
C
        END DO
      END DO
C
C----------------------------------------------------------------------
C     JACOBIAN MATERIAL MATRIX
C----------------------------------------------------------------------
      AUX0=(FOUR/THREE)*SCALE1*CBARI1
      AUX1=TWO*(FOUR/THREE)*SCALE1*CBARI2
      AUX2=(FOUR/THREE)*SCALE1*CBARI4
      AUX3=(FOUR/THREE)*SCALE1*CBARI6
C     SPATIAL TANGENT MODULI BASED ON THE JAUMANN
C                        RATE OF THE KIRCHHOFF STRESS TENSOR
      DO I=1,NDI
         DO J=1,NDI
            DO K=1,NDI
               DO L=1,NDI
               UNIT4=(ONE/TWO)*(UNIT(I,K)*UNIT(J,L)+UNIT(I,L)*UNIT(J,K))
C         ------JAUMMAN RATE PART------
                   DDSDDEJR(I,J,K,L)=
     +             (ONE/TWO)*(UNIT(I,K)*SIGMA(J,L)
     +             +SIGMA(I,K)*UNIT(J,L)+UNIT(I,L)*SIGMA(J,K)
     +             +SIGMA(I,L)*UNIT(J,K))
C         ------ISOTROPIC PART------
                   DDSDDEISO(I,J,K,L)=
     +             -(TWO/THREE)*(S2(I,J)*UNIT(K,L)+UNIT(I,J)*S2(K,L))+
     +             AUX0*DUDI1*(UNIT4-(ONE/THREE)*(UNIT(I,J)*UNIT(K,L)))-
     +             (FOUR/THREE)*(S3(I,J)*UNIT(K,L)+UNIT(I,J)*S3(K,L))+
     +             AUX1*DUDI2*(UNIT4-(TWO/THREE)*(UNIT(I,J)*UNIT(K,L)))-
     +             FOUR*SCALE1*DUDI2*
     +                  (BBAR(I,K)*BBAR(L,J)
     +                       -BBAR(I,J)*BBAR(K,L))+
     +             DET*(
     +             D2UD2I1*SS2(I,J)*SS2(K,L)+
     +             D2UD2I2*SS3(I,J)*SS3(K,L))
C         ------ANIISOTROPIC PART------
                   DDSDDEFIB1(I,J,K,L)=
C                  FAMILY 1 CONTRIBUTION
     +              -(TWO/THREE)*(S4(I,J)*UNIT(K,L)
     +                                          +UNIT(I,J)*S4(K,L))+
     +             AUX2*DUDI4*(UNIT4
     +                             -(ONE/THREE)*UNIT(I,J)*UNIT(K,L))+
     +             DET*(
     +            (-ONE/LAMBDA)*(S5(I,J)*SS5(K,L))+
     +             D2UD2I4*SS4(I,J)*SS4(K,L)+
     +             D2UD2LAMBDA*SS5(I,J)*SS5(K,L)+
     +             D2DUDI1DI4*(SS2(I,J)*SS4(K,L)+SS4(I,J)*SS2(K,L))+
     +             D2DUDI1DLAMBDA*(SS2(I,J)*SS5(K,L)+SS5(I,J)*SS2(K,L))+
     +             D2DUDI2DI4*(SS3(I,J)*SS4(K,L)+SS4(I,J)*SS3(K,L))+
     +             D2DUDI2DLAMBDA*(SS3(I,J)*SS5(K,L)+SS5(I,J)*SS3(K,L))+
     +             D2DUDI4DLAMBDA*(SS4(I,J)*SS5(K,L)+SS5(I,J)*SS4(K,L))
     +                )
                   DDSDDEFIB2(I,J,K,L)=
C                  FAMILY 2 CONTRIBUTION
     +              -(TWO/THREE)*(S6(I,J)*UNIT(K,L)
     +                                          +UNIT(I,J)*S6(K,L))+
     +             AUX3*DUDI6*(UNIT4
     +                             -(ONE/THREE)*UNIT(I,J)*UNIT(K,L))+
     +             DET*(
     +            (-ONE/LAMBDA2)*(S7(I,J)*SS7(K,L))+
     +            D2UD2I6*SS6(I,J)*SS6(K,L)+
     +            D2UD2LAMBDA2*SS7(I,J)*SS7(K,L)+
     +            D2DUDI1DI6*(SS2(I,J)*SS6(K,L)+SS6(I,J)*SS2(K,L))+
     +            D2DUDI1DLAMBDA2*(SS2(I,J)*SS7(K,L)+SS7(I,J)*SS2(K,L))+
     +            D2DUDI2DI6*(SS3(I,J)*SS6(K,L)+SS6(I,J)*SS3(K,L))+
     +            D2DUDI2DLAMBDA2*(SS3(I,J)*SS7(K,L)+SS7(I,J)*SS3(K,L))+
     +            D2DUDI6DLAMBDA2*(SS6(I,J)*SS7(K,L)+SS7(I,J)*SS6(K,L))
     +                )
                   DDSDDEFIB3(I,J,K,L)=
C                  INTERACTION BETWEEN FIBERS CONTRIBUTION
     +              -(TWO/THREE)*(S8(I,J)*UNIT(K,L)
     +                                          +UNIT(I,J)*S8(K,L))+
     +             AUX2*DUDI8*(UNIT4-(ONE/THREE)
     +                                          *UNIT(I,J)*UNIT(K,L))+
     +             DET*(
     +             D2UD2I8*SS8(I,J)*SS8(K,L)+
     +             D2DUDI1DI8*(SS2(I,J)*SS8(K,L)+SS8(I,J)*SS2(K,L))+
     +             D2DUDI2DI8*(SS3(I,J)*SS8(K,L)+SS8(I,J)*SS3(K,L))
     +                )
C        ------VOLUMETRIC PART------
                   DDSDDEVOL(I,J,K,L)=
C                  PURE VOLUMETRIC CONTRIBUTION
     +             UNIT(I,J)*S1(K,L)+S1(I,J)*UNIT(K,L)-
     +             DUDJ*(UNIT(I,J)*UNIT(K,L))-TWO*DUDJ*UNIT4+
     +             DET*(
     +             D2UDJ2*SS1(I,J)*SS1(K,L)+
C                  VOLUMETRIC-ISOTROPIC CONTRIBUTION
     +             D2DUDJDI1*(SS1(I,J)*SS2(K,L)+SS2(I,J)*SS1(K,L))+
     +             D2DUDJDI2*(SS1(I,J)*SS3(K,L)+SS3(I,J)*SS1(K,L))+
     +             D2DUDI1DI2*(SS2(I,J)*SS3(K,L)+SS3(I,J)*SS2(K,L))+
C                  VOLUMETRIC-ANISOTRPIC CONTRIBUTION
     +             D2DUDJDI4*(SS1(I,J)*SS4(K,L)+SS4(I,J)*SS1(K,L))+
     +             D2DUDJDLAMBDA*(SS1(I,J)*SS5(K,L)+SS5(I,J)*SS1(K,L))+
     +             D2DUDJDI6*(SS1(I,J)*SS6(K,L)+SS6(I,J)*SS1(K,L))+
     +             D2DUDJDLAMBDA2*(SS1(I,J)*SS7(K,L)+SS7(I,J)*SS1(K,L))+
     +             D2DUDJDI8*(SS1(I,J)*SS8(K,L)+SS8(I,J)*SS1(K,L))
     +                )
C        -----SPATIAL TANGENT MODULI------
                   DDSIGDDE(I,J,K,L)=
     +                                DDSDDEJR(I,J,K,L)+
     +                                DDSDDEISO(I,J,K,L)+
     +                                DDSDDEFIB1(I,J,K,L)+
     +                                DDSDDEFIB2(I,J,K,L)+
     +                                DDSDDEFIB3(I,J,K,L)+
     +                                DDSDDEVOL(I,J,K,L)
               END DO
             END DO
          END DO
       END DO
C----------------------------------------------------------------------
C     STRESS AND JACOBIAN MATRIX STORAGE
C----------------------------------------------------------------------
      DO I1=1,NTENS
C       STRESS VECTOR
         STRESS(I1)=SIGMA(II1(I1),II2(I1))
         DO J1=1,NTENS
C       DDSDDE - FULLY SIMMETRY IMPOSED
            PP1=DDSIGDDE(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=DDSIGDDE(II1(I1),II2(I1),II2(J1),II1(J1))
            DDSDDE(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     PRINTS: DEBUGGING PURPOSES
C      IF ((NOEL.EQ.1).AND.(NPT.EQ.1))THEN
C        WRITE(7,*)   SIGMA
C        WRITE(7,*) ""
C      ENDIF
C----------------------------------------------------------------------
C     STATE VARIABLES
C----------------------------------------------------------------------
C      DO K1 = 1, NTENS
        STATEV(1) = DET
        STATEV(2) = A
        STATEV(3) = E1
        STATEV(4) = EE1
C      END DO
C----------------------------------------------------------------------
C----------------------------------------------------------------------

C----------------------------------------------------------------------
      RETURN
      END
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C END OF MAIN UMAT ROUTINE
C****************************************************************************
C     UTILITY SUBROUTINES
C****************************************************************************
C
      SUBROUTINE MATINV3D(A,A_INV,DET_A,ISTAT)
C
C      RETURNS A_INV, THE INVERSE AND DET_A, THE DETERMINANT
C     DET OF THE ORIGINAL MATRIX, NOT OF THE INVERSE
C      RECURSIVE EXPRESSION OBTAINED FROM MATHEMATICA
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(3,3),A_INV(3,3),DET_A,DET_A_INV
C
C
      ISTAT = 1
C
      DET_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV3D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
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
C
      RETURN
      END SUBROUTINE MATINV3D
!****************************************************************************
C
      SUBROUTINE MATINV2D(A,A_INV,DET_A,ISTAT)
C
C     RETURNS A_INV, THE INVERSE, AND DET_A, THE DETERMINANT
C     NOTE THAT THE DET IS OF THE ORIGINAL MATRIX, NOT THE
C     INVERSE
C
      IMPLICIT NONE
C
      INTEGER ISTAT
C
      REAL*8 A(2,2),A_INV(2,2),DET_A,DET_A_INV
C
C
      ISTAT = 1.D0
C
      DET_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
      IF (DET_A .LE. 0.D0) THEN
        WRITE(*,*) 'WARNING: SUBROUTINE MATINV2D:'
        WRITE(*,*) 'WARNING: DET OF MAT=',DET_A
        ISTAT = 0
        RETURN
      END IF
C
      DET_A_INV = 1.D0/DET_A
C
      A_INV(1,1) =  DET_A_INV*A(2,2)
      A_INV(1,2) = -DET_A_INV*A(1,2)
      A_INV(2,1) = -DET_A_INV*A(2,1)
      A_INV(2,2) =  DET_A_INV*A(1,1)
C
C
      RETURN
      END SUBROUTINE MATINV2D
C
!****************************************************************************
C
      SUBROUTINE MDET(A,DET)
C
C      THIS SUBROUTINE CALCULATES THE DETERMINANT
C      OF A 3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      REAL*8  A(3,3),DET
C
C
      DET = A(1,1)*A(2,2)*A(3,3)
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)
C
C
      RETURN
      END SUBROUTINE MDET
C
!****************************************************************************
C
      SUBROUTINE ONEM(A)
C
C      THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE
C      3 BY 3 MATRIX [A]
C
      IMPLICIT NONE
C
      INTEGER I,J
C
      DOUBLE PRECISION A(3,3)
C
C
      DO I=1,3
         DO J=1,3
	           IF (I .EQ. J) THEN
              A(I,J) = 1.D0
            ELSE
              A(I,J) = 0.D0
            END IF
         END DO
      END DO
C
C
      RETURN
      END SUBROUTINE ONEM
C
!****************************************************************************
