C********************************************************************
C Record of revisions:                                              |
C        Date        Programmer        Description of change        |
C        ====        ==========        =====================        |
C     15/02/2015    Joao Ferreira      adapted from Jabaren eqs     |                            
C                                                                   |                                          
C--------------------------------------------------------------------
C Description:
C This script contains a subroutine for: 
C OGDEN isotropic materials
C   
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'aba_param.inc'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C---------------------------------------------------------------------
C     LOCAL ARRAYS (CHANGE THIS LIST LATER!!)
C---------------------------------------------------------------------
C     EELAS  - LOGARITHMIC ELASTIC STRAINS
C     EELASP - PRINCIPAL ELASTIC STRAINS
C     BBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     BBARP  - PRINCIPAL VALUES OF BBAR
C     BBARN  - PRINCIPAL DIRECTION OF BBAR (AND EELAS)
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C     C      -  RIGHT CAUCHY-GREEN TENSOR
C     CBAR   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     CBAR2  - SQUARE OF THE DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     DI1DE  - 1ST INVARIANT DERIVATIVE IN ORDER TO GREEN-LAGRANGE TENSOR
C     DI2DE  - 2ND INVARIANT DERIVATIVE IN ORDER TO GREEN-LAGRANGE TENSOR
C     EGL    - GREEN-LAGRANGE STRAIN TENSOR
C     CBARI1 - CBAR 1ST INVARIANT
C     CBARI2 - CBAR 1ND INVARIANT
C     DWDI1  - STRAIN-ENERGY DERIVATIVE IN ORDER TO CBARI1
C     DWDI2  - STRAIN-ENERGY DERIVATIVE IN ORDER TO CBARI2
C     S      - 2ND PIOLA-KIRCHHOFF STRESS TENSOR
C     SIGMA  - CAUCHY STRESS TENSOR
C----------------------------------------------------------------------
C
      PARAMETER(ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0,
     1          FOUR = 4.D0)
      INTEGER   II1(6),II2(6),NN

      REAL*8    UNIT(3,3), C(3,3),B(3,3),CBAR(3,3),BBAR(3,3),
     2          C2(3,3),B2BAR(3,3),PS(3),
     3          SS1(3,3), SS2(3,3), SS3(3,3),  
     4          S1(3,3), S2(3,3), S3(3,3), SIGMA(3,3),
     5          DDSIGDDE(3,3,3,3),DD(3),AN(3,3),VCBAR(6)
C----------------------------------------------------------------------
C     UMAT FOR COMPRESSIBLE OGDEN HYPERELASTICITY 
C     CANNOT BE USED FOR PLANE STRESS
C----------------------------------------------------------------------
C     PROPS(1)  - C10
C     PROPS(2)  - C01
C----------------------------------------------------------------------
C
C     ELASTIC PROPERTIES
C
C      EMOD = PROPS(1)
C      ENU  = PROPS(2)
C      C10  = EMOD / (FOUR * (ONE + ENU))
C      D1   = 6.0 * (ONE - TWO * ENU) / EMOD
      D1   = PROPS(NPROPS)
C----------------------------------------------------------------------
C     INITIALIZATIONS
C----------------------------------------------------------------------
C     AUXILIAR VARIABLES
      ISTAT=ONE
      DET=ZERO
      SCALE=ZERO
      SCALE0=ZERO
      SCALE1=ZERO
      SCALE2=ZERO
      SCALE3=ZERO
      SUM=ZERO
      SUMB=ZERO
      CI1=ZERO
      C2I1=ZERO
      CBARI1=ZERO
      C2BARI1=ZERO
      CBARI2=ZERO
      DUDI1=ZERO
      D2UD2I1=ZERO
      DUDI2=ZERO
      D2UD2I2=ZERO
      D2DUDI1DI2=ZERO
      D2DUDIJDI1=ZERO
      D2DUDIJDI2=ZERO
      DUDJ=ZERO
      D2UDJ2=ZERO
      AUX0=ZERO
      AUX1=ZERO
      UNIT4=ZERO
C----------------------------------------------------------------------
c      DEBUGGING!!!!!!!!
c      WRITE(6,*) "F(1,2)",DFGRD1(1,2)
c      WRITE(6,*) "gp", NPT
C----------------------------------------------------------------------
C----------------------------------------------------------------------
c * IDENTITY TENSOR DEFINITION                                       *
C----------------------------------------------------------------------
C     
         DO I=1,3
         DO J=1,3
	    IF (I .EQ. J) THEN
              UNIT(I,J) = ONE
            ELSE
              UNIT(I,J) = ZERO
            END IF
         END DO
      END DO
C----------------------------------------------------------------------
c * POINTERS DEFINITION                                              *
C----------------------------------------------------------------------
C     POINTERS TO STORAGE INTO "STRESS" AND "DDSDDE" ARRAYS
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
C----------------------------------------------------------------------
C     CALCULATE RIGHT AND LEFT CAUCHY-GREEN TENSOR
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
      C2=MATMUL(C,C)
      B2BAR=MATMUL(BBAR,BBAR)
      DO I1=1,NDI
        C2I1=C2I1+C2(I1,I1)
      END DO
C    SECOND INVARIANT OF CBAR
      CI2=(ONE/TWO)*(CI1*CI1-C2I1)
      CBARI2=SCALE3*CI2
C----------------------------------------------------------------------
C     STRAIN-ENERGY DERIVATIVES WITH RESPECT TO INVARIANTS
C----------------------------------------------------------------------

C     EIGENVALUES AND EIGENVECTORS OF CBAR
C      CALL SPRIND(VCBAR, PS, AN, 2, NDI, NSHR)
      CALL SPECTRAL(CBAR,PS,AN,1)
C     
C        
            CHECK12=ABS(PS(1)-PS(2))
            CHECK13=ABS(PS(1)-PS(3))
            CHECK23=ABS(PS(2)-PS(3))
C
C           NEIG IS THE NUMBER OF EQUAL EIGENVALUES
          TOL=1.D-5
          MINCHECK=CHECK12

            IF((CHECK12.LT.TOL).OR.
     &         (CHECK13.LT.TOL).OR.
     &         (CHECK23.LT.TOL)) THEN
C          CHECK FOR THREE EQUAL EIGENVALUES
            IF ((CHECK12.LT.TOL).AND.
     &              (CHECK13.LT.TOL).AND.
     &                  (CHECK23.LT.TOL)) THEN
                NEIG=3
             ELSE
C          CHECK FOR TWO EQUAL EIGENVALUES
                NEIG=2
                 AA=(MINCHECK/TWO)**TWO
                 BB=(PS(1)+PS(2))/TWO
                IF (CHECK13 .LT. MINCHECK) THEN
                MINCHECK=CHECK13
                AA=(MINCHECK/TWO)**TWO
                BB=(PS(1)+PS(3))/TWO
                ELSEIF (CHECK23 .LT. MINCHECK) THEN
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
      NN=(NPROPS-1)/2
C
       IF (NEIG==1) THEN
C---------------------------------
C      THREE DIFFERENT EIGENVALUES
C---------------------------------
      DD(1)=(PS(1)-PS(2))*(PS(1)-PS(3))
      DD(2)=(PS(2)-PS(1))*(PS(2)-PS(3))
      DD(3)=(PS(3)-PS(1))*(PS(3)-PS(2))
      DO I1=1,NN
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
C      TWO DIFFERENT EIGENVALUES
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
       DO I1=1,NN
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
      DO I1=1,NN
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
C
C
C  FIRST DERIVATIVE UJ  WITH RESPECT TO J 
       DUDJ=(TWO/D1)*(DET-ONE)
C  SECOND DERIVATIVE UJ  WITH RESPECT TO J 
       D2UDJ2=(TWO/D1)
       D2DUDJDI1=ZERO
       D2DUDJDI2=ZERO
C  

C----------------------------------------------------------------------
C     CALCULATE THE STRESS TENSOR
C----------------------------------------------------------------------
C     
C     PUSH FORWARD OPERATION
C     CAUCHY TENSOR SIGMA + AUXILIAR SS TENSORS
      DO I1 = 1, NDI
        DO J1=1,NDI
               SS1(I1,J1)=UNIT(I1,J1)
               SS2(I1,J1)=TWO*SCALE1*(BBAR(I1,J1)-
     +                           (ONE/THREE)*CBARI1*UNIT(I1,J1))
               SS3(I1,J1)=TWO*SCALE1*(CBARI1*BBAR(I1,J1)-
     +                B2BAR(I1,J1)
     +                  -(TWO/THREE)*CBARI2*UNIT(I1,J1))   
               S1(I1,J1)=DUDJ*SS1(I1,J1)
               S2(I1,J1)=DUDI1*SS2(I1,J1)
               S3(I1,J1)=DUDI2*SS3(I1,J1)
               SIGMA(I1,J1) =S1(I1,J1)+S2(I1,J1)+S3(I1,J1)    
        END DO
      END DO
C
C      CALL ROTSIG(SIGMA,R,SIGMA,1,NDI,NSHR)
C        SIGMA=MATMUL(MATMUL(TRANSPOSE(AN),SIGMA),AN)
C----------------------------------------------------------------------
C     CALCULATE THE JACOBIAN MATERIAL MATRIX
C----------------------------------------------------------------------
      AUX0=(FOUR/THREE)*SCALE1*CBARI1
      AUX1=TWO*(FOUR/THREE)*SCALE1*CBARI2
C     SPATIAL TANGENT MODULI BASED ON THE JAUMANN 
C                        RATE OF THE KIRCHHOFF STRESS TENSOR
      DO I=1,NDI
         DO J=1,NDI
            DO K=1,NDI
               DO L=1,NDI
               UNIT4=(ONE/TWO)*(UNIT(I,K)*UNIT(J,L)+UNIT(I,L)*UNIT(J,K))
C   
                  DDSIGDDE(I,J,K,L)=
     +             (ONE/TWO)*(UNIT(I,K)*SIGMA(J,L)             
     +             +SIGMA(I,K)*UNIT(J,L)+UNIT(I,L)*SIGMA(J,K)
     +             +SIGMA(I,L)*UNIT(J,K))+
     +             UNIT(I,J)*S1(K,L)+S1(I,J)*UNIT(K,L)-
     +             DUDJ*(UNIT(I,J)*UNIT(K,L))-TWO*DUDJ*UNIT4-
     +             (TWO/THREE)*(S2(I,J)*UNIT(K,L)+UNIT(I,J)*S2(K,L))+
     +             AUX0*DUDI1*(UNIT4-(ONE/THREE)*(UNIT(I,J)*UNIT(K,L)))-
     +             (FOUR/THREE)*(S3(I,J)*UNIT(K,L)+UNIT(I,J)*S3(K,L))+
     +             AUX1*DUDI2*(UNIT4-(TWO/THREE)*(UNIT(I,J)*UNIT(K,L)))-
     +             FOUR*SCALE1*DUDI2*
     +                  (BBAR(I,K)*BBAR(L,J)
     +                       -BBAR(I,J)*BBAR(K,L))+
     +             DET*(
     +             D2UDJ2*SS1(I,J)*SS1(K,L)+
     +             D2UD2I1*SS2(I,J)*SS2(K,L)+
     +             D2UD2I2*SS3(I,J)*SS3(K,L)+
     +             D2DUDJDI1*(SS1(I,J)*SS2(K,L)+SS2(I,J)*SS1(K,L))+
     +             D2DUDJDI2*(SS1(I,J)*SS3(K,L)+SS3(I,J)*SS1(K,L))+
     +             D2DUDI1DI2*(SS2(I,J)*SS3(K,L)
     +                                  +SS3(I,J)*SS2(K,L)))           
               END DO
             END DO
          END DO
       END DO        
C----------------------------------------------------------------------
C     STRESS AND JACOBIAN MATRIX STORAGE  
C----------------------------------------------------------------------
      DO I1=1,6
C       STRESS VECTOR
         STRESS(I1)=SIGMA(II1(I1),II2(I1))
         DO J1=1,6
C       DDSDDE 
            PP1=DDSIGDDE(II1(I1),II2(I1),II1(J1),II2(J1))
            PP2=DDSIGDDE(II1(I1),II2(I1),II2(J1),II1(J1))
            DDSDDE(I1,J1)=(ONE/TWO)*(PP1+PP2)
         END DO
      END DO
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
C****************************************************************************
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
C****************************************************************************
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
C****************************************************************************
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
      REAL*8 A(3,3)
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
C****************************************************************************
      SUBROUTINE SPECTRAL(A,D,V,ISTAT)
C     
C      THIS SUBROUTINE CALCULATES THE EIGENVALUES AND EIGENVECTORS OF
C      A SYMMETRIC 3 BY 3 MATRIX A.
C     
C      THE OUTPUT CONSISTS OF A VECTOR D CONTAINING THE THREE
C      EIGENVALUES IN ASCENDING ORDER, AND A MATRIX V WHOSE
C      COLUMNS CONTAIN THE CORRESPONDING EIGENVECTORS.
C
      IMPLICIT NONE
C
      INTEGER NP,NROT,I,J,ISTAT
      PARAMETER(NP=3)
C
      REAL*8 D(3),V(3,3),A(3,3),E(3,3)
C
      E = A
C
      CALL JACOBI(E,3,NP,D,V,NROT,ISTAT)
      CALL EIGSRT(D,V,3,NP)
C
      RETURN
      END SUBROUTINE SPECTRAL
C	
C****************************************************************************
C
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT,ISTAT)
C
C     COMPUTES ALL EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX A, WHICH IS OF SIZE N BY N, STORED IN A PHYSICAL
C      NP BY NP ARRAY.  ON OUTPUT, ELEMENTS OF A ABOVE THE DIAGONAL
C     ARE DESTROYED, BUT THE DIAGONAL AND SUB-DIAGONAL ARE UNCHANGED
C     AND GIVE FULL INFORMATION ABOUT THE ORIGINAL SYMMETRIC MATRIX.
C      VECTOR D RETURNS THE EIGENVALUES OF A IN ITS FIRST N ELEMENTS.
C      V IS A MATRIX WITH THE SAME LOGICAL AND PHYSICAL DIMENSIONS AS
C     A WHOSE COLUMNS CONTAIN, UPON OUTPUT, THE NORMALIZED
C      EIGENVECTORS OF A.  NROT RETURNS THE NUMBER OF JACOBI ROTATION
C     WHICH WERE REQUIRED.
C
C      THIS SUBROUTINE IS AVAILABLE AT 'NUMERICAL RECIPES IN FORTRAN.'
C
      IMPLICIT NONE
      !
      INTEGER IP,IQ,N,NMAX,NP,NROT,I,J,ISTAT
      PARAMETER (NMAX=100)
      !
      REAL*8 A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX),
     +  SM,TRESH,G,T,H,THETA,S,C,TAU
C
C      INITIALIZE V TO THE IDENTITY MATRIX
C
      CALL ONEM(V)
C
      
C     INITIALIZE B AND D TO THE DIAGONAL OF A, AND Z TO ZERO.
C     THE VECTOR Z WILL ACCUMULATE TERMS OF THE FORM T*A_PQ AS
C     IN EQUATION (11.1.14)
C     
      DO IP = 1,N
	B(IP) = A(IP,IP)
	D(IP) = B(IP)
	Z(IP) = 0.D0
      END DO
C
C     BEGIN ITERATION
C
      NROT = 0
      DO I=1,50
C
C     SUM OFF-DIAGONAL ELEMENTS
C
          SM = 0.D0
          DO IP=1,N-1
            DO IQ=IP+1,N
	      SM = SM + DABS(A(IP,IQ))
            END DO
          END DO
C
C        IF SM = 0., THEN RETURN.  THIS IS THE NORMAL RETURN,
C       WHICH RELIES ON QUADRATIC CONVERGENCE TO MACHINE
C      UNDERFLOW.
C
          IF (SM.EQ.0.D0) RETURN
C
C     IN THE FIRST THREE SWEEPS CARRY OUT THE PQ ROTATION ONLY IF
C     |A_PQ| > TRESH, WHERE TRESH IS SOME THRESHOLD VALUE,
C      SEE EQUATION (11.1.25).  THEREAFTER TRESH = 0.
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
C     AFTER FOUR SWEEPS, SKIP THE ROTATION IF THE 
C     OFF-DIAGONAL ELEMENT IS SMALL.
C
	      IF ((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP)))
     +            .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ)))) THEN
                A(IP,IQ) = 0.D0
              ELSE IF (DABS(A(IP,IQ)).GT.TRESH) THEN
                H = D(IQ) - D(IP)
                IF (DABS(H)+G.EQ.DABS(H)) THEN
C
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
C       CASE OF ROTATIONS 1 <= J < P
C	
	        DO J=1,IP-1
	          G = A(J,IP)
	          H = A(J,IQ)
	          A(J,IP) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
	        END DO
C
C       CASE OF ROTATIONS P < J < Q
C
	        DO J=IP+1,IQ-1
	          G = A(IP,J)
	          H = A(J,IQ)
	          A(IP,J) = G - S*(H + G*TAU)
	          A(J,IQ) = H + S*(G - H*TAU)
	        END DO
C
C      CASE OF ROTATIONS Q < J <= N
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
C      UPDATE D WITH THE SUM OF T*A_PQ, AND REINITIALIZE Z
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
      ISTAT = 0
C
      RETURN
      END SUBROUTINE JACOBI
C
C****************************************************************************
C
      SUBROUTINE EIGSRT(D,V,N,NP)
C
C     GIVEN THE EIGENVALUES D AND EIGENVECTORS V AS OUTPUT FROM
C     JACOBI, THIS SUBROUTINE SORTS THE EIGENVALUES INTO ASCENDING
C      ORDER AND REARRANGES THE COLMNS OF V ACCORDINGLY.
C
C THE SUBROUTINE WAS TAKEN FROM 'NUMERICAL RECIPES.'
C
      IMPLICIT NONE
C
      INTEGER N,NP,I,J,K
C
      REAL*8 D(NP),V(NP,NP),P
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