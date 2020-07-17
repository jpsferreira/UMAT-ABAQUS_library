C********************************************************************
C Record of revisions:                                              |
C        Date        Programmer        Description of change        |
C        ====        ==========        =====================        |
C     15/02/2015    Joao Ferreira      adapted from Jabaren eqs     |                            
C                                                                   |                                          
C--------------------------------------------------------------------
C Description:
C This script contains a subroutine for: 
C NEOHOOKE isotropic materials
C   
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
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
      INTEGER   II1(6),II2(6)

      REAL*8    UNIT(3,3), C(3,3),B(3,3),CBAR(3,3),BBAR(3,3),
     2          C2BAR(3,3),B2BAR(3,3),
     3          SS1(3,3), SS2(3,3), SS3(3,3),  
     4          S1(3,3), S2(3,3), S3(3,3), SIGMA(3,3),
     5          DDSIGDDE(3,3,3,3)
C SEE LATER WHICH VARIABLES ARE INITIATED AT "REAL*8" AND "DIMENSION"!!
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
      C10  = PROPS(1)
      D1   = PROPS(2)
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
      CALL ONEM(UNIT)
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
      C2BAR=MATMUL(CBAR,CBAR)
      B2BAR=MATMUL(BBAR,BBAR)
      DO I1=1,NDI
        C2BARI1=C2BARI1+C2BAR(I1,I1)
      END DO
      CBARI2=CBARI1*CBARI1-C2BARI1
C----------------------------------------------------------------------
C     STRAIN-ENERGY DERIVATIVES WITH RESPECT TO INVARIANTS
C----------------------------------------------------------------------
C     FIRST AND SECOND DERIVATIVES UI 
C
C  FIRST DERIVATIVE UI WITH RESPECT TO I1    
       DUDI1=C10
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1    
       D2UD2I1=ZERO
C  FIRST DERIVATIVE UI WITH RESPECT TO I2    
       DUDI2=ZERO
C  SECOND DERIVATIVE UI  WITH RESPECT TO I2    
       D2UD2I2=ZERO
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1 AND I2
       D2DUDI1DI2=ZERO
C  SECOND DERIVATIVE UI  WITH RESPECT TO I1 AND J
       D2DUDJDI1=ZERO
C  SECOND DERIVATIVE UI  WITH RESPECT TO I2 AND J
       D2DUDJDI2=ZERO
C  FIRST DERIVATIVE UJ  WITH RESPECT TO J 
       DUDJ=(TWO/D1)*(DET-ONE)
C  SECOND DERIVATIVE UJ  WITH RESPECT TO J 
       D2UDJ2=(TWO/D1)
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
     +                B2BAR(I1,J1)-(TWO/THREE)*CBARI2*UNIT(I1,J1))     
               S1(I1,J1)=DUDJ*SS1(I1,J1)
               S2(I1,J1)=DUDI1*SS2(I1,J1)
               S3(I1,J1)=DUDI2*SS3(I1,J1)
               SIGMA(I1,J1) =S1(I1,J1)+S2(I1,J1)+S3(I1,J1)        
        END DO
      END DO
C
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
C     +            UNIT4*SIGMA(K,L)+SIGMA(I,J)*UNIT4+
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
     +                  (BBAR(I,J)*UNIT4*BBAR(K,L)-BBAR(I,J)*BBAR(K,L))+
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
!****************************************************************************