      SUBROUTINE CMATISOMATFIC(CMISOMATFIC,CBAR,CBARI1,CBARI2,
     1                          DISO,UNIT2,UNIT4,DET,NDI)
C>    ISOTROPIC MATRIX: MATERIAL 'FICTICIOUS' ELASTICITY TENSOR
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
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
