      SUBROUTINE BANGLE(ANG,F,MF,NOEL,NDI)
C>    ANGLE BETWEEN FILAMENT AND PREFERED DIRECTION
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
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

