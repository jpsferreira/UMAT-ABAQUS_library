       SUBROUTINE CHEMICALSTAT(FRAC,FRAC0,K,DTIME)
C
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION FRAC(4),FRAC0(4),STIFF(4,4),K(7)
C
      DOUBLE PRECISION DTIME,AUX
      INTEGER I1,J1
C
        FRAC=ZERO
        STIFF=ZERO
        STIFF(1,1)=-K(1)
        STIFF(1,2)=K(2)
        STIFF(1,4)=K(7)
        STIFF(2,1)=K(1)
        STIFF(2,2)=-K(2)-K(3)
        STIFF(2,3)=K(4)
        STIFF(3,2)=K(3)
        STIFF(3,3)=-K(4)-K(5)
        STIFF(3,4)=K(6)
        STIFF(4,3)=K(5)
        STIFF(4,4)=-K(6)-K(7)
C        
      DO I1=1,4
        AUX=ZERO
       DO J1=1,4
          AUX=AUX+STIFF(I1,J1)*FRAC0(J1)
       ENDDO
       FRAC(I1)=AUX*DTIME+FRAC0(I1)
      ENDDO        
C      
C      FRAC0=FRAC
C
      RETURN
      END SUBROUTINE CHEMICALSTAT
