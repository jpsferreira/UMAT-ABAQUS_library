      SUBROUTINE FIBDIR(FIB,ST0,ST,NE,NOEL,NDI,VORIF,VD,DISTGR,DFGRD1)
C
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      INTEGER NDI, NE, NOEL,INOEL,I,J,I1,J1
      DOUBLE PRECISION SUM1, DFGRD1(3,3), DNORM
      DOUBLE PRECISION VORIF(3),ST(3,3),VD(3),ST0(3,3),DISTGR(3,3)
      DOUBLE PRECISION FIB(NE,4)
C
        INOEL=0
        I=0
        DO I=1,NE
C               ELEMENT IDENTIFICATION
            IF(NOEL.EQ.INT(FIB(I,1))) THEN
                INOEL=I
            ENDIF
        ENDDO
C
C     FIB - FIBER ORIENTATION
             DNORM=DSQRT(FIB(INOEL,2)*FIB(INOEL,2)+
     1                   FIB(INOEL,3)*FIB(INOEL,3)+
     2                   FIB(INOEL,4)*FIB(INOEL,4))
C
C       UNDERFORMED FIBER ORIENTATION TENSOR
C
        DO I1=1,NDI
        J1=I1+1
C       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 1
        VORIF(I1)=FIB(INOEL,J1)/DNORM
        END DO
C

      DO I=1,NDI
         SUM1=ZERO
         DO J=1,NDI
          SUM1=SUM1+DFGRD1(I,J)*VORIF(J)
         ENDDO
C     FIBER DIRECTIONS IN THE DEFORMED CONFIGURATION
C               -FAMILY 1
         VD(I)=SUM1
      ENDDO
      DNORM=DSQRT(VD(1)*VD(1)+
     1             VD(2)*VD(2)+
     2             VD(3)*VD(3))
C           COSINE OF THE ANGLE BETWEEN FIBERS
C
C
C--------------------------------------------------------------------------
      DO I=1,NDI
       DO J=1,NDI
C       STRUCTURAL TENSOR - FAMILY 1
       ST0(I,J)=VORIF(I)*VORIF(J)
       END DO
      END DO
C
C       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 1
      ST=MATMUL(ST0,TRANSPOSE(DISTGR))
      ST=MATMUL(DISTGR,ST)
C
C
      RETURN
      END SUBROUTINE FIBDIR
