      SUBROUTINE ISOMAT(SSEISO,DISO,C10,C01,CBARI1,CBARI2)
C>     ISOTROPIC MATRIX : ISOCHORIC SEF AND DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION SSEISO,DISO(5)
      DOUBLE PRECISION C10,C01,CBARI1,CBARI2
C
      SSEISO=C10*(DEXP(C01*(CBARI1-THREE))-ONE)
C
      !FIRST DERIVATIVE OF SSEISO IN ORDER TO I1
      DISO(1)=C10*C01*DEXP(C01*(CBARI1-THREE))
      !FIRST DERIVATIVE OF SSEISO IN ORDER TO I2
      DISO(2)=ZERO
      !SECOND DERIVATIVE OF SSEISO IN ORDER TO I1
      DISO(3)=C01*C01*C10*DEXP(C01*(CBARI1-THREE))
      !SECOND DERIVATIVE OF SSEISO IN ORDER TO I2
      DISO(4)=ZERO
      !SECOND DERIVATIVE OF SSEISO IN ORDER TO I1 AND I2
      DISO(5)=ZERO
C
      RETURN
      END SUBROUTINE ISOMAT
