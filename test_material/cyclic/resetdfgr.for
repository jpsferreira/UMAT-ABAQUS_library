      SUBROUTINE RESETDFGRD(DFGRD,NDI)
C
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
      
      INTEGER NDI
      DOUBLE PRECISION DFGRD(NDI,NDI)

        DFGRD(1,1)=  ONE
        DFGRD(1,2)=  ZERO
        DFGRD(1,3)=  ZERO
        DFGRD(2,1)=  ZERO
        DFGRD(2,2)=  ONE
        DFGRD(2,3)=  ZERO
        DFGRD(3,1)=  ZERO
        DFGRD(3,2)=  ZERO
        DFGRD(3,3)=  ONE

      END