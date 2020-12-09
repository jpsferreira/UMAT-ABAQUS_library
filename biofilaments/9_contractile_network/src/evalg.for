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