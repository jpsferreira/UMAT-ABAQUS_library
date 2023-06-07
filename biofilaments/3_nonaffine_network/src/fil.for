      SUBROUTINE FIL(F,FF,DW,DDW,LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0)
C>    SINGLE FILAMENT: STRAIN ENERGY DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION DW,DDW
      DOUBLE PRECISION LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0
      DOUBLE PRECISION A,B,MACHEP,T
      DOUBLE PRECISION AUX,F,FF,PI,ALPHA
      DOUBLE PRECISION AUX0,AUX1,AUX2,AUX3,AUX4,AUX5,AUX6,Y
C
      A=ZERO
      B=1.0E09
      MACHEP=2.2204e-16
      T=1.0E-6
      F=ZERO
C
      CALL PULLFORCE(F, A, B, MACHEP, T,
     1                   LAMBDA,LAMBDA0,LL,R0,MU0,BETA,B0)
C
      PI=FOUR*ATAN(ONE)
      FF=F*LL*(PI*PI*B0)**(-ONE)
      ALPHA=PI*PI*B0*(LL*LL*MU0)**(-ONE)
C
      AUX0=BETA/ALPHA
      AUX=ALPHA*FF
      AUX1=ONE+FF+AUX*FF
      AUX2=ONE+TWO*AUX
      AUX3=ONE+AUX
      AUX4=LAMBDA0*R0*R0*MU0*(LL**(-ONE))
      AUX5=((ONE+AUX)*(AUX1**(-ONE)))**BETA
      AUX6=ONE-R0*((LL)**(-ONE))
C      
      Y=AUX0*(AUX2*AUX2*(AUX1**(-ONE)))-BETA*(AUX2*(AUX3**(-ONE)))-TWO
C      
      DW=LAMBDA0*R0*F
      DDW=AUX4*((ONE+Y*AUX5*AUX6)**(-ONE))
C
      RETURN
      END SUBROUTINE FIL
