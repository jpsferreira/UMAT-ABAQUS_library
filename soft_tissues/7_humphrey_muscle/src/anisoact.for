      SUBROUTINE ANISOACT(SSEANISO,DANISO,ACT,CBARI4,PREDEF,DPRED,T0M)
C>     ANISOTROPIC PART : ISOCHORIC SEF AND DERIVATIVES WITH ACTIVATION
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C
      DOUBLE PRECISION SSEANISO,DANISO(4),DISO(5)
      DOUBLE PRECISION PREDEF(1),DPRED(1)
      DOUBLE PRECISION T0M,CBARI4,E1,ACT,BARLAMBDA
      DOUBLE PRECISION DUACT,D2UACT,DUDI4,D2UD2I4
C      
C
      DUDI4 = DANISO(1)
      D2UD2I4 = DANISO(2)
      BARLAMBDA = DSQRT(CBARI4)
      E1 = BARLAMBDA - ONE
      ACT=PREDEF(1)+DPRED(1)
C      
      SSEANISO=SSEANISO

      IF((BARLAMBDA.GT.0.5d0).AND.(BARLAMBDA.LT.1.5d0)) THEN
C
         DUACT=ACT*T0M*(-FOUR*E1**TWO+ONE)
         D2UACT=ACT*T0M*(-8.0d0*E1)
         SSEANISO=SSEANISO+ACT*T0M*(-(FOUR/THREE)*E1**THREE+BARLAMBDA)
      ELSE
         DUACT=ZERO
         D2UACT=ZERO
      END IF
C
c      DUDI4=DUDI4+DUACT
       DUDI4=DUDI4+DUACT*(ONE/TWO)*CBARI4**(-ONE/TWO)
C      D2UD2I4=D2UD2I4+D2UACT
       D2UD2I4=D2UD2I4+D2UACT*(-ONE/FOUR)*(CBARI4**(-THREE/TWO))+D2UACT*((ONE/TWO)*(CBARI4**(-ONE/TWO))**TWO)
C
      !FIRST DERIVATIVE OF SSEANISO IN ORDER TO I1
      DANISO(1)=DUDI4
      !FIRST DERIVATIVE OF SSEANISO IN ORDER TO I2
      DANISO(2)=D2UD2I4
C
      RETURN
      END SUBROUTINE ANISOACT
