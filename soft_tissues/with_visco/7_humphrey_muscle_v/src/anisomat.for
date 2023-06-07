      SUBROUTINE ANISOMAT(SSEANISO,DANISO,DISO,K1,K2,KDISP,I4,I1)
C>     ANISOTROPIC PART : ISOCHORIC SEF AND DERIVATIVES
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION DANISO(4),DISO(5)
      DOUBLE PRECISION K1,K2,KDISP,I4,I1,LAMBDA
      DOUBLE PRECISION DUDI1,D2UD2I1,SSEANISO
      DOUBLE PRECISION E1,EE2,EE3,DUDI4,D2UD2I4,D2DUDI1DI4,D2DUDI2DI4
      DOUBLE PRECISION E2,E3,E4,E5
C
      DUDI1=DISO(1)
      D2UD2I1=DISO(3)
      LAMBDA=DSQRT(I4)
C      
C
      SSEANISO=K1*(DEXP(K2*(LAMBDA-ONE)*(LAMBDA-ONE))-ONE)
C
C      E1=I4*(ONE-THREE*KDISP)+I1*KDISP-ONE
      E1=LAMBDA-ONE
      IF(E1.GT.ZERO) THEN
C
      E3=DEXP(K2*E1*E1)
      E2=TWO*K1*K2
      E4=ONE+TWO*K2*E1*E1*LAMBDA
      E5=FOUR*(I4**(THREE/TWO))
C
      DUDI4=E2*E1*E3/(TWO*LAMBDA)
      D2UD2I4=E2*E3*E4/E5
C
      ELSE
      DUDI4=ZERO
      D2UD2I4=ZERO
C
      ENDIF
C     this commented portion of the code is to change later if one wants to implement the model with dispersion
!       IF(E1.GT.ZERO) THEN
! C
!       EE2=DEXP(K2*E1*E1)
!       EE3=(ONE+TWO*K2*E1*E1)
! C
!       DUDI1=DUDI1+K1*KDISP*E1*EE2
!       D2UD2I1=D2UD2I1+K1*KDISP*KDISP*EE3*EE2
! C      
!       DUDI4=K1*(ONE-THREE*KDISP)*E1*EE2
! C
!       D2UD2I4=K1*((ONE-THREE*KDISP)**TWO)*EE3*EE2
      
!       D2DUDI1DI4=K1*(ONE-THREE*KDISP)*KDISP*EE3*EE2
!       D2DUDI2DI4=ZERO
! C
!       ELSE
!       DUDI4=ZERO
!       D2UD2I4=ZERO
!       D2DUDI1DI4=ZERO
!       D2DUDI2DI4=ZERO
!       D2UD2I1=ZERO
! C
!       END IF

      ! DISO(1)=DUDI1
      ! DISO(3)=D2UD2I1

C
      !FIRST DERIVATIVE OF SSEANISO IN ORDER TO I1
      DANISO(1)=DUDI4
      !FIRST DERIVATIVE OF SSEANISO IN ORDER TO I2
      DANISO(2)=D2UD2I4
      !2ND DERIVATIVE OF SSEANISO IN ORDER TO I1
      DANISO(3)=D2DUDI1DI4
      !2ND DERIVATIVE OF SSEANISO IN ORDER TO I2
      DANISO(4)=D2DUDI2DI4
C
      RETURN
      END SUBROUTINE ANISOMAT
