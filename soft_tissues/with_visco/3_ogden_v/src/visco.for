      SUBROUTINE VISCO(PK,CMAT,VV,PKVOL,PKISO,CMATVOL,CMATISO,DTIME,
     1                                              VSCPROPS,STATEV,NDI)
C>    VISCOUS DISSIPATION: MAXWELL SPRINGS AND DASHPOTS SCHEME
      IMPLICIT NONE
      INCLUDE 'param_umat.inc'
C
      DOUBLE PRECISION PK(NDI,NDI),PKVOL(NDI,NDI),PKISO(NDI,NDI),
     1                  CMAT(NDI,NDI,NDI,NDI),CMATVOL(NDI,NDI,NDI,NDI),
     2                  CMATISO(NDI,NDI,NDI,NDI),VSCPROPS(6)
      DOUBLE PRECISION Q(NDI,NDI),QV(NDI,NDI),HV(NDI,NDI),
     1                  HV0(NDI,NDI),STATEV(NSDV)
      DOUBLE PRECISION DTIME,TETA,TAU,AUX,AUXC
      INTEGER I1,J1,K1,L1,NDI,VV,V1
C      
      Q=ZERO
      QV=ZERO
      HV=ZERO
      AUXC=ZERO
C      
C     ( GENERAL MAXWELL DASHPOTS)
      DO V1=1,VV 
C      
      TAU=VSCPROPS(2*V1-1)
      TETA=VSCPROPS(2*V1)
C
C      READ STATE VARIABLES
      CALL HVREAD(HV,STATEV,V1,NDI)
      HV0=HV
C        RALAXATION TENSORS      
      CALL RELAX(QV,HV,AUX,HV0,PKISO,DTIME,TAU,TETA,NDI)
      AUXC=AUXC+AUX     
C        WRITE STATE VARIABLES      
      CALL HVWRITE(STATEV,HV,V1,NDI)
C
      Q=Q+QV
C
      END DO
C              
      AUXC=ONE+AUXC
      PK=PKVOL+PKISO
C 
      DO I1=1,NDI
       DO J1=1,NDI
        PK(I1,J1)=PK(I1,J1)+Q(I1,J1)
        DO K1=1,NDI
         DO L1=1,NDI
          CMAT(I1,J1,K1,L1)= CMATVOL(I1,J1,K1,L1)+
     1                        AUXC*CMATISO(I1,J1,K1,L1)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
C
C
      RETURN
      END SUBROUTINE VISCO
