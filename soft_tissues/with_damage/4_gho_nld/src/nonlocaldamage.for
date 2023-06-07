      SUBROUTINE NONLOCALDAMAGE(LOP)
C
      INCLUDE 'aba_param.inc'
      INCLUDE 'param_umat.inc'
C
      COMMON /KPASSAGE/CMNV
C
      DOUBLE PRECISION SEF,SEF0,DMG,DMGR,DMGDIFF,BETA,PSIHALF
      DOUBLE PRECISION DMGNL,DMGNLR,DMGDIFFNL
      DOUBLE PRECISION CMNV(NUEL,NIPP,6)

C     NONLOCAL AVERAGING AT THE END OF EACH LOCAL INCREMENT (LOP=2)
      IF(LOP.EQ.2) THEN
C
      DO NOEL=1,NUEL ! select element
       DO KK=1,NIPP ! select GP
        ALPHAW=ZERO
        RVOL=ZERO
        DMGNL=ZERO
        DO  II=1,NUEL ! search in element
          DO JJ=1,NIPP ! cycle trough GP of search element
             DIST=(CMNV(NOEL,KK,1)-CMNV(II,JJ,1))**2+
     1                    (CMNV(NOEL,KK,2)-CMNV(II,JJ,2))**2+
     2                    (CMNV(NOEL,KK,3)-CMNV(II,JJ,3))**2
             DISTS=SQRT(DIST)
C
       IF((DISTS.LE.DISTANCE))THEN
          ALPHAW=ONE-DIST/(DISTANCE**TWO)
          ALPHAW=ALPHAW*ALPHAW
          DMGNL=DMGNL+ALPHAW*CMNV(II,JJ,5)*CMNV(II,JJ,4)
          RVOL=RVOL+CMNV(II,JJ,4)*ALPHAW
       ENDIF
        END DO
      END DO
C     CALC NONLOCAL DMG
          DMGNL=DMGNL/RVOL
          CMNV(NOEL,KK,6)=DMGNL
        ENDDO
       ENDDO
C
      END IF

C
      RETURN
C
      END SUBROUTINE NONLOCALDAMAGE
