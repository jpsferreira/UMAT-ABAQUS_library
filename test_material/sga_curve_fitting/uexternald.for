      SUBROUTINE UEXTERNALD(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'aba_param.inc'
      INCLUDE 'param_umat.inc'
C       this subroutine get the fibers directions resorting the
C      .inc files and read the fibers direction
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KPAR/MATPAR
      DIMENSION TIME(2)
      DOUBLE PRECISION MATPAR(NNPAR)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
        CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIRC
          OPEN(12,FILE=FILENAME)
C          DO I=1,4
             READ(12,*) (MATPAR(I),I=1,NNPAR)
C          END DO
           CLOSE(12)
      END IF
C
      RETURN
C
      END
