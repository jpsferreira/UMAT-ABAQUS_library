      SUBROUTINE UEXTERNALDBB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'aba_param.inc'
      INCLUDE 'param_umat.inc'
C       this subroutine get the fibers directions resorting the
C      .inc files and read the fibers direction
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KDATA/PTS
      DIMENSION TIME(2)
      DOUBLE PRECISION PTS(NPTS,2)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
        CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
          FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIRA
          OPEN(11,FILE=FILENAME)
          DO I=1,NPTS
             READ(11,*) (PTS(I,J),J=1,2)
          END DO
           CLOSE(11)
      END IF
C
      RETURN
C
      END
