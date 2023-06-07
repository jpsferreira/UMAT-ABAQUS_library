      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C>    READ MESH DATA
      INCLUDE 'aba_param.inc'
      INCLUDE 'param_umat.inc'
C
C     UEXTERNAL just called once; work in parallel computing
C     ADD COMMON BLOCKS HERE IF NEEDED (and in UMAT)
C      COMMON /KBLOCK/KBLOCK
      COMMON /KFIB/FIBORI
C           
      REAL*8 DTIME
      DIMENSION TIME(2)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
      INTEGER  LENJOBDIR
      
      REAL*8 FIBORI(NELEM,4)
      
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
       CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C        DIR1 DEFNIED IN param_umat.inc
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
C
         OPEN(15,FILE=FILENAME)
         DO I=1,NELEM
            READ(15,*) (FIBORI(I,J),J=1,4)
         END DO
          CLOSE(15)
!C         
      END IF    
C
      RETURN
C
      END SUBROUTINE UEXTERNALDB
