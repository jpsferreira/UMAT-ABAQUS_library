      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C>    READ FILAMENTS ORIENTATION AND PREFERED DIRECTIONS
      INCLUDE 'aba_param.inc'
      INCLUDE 'param_umat.inc'
C       this subroutine get the directions and weights for
C      the numerical integration
C
C     UEXTERNAL just called once; work in parallel computing
      COMMON /KFIL/MF0
      COMMON /KFILR/RW
      COMMON /KFILP/PREFDIR
      COMMON /KFILK/KCH
      COMMON /KFILF/FRAC0
C      
      DIMENSION TIME(2)
      DOUBLE PRECISION  SPHERE(NWP,5),MF0(NWP,3),RW(NWP),
     1                  PREFDIR(NELEM,4),KCH(7),FRAC0(4),CHEM(2,7)
      CHARACTER(256) FILENAME
      CHARACTER(256) JOBDIR
C     LOP=0 --> START OF THE ANALYSIS
      IF(LOP.EQ.0.OR.LOP.EQ.4) THEN
C
       CALL GETOUTDIR(JOBDIR,LENJOBDIR)
C
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR1
         OPEN(15,FILE=FILENAME)
         DO I=1,NWP
            READ(15,*) (SPHERE(I,J),J=1,5)
         END DO
          CLOSE(15)
C
         DO I=1,NWP
           DO J=1,3
            K=J+1
           MF0(I,J)=SPHERE(I,K)
          END DO
          RW(I)=TWO*SPHERE(I,5)
         END DO
C
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR2
         OPEN(16,FILE=FILENAME)
         DO I=1,NELEM
            READ(16,*) (PREFDIR(I,J),J=1,4)
         END DO
         CLOSE(16)
C
         FILENAME=JOBDIR(:LENJOBDIR)//'/'//DIR3
         OPEN(17,FILE=FILENAME)
         DO I=1,2
            READ(17,*) (CHEM(I,J),J=1,7)
         END DO   
         DO J=1,4
         FRAC0(J)=CHEM(1,J)
         END DO
         DO J=1,7
         KCH(J)=CHEM(2,J)
         END DO
         CLOSE(17)
C         
      END IF    
C
      RETURN
C
      END SUBROUTINE UEXTERNALDB
