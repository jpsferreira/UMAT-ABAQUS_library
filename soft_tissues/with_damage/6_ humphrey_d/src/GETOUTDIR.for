      SUBROUTINE GETOUTDIR(OUTDIR, LENOUTDIR)
C>     GET CURRENT WORKING DIRECTORY
      INCLUDE 'aba_param.inc'
C 
      CHARACTER*256 OUTDIR
      INTEGER LENOUTDIR
C
      CALL GETCWD(OUTDIR)
c        OUTDIR=OUTDIR(1:SCAN(OUTDIR,'\',BACK=.TRUE.)-1)
      LENOUTDIR=LEN_TRIM(OUTDIR)
C
      RETURN
      END SUBROUTINE GETOUTDIR