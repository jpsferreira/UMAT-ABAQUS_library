SUBROUTINE getoutdir(outdir, lenoutdir)



!>     GET CURRENT WORKING DIRECTORY
INCLUDE 'aba_param.inc'


CHARACTER (LEN=256), INTENT(IN OUT)      :: outdir
INTEGER, INTENT(OUT)                     :: lenoutdir



CALL getcwd(outdir)
!        OUTDIR=OUTDIR(1:SCAN(OUTDIR,'\',BACK=.TRUE.)-1)
lenoutdir=len_trim(outdir)

RETURN
END SUBROUTINE getoutdir
