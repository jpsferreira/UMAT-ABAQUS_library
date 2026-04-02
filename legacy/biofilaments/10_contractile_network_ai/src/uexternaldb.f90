SUBROUTINE uexternaldb(lop,lrestart,time,dtime,kstep,kinc)



!>    READ FILAMENTS ORIENTATION AND PREFERED DIRECTIONS
use global
INCLUDE 'aba_param.inc'
!       this subroutine get the directions and weights for
!      the numerical integration

!     UEXTERNAL just called once; work in parallel computing

INTEGER, INTENT(IN OUT)                  :: lop
INTEGER, INTENT(IN OUT)                  :: lrestart
REAL, INTENT(IN OUT)                     :: time(2)
real(8), INTENT(IN OUT)                   :: dtime
INTEGER, INTENT(IN OUT)                  :: kstep
INTEGER, INTENT(IN OUT)                  :: kinc

COMMON /kfil/mf0
COMMON /kfilr/rw
COMMON /kfilp/prefdir
COMMON /kinit/init


DOUBLE PRECISION :: sphere(nwp,5),mf0(nwp,3),rw(nwp),  &
                    prefdir(nelem,4), init(2)
CHARACTER (LEN=256) ::  filename, jobdir
INTEGER :: lenjobdir,i,j,k

!     LOP=0 --> START OF THE ANALYSIS
IF(lop == 0.OR.lop == 4) THEN
  
  CALL getoutdir(jobdir,lenjobdir)
  
  !preferential direction
  filename=jobdir(:lenjobdir)//'/'//dir2
  OPEN(16,FILE=filename)
  DO i=1,nelem
    READ(16,*) (prefdir(i,j),j=1,4)
  END DO
  CLOSE(16)
  !initial conditions
  filename=jobdir(:lenjobdir)//'/'//dir3
  OPEN(18,FILE=filename)
  READ(18,*) (init(j),j=1,2)
  CLOSE(18)
  
END IF

RETURN

END SUBROUTINE uexternaldb
