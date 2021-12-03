SUBROUTINE fibdir(fib,st0,st,noel,ndi,vorif,vd,distgr,dfgrd1)


DOUBLE PRECISION, INTENT(IN)             :: fib(NE,4)
DOUBLE PRECISION, INTENT(OUT)            :: st0(3,3)
DOUBLE PRECISION, INTENT(OUT)            :: st(3,3)
INTEGER, INTENT(IN OUT)                  :: noel
INTEGER, INTENT(IN)                      :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: vorif(3)
DOUBLE PRECISION, INTENT(OUT)            :: vd(3)
DOUBLE PRECISION, INTENT(IN OUT)         :: distgr(3,3)
DOUBLE PRECISION, INTENT(IN)             :: dfgrd1(3,3)


INTEGER :: inoel,i,j,i1,j1
DOUBLE PRECISION :: sum1, dnorm



inoel=0
i=0
DO i=1,nelem
!               ELEMENT IDENTIFICATION
  IF(noel == INT(fib(i,1))) THEN
    inoel=i
  END IF
END DO

!     FIB - FIBER ORIENTATION
dnorm=DSQRT(fib(inoel,2)*fib(inoel,2)+ fib(inoel,3)*fib(inoel,3)+  &
    fib(inoel,4)*fib(inoel,4))

!       UNDERFORMED FIBER ORIENTATION TENSOR

DO i1=1,ndi
  j1=i1+1
!       FIBER ORIENTATION NORMALIZED VECTOR - FAMILY 1
  vorif(i1)=fib(inoel,j1)/dnorm
END DO


DO i=1,ndi
  sum1=zero
  DO j=1,ndi
    sum1=sum1+dfgrd1(i,j)*vorif(j)
  END DO
!     FIBER DIRECTIONS IN THE DEFORMED CONFIGURATION
!               -FAMILY 1
  vd(i)=sum1
END DO
dnorm=DSQRT(vd(1)*vd(1)+ vd(2)*vd(2)+  &
    vd(3)*vd(3))
!           COSINE OF THE ANGLE BETWEEN FIBERS


!--------------------------------------------------------------------------
DO i=1,ndi
  DO j=1,ndi
!       STRUCTURAL TENSOR - FAMILY 1
    st0(i,j)=vorif(i)*vorif(j)
  END DO
END DO

!       STRUCTURE TENSOR IN THE DEFORMED CONFIGURATION - FAMILY 1
st=matmul(st0,transpose(distgr))
st=matmul(distgr,st)


RETURN
END SUBROUTINE fibdir
