SUBROUTINE solvephi(root,args,nargs,rootold)



! This subroutine will numerically solve for the polymer
!  volume fraction based on the current osmotic pressure
!  and the previous state.  See numerical recipies RTSAFE.

use global

real(8), INTENT(OUT)                      :: root
real(8), INTENT(IN)                       :: args(nargs)
INTEGER, INTENT(IN OUT)                  :: nargs
real(8), INTENT(IN OUT)                   :: rootold

INTEGER :: j

real(8) f,df,fl,fh,xl,xh,x1,x2,swap, dxold
real(8) dx,  temp,rootmax,rootmin

INTEGER, PARAMETER :: maxit=50
real(8), PARAMETER :: xacc=1.d-6

rootmax = 0.9999D0 ! corresponds to nearly 100% dry polymer
rootmin = 0.05D0   ! corresponds to nearly 100% fluid

x1 = rootmin
x2 = rootmax
CALL phifunc(x1,fl,df,args,nargs)
CALL phifunc(x2,fh,df,args,nargs)

IF(fl*fh >= zero) THEN
  root = rootold
  WRITE(*,*) 'FYI, root not bracketed on phi'
  WRITE(*,*) 'fl=',fl
  WRITE(*,*) 'fh=',fh
  WRITE(*,*) 'rootOld=',rootold
!write(80,*) 'FYI, the root is not bracketed on phi'
!write(80,*) 'fl=',fl
!write(80,*) 'fh=',fh
!write(80,*) 'rootOld=',rootOld
  
  WRITE(*,*) 'mu =',args(1)
  WRITE(*,*) 'mu0=',args(2)
  WRITE(*,*) 'Rgas=',args(3)
  WRITE(*,*) 'theta=',args(4)
  WRITE(*,*) 'chi=',args(5)
  WRITE(*,*) 'Vmol=',args(6)
  WRITE(*,*) 'Kbulk=',args(7)
  WRITE(*,*) 'detF=',args(8)
  
  CALL xit()
  RETURN
END IF


!  ORIENT THE SEARCH SO THAT F(XL) < 0.

IF( fl < 0.d0 ) THEN
  xl = x1
  xh = x2
ELSE
  xh = x1
  xl = x2
  swap = fl
  fl = fh
  fh = swap
END IF

!  INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
!  BEFORE LAST'', AND THE LAST STEP

IF(rootold < rootmin) rootold = rootmin
IF(rootold > rootmax) rootold = rootmax
root = rootold !0.5D0 *( X1 + X2)
dxold = DABS(x2 - x1)
dx = dxold

CALL phifunc(root,f,df,args,nargs)


!   LOOP OVER ALLOWED ITERATIONS

DO  j = 1,maxit
  
!   BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
!   FAST ENOUGH.
  
  IF( ((root-xh)*df - f)*((root - xl)*df -f) >= 0.d0  &
        .OR. DABS(2.d0*f) > DABS(dxold*df) ) THEN
    
    dxold = dx
    dx = 0.5D0*(xh-xl)
    root = xl + dx
    IF( xl == root ) THEN
      
!   CHANGE IN ROOT IS NEGLIGIBLE
      
      RETURN
    END IF
    
  ELSE
    
!   NEWTON STEP IS ACCEPTABLE. TAKE IT.
    
    dxold = dx
    dx = f/df
    temp = root
    root = root - dx
    IF( temp == root) THEN
      
!    CHANGE IN ROOT IS NEGLIGIBLE
      
      RETURN
    END IF
    
  END IF
  
!  CONVERVEGENCE CRITERION
  
  IF( DABS(dx) < xacc) RETURN
  
  
!   THE ONE NEW FUNCTION EVALUATION PER ITERATION
  
  CALL phifunc(root,f,df,args,nargs)
  
  
!  MAINTAIN THE BRACKET ON THE ROOT
  
  IF( f < 0.d0) THEN
    xl = root
    fl = f
  ELSE
    xh = root
    fh = f
  END IF
  
END DO

WRITE(*,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'
WRITE(80,'(/1X,A)') 'solvePhi EXCEEDING MAXIMUM ITERATIONS'

RETURN
END SUBROUTINE solvephi

!***************************************************************************
