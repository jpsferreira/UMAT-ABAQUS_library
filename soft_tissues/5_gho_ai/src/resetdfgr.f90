SUBROUTINE resetdfgrd(dfgrd,ndi)

use global
implicit none
INTEGER, INTENT(IN OUT)                  :: ndi
DOUBLE PRECISION, INTENT(OUT)            :: dfgrd(ndi,ndi)

dfgrd(1,1)=  one
dfgrd(1,2)=  zero
dfgrd(1,3)=  zero
dfgrd(2,1)=  zero
dfgrd(2,2)=  one
dfgrd(2,3)=  zero
dfgrd(3,1)=  zero
dfgrd(3,2)=  zero
dfgrd(3,3)=  one

END SUBROUTINE resetdfgrd
