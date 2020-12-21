module global

! This module is used to transfer SDV's from the UEL
!  to the UVARM so that SDV's can be visualized on a
!  dummy mesh
!
!  globalSdv(X,Y,Z)
!   X - element pointer
!   Y - integration point pointer
!   Z - SDV pointer
!
!  numElem
!   Total number of elements in the real mesh, the dummy
!   mesh needs to have the same number of elements, and
!   the dummy mesh needs to have the same number of integ
!   points.  You must set that parameter value here.
!
!  ElemOffset
!   Offset between element numbers on the real mesh and
!    dummy mesh.  That is set in the input file, and
!    that value must be set here the same.

INTEGER :: numelem,elemoffset,ERR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the number of UEL elements used here
PARAMETER(numelem=1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the offset here for UVARM plotting, must match input file!
PARAMETER(elemoffset=1000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8), allocatable :: globalsdv(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the control parameters to run the material-related routines

INTEGER NWP,NELEM,NCH,NCNT,NSDV,NVSC,NMECH,NQUEM
PARAMETER (NELEM=1)
PARAMETER (NWP=21,NVSC=27,NMECH=1,NCH=4,NCNT=4,NQUEM=5)
!      PARAMETER (NSDV=NWP+NCH+NVSC+NMECH+NCNT+NQUEM)
PARAMETER (NSDV=3)
DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
PARAMETER (ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
PARAMETER (THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
DOUBLE PRECISION HALF,THIRD
PARAMETER (HALF=0.5d0,THIRD=1.d0/3.d0)
CHARACTER(256) DIR1,DIR2,DIR3
PARAMETER (DIR1='sphere_int21c.inp',DIR2='prefdir.inp')
PARAMETER (DIR3='initial_cond.inp')


END module global
