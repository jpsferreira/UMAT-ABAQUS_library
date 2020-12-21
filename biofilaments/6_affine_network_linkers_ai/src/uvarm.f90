subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname, &
        nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord, &
        jmac,jmatyp,matlayo,laccfla)

! this subroutine is used to transfer sdv's from the uel
!  onto the dummy mesh for viewing.  note that an offset of
!  elemoffset is used between the real mesh and the dummy mesh.
!  if your model has more than elemoffset uel elements, then
!  this will need to be modified.

use global

include 'aba_param.inc'

character*80 cmname,orname
character*3 flgray(15)
dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)
integer i1

!     the dimensions of the variables flgray, array and jarray
!     must be set equal to or greater than 15.

!      uvar(1) = globalsdv(noel-elemoffset,npt,1)
!      for example
!      uvar(2) = globalsdv(noel-elemoffset,npt,2)
do i1=1,nsdv
uvar(i1) = globalsdv(noel-elemoffset,npt,i1)
enddo
!      uvar(3) = globalsdv(noel-elemoffset,npt,3)
!      uvar(4) = globalsdv(noel-elemoffset,npt,4)

return
end subroutine uvarm
 
