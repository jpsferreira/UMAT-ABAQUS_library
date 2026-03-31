!> @brief ABAQUS UEXTERNALDB callback for loading external data files.
!>
!> Called by ABAQUS at the start of an analysis to load quadrature directions
!> (for network models) and preferred fiber directions.
!>
!> Required data files (only for network models):
!>   - sphere_intXXc.inp : NWP lines of [x, y, z, weight] quadrature on unit sphere
!>   - prefdir.inp       : NELEM lines of [x, y, z, angle] preferred directions
!>
!> Data is stored in COMMON blocks accessible by the UMAT.

subroutine uexternaldb(lop, lrestart, time, dtime, kstep, kinc)

  implicit none
  include 'aba_param.inc'

  integer, intent(in out) :: lop, lrestart, kstep, kinc
  real,    intent(in out) :: time(2)
  real(8), intent(in out) :: dtime

  ! Maximum quadrature points (must match umat_builder network_contribution)
  integer, parameter :: MAX_NWP  = 720
  integer, parameter :: MAX_ELEM = 1

  double precision :: mf0(MAX_NWP, 3), rw(MAX_NWP)
  double precision :: prefdir(MAX_ELEM, 4)
  integer          :: nwp_active

  common /kfil/  mf0
  common /kfilr/ rw
  common /kfilp/ prefdir
  common /knwp/  nwp_active

  character(len=256) :: filename, outdir
  integer :: lenoutdir, njob, ios, i

  ! LOP == 0: called at start of analysis
  if (lop == 0) then

    ! Get ABAQUS job output directory
    call getoutdir(outdir, lenoutdir)

    ! ---------------------------------------------------------------
    ! Load sphere integration quadrature points
    ! Try common sizes: 720, 300, 60, 21
    ! ---------------------------------------------------------------
    nwp_active = 0
    mf0 = 0.0d0
    rw  = 0.0d0

    ! Try sphere_int720c.inp first, then smaller
    filename = trim(outdir) // '/sphere_int720c.inp'
    open(unit=80, file=trim(filename), status='old', iostat=ios)
    if (ios /= 0) then
      filename = trim(outdir) // '/sphere_int300c.inp'
      open(unit=80, file=trim(filename), status='old', iostat=ios)
    end if
    if (ios /= 0) then
      filename = trim(outdir) // '/sphere_int60c.inp'
      open(unit=80, file=trim(filename), status='old', iostat=ios)
    end if
    if (ios /= 0) then
      filename = trim(outdir) // '/sphere_int21c.inp'
      open(unit=80, file=trim(filename), status='old', iostat=ios)
    end if

    if (ios == 0) then
      ! Read quadrature data: x, y, z, weight per line
      do i = 1, MAX_NWP
        read(80, *, iostat=ios) mf0(i,1), mf0(i,2), mf0(i,3), rw(i)
        if (ios /= 0) exit
        nwp_active = i
      end do
      close(80)
    end if

    ! ---------------------------------------------------------------
    ! Load preferred direction (for orientation density in networks)
    ! ---------------------------------------------------------------
    prefdir = 0.0d0
    filename = trim(outdir) // '/prefdir.inp'
    open(unit=81, file=trim(filename), status='old', iostat=ios)
    if (ios == 0) then
      do i = 1, MAX_ELEM
        read(81, *, iostat=ios) prefdir(i,1), prefdir(i,2), prefdir(i,3), prefdir(i,4)
        if (ios /= 0) exit
      end do
      close(81)
    end if

  end if

  return
end subroutine uexternaldb
