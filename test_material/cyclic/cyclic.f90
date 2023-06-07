!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'aba_param.inc'
INCLUDE 'param_umat.inc'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 13, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),STATEVP(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

CHARACTER*8 filename
CHARACTER(*), PARAMETER :: uni_f_path = './stress_curves/freq_sweep/uniaxial/'
CHARACTER(*), PARAMETER :: bi_f_path = './stress_curves/freq_sweep/equibiaxial/'
CHARACTER(*), PARAMETER :: sh_f_path = './stress_curves/freq_sweep/shear/'
CHARACTER(*), PARAMETER :: ssh_f_path = './stress_curves/freq_sweep/sshear/'
!
CHARACTER(*), PARAMETER :: uni_a_path = './stress_curves/amp_sweep/uniaxial/'
CHARACTER(*), PARAMETER :: bi_a_path = './stress_curves/amp_sweep/equibiaxial/'
CHARACTER(*), PARAMETER :: sh_a_path = './stress_curves/amp_sweep/shear/'
CHARACTER(*), PARAMETER :: ssh_a_path = './stress_curves/amp_sweep/sshear/'
!
integer un
!
real*8,ALLOCATABLE::curvepoints(:,:)
real*8::maxstresstime,maxstress,maxgammatime,maxgamma,minstresstime,minstress
real*8:: maxstretch, maxstretchtime
real*8:: delta,lossmodulus

i=1.0d0
j=1.0d0
DO i=1,NTENS
    DO j=1,NTENS
        DDSDDE(i,j)=0.0D0
    ENDDO
    STRESS(i)=0.0D0
ENDDO
!!!
time(1)=ZERO
time(2)=ZERO
call UEXTERNALDB(0,0,time,ZERO,0,0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=2.d0/1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELASTIC PARAM
PROPS(2)=1.00d0  ! C10
PROPS(3)=0.00d0 ! C01
PROPS(4)=1.00d0 !K1
PROPS(5)=1.d0 !K2
PROPS(6)=0.1d0 !kdisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(7)=1
!tau1 %
PROPS(8)=2.0d0
!teta1
PROPS(9)=0.835d0
!tau2 %
PROPS(10)=1.2d0
!teta2
PROPS(11)=7.0d0
!tau3 %
PROPS(12)=12.d0
!teta3
PROPS(13)=2.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STATEV=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
erf=ZERO
RHO=ZERO
PI=FOUR*ATAN(ONE)
!
!################################################################################################!
NSTEPS=100 !NUMBER OF POINTS PER CYCLE
!
!################################################################################################!
!CYLCIC RELATED VARIABLES
NTESTSF=61 !NUMBER OF TESTS
NTESTSA=31 !NUMBER OF TESTS
NTESTSAS=21 !NUMBER OF TESTS
!
NCYCLES=3 !NUMBER OF CYCLES PER TEST
!
ALLOCATE(CURVEPOINTS(NTESTSF,3))
!
!FREQUENCY SWEEP
FREQMIN=0.003d0
DFREQ=0.01d0
!
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - FREQ VALUES
DO I1=1,NTESTSF
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,1)=FREQMIN*10.0d0**COUNTER 
ENDDO
!
!AMPLITUDE SWEEP
FREQ0 = 1.d0
AMP_STRETCH=0.05D0
PRE_STRETCH=1.1D0

AMP_MIN = 0.01D0
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - AMP STRETCH VALUES
DO I1=1,NTESTSA
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,2)=AMP_MIN*10.0d0**COUNTER 
ENDDO
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - AMP GAMMA VALUES
DO I1=1,NTESTSAS
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,3)=AMP_MIN*10.0d0**COUNTER 
ENDDO

AMP_GAMMA=0.3D0
PRE_GAMMA=0.1D0
!########################################################################################!
!############################## |||   FREQUENCY SWEEP   ||| #############################!
!########################################################################################!
!
!################################## |||   UNIAXIAL   ||| ################################!
 CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=150, FILE='./stress_curves/freq_sweep/uniaxial.out', STATUS='UNKNOWN')
        rewind(150)
!
DO I1=1,NTESTSF   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=uni_f_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ
 DTIME=CYCLETIME/NSTEPS
 maxstress=0.0d0
 maxstresstime=0.0d0
 minstress=0.0d0
 minstresstime=0.0d0
!
!
MAXSTRETCH = 0.D0
MAXSTRETCHTIME= 0.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    STRETCH=AMP_STRETCH*DSIN(F*TIME(1)) + PRE_STRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=ONE/SQRT(STRETCH)
    DFGRD1(3,3)=ONE/SQRT(STRETCH)

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    
    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
     
      if(stress(1).gt.maxstress) THEN
      maxstress=stress(1)
      maxstresstime=time(1)
      endif
      
      if(stress(1).lt.minstress) THEN
      minstress=stress(1)
      minstresstime=time(1)
      endif
      
      if(STRETCH.gt.MAXSTRETCH) THEN
      MAXSTRETCH=STRETCH
      MAXSTRETCHTIME=time(1)
      endif
      
    endif
    
    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),STRETCH,STRESS(1)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=MAXSTRETCHTIME-maxstresstime !input - output shift

 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_stretch)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_stretch)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  ! write(*,*) props
  ! write(*,*)
  ! write(*,*) freq,storagemodulus,lossmodulus,tandelta
  ! write(*,*)
  write(150,*) freq,storagemodulus,lossmodulus,tandelta
  close(unit=un)
ENDDO !
close(150)


!################################## |||   BIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=151, FILE='stress_curves/freq_sweep/biaxial.out', STATUS='UNKNOWN')
        rewind(151)
TIME(1)     = 0.d0
!
DO I1=1,NTESTSF   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=bi_f_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ
 DTIME=CYCLETIME/NSTEPS
 maxstress=-1000.0d0
 maxstresstime=-10000.0d0
 minstress=1000.0d0
 minstresstime=10000.0d0
!
!
MAXSTRETCH = -1000.D0
MAXSTRETCHTIME=-1000.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    STRETCH=AMP_STRETCH*DSIN(F*TIME(1)) + PRE_STRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=STRETCH
    DFGRD1(3,3)=ONE/(STRETCH*STRETCH)
!

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    

    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
      if(stress(1).gt.maxstress) THEN
      maxstress=stress(1)
      maxstresstime=time(1)
      endif
      
      if(stress(1).lt.minstress) THEN
      minstress=stress(1)
      minstresstime=time(1)
      endif
      
      if(STRETCH.gt.MAXSTRETCH) THEN
      MAXSTRETCH=STRETCH
      MAXSTRETCHTIME=time(1)
      endif
    endif

    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),STRETCH,STRESS(1)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=-maxstresstime+MAXSTRETCHTIME
        
 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_stretch)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_stretch)*dsin(delta)!!pi/two)
 tandelta=lossmodulus/storagemodulus
  write(151,*) freq,storagemodulus,tandelta,lossmodulus
  close(unit=un)
ENDDO !
close(151)
! 
!################################## |||   SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=152, FILE='stress_curves/freq_sweep/shear.out', STATUS='UNKNOWN')
        rewind(152)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTSF   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=sh_f_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ
 DTIME=CYCLETIME/NSTEPS
 maxstress=-1000.0d0
 maxstresstime=-10000.0d0
 minstress=1000.0d0
 minstresstime=10000.0d0
!
!
MAXGAMMA = -1000.D0
MAXGAMMATIME=-1000.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    GAMMA=AMP_GAMMA*DSIN(F*TIME(1))
    DFGRD1(1,2)=GAMMA
    DFGRD1(2,1)=GAMMA
!

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    

    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
      if(stress(4).gt.maxstress) THEN
      maxstress=stress(4)
      maxstresstime=time(1)
      endif
      
      if(stress(4).lt.minstress) THEN
      minstress=stress(4)
      minstresstime=time(1)
      endif
      
      if(GAMMA.gt.MAXGAMMA) THEN
      MAXGAMMA=GAMMA
      MAXGAMMATIME=time(1)
      endif
    endif

    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),GAMMA,STRESS(4)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=-maxstresstime+maxgammatime
        
 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_gamma)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_gamma)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  write(152,*) freq,storagemodulus,tandelta,lossmodulus
  close(unit=un)
ENDDO !
close(152)
!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=153, FILE='stress_curves/freq_sweep/sshear.out', STATUS='UNKNOWN')
        rewind(153)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTSF   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=ssh_f_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ
 DTIME=CYCLETIME/NSTEPS
 maxstress=-1000.0d0
 maxstresstime=-10000.0d0
 minstress=1000.0d0
 minstresstime=10000.0d0
!
!
MAXGAMMA = -1000.D0
MAXGAMMATIME=-1000.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    GAMMA=AMP_GAMMA*DSIN(F*TIME(1))+PRE_GAMMA
    DFGRD1(1,2)=GAMMA
!

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    

    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
      if(stress(4).gt.maxstress) THEN
      maxstress=stress(4)
      maxstresstime=time(1)
      endif
      
      if(stress(4).lt.minstress) THEN
      minstress=stress(4)
      minstresstime=time(1)
      endif
      
      if(GAMMA.gt.MAXGAMMA) THEN
      MAXGAMMA=GAMMA
      MAXGAMMATIME=time(1)
      endif
    endif

    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),GAMMA,STRESS(4)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=-maxstresstime+maxgammatime
        
 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_gamma)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_gamma)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  write(153,*) freq,storagemodulus,tandelta,lossmodulus
  close(unit=un)
ENDDO !
close(153)


!########################################################################################!
!############################## |||   AMPLITUDE SWEEP   ||| #############################!
!########################################################################################!

!################################## |||   UNIAXIAL   ||| ################################!
 CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=150, FILE='./stress_curves/amp_sweep/uniaxial.out', STATUS='UNKNOWN')
        rewind(150)
!
!
DO I1=1,NTESTSA   !TESTS LOOP
!
 AMP_STRETCH=CURVEPOINTS(I1,2)
 TIME(1)=ZERO
 F=FREQ0*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=uni_a_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ0
 DTIME=CYCLETIME/NSTEPS
 maxstress=0.0d0
 maxstresstime=0.0d0
 minstress=0.0d0
 minstresstime=0.0d0
!
!
MAXSTRETCH = 0.D0
MAXSTRETCHTIME= 0.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    STRETCH=AMP_STRETCH*DSIN(F*TIME(1)) + PRE_STRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=ONE/SQRT(STRETCH)
    DFGRD1(3,3)=ONE/SQRT(STRETCH)


    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    
 
    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
     
      if(stress(1).gt.maxstress) THEN
      maxstress=stress(1)
      maxstresstime=time(1)
      endif
      
      if(stress(1).lt.minstress) THEN
      minstress=stress(1)
      minstresstime=time(1)
      endif
      
      if(STRETCH.gt.MAXSTRETCH) THEN
      MAXSTRETCH=STRETCH
      MAXSTRETCHTIME=time(1)
      endif
      
    endif
    
    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),STRETCH,STRESS(1)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=MAXSTRETCHTIME-maxstresstime !input - output shift

 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_stretch)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_stretch)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  ! write(*,*) props
  ! write(*,*)
  ! write(*,*) freq,storagemodulus,lossmodulus,tandelta
  ! write(*,*)
  write(150,*) AMP_STRETCH,storagemodulus,lossmodulus,tandelta
  close(unit=un)
ENDDO !
close(150)

!################################## |||   BIAXIAL   ||| ################################!
 CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=150, FILE='./stress_curves/amp_sweep/biaxial.out', STATUS='UNKNOWN')
        rewind(150)
!
!
DO I1=1,NTESTSA   !TESTS LOOP
!
 AMP_STRETCH=CURVEPOINTS(I1,2)
 TIME(1)=ZERO
 F=FREQ0*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=bi_a_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ0
 DTIME=CYCLETIME/NSTEPS
 maxstress=0.0d0
 maxstresstime=0.0d0
 minstress=0.0d0
 minstresstime=0.0d0
!
!
MAXSTRETCH = 0.D0
MAXSTRETCHTIME= 0.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    STRETCH=AMP_STRETCH*DSIN(F*TIME(1)) + PRE_STRETCH
    DFGRD1(1,1)=STRETCH
    DFGRD1(2,2)=STRETCH
    DFGRD1(3,3)=ONE/(STRETCH*STRETCH)


    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    
 
    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
     
      if(stress(1).gt.maxstress) THEN
      maxstress=stress(1)
      maxstresstime=time(1)
      endif
      
      if(stress(1).lt.minstress) THEN
      minstress=stress(1)
      minstresstime=time(1)
      endif
      
      if(STRETCH.gt.MAXSTRETCH) THEN
      MAXSTRETCH=STRETCH
      MAXSTRETCHTIME=time(1)
      endif
      
    endif
    
    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),STRETCH,STRESS(1)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=MAXSTRETCHTIME-maxstresstime !input - output shift

 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_stretch)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_stretch)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  ! write(*,*) props
  ! write(*,*)
  ! write(*,*) freq,storagemodulus,lossmodulus,tandelta
  ! write(*,*)
  write(150,*) AMP_STRETCH,storagemodulus,lossmodulus,tandelta
  close(unit=un)
ENDDO !
close(150)

! 
!################################## |||   SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=152, FILE='stress_curves/amp_sweep/shear.out', STATUS='UNKNOWN')
        rewind(152)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTSAS   !TESTS LOOP
!
 AMP_GAMMA=CURVEPOINTS(I1,3)
 TIME(1)=ZERO
 F=FREQ0*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=sh_a_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ0
 DTIME=CYCLETIME/NSTEPS
 maxstress=-1000.0d0
 maxstresstime=-10000.0d0
 minstress=1000.0d0
 minstresstime=10000.0d0
!
!
MAXGAMMA = -1000.D0
MAXGAMMATIME=-1000.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    GAMMA=AMP_GAMMA*DSIN(F*TIME(1))
    DFGRD1(1,2)=GAMMA
    DFGRD1(2,1)=GAMMA
!

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    

    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
      if(stress(4).gt.maxstress) THEN
      maxstress=stress(4)
      maxstresstime=time(1)
      endif
      
      if(stress(4).lt.minstress) THEN
      minstress=stress(4)
      minstresstime=time(1)
      endif
      
      if(GAMMA.gt.MAXGAMMA) THEN
      MAXGAMMA=GAMMA
      MAXGAMMATIME=time(1)
      endif
    endif

    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),GAMMA,STRESS(4)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=-maxstresstime+maxgammatime
        
 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_gamma)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_gamma)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  write(152,*) AMP_GAMMA,storagemodulus,tandelta,lossmodulus
  close(unit=un)
ENDDO !
close(152)
!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=153, FILE='stress_curves/amp_sweep/sshear.out', STATUS='UNKNOWN')
        rewind(153)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTSAS   !TESTS LOOP
!
 AMP_GAMMA=CURVEPOINTS(I1,3)
 TIME(1)=ZERO
 F=FREQ0*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=ssh_a_path//filename,status='UNKNOWN')
!
 CYCLETIME=ONE/FREQ0
 DTIME=CYCLETIME/NSTEPS
 maxstress=-1000.0d0
 maxstresstime=-10000.0d0
 minstress=1000.0d0
 minstresstime=10000.0d0
!
!
MAXGAMMA = -1000.D0
MAXGAMMATIME=-1000.D0
!
 DO KK=1,NCYCLES    !CYCLES LOOP

   DO KSTEP=1,NSTEPS   !DEFORMATION LOOP
!
    GAMMA=AMP_GAMMA*DSIN(F*TIME(1))+PRE_GAMMA
    DFGRD1(1,2)=GAMMA
!

    CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
    DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!    

    if(kk.eq.ncycles)THEN  !LAST CYCLE VERIFICATIONS 
      if(stress(4).gt.maxstress) THEN
      maxstress=stress(4)
      maxstresstime=time(1)
      endif
      
      if(stress(4).lt.minstress) THEN
      minstress=stress(4)
      minstresstime=time(1)
      endif
      
      if(GAMMA.gt.MAXGAMMA) THEN
      MAXGAMMA=GAMMA
      MAXGAMMATIME=time(1)
      endif
    endif

    TIME(1)   = TIME(1)+DTIME
    write(un,*) TIME(1),GAMMA,STRESS(4)
      !
   ENDDO !nsteps
 ENDDO !ncycles

 stressamp=maxstress-minstress
 delta=-maxstresstime+maxgammatime
        
 delta=delta*two*PI/cycletime
        
!   storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
 storagemodulus=(maxstress/amp_gamma)*dcos(delta)!-pi/two)
!  lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
 lossmodulus=(maxstress/amp_gamma)*dsin(delta)!-pi/two)
 tandelta=lossmodulus/storagemodulus
  write(153,*) AMP_GAMMA,storagemodulus,tandelta,lossmodulus
  close(unit=un)
ENDDO !
close(153)
!
! !################################################################################################!
CALL SYSTEM('gnuplot -p data_cyclicfreq.plt')
! !################################################################################################!
!
! !################################################################################################!
CALL SYSTEM('gnuplot -p data_cyclicamp.plt')
! !################################################################################################!
END PROGRAM
