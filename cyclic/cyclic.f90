!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'ABA_PARAM.INC'
INCLUDE 'PARAM_UMAT.INC'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 13, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),STATEVP(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

CHARACTER*8 filename
CHARACTER(*), PARAMETER :: uni_path = './stress_curves/freq_sweep/uniaxial/'
CHARACTER(*), PARAMETER :: bi_path = './stress_curves/freq_sweep/equibiaxial/'
CHARACTER(*), PARAMETER :: sh_path = './stress_curves/freq_sweep/shear/'
CHARACTER(*), PARAMETER :: ssh_path = './stress_curves/freq_sweep/sshear/'

integer un

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
PROPS(4)=0.00d0 !K1
PROPS(5)=0.1d0 !K2
PROPS(6)=0.1d0 !kdisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(7)=3
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
NSTEPS=400 !NUMBER OF POINTS PER CYCLE
!
STRETCH_MAX = 2.D0
STRETCH_INITIAL = 0.5d0
DSTRETCH    = (STRETCH_MAX-STRETCH_INITIAL)/NSTEPS
!
GAMMA_MAX = 0.6D0
GAMMA_INITIAL = -.6D0
DGAMMA   = (GAMMA_MAX-GAMMA_INITIAL)/NSTEPS
!################################################################################################!
!CYLCIC RELATED VARIABLES
NTESTS=61 !NUMBER OF TESTS
NCYCLES=3 !NUMBER OF CYCLES PER TEST
!
ALLOCATE(CURVEPOINTS(NTESTS,3))
!
!FREQUENCY SWEEP
FREQMIN=0.003d0
DFREQ=0.01d0
!
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS
DO I1=1,NTESTS
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,1)=FREQMIN*10.0d0**COUNTER 
ENDDO
!
!AMPLITUDE SWEEP
AMP_STRETCH=0.05D0
PRE_STRETCH=1.1D0

AMP_GAMMA=0.3D0
! 
!################################## |||   UNIAXIAL   ||| ################################!
 CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=150, FILE='./stress_curves/freq_sweep/uniaxial.out', STATUS='UNKNOWN')
        rewind(150)
!

!
DO I1=1,NTESTS   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=uni_path//filename,status='UNKNOWN')
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
ENDDO !
close(150)


!################################## |||   BIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=151, FILE='stress_curves/freq_sweep/biaxial.out', STATUS='UNKNOWN')
        rewind(151)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTS   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=bi_path//filename,status='UNKNOWN')
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
DO I1=1,NTESTS   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=sh_path//filename,status='UNKNOWN')
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
ENDDO !
close(152)
!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
 OPEN (UNIT=153, FILE='stress_curves/freq_sweep/sshear.out', STATUS='UNKNOWN')
        rewind(153)
TIME(1)     = 0.d0
!
!
DO I1=1,NTESTS   !TESTS LOOP
!
 FREQ=CURVEPOINTS(I1,1)
 TIME(1)=ZERO
 F=FREQ*TWO*PI
 WRITE(filename,fmt='(i0,a)')I1,'.out'
 un=un+2;
 OPEN(unit=un,file=ssh_path//filename,status='UNKNOWN')
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
ENDDO !
close(153)
! !################################################################################################!
CALL SYSTEM('gnuplot -p data_cyclic.plt')
! !################################################################################################!
! AMPLITUDE SWWEEP
  AMPMIN=0.00d0
  DAMP=0.01d0
!  DAMP=0.125d0
  FREQ=0.05d0
  NTEST=18
!  NTEST=1
  NCYCLES=3
  NSTEPS=400
!  CALL AMPSWEEP(AMPMIN,DAMP,FREQ,NTEST,NCYCLES,NSTEPS,PROPS)
!################################################################################################!
!################################################################################################!

!  ALLOCATE(curvepoints(NPOINTS,3))
!  do i=1,npoints
!     curvepoints(i,1)=AMPMIN+I*ONE*DAMP
!  enddo

!      DO I1=1,NPOINTS
! !C
!        AMP=curvepoints(I1,1)
!        TIME(1)=ZERO
!        PI=FOUR*ATAN(ONE)
!        F=FREQ*TWO*PI
! !
!        write(filename,fmt='(i0,a)')I1,'.out'
!        un=un+2;
!        open(unit=un,file=filename,status='UNKNOWN')
! !
!        cycletime=one/freq
! !
!        dtime=cycletime/nsteps
! !
!        maxstress=0.0d0
!        maxstresstime=0.0d0
!        minstress=0.0d0
!        minstresstime=0.0d0
!        maxsgama=0.d0
!        maxgamatime=0.d0

!        DO KK=1,NCYCLES

!           DO KSTEP=1,NSTEPS

!              GAMMA=AMP*DSIN(F*TIME(1))
!              DFGRD1(1,2)=GAMMA
! !C>---------------------------------------------------------------------
!      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

!      if(kk.eq.ncycles)THEN  !pesquisa apenas no ultimo ciclo
! !      write(50,*) time(1),gamma,stress(4)
!      if(stress(4).gt.maxstress) THEN
!      maxstress=stress(4)
!      maxstresstime=time(1)
!      endif

!      if(stress(4).lt.minstress) THEN
!      minstress=stress(4)
!      minstresstime=time(1)
!      endif

!      if(gamma.gt.maxgamma) THEN
!      maxsgama=gamma
!      maxgamatime=time(1)
!      endif
!      endif

!      time(1)=time(1)+dtime
!      write(un,*) time(1),gamma,stress(4)
!           ENDDO !KSTEP

!        !   write(*,*) KK,time(1)
!        ENDDO !NCYCLES

!        stressamp=maxstress-minstress
!        delta=-maxstresstime+maxgamatime

!        delta=delta*TWO*PI/cycletime

! !        storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
!        storagemodulus=(maxstress/amp)*dcos(delta-pi/two)
! !           lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
!           lossmodulus=(maxstress/amp)*dsin(delta-pi/two)

! !              write(*,*) time(1),gamma,stress(4)
!       ! close(unit=un)
! !       write(*,*) KK,maxstress,minstress,delta
!       write(*,*) AMP,storagemodulus,lossmodulus
!       write(150,*) AMP,storagemodulus,lossmodulus

!      ENDDO

!  close(150)
!  CALL SYSTEM('gnuplot -p data_xyg.plt')
! 





! !################################################################################################!
! !################################################################################################!
! !################################################################################################!
! !################################################################################################!
! ! FREQ SWEEP
! !  FREQMIN=0.00d0
! !  DFREQ=0.01d0
! !  AMP=0.005d0
! !  NTEST=61
! !  NCYCLES=3
! !  NSTEPS=200
! !!  CALL FREQSWEEP(FREQMIN,DFREQ,AMP,NTEST,NCYCLES,NSTEPS,PROPS)
! !
! !!################################################################################################!
! !!################################################################################################!
! !! FREQ SWEEP WORKING
! !  FREQMIN=1.00d0
! !  DFREQ=0.0d0
! !  AMP=0.005d0
! !  NTEST=1
! !  NCYCLES=3
! !  NSTEPS=400
! !  TIME=0.0d0
! !  FREQMIN=0.00d0
! !  DFREQ=0.01d0
! !  AMP=0.03d0
! !!  NTEST=30
! !  NCYCLES=10
! !  NSTEPS=400
! !  NPOINTS=71
! !  ALLOCATE(curvepoints(NPOINTS,3))
! !  counter=-0.05d0
! !
! !  ! o ciclo seguinte vai gerar 61 pontos uniformemente
! !  ! distribuidos numa escala logaritmica.
! !  do i=1,npoints
! !  counter=counter+0.05d0
! !  curvepoints(i,1)=0.003d0*10.0d0**counter
! !  enddo
! !
! !      DO I1=1,NPOINTS
! !!C
! !        freq=curvepoints(I1,1)
! !        TIME(1)=ZERO
! !        PI=FOUR*ATAN(ONE)
! !        F=FREQ*two*PI
! !!C        dtime=0.02d0
! !         write(filename,fmt='(i0,a)')I1,'.out'
! !        un=un+2;
! !        open(unit=un,file=filename,status='UNKNOWN')
! !
! !        cycletime=one/freq
! !
! !        dtime=cycletime/nsteps
! !!C
! !        maxstress=0.0d0
! !        maxstresstime=0.0d0
! !        minstress=0.0d0
! !        minstresstime=0.0d0
! !
! !        DO KK=1,NCYCLES
! !
! !           DO KSTEP=1,NSTEPS
! !
! !              GAMMA=AMP*DSIN(F*TIME(1))
! !              DFGRD1(1,2)=GAMMA
! !!C>---------------------------------------------------------------------
! !      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
! !      RPL,DDSDDT,DRPLDE,DRPLDT,&
! !      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
! !      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
! !      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !
! !      if(kk.eq.ncycles)THEN  !pesquisa apenas no ultimo ciclo
! !!      write(50,*) time(1),gamma,stress(4)
! !      if(stress(4).gt.maxstress) THEN
! !      maxstress=stress(4)
! !      maxstresstime=time(1)
! !      endif
! !
! !      if(stress(4).lt.minstress) THEN
! !      minstress=stress(4)
! !      minstresstime=time(1)
! !      endif
! !
! !      if(gamma.gt.maxgamma) THEN
! !      maxsgama=gamma
! !      maxgamatime=time(1)
! !      endif
! !      endif
! !
! !      time(1)=time(1)+dtime
! !      write(un,*) time(1),gamma,stress(4)
! !           ENDDO !KSTEP
! !
! !        !   write(*,*) KK,time(1)
! !        ENDDO !NCYCLES
! !
! !        stressamp=maxstress-minstress
! !        delta=-maxstresstime+maxgamatime
! !
! !        delta=delta*two*PI/cycletime
! !
! !!        storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
! !        storagemodulus=(maxstress/amp)*dcos(delta-pi/two)
! !!           lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
! !           lossmodulus=(maxstress/amp)*dsin(delta-pi/two)
! !           tandelta=lossmodulus/storagemodulus
! !
! !!              write(*,*) time(1),gamma,stress(4)
! !       ! close(unit=un)
! !       write(*,*) freq*(two*pi),storagemodulus,lossmodulus,tandelta
! !       write(150,*) freq*(two*pi),storagemodulus,tandelta,lossmodulus
! !
! !      ENDDO
! !
! !  close(150)
! !  CALL SYSTEM('gnuplot -p data_xyg.plt')
! !################################################################################################!
! !################################################################################################!
! ! SHEAR CYCLING TEST

!  DFGRD1(1,1)=  1.0D0
!  DFGRD1(1,2)=  0.0D0
!  DFGRD1(1,3)=  0.0d0
!  DFGRD1(2,1)=  0.0d0
!  DFGRD1(2,2)=  1.0D0
!  DFGRD1(2,3)=  0.0d0
!  DFGRD1(3,1)=  0.0d0
!  DFGRD1(3,2)=  0.0d0
!  DFGRD1(3,3)=  1.0D0
! !statev=0.d0
! !gamma=0.d0
! !time(1)=0.d0
! !dtime=0.02d0
! !kstep=1
! !!
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! ! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! ! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !freq=0.05d0
! !dtime=0.2
! !period=1.d0/freq
! !amp=0.12d0
! !angfreq=2.d0*4.d0*atan(1.d0)*freq
! !!gamma=amp*sin(angfreq*time(1))
! ! write(150,*) time(1),gamma,stress(4)
! !do k1=1,20 !cycles loop
! !write(*,*) k1,time(1)
! !do kstep=1,25
! !!gamma=gamma+0.0012d0
! !time(1)=time(1)+dtime
! !gamma=amp*sin(angfreq*time(1))
! !DFGRD1(1,2)=  gamma
! !!write(*,*) kstep,gamma
! !!write(*,*)
! ! CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! ! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! ! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! ! !write(*,*) statev(10:13)
! ! write(150,*) time(1),gamma,stress(4)
! !enddo
! !do kstep=26,75
! !time(1)=time(1)+dtime
! !gamma=amp*sin(angfreq*time(1))
! !DFGRD1(1,2)=  gamma
! !!write(*,*) kstep,gamma
! !!write(*,*)
! ! CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! ! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! ! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! ! ! write(*,*) statev(10:13)
! ! write(150,*) time(1),gamma,stress(4)
! !enddo
! !do kstep=76,100
! !time(1)=time(1)+dtime
! !gamma=amp*sin(angfreq*time(1))
! !DFGRD1(1,2)=  gamma
! !!write(*,*) kstep,gamma
! !!write(*,*)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !! write(*,*) statev(10:13)
! !write(150,*) time(1),gamma,stress(4)
! !enddo
! !write(*,*) k1,time(1)
! !enddo
! !
! !
! !close(150)
! !  CALL SYSTEM('gnuplot -p data_xy.plt')
! !################################################################################################!
! !0.740469148776867d0, 0.672090350852965d0, 0.d0
! !! CREEP RELAXATION TEST
! time(1)=0.d0
! dtime=0.00025d0
! !
! !dtime=0.1d0
! amp=0.125d0
! nsteps=100
! !
! !time(1)=0.d0
! !dtime=0.02d0
! !freq=2.5d0
! !dtime=1/(freq*4*nsteps)
! !!period=1.d0/freq
! !amp=0.125d0
! !angfreq=2.d0*4.d0*atan(1.d0)*freq
! !!gamma=0.125d0
! !do kstep=1,nsteps
! !gamma=amp*sin(angfreq*time(1))
! !DFGRD1(1,2)=gamma
! !CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !s1=stress(4)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !write(150,*) time(1),gamma,s1,stress(4)
! !time(1)=time(1)+dtime
! !enddo
! !
! !gamma=1.1d0
! !write(*,*)
! !DFGRD1(2,2)=1.d0/sqrt(gamma)
! !DFGRD1(1,1)=gamma
! !DFGRD1(3,3)=1.d0/sqrt(gamma)
! !gamma=0.1d0
! !DFGRD1(1,2)=gamma
! !dtime=0.1d0
! !do kstep=1,1000
! !!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !
! !write(150,*) gamma,time(1),stress(1)
! !time(1)=time(1)+dtime
! !enddo
! !close(150)
! !  CALL SYSTEM('gnuplot -p data_xy.plt')
! !################################################################################################!
! !!     TENSILE MONOTONIC LOAD TEST

!  DFGRD1(1,1)=  1.2D0
!  DFGRD1(1,2)=  0.0D0
!  DFGRD1(1,3)=  0.0d0
!  DFGRD1(2,1)=  0.0d0
!  DFGRD1(2,2)=  1/sqrt(DFGRD1(1,1))
!  DFGRD1(2,3)=  0.0d0
!  DFGRD1(3,1)=  0.0d0
!  DFGRD1(3,2)=  0.0d0
!  DFGRD1(3,3)=  1/sqrt(DFGRD1(1,1))


!  CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

!  write(*,*) stress
!   write(*,*)
!    write(*,*) DDSDDE

! ! DFGRD1(1,1)=  1.0D0
! ! DFGRD1(1,2)=  0.0D0
! ! DFGRD1(1,3)=  0.0d0
! ! DFGRD1(2,1)=  0.0d0
! ! DFGRD1(2,2)=  1.0D0
! ! DFGRD1(2,3)=  0.0d0
! ! DFGRD1(3,1)=  0.0d0
! ! DFGRD1(3,2)=  0.0d0
! ! DFGRD1(3,3)=  1.0D0
! !
! !gamma=0.d0
! !time(1)=0.d0
! !dtime=0.01d0
! !do kstep=1,100
! !!gamma=0.1942d0
! !gamma=gamma+0.03d0
! !!gamma=gamma+0.00025d0
! !DFGRD1(1,1)=1.d0+gamma
! !DFGRD1(2,2)=  1.0D0/sqrt(DFGRD1(1,1))
! !DFGRD1(3,3)=  1.0D0/sqrt(DFGRD1(1,1))
! !!!
! !!!s1=stress(4)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!!
! !write(150,*) time(1),gamma,stress(1)
! !!gamma=gamma+0.0002d0
! !time(1)=time(1)+dtime
! !
! !enddo
! ! CALL SYSTEM('gnuplot -p data_xy.plt')
! !################################################################################################!
! !write(50,*) time(1),gamma,s1,stress(4)
! !! SHEAR MONOTONIC LOAD TEST

! ! DFGRD1(1,1)=  1.0D0
! ! DFGRD1(1,2)=  0.0D0
! ! DFGRD1(1,3)=  0.0d0
! ! DFGRD1(2,1)=  0.0d0
! ! DFGRD1(2,2)=  1.0D0
! ! DFGRD1(2,3)=  0.0d0
! ! DFGRD1(3,1)=  0.0d0
! ! DFGRD1(3,2)=  0.0d0
! ! DFGRD1(3,3)=  1.0D0
! !
! !gamma=0.d0
! !time(1)=0.d0
! !dtime=0.01d0
! !do kstep=1,100
! !!gamma=0.1942d0
! !gamma=gamma+0.003d0
! !!gamma=gamma+0.00025d0
! !DFGRD1(1,2)=gamma
! !!!
! !!!s1=stress(4)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!!
! !write(150,*) time(1),gamma,stress(4)
! !!gamma=gamma+0.0002d0
! !time(1)=time(1)+dtime
! !
! !enddo


! !
! !time(1)=0.d0
! !dtime=0.2d0
! !gamma=0.19428d0
! !!gamma=0.0d0
! !!time(1)=time(1)+dtime
! !DFGRD1(1,2)=gamma
! !do kstep=1,100
! !!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !!s1=stress(4)
! !!
! !!!!!!!!!!!!!!!!!!
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !s2=stress(4)
! !write(150,*) time(1),gamma,STRESS(4)
! !!gamma=gamma+0.019428d0
! !!gamma=0.19428d0
! !time(1)=time(1)+dtime
! !!DFGRD1(1,2)=gamma
! !enddo
!  ! CALL SYSTEM('gnuplot -p data_xy.plt')
! !################################################################################################!
! !
! ! !DFGRD1(1,1)=  1.0D0
! ! DFGRD1(1,2)=  0.0D0
! ! DFGRD1(1,3)=  0.0d0
! ! DFGRD1(2,1)=  0.0d0
! ! DFGRD1(2,2)=  1.0D0
! ! DFGRD1(2,3)=  0.0d0
! ! DFGRD1(3,1)=  0.0d0
! ! DFGRD1(3,2)=  0.0d0
! ! DFGRD1(3,3)=  1.0D0
! !
! !gamma=0.d0
! !dtime=0.1
! !
! !!ISOMETRIC ACTIVATION
! !do kstep=1,1000
! !!gamma=0.1942d0
! !!gamma=gamma+0.0001942d0
! !DFGRD1(1,1)=1.d0
! !CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !s1=stress(1)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !write(50,*) time(1),gamma,s1,stress(1)
! !!gamma=gamma+0.0002d0
! !!
! !time(1)=time(1)+dtime
! ! WRITE(*,*) KSTEP
! !enddo
! !
! !!ACTIVE CONTRACTION
! !do KSTEP=1001,2000
! !!gamma=0.1942d0
! !CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !s1=stress(1)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !write(50,*) time(1),gamma,s1,stress(1)
! !!gamma=gamma+0.0002d0
! !!
! !time(1)=time(1)+dtime
! ! WRITE(*,*) KSTEP
! !enddo
! !!STRETCHING WITH ACTIVE CONTRACTION
! !do KSTEP=1001,2000
! !!gamma=0.1942d0
! !gamma=gamma+0.00001d0
! !DFGRD1(1,1)=1.d0+gamma
! !DFGRD1(2,2)=  1.0D0/sqrt(DFGRD1(1,1))
! !DFGRD1(3,3)=DFGRD1(2,2)
! !CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !s1=stress(1)
! !CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! !DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! !NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !!
! !write(50,*) time(1),gamma,s1,stress(1)
! !!gamma=gamma+0.0002d0
! !!
! !time(1)=time(1)+dtime
! ! WRITE(*,*) KSTEP
! !enddo
! !################################################################################################!
close(150)
!
END PROGRAM
