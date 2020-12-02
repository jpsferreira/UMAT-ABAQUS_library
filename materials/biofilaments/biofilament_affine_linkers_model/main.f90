!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'ABA_PARAM.INC'
INCLUDE 'PARAM_UMAT.INC'

      COMMON /KFIL/MF0
      COMMON /KFILR/RW
      COMMON /KFILP/PREFDIR
      COMMON /KFILC/CHEM
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),PREFDIR(NELEM,4),CHEM(2,7)
PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 21, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),STATEVP(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

real*8,ALLOCATABLE::curvepoints(:,:)
real*8::maxstresstime,maxstress,maxgamatime,maxgama,minstresstime,minstress
real*8:: delta,lossmodulus

 CHARACTER*8 filename
integer un


i=1.0d0
j=1.0d0
DO i=1,NTENS
    DO j=1,NTENS
        DDSDDE(i,j)=0.0D0
    ENDDO
    STRESS(i)=0.0D0
ENDDO



!open(unit=20,file='output.txt',status='replace',action="write", &
!            position="append",iostat=ierr,access='sequential')

!
! DEFORMATION GRADIENT
 !PURE SHEAR
! DFGRD1(1,1)=   1.2D0
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.D0
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.D0/DFGRD1(1,1)
!
! 
!  !ARBITRARY
! DFGRD1(1,1)=   1.10000000000000
! DFGRD1(1,2)=  0.200000000000000E+000
! DFGRD1(1,3)=  0.300000000000000E+000
! DFGRD1(2,1)=  0.400000000000000E+000
! DFGRD1(2,2)=   2.50000000000000
! DFGRD1(2,3)=  0.600000000000000E+000
! DFGRD1(3,1)=  0.700000000000000E+000
! DFGRD1(3,2)=  0.800000000000000E+000
! DFGRD1(3,3)=   3.90000000000000

!
!      !DILATATION
! DFGRD1(1,1)=   1.20000000000000
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.1d0
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.1d0
!
!   !EXTENSION
! DFGRD1(1,1)=   1.1d0
! DFGRD1(1,2)=  0.0E+000
! DFGRD1(1,3)=  0.0E+000
! DFGRD1(2,1)=  0.0E+000
! DFGRD1(2,2)=  1.d0/SQRT(DFGRD1(1,1))
! DFGRD1(2,3)=  0.0E+000
! DFGRD1(3,1)=  0.0E+000
! DFGRD1(3,2)=  0.0E+000
! DFGRD1(3,3)=   1.d0/SQRT(DFGRD1(1,1))
 
  DFGRD1(1,1)=   1.0D0
 !DFGRD1(1,2)=  0.3d0
 DFGRD1(1,3)=  0.0E+000
 DFGRD1(2,1)=  0.0D0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0E+000
 DFGRD1(3,1)=  0.0E+000
 DFGRD1(3,2)=  0.0E+000
 DFGRD1(3,3)=   1.0D0 
!!! 
call UEXTERNALDB(0,0,0.0,0.0,0,0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ISOTROPIC MATRIX PARAMS
! C10=
PROPS(2)=0.00d0
! C01
PROPS(3)=0.00d0
!PHI....
PROPS(4)=1.0000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SINGLE FILAMENT PARAMS
!L
PROPS(5)=1.96d0
!R0F
PROPS(6)=1.63d0
!R0C
PROPS(7)=0.014d0
!ETAC
PROPS(8)=0.3333d0
!mu0
!PROPS(7)=38600.0d0
PROPS(9)=111111111111111111111111111138600.0d0
!beta
PROPS(10)=0.5d0
!PROPS(8)=0.5d0
!B0 = tk*lp*k0
PROPS(11)=294.d0*16.d0*1.38d-5
!lambda0.
PROPS(12)=1.00d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!AFFINE NETWORK PARAMS
!n - isotropic filaments per unit volume
PROPS(13)=7.66D0
!B....
PROPS(14)=0.001d0
!PROPS(16)=1.5e1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(15)=0
!tau1 %
PROPS(16)=0.2d0
!teta1
PROPS(17)=10.0d0
!tau2 %
PROPS(18)=1.2d0
!teta2
PROPS(19)=7.0d0
!tau3 %
PROPS(20)=12.d0
!teta3
PROPS(21)=2.0d0
!
!!tau1  % for low freq
!PROPS(18)=12.0d0  
!!teta1
!PROPS(19)=0.15d0
!!tau2  %for med freq 
!PROPS(20)=0.9d0
!!teta2
!PROPS(21)=1.6d0
!!tau3 %for high freq 
!PROPS(22)=0.05d0
!!teta3
!PROPS(23)=3.6d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STATEV=0.D0       
STATEVP=0.D0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
erf=0.d0
RHO=0.D0
!
OPEN (UNIT=150, FILE='data.out', STATUS='UNKNOWN')
        rewind(150)
!   
!    
 DFGRD1(1,1)=  1.0D0
 DFGRD1(1,2)=  0.0D0
 DFGRD1(1,3)=  0.0d0
 DFGRD1(2,1)=  0.0d0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0d0
 DFGRD1(3,1)=  0.0d0
 DFGRD1(3,2)=  0.0d0
 DFGRD1(3,3)=  1.0D0
 
!################################################################################################!
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
! AMP SWEEP WORKING
!  NCYCLES=20
!  NSTEPS=200
!  NPOINTS=18
!  ALLOCATE(curvepoints(NPOINTS,3))
!  do i=1,npoints
!     curvepoints(i,1)=AMPMIN+I*ONE*DAMP
!  enddo
!
!      DO I1=1,NPOINTS
!!C
!        AMP=curvepoints(I1,1)
!        TIME(1)=ZERO
!        PI=FOUR*ATAN(ONE)
!        F=FREQ*TWO*PI
!!        
!         write(filename,fmt='(i0,a)')I1,'.out'
!        un=un+2;
!        open(unit=un,file=filename,status='UNKNOWN')
!!        
!        cycletime=one/freq
!!       
!        dtime=cycletime/nsteps
!!        
!        maxstress=0.0d0
!        maxstresstime=0.0d0
!        minstress=0.0d0
!        minstresstime=0.0d0
!        maxsgama=0.d0
!        maxgamatime=0.d0
!
!        DO KK=1,NCYCLES
!                
!           DO KSTEP=1,NSTEPS
!
!              GAMMA=AMP*DSIN(F*TIME(1))
!              DFGRD1(1,2)=GAMMA
!!C>---------------------------------------------------------------------
!      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!       
!      if(kk.eq.ncycles)THEN  !pesquisa apenas no ultimo ciclo    
!!      write(50,*) time(1),gamma,stress(4)
!      if(stress(4).gt.maxstress) THEN
!      maxstress=stress(4)
!      maxstresstime=time(1)
!      endif
!      
!      if(stress(4).lt.minstress) THEN
!      minstress=stress(4)
!      minstresstime=time(1)
!      endif
!      
!      if(gamma.gt.maxgamma) THEN
!      maxsgama=gamma
!      maxgamatime=time(1)
!      endif
!      endif
!
!      time(1)=time(1)+dtime
!      write(un,*) time(1),gamma,stress(4)
!           ENDDO !KSTEP
!           
!        !   write(*,*) KK,time(1)
!        ENDDO !NCYCLES
!        
!        stressamp=maxstress-minstress
!        delta=-maxstresstime+maxgamatime
!        
!        delta=delta*TWO*PI/cycletime
!        
!!        storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
!        storagemodulus=(maxstress/amp)*dcos(delta-pi/two)
!!           lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
!           lossmodulus=(maxstress/amp)*dsin(delta-pi/two)
!
!!              write(*,*) time(1),gamma,stress(4)
!       ! close(unit=un)
!!       write(*,*) KK,maxstress,minstress,delta
!       write(*,*) AMP,storagemodulus,lossmodulus
!       write(150,*) AMP,storagemodulus,lossmodulus
!
!      ENDDO
!      
!  close(150)
!  CALL SYSTEM('gnuplot -p data_xyg.plt')        
!stop              
!################################################################################################!
!################################################################################################!      
!################################################################################################!
!################################################################################################!
! FREQ SWEEP
!  FREQMIN=0.00d0
!  DFREQ=0.01d0
!  AMP=0.005d0
!  NTEST=61
!  NCYCLES=3
!  NSTEPS=200
!!  CALL FREQSWEEP(FREQMIN,DFREQ,AMP,NTEST,NCYCLES,NSTEPS,PROPS)
!
!!################################################################################################!
!!################################################################################################!
!! FREQ SWEEP WORKING
!  FREQMIN=1.00d0
!  DFREQ=0.0d0
!  AMP=0.005d0
!  NTEST=1
!  NCYCLES=3
!  NSTEPS=400
!  TIME=0.0d0
!  FREQMIN=0.00d0
!  DFREQ=0.01d0
!  AMP=0.03d0
!!  NTEST=30
!  NCYCLES=10
!  NSTEPS=400
!  NPOINTS=71
!  ALLOCATE(curvepoints(NPOINTS,3))
!  counter=-0.05d0
!  
!  ! o ciclo seguinte vai gerar 61 pontos uniformemente
!  ! distribuidos numa escala logaritmica.
!  do i=1,npoints
!  counter=counter+0.05d0
!  curvepoints(i,1)=0.003d0*10.0d0**counter
!  enddo
!
!      DO I1=1,NPOINTS
!!C
!        freq=curvepoints(I1,1)
!        TIME(1)=ZERO
!        PI=FOUR*ATAN(ONE)
!        F=FREQ*two*PI
!!C        dtime=0.02d0
!         write(filename,fmt='(i0,a)')I1,'.out'
!        un=un+2;
!        open(unit=un,file=filename,status='UNKNOWN')
!        
!        cycletime=one/freq
!        
!        dtime=cycletime/nsteps
!!C
!        maxstress=0.0d0
!        maxstresstime=0.0d0
!        minstress=0.0d0
!        minstresstime=0.0d0
!
!        DO KK=1,NCYCLES
!                
!           DO KSTEP=1,NSTEPS
!
!              GAMMA=AMP*DSIN(F*TIME(1))
!              DFGRD1(1,2)=GAMMA
!!C>---------------------------------------------------------------------
!      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!      RPL,DDSDDT,DRPLDE,DRPLDT,&
!      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!       
!      if(kk.eq.ncycles)THEN  !pesquisa apenas no ultimo ciclo    
!!      write(50,*) time(1),gamma,stress(4)
!      if(stress(4).gt.maxstress) THEN
!      maxstress=stress(4)
!      maxstresstime=time(1)
!      endif
!      
!      if(stress(4).lt.minstress) THEN
!      minstress=stress(4)
!      minstresstime=time(1)
!      endif
!      
!      if(gamma.gt.maxgamma) THEN
!      maxsgama=gamma
!      maxgamatime=time(1)
!      endif
!      endif
!
!      time(1)=time(1)+dtime
!      write(un,*) time(1),gamma,stress(4)
!           ENDDO !KSTEP
!           
!        !   write(*,*) KK,time(1)
!        ENDDO !NCYCLES
!        
!        stressamp=maxstress-minstress
!        delta=-maxstresstime+maxgamatime
!        
!        delta=delta*two*PI/cycletime
!        
!!        storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
!        storagemodulus=(maxstress/amp)*dcos(delta-pi/two)
!!           lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
!           lossmodulus=(maxstress/amp)*dsin(delta-pi/two)
!           tandelta=lossmodulus/storagemodulus
!
!!              write(*,*) time(1),gamma,stress(4)
!       ! close(unit=un)
!       write(*,*) freq*(two*pi),storagemodulus,lossmodulus,tandelta
!       write(150,*) freq*(two*pi),storagemodulus,tandelta,lossmodulus
!
!      ENDDO
!      
!  close(150)
!  CALL SYSTEM('gnuplot -p data_xyg.plt')      
!################################################################################################!
!################################################################################################!      
! SHEAR CYCLING TEST
!statev=0.d0
!gamma=0.d0
!time(1)=0.d0
!dtime=0.02d0
!kstep=1
!!
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!freq=0.05d0
!dtime=0.2
!period=1.d0/freq
!amp=0.1d0
!angfreq=2.d0*4.d0*atan(1.d0)*freq
!!gamma=amp*sin(angfreq*time(1))
! write(50,*) time(1),gamma,stress(4)
!do k1=1,20 !cycles loop
!write(*,*) k1,time(1)
!do kstep=1,25
!!gamma=gamma+0.0012d0
!time(1)=time(1)+dtime
!gamma=amp*sin(angfreq*time(1))
!DFGRD1(1,2)=  gamma
!!write(*,*) kstep,gamma
!!write(*,*)
! CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! !write(*,*) statev(10:13)
! write(50,*) time(1),gamma,stress(4)
!enddo
!do kstep=26,75
!time(1)=time(1)+dtime
!gamma=amp*sin(angfreq*time(1))
!DFGRD1(1,2)=  gamma
!!write(*,*) kstep,gamma
!!write(*,*)
! CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
! DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
! NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
! ! write(*,*) statev(10:13)
! write(50,*) time(1),gamma,stress(4)
!enddo
!do kstep=76,100
!time(1)=time(1)+dtime
!gamma=amp*sin(angfreq*time(1))
!DFGRD1(1,2)=  gamma
!!write(*,*) kstep,gamma
!!write(*,*)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!! write(*,*) statev(10:13)
!write(50,*) time(1),gamma,stress(4)
!enddo
!write(*,*) k1,time(1)
!enddo
!
!
!close(50)
!  CALL SYSTEM('gnuplot -p data_xy.plt')
!################################################################################################!
!0.740469148776867d0, 0.672090350852965d0, 0.d0
!! CREEP RELAXATION TEST
!time(1)=0.d0
!dtime=0.00025d0
!!
!!dtime=0.1d0
!amp=0.125d0
!nsteps=100
!!
!!time(1)=0.d0
!!dtime=0.02d0
!freq=2.5d0
!dtime=1/(freq*4*nsteps)
!!period=1.d0/freq
!amp=0.125d0
!angfreq=2.d0*4.d0*atan(1.d0)*freq
!!gamma=0.125d0
!do kstep=1,nsteps
!gamma=amp*sin(angfreq*time(1))
!DFGRD1(1,2)=gamma
!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s1=stress(4)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!write(50,*) time(1),gamma,s1,stress(4)
!time(1)=time(1)+dtime
!enddo
!!
!write(*,*) 
!DFGRD1(1,2)=gamma
!dtime=0.1d0
!do kstep=nsteps+1,10*nsteps+nsteps-5
!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s1=stress(4)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
!write(50,*) time(1),gamma,s1,stress(4)
!time(1)=time(1)+dtime
!enddo
!close(50)
!  CALL SYSTEM('gnuplot -p data_xy.plt')

!################################################################################################!
!write(50,*) time(1),gamma,s1,stress(4)
!! SHEAR MONOTONIC LOAD TEST

 DFGRD1(1,1)=  1.0D0
 DFGRD1(1,2)=  0.0D0
 DFGRD1(1,3)=  0.0d0
 DFGRD1(2,1)=  0.0d0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0d0
 DFGRD1(3,1)=  0.0d0
 DFGRD1(3,2)=  0.0d0
 DFGRD1(3,3)=  1.0D0
 
gamma=0.d0
time(1)=0.d0
dtime=0.01d0
do kstep=1,100
!gamma=0.1942d0
gamma0=gamma
stress0=stress(4)
gamma=gamma+0.003d0
!gamma=gamma+0.00025d0
DFGRD1(1,2)=gamma
!!
!!s1=stress(4)
CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
write(150,*) time(1),gamma,stress(4)!,(stress(4)-stress0)/(gamma-gamma0)
!gamma=gamma+0.0002d0
time(1)=time(1)+dtime

enddo
!
!time(1)=0.d0
!dtime=0.2d0
!gamma=0.19428d0
!!gamma=0.0d0
!!time(1)=time(1)+dtime
!DFGRD1(1,2)=gamma
!do kstep=1,100
!!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!!s1=stress(4)
!!
!!!!!!!!!!!!!!!!!!
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s2=stress(4)
!write(150,*) time(1),gamma,STRESS(4)
!!gamma=gamma+0.019428d0
!!gamma=0.19428d0
!time(1)=time(1)+dtime
!!DFGRD1(1,2)=gamma
!enddo
  CALL SYSTEM('gnuplot -p data_xy.plt')
!################################################################################################!
!
! !DFGRD1(1,1)=  1.0D0
! DFGRD1(1,2)=  0.0D0
! DFGRD1(1,3)=  0.0d0
! DFGRD1(2,1)=  0.0d0
! DFGRD1(2,2)=  1.0D0
! DFGRD1(2,3)=  0.0d0
! DFGRD1(3,1)=  0.0d0
! DFGRD1(3,2)=  0.0d0
! DFGRD1(3,3)=  1.0D0
! 
!gamma=0.d0
!dtime=0.1
!
!!ISOMETRIC ACTIVATION 
!do kstep=1,1000
!!gamma=0.1942d0
!!gamma=gamma+0.0001942d0
!DFGRD1(1,1)=1.d0
!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s1=stress(1)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!write(50,*) time(1),gamma,s1,stress(1)
!!gamma=gamma+0.0002d0
!!
!time(1)=time(1)+dtime
! WRITE(*,*) KSTEP
!enddo
!
!!ACTIVE CONTRACTION
!do KSTEP=1001,2000
!!gamma=0.1942d0
!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s1=stress(1)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!write(50,*) time(1),gamma,s1,stress(1)
!!gamma=gamma+0.0002d0
!!
!time(1)=time(1)+dtime
! WRITE(*,*) KSTEP
!enddo
!!STRETCHING WITH ACTIVE CONTRACTION
!do KSTEP=1001,2000
!!gamma=0.1942d0
!gamma=gamma+0.00001d0
!DFGRD1(1,1)=1.d0+gamma
!DFGRD1(2,2)=  1.0D0/sqrt(DFGRD1(1,1))
!DFGRD1(3,3)=DFGRD1(2,2)
!CALL UMATP(STRESS,STATEVP,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!s1=stress(1)
!CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
!DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
!NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
!write(50,*) time(1),gamma,s1,stress(1)
!!gamma=gamma+0.0002d0
!!
!time(1)=time(1)+dtime
! WRITE(*,*) KSTEP
!enddo
!################################################################################################! 
close(150)
!
END PROGRAM
