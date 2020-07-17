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
PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 19, NDI=3, NSHR=3)
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
PROPS(5)=0.1161d0
!R0
PROPS(6)=0.0839d0
!mu0
PROPS(7)=38600.0d0
!PROPS(7)=111111111111111111111111111138600.0d0
!beta
PROPS(8)=0.438d0
!PROPS(8)=0.5d0
!B0 = tk*lp*k0
PROPS(9)=294.d0*16.d0*1.38d-5
!lambda0.
PROPS(10)=1.0039d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!NON-AFFINE NETWORK PARAMS
!n - isotropic filaments per unit volume
PROPS(11)=129.19d0
!P....
PROPS(12)=13.2258d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PROPS(16)=1.5e1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(13)=0
!tau1 %
PROPS(14)=2.0d0
!teta1
PROPS(15)=0.835d0
!tau2 %
PROPS(16)=0.007d0
!teta2
PROPS(17)=1.6d0
!tau3 %
PROPS(18)=15.d0
!teta3
PROPS(19)=1.4d0
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
! !    
!  DFGRD1(1,1)=  1.0D0
!  DFGRD1(1,2)=  0.0D0
!  DFGRD1(1,3)=  0.0d0
!  DFGRD1(2,1)=  0.0d0
!  DFGRD1(2,2)=  1.0D0
!  DFGRD1(2,3)=  0.0d0
!  DFGRD1(3,1)=  0.0d0
!  DFGRD1(3,2)=  0.0d0
!  DFGRD1(3,3)=  1.0D0
 
! !################################################################################################!
! ! AMPLITUDE SWWEEP
!   AMPMIN=0.00d0
!   DAMP=0.01d0
!   AMP=0.125
! !  DAMP=0.125d0
!   FREQ=0.5d0
!   NTEST=18
! !  NTEST=1
!   NCYCLES=1
!   NSTEPS=100      
! !  CALL AMPSWEEP(AMPMIN,DAMP,FREQ,NTEST,NCYCLES,NSTEPS,PROPS)
! !################################################################################################!
! !################################################################################################!
! ! AMP SWEEP WORKING
! !  NCYCLES=20
! !  NSTEPS=200
! !  NPOINTS=18
! !  ALLOCATE(curvepoints(NPOINTS,3))
! !  do i=1,npoints
! !     curvepoints(i,1)=AMPMIN+I*ONE*DAMP
! !  enddo
! !
! !      DO I1=1,NPOINTS
! !!C
! !        AMP=curvepoints(I1,1)
!         TIME(1)=ZERO
!         PI=FOUR*ATAN(ONE)
!         F=FREQ*TWO*PI
! !!        
! !         write(filename,fmt='(i0,a)')I1,'.out'
! !        un=un+2;
! !        open(unit=un,file=filename,status='UNKNOWN')
! !!        
!         cycletime=one/freq
! !       
!         dtime=cycletime/nsteps
! !        
!         maxstress=0.0d0
!         maxstresstime=0.0d0
!         minstress=0.0d0
!         minstresstime=0.0d0
!         maxsgama=0.d0
!         maxgamatime=0.d0

!         DO KK=1,NCYCLES
                
!            DO KSTEP=1,NSTEPS

!               GAMMA=AMP*DSIN(F*TIME(1))
!               DFGRD1(1,2)=GAMMA
! !C>---------------------------------------------------------------------
!       CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!       RPL,DDSDDT,DRPLDE,DRPLDT,&
!       STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
!       NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
!       CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
       
!       if(kk.eq.ncycles)THEN  !pesquisa apenas no ultimo ciclo    
! !      write(50,*) time(1),gamma,stress(4)
!       if(stress(4).gt.maxstress) THEN
!       maxstress=stress(4)
!       maxstresstime=time(1)
!       endif
      
!       if(stress(4).lt.minstress) THEN
!       minstress=stress(4)
!       minstresstime=time(1)
!       endif
      
!       if(gamma.gt.maxgamma) THEN
!       maxsgama=gamma
!       maxgamatime=time(1)
!       endif
!       endif

!       time(1)=time(1)+dtime
!      write(*,*) time(1),gamma, stress(4)
!      write(150,*) time(1),gamma,stress(4)
!            ENDDO !KSTEP
           
!         !   write(*,*) KK,time(1)
!         ENDDO !NCYCLES
        
! !        stressamp=maxstress-minstress
! !        delta=-maxstresstime+maxgamatime
! !        
! !        delta=delta*TWO*PI/cycletime
! !        
! !!        storagemodulus=(-maxstress+minstress)*0.5d0/AMP*dcos(delta)
! !        storagemodulus=(maxstress/amp)*dcos(delta-pi/two)
! !!           lossmodulus=(-maxstress+minstress)*0.5d0/AMP*dsin(delta)
! !           lossmodulus=(maxstress/amp)*dsin(delta-pi/two)
! !
! !!              write(*,*) time(1),gamma,stress(4)
! !       ! close(unit=un)
! !!       write(*,*) KK,maxstress,minstress,delta
! !       write(*,*) AMP,storagemodulus,lossmodulus
! !       write(150,*) AMP,storagemodulus,lossmodulus
! !
! !      ENDDO
      
!   close(150)
!   CALL SYSTEM('gnuplot -p data_xy.plt')        
! stop              
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
!  AMP=0.005d0
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

! DFGRD1(1,1)=  1.0D0
! DFGRD1(1,2)=  0.0D0
! DFGRD1(1,3)=  0.0d0
! DFGRD1(2,1)=  0.0d0
! DFGRD1(2,2)=  1.0D0
! DFGRD1(2,3)=  0.0d0
! DFGRD1(3,1)=  0.0d0
! DFGRD1(3,2)=  0.0d0
! DFGRD1(3,3)=  1.0D0
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
!period=1.d0/freq
!nsteps=100
!dtime=period/nsteps
!!dtime=0.2
!amp=0.12d0
!angfreq=2.d0*4.d0*atan(1.d0)*freq
!!gamma=amp*sin(angfreq*time(1))
! write(150,*) time(1),gamma,stress(4)
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
! write(150,*) time(1),gamma,stress(4)
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
! write(150,*) time(1),gamma,stress(4)
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
!write(150,*) time(1),gamma,stress(4)
!enddo
!write(*,*) k1,time(1)
!enddo
!
!
!close(150)
!  CALL SYSTEM('gnuplot -p data_xy.plt')
!################################################################################################!
!0.740469148776867d0, 0.672090350852965d0, 0.d0

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
gamma=gamma+0.001942d0
!gamma=gamma+0.00025d0
DFGRD1(1,2)=gamma
!!
!!s1=stress(4)
CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!!
write(150,*) time(1),gamma,stress(4)
!gamma=gamma+0.0002d0
time(1)=time(1)+dtime

enddo

CALL SYSTEM('gnuplot -p data_xy.plt')
!
!################################################################################################! 
close(150)
!
END PROGRAM
