      SUBROUTINE AMPSWEEP(AMPMIN,DAMP,FREQ,NTEST,NCYCLES,NSTEPS,PROPS)
C      
      IMPLICIT NONE
      INCLUDE 'PARAM_UMAT.INC'
C      
C     FILAMENTS DIRECTION
      COMMON /KFIL/MF0
C     FILAMENTS WEIGHT      
      COMMON /KFILR/RW
C     CHEMICAL DYNAMICS MATRIX   
      COMMON /KFILF/FRAC0     
      COMMON /KFILK/KCH
C     PREFERED DIRETION
      COMMON /KFILP/PREFDIR
C
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),PREFDIR(NELEM,4),FRAC0(NCH),
     1                        KCH(7)
      PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 23, NDI=3)
      PARAMETER(NOEL = 1, NPT = 8,NSHR=3)
C
      CHARACTER*8 CMNAME
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      DOUBLE PRECISION SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,
     1                 DTEMP,PNEWDT,CELENT
C
      INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT,
     1        LAYER, KSPT, KSTEP, KINC   
      INTEGER NTEST,NCYCLES,NSTEPS,I1,KK
      DOUBLE PRECISION AMPMIN,AMPMAX,FREQ,F,AMP,GAMMA,DAMP,TT
      CHARACTER*8 filename
      integer n,un
C      ###########################################################
      F=FREQ*TWO*FOUR*ATAN(ONE)
      dtime=0.02d0
      dtime=one/(freq*(nsteps))
C
          DFGRD1(1,1)=   1.0D0
         !DFGRD1(1,2)=  0.3d0
         DFGRD1(1,3)=  0.D0
         DFGRD1(2,1)=  0.0D0
         DFGRD1(2,2)=  1.0D0
         DFGRD1(2,3)=  0.D0
         DFGRD1(3,1)=  0.D0
         DFGRD1(3,2)=  0.D0
         DFGRD1(3,3)=  1.0D0 
C
      DO I1=1,NTEST
C      
        AMP=AMPMIN+(I1)*ONE*DAMP
        TIME(1)=ZERO
        tt=0.d0
        write(filename,fmt='(i0,a)')I1,'.out'
        un=un+2;
        open(unit=un,file=filename,status='UNKNOWN')
        write(*,*) 'amp=',amp
!        write(un,*) 'freq=',freq
        OPEN (UNIT=50, FILE='data.out', STATUS='UNKNOWN')
        rewind(50)

C        
        DO KK=1,NCYCLES
        
           DO KSTEP=1,NSTEPS
              GAMMA=AMP*DSIN(F*tt)
              DFGRD1(1,2)=GAMMA
              
C>---------------------------------------------------------------------
      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C       write(*,*) stress
      write(50,*) time(1),gamma,stress(4)
      time(1)=time(1)+dtime
      tt=tt+dtime
      write(un,*) time(1),gamma,stress(4)
           ENDDO
           write(*,*) KK,time(1)
        ENDDO
        close(unit=un)
      ENDDO
      close(50)
        CALL SYSTEM('gnuplot -p data_xy.plt')
C
      RETURN
C
      END SUBROUTINE
