!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST

use,intrinsic :: ISO_Fortran_env
 
INCLUDE 'ABA_PARAM.INC'
INCLUDE 'PARAM_UMAT.INC'

      COMMON /KFIL/MF0
      COMMON /KFIL/RW
      COMMON /KDATA/PTS
      COMMON /KPAR/MATPAR
      DOUBLE PRECISION MF0(NWP,3),RW(NWP),PTS(NPTS,2),MATPAR(NNPAR)
      double precision aux2(npts),aux3(npts)
      PARAMETER(NTENS = 6, NSTATEV = 2, NPROPS = 6, NDI=3, NSHR=3)

CHARACTER*8 CMNAME,frmt
DIMENSION STRESS(NTENS),STATEV(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
       DOUBLE PRECISION GAMMA,mean,meandistsqd,errorsqd,r2
       logical :: exist
!
!read data points
CALL UEXTERNALDBB(0,0,0.d0,0.d0,0,0)
!
 DFGRD1(1,1)=   1.0D0
 DFGRD1(1,3)=  0.0E+000
! !DFGRD1(1,2)=  GAMMA
 DFGRD1(2,1)=  0.0D0
 DFGRD1(2,2)=  1.0D0
 DFGRD1(2,3)=  0.0E+000
 DFGRD1(3,1)=  0.0E+000
 DFGRD1(3,2)=  0.0E+000
 DFGRD1(3,3)=   1.0D0
 !
 !read fibres  directions
call UEXTERNALDB(0,0,0.0,0.0,0,0)
! run genetic algorithm for fitting
call gafortran
!write best guess of mat parameters
call UEXTERNALD(0,0,0.0,0.0,0,0)
!
        !  DFGRD1(1,1)=   1.0D0
          DFGRD1(1,3)=  0.D0
          DFGRD1(1,2)=  0.D0
          DFGRD1(2,1)=  0.0D0
          !DFGRD1(2,2)=  1.0D0
          DFGRD1(2,3)=  0.D0
          DFGRD1(3,1)=  0.D0
          DFGRD1(3,2)=  0.D0
         ! DFGRD1(3,3)=  1.0D0
          NOEL=1
          NPT=8
     ! FIXED MATERIAL PROPERTIES
     
!      C10 = PROPS(2)
!      C01 = PROPS(3)
!      D1  = PROPS(1)
!      C3  = PROPS(4)
!      C4  = PROPS(5)
!      C5  = PROPS(6)
!      C6  = PROPS(7)
!      KAPPA = PROPS(8)
!C
!      BETAM  = PROPS(9)
!      SEMMIN = PROPS(10)
!      SEMMAX = PROPS(11)
!C
!      BETAF  = PROPS(12)
!      SEFMIN = PROPS(13)
!      SEFMAX = PROPS(14)
         ! KBULK
         PROPS(1)=0.000001d0
         ! C10=
         PROPS(2)=1.00d0
         ! C01
         PROPS(3)=0.00d0
         !k1
         PROPS(4)=0.0d0
         !k2
         PROPS(5)=.0001d0
         !kdisp
         PROPS(6)=0.0d0     
!
         ! FREE MATERIAL PROPERTIES
         PROPS(2)=matpar(1)
         PROPS(4)=matpar(2)
         PROPS(5)=matpar(3)
         PROPS(6)=matpar(4)
         !  11   continue
          gamma=0.d0
          OPEN (UNIT=50, FILE=DIRCC, STATUS='UNKNOWN')
          rewind 50
          aux=0.d0
        do  j1=1,npts
           gamma=pts(j1,1)
           sexp=pts(j1,2)
           aux=aux+sexp
           DFGRD1(1,1)=gamma
           DFGRD1(2,2)=ONE/SQRT(gamma)
           DFGRD1(3,3)=ONE/SQRT(gamma)
         CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
         DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
         NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
          write(50,*) gamma,sexp,stress(1)
          aux2(j1)=sexp-stress(1)
          aux2(j1)=aux2(j1)*aux2(j1)
        enddo
        close(50)
        mean=aux/npts
        errorsqd=0.d0
        meandistsqd=0.d0
        do j1=1,npts
        errorsqd=errorsqd+aux2(j1)
          aux3(j1)=pts(j1,2)-mean
          aux3(j1)=aux3(j1)*aux3(j1)
        meandistsqd=meandistsqd+aux3(j1)
        enddo
        r2=1.d0-(errorsqd/meandistsqd)
!        
        write (6,1191) 
!        
        write(*,*) 'R^2=',r2
        OPEN (UNIT=59, FILE='plot_r2.out', STATUS='UNKNOWN')
        rewind 59
        write(59,*) r2
        write(59,*) matpar(1), matpar(2), matpar(3), matpar(4), matpar(3)/(4.d0*matpar(1)*matpar(4))
        CLOSE(59)
!       
        ! CALL SYSTEM('gnuplot -p'// ' ' // DIRDD)
        ! CALL SYSTEM('gnuplot -p'// ' ' // DIRD)
        write(*,*) 'PLOTS generated at:'
        write(*,*)
        write(*,*) dirdd
        write(*,*) dird
!
 1191 format(//'#################  FITTING QUALITY  #################')
 990 format (' ',F5.3)
 991 format (' ',4E12.5E2)
END PROGRAM
