      INTEGER KNODES, NELEM, KNPE, KNPG, NUEL, NIPP, NWP
      PARAMETER (KNODES=8,NELEM=1,KNPE=8, KNPG=8,NUEL=2,NIPP=8,NWP=21)
      DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
      PARAMETER (ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
      PARAMETER (THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
      CHARACTER(256) DIR1,DIR2, DIRA, DIRC, DIRCC, DIRD, DIRDD
      INTEGER NPTS,NNPAR
      PARAMETER (NPTS=21,NNPAR=4) 
      PARAMETER (DIR1='fibers.inp',DIR2='fibers.inp')
      PARAMETER (DIRA='uterus_parous.csv')
      PARAMETER (DIRC='par_uterus_parous.out')
      PARAMETER (DIRCC='plot_uterus_parous.out')
      PARAMETER (DIRD='fit_uterus_parous.plt')
      PARAMETER (DIRDD= 'ga_uterus_parous.plt') 
