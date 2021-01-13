c
c  Any users new to the GA world are encouraged to read David Goldberg's
c  "Genetic Algorithms in Search, Optimization and Machine Learning,"
c  Addison-Wesley, 1989.
c
c  Other associated files are:  ga.inp
c                               ga.out
c                               ga.restart
c                               params.inc
c                               ReadMe
c                               ga2.inp (w/ different namelist identifier)
c
c  I have provided a sample subroutine "func", but ultimately
c  the user must supply this subroutine "func" which should be your
c  cost function.  You should be able to run the code with the
c  sample subroutine "func" and the provided ga.inp file and obtain
c  the optimal function value of 1.0000 at generation 187 with the
c  uniform crossover micro-GA enabled (this is 935 function evaluations).
c
c  The code is presently set for a maximum population size of 200,
c  30 chromosomes (binary bits) and 8 parameters.  These values can be
c  changed in params.inc as appropriate for your problem.  Correspondingly
c  you will have to change a few 'write' and 'format' statements if you
c  change nchrome and/or nparam.  In particular, if you change nchrome
c  and/or nparam, then you should change the 'format' statement numbers
c  1050, 1075, 1275, and 1500 (see ReadMe file).
c
c  Please feel free to contact me with questions, comments, or errors
c  (hopefully none of latter).
c
c  Disclaimer:  this program is not guaranteed to be free of error
c  (although it is believed to be free of error), therefore it should
c  not be relied on for solving problems where an error could result in
c  injury or loss.  If this code is used for such solutions, it is
c  entirely at the user's risk and the author disclaims all liability.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'params.inc'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension ibest(nchrmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
      dimension geni(1000000),genavg(1000000),genmax(1000000)
c      real*4 cpu,cpu0,cpu1,tarray(2)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Input variable definitions:
c
c  icreep   = 0 for no creep mutations
c           = 1 for creep mutations; creep mutations are recommended.
c  idum     The initial random number seed for the GA run.  Must equal
c           a negative integer, e.g. idum=-1000.
c  ielite   = 0 for no elitism (best individual not necessarily
c               replicated from one generation to the next).
c           = 1 for elitism to be invoked (best individual replicated
c               into next generation); elitism is recommended.
c  iend         = 0 for normal GA run (this is standard).
c           = number of last population member to be looked at in a set
c             of individuals.  Setting iend-0 is only used for debugging
c             purposes and is commonly used in conjunction with iskip.
c  iniche   = 0 for no niching
c           = 1 for niching; niching is recommended.
c  irestrt  = 0 for a new GA run, or for a single function evaluation
c           = 1 for a restart continuation of a GA run.
c  iskip    = 0 for normal GA run (this is standard).
c           = number in population to look at a specific individual or
c             set of individuals.  Setting iskip-0 is only used for
c             debugging purposes.
c  itourny  No longer used.  The GA is presently set up for only
c           tournament selection.
c  itwopoint= 1 for two-point crossover, only if iunifrm=0
c  iunifrm  = 0 for single-point crossover
c           = 1 for uniform crossover; uniform crossover is recommended.
c  kountmx  = the maximum value of kount before a new restart file is
c             written; presently set to write every fifth generation.
c             Increasing this value will reduce I/O time requirements
c             and reduce wear and tear on your storage device
c  maxgen   The maximum number of generations to run by the GA.
c           For a single function evaluation, set equal to 1.
c  microga  = 0 for normal conventional GA operation (not activated right now)
c           = 1 for micro-GA operation (this will automatically reset
c             some of the other input flags).  I recommend using
c             npopsiz=5 when microga=1.
c  nchild   = 1 for one child per pair of parents (this is what I
c               typically use).
c           = 2 for two children per pair of parents (2 is more common
c               in GA work).
c  nichflg  = array of 1/0 flags for whether or not niching occurs on
c             a particular parameter.  Set to 0 for no niching on
c             a parameter, set to 1 for niching to operate on parameter.
c             The default value is 1, but the implementation of niching
c             is still controlled by the flag iniche.
c  nowrite  = 0 to write detailed mutation and parameter adjustments
c           = 1 to not write detailed mutation and parameter adjustments
c  nparam   Number of parameters (groups of bits) of each individual.
c           Make sure that nparam matches the number of values in the
c           parmin, parmax and nposibl input arrays.
c  npopsiz  The population size of a GA run (typically 100 works well).
c           For a single calculation, set equal to 1.
c  nposibl  = array of integer number of possibilities per parameter.
c             For optimal code efficiency set nposibl=2**n, i.e. 2, 4,
c             8, 16, 32, 64, etc.
c  parmax   = array of the maximum allowed values of the parameters
c  parmin   = array of the minimum allowed values of the parameters
c  pcreep   The creep mutation probability.  Typically set this
c           = (nchrome/nparam)/npopsiz.
c  pcross   The crossover probability.  For single-point crossover, a
c           value of 0.6 or 0.7 is recommended.  For uniform crossover,
c           a value of 0.5 is suggested.
c  pmutate  The jump mutation probability.  Typically set = 1/npopsiz.
c
c
c  For single function evaluations, set npopsiz=1, maxgen=1, & irestrt=0.
c
c  My favorite initial choices of GA parameters are:
c     microga=1, npopsiz=5, iunifrm=1, maxgen=200
c     microga=1, npopsiz=5, iunifrm=0, maxgen=200
c  I generally get good performance with both the uniform and single-
c  point crossover micro-GA.
c
c  For those wishing to use the more conventional GA techniques,
c  my old favorite choice of GA parameters was:
c     iunifrm=1, iniche=1, ielite=1, itourny=1, nchild=1
c  For most problems I have dealt with, I get good performance using
c     npopsiz=100, pcross=0.5, pmutate=0.01, pcreep=0.02, maxgen=26
c  or
c     npopsiz= 50, pcross=0.5, pmutate=0.02, pcreep=0.04, maxgen=51
c
c  Any negative integer for idum should work.  I typically arbitrarily
c  choose idum=-10000 or -20000.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Code variable definitions (those not defined above):
c
c  best     = the best fitness of the generation
c  child    = the floating point parameter array of the children
c  cpu      = cpu time of the calculation
c  cpu0,cpu1= cpu times associated with 'etime' timing function
c  creep    = +1 or -1, indicates which direction parameter creeps
c  delta    = del/nparam
c  diffrac  = fraction of total number of bits which are different
c             between the best and the rest of the micro-GA population.
c             Population convergence arbitrarily set as diffrac<0.05.
c  evals    = number of function evaluations
c  fbar     = average fitness of population
c  fitness  = array of fitnesses of the parents
c  fitsum   = sum of the fitnesses of the parents
c  genavg   = array of average fitness values for each generation
c  geni     = generation array
c  genmax   = array of maximum fitness values for each generation
c  g0       = lower bound values of the parameter array to be optimized.
c             The number of parameters in the array should match the
c             dimension set in the above parameter statement.
c  g1       = the increment by which the parameter array is increased
c             from the lower bound values in the g0 array.  The minimum
c             parameter value is g0 and the maximum parameter value
c             equals g0+g1*(2**g2-1), i.e. g1 is the incremental value
c             between min and max.
c  ig2      = array of the number of bits per parameter, i.e. the number
c             of possible values per parameter.  For example, ig2=2 is
c             equivalent to 4 (=2**2) possibilities, ig2=4 is equivalent
c             to 16 (=2**4) possibilities.
c  ig2sum   = sum of the number of possibilities of ig2 array
c  ibest    = binary array of chromosomes of the best individual
c  ichild   = binary array of chromosomes of the children
c  icount   = counter of number of different bits between best
c             individual and other members of micro-GA population
c  icross   = the crossover point in single-point crossover
c  indmax   = maximum # of individuals allowed, i.e. max population size
c  iparent  = binary array of chromosomes of the parents
c  istart   = the generation to be started from
c  jbest    = the member in the population with the best fitness
c  jelite   = a counter which tracks the number of bits of an individual
c             which match those of the best individual
c  jend     = used in conjunction with iend for debugging
c  jstart   = used in conjunction with iskip for debugging
c  kount    = a counter which controls how frequently the restart
c             file is written
c  kelite   = kelite set to unity when jelite=nchrome, indicates that
c             the best parent was replicated amongst the children
c  mate1    = the number of the population member chosen as mate1
c  mate2    = the number of the population member chosen as mate2
c  nchrmax  = maximum # of chromosomes (binary bits) per individual
c  nchrome  = number of chromosomes (binary bits) of each individual
c  ncreep   = # of creep mutations which occurred during reproduction
c  nmutate  = # of jump mutations which occurred during reproduction
c  nparmax  = maximum # of parameters which the chromosomes make up
c  paramav  = the average of each parameter in the population
c  paramsm  = the sum of each parameter in the population
c  parent   = the floating point parameter array of the parents
c  pardel   = array of the difference between parmax and parmin
c  rand     = the value of the current random number
c  npossum  = sum of the number of possible values of all parameters
c  tarray   = time array used with 'etime' timing function
c  time0    = clock time at start of run
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutines:
c  ____________
c
c  code     = Codes floating point value to binary string.
c  crosovr  = Performs crossover (single-point or uniform).
c  decode   = Decodes binary string to floating point value.
c  evalout  = Evaluates the fitness of each individual and outputs
c             generational information to the 'ga.out' file.
c  func     = The function which is being evaluated.
c  gamicro  = Implements the micro-GA technique.
c  input    = Inputs information from the 'ga.inp' file.
c  initial  = Program initialization and inputs information from the
c             'ga.restart' file.
c  mutate   = Performs mutation (jump and/or creep).
c  newgen   = Writes child array back into parent array for new
c             generation; also checks to see if best individual was
c             replicated (elitism).
c  niche    = Performs niching (sharing) on population.
c  possibl  = Checks to see if decoded binary string falls within
c             specified range of parmin and parmax.
c  ran3     = The random number generator.
c  restart  = Writes the 'ga.restart' file.
c  select   = A subroutine of 'selectn'.
c  selectn  = Performs selection; tournament selection is the only
c             option in this version of the code.
c  shuffle  = Shuffles the population randomly for selection.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
