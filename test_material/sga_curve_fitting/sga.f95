!#######################################################################
!
      subroutine gafortran
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension ibest(nchrmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
      dimension geni(1000000),genavg(1000000),genmax(1000000)
      real*4 cpu,cpu0,cpu1,tarray(2),t_results
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt, &
         &             itourny,ielite,icreep,itwopoint,iunifrm,iniche,     &
         &             iskip,iend,nchild,microga,kountmx
!
!
!
      call etime(tarray,t_results)
      write(6,*) tarray(1),tarray(2)
      cpu0=tarray(1)
!
!  Call the input subroutine.
      TIME0=SECNDS(0.0)
      call input
!
!  Perform necessary initialization and read the ga.restart file.
      call initial(istart,npossum,ig2sum)
!
!  $$$$$ Main generational processing loop. $$$$$
      kount=0
      do 20 i=istart,maxgen+istart-1
         write (6,1111) i
         write (24,1111) i
         write(24,1050)
!
!  Evaluate the population, assign fitness, establish the best
!  individual, and write output information.
         call evalout(iskip,iend,ibest,fbar,best)
         geni(i)=float(i)
         genavg(i)=fbar
         genmax(i)=best
         write(48,*) i,fbar, best
         if(npopsiz.eq.1 .or. iskip.ne.0) then
            close(24)
            close(48)
            close(49)
            stop
         endif
!
!  Implement "niching".
         if (iniche.ne.0) call niche
!
!  Enter selection, crossover and mutation loop.
         ncross=0
         ipick=npopsiz
         do 45 j=1,npopsiz,nchild
!
!  Perform selection.
            call selectn(ipick,j,mate1,mate2)
!
!  Now perform crossover between the randomly selected pair.
            call crosovr(ncross,j,mate1,mate2)
 45      continue
         write(6,1225) ncross
         write(24,1225) ncross
!
!  Now perform random mutations.  If running micro-GA, skip mutation.
         if (microga.eq.0) call mutate
!
!  Write child array back into parent array for new generation.  Check
!  to see if the best parent was replicated.
         call newgen(ielite,npossum,ig2sum,ibest)
!
!  Implement micro-GA if enabled.
         if (microga.ne.0) call gamicro(i,npossum,ig2sum,ibest)
!
!  Write to restart file.
         call restart(i,istart,kount)
 20   continue
!  $$$$$ End of main generational processing loop. $$$$$
! 999  continue
      write(24,3000)
      do 100 i=istart,maxgen+istart-1
         evals=float(npopsiz)*geni(i)
         write(24,3100) geni(i),evals,genavg(i),genmax(i)
 100  continue
      call etime(tarray,t_results)
!      write(6,*) tarray(1),tarray(2)
      cpu1=tarray(1)
      cpu=(cpu1-cpu0)
!      write(6,1400) cpu,cpu/60.0
      write(24,1400) cpu,cpu/60.0
      CLOSE (24)
      CLOSE (48)
      CLOSE (49)
!
 1400 format(2X,'CPU time for all generations=',E12.6,' sec'/  &
          2X,'                             ',E12.6,' min') 
! 1050 format(1x,' #      Binary Code',48x,'  Param1  Param2  Param3  Param4  Fitness')
 1050 format(1x,' #      Binary Code',32x,'  Param1  Param2  Param3  Param4  Fitness')
  !1050 format(1x,' #      Binary Code',4x,'  Param1  Param2  Param3  Param4  Fitness')
 !1050 format(1x,' #      Binary Code',2x,'  Param1  Param2  Param3  Param4  Fitness')
 1111 format(//'#################  Generation',i5,'  #################')
 1225 format(/'  Number of Crossovers      =',i5)
 3000 format(2x//'Summary of Output'/&
    &       2x,'Generation   Evaluations   Avg.Fitness   Best Fitness')
 3100 format(2x,3(e10.4,4x),e11.5)
!
!      CALL SYSTEM('gnuplot -p data_plot.plt')
!      stop
      end
!
!#######################################################################
      subroutine input
!
!  This subroutine inputs information from the ga.inp (gafort.in) file.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      include 'param_umat.inc'
      dimension nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,        &
                      itourny,ielite,icreep,itwopoint,iunifrm,iniche,        &
                      iskip,iend,nchild,microga,kountmx
!
      namelist / ga   / irestrt,npopsiz,pmutate,maxgen,idum,pcross,        &
                       itourny,ielite,icreep,pcreep,itwopoint,iunifrm,iniche,   &
                       iskip,iend,nchild,nparam,parmin,parmax,nposibl,        &
                       nowrite,nichflg,microga,kountmx
!
      kountmx=5
      irestrt=0
      itourny=0
      ielite=0
      itwopoint=0
      iunifrm=0
      iniche=0
      iskip=0
      iend=0
      nchild=1
      do 2 i=1,nparmax
         nichflg(i)=1
 2    continue
      microga=0
!
      OPEN (UNIT=24, FILE='ga.out', STATUS='UNKNOWN')
      rewind 24
      OPEN (UNIT=48, FILE='plot_ga.out', STATUS='UNKNOWN')
      rewind 48
      OPEN (UNIT=49, FILE=DIRC, STATUS='UNKNOWN')
      rewind 49
      OPEN (UNIT=23, FILE='ga.inp', STATUS='OLD')
      READ (23, NML = ga)
      CLOSE (23)
      itourny=1
!      if (itourny.eq.0) nchild=2
!
!  Check for array sizing errors.
      if (npopsiz.gt.indmax) then
         write(6,1600) npopsiz
         write(24,1600) npopsiz
         close(24)
         stop
      endif
      if (nparam.gt.nparmax) then
         write(6,1700) nparam
         write(24,1700) nparam
         close(24)
         stop
      endif
!
!  If using the microga option, reset some input variables
      if (microga.ne.0) then
         pmutate=0.0d0
         pcreep=0.0d0
         itourny=1
         ielite=1
         iniche=0
         nchild=1
         if (iunifrm.eq.0) then
            pcross=1.0d0
         else
            pcross=0.5d0
         endif
      endif
!
 1600 format(1x,'ERROR: npopsiz > indmax.  Set indmax = ',i6)
 1700 format(1x,'ERROR: nparam > nparmax.  Set nparmax = ',i6)
!
      return
      end
!
!#######################################################################
      subroutine initial(istart,npossum,ig2sum)
!
!  This subroutine sets up the program by generating the g0, g1 and
!  ig2 arrays, and counting the number of chromosomes required for the
!  specified input.  The subroutine also initializes the random number
!  generator, parent and iparent arrays (reads the ga.restart file).
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension nposibl(nparmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,        &
                      itourny,ielite,icreep,itwopoint,iunifrm,iniche,        &
                      iskip,iend,nchild,microga,kountmx
!
!
      do 3 i=1,nparam
         g0(i)=parmin(i)
         pardel(i)=parmax(i)-parmin(i)
         g1(i)=pardel(i)/dble(nposibl(i)-1)
 3    continue
      do 6 i=1,nparam
         do 7 j=1,30
            n2j=2**j
            if (n2j.ge.nposibl(i)) then
               ig2(i)=j
               goto 8
            endif
            if (j.ge.30) then
               write(6,2000)
               write(24,2000)
               close(24)
               stop
            endif
 7       continue
 8       continue
 6    continue
!
!  Count the total number of chromosomes (bits) required
      nchrome=0
      npossum=0
      ig2sum=0
      do 9 i=1,nparam
         nchrome=nchrome+ig2(i)
         npossum=npossum+nposibl(i)
         ig2sum=ig2sum+(2**ig2(i))
 9    continue
      if (nchrome.gt.nchrmax) then
         write(6,1800) nchrome
         write(24,1800) nchrome
         close(24)
         stop
      endif
!
      if (npossum.lt.ig2sum .and. microga.ne.0) then
         write(6,2100)
         write(24,2100)
      endif
!
!  Initialize random number generator
      call ran3(idum,rand)
!
      if(irestrt.eq.0) then
!  Initialize the random distribution of parameters in the individual
!  parents when irestrt=0.
         istart=1
         do 10 i=1,npopsiz
            do 15 j=1,nchrome
               call ran3(1,rand)
               iparent(j,i)=1
               if(rand.lt.0.5d0) iparent(j,i)=0
 15         continue
 10      continue
         if (npossum.lt.ig2sum) call possibl(parent,iparent)
      else
!  If irestrt.ne.0, read from restart file.
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         read(25,*) istart,npopsiz
         do 1 j=1,npopsiz
            read(25,*) k,(iparent(l,j),l=1,nchrome)
 1       continue
         CLOSE (25)
      endif
!
      if(irestrt.ne.0) call ran3(idum-istart,rand)
!
 1800 format(1x,'ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
 2000 format(1x,'ERROR: You have a parameter with a number of '/ &
             1x,'   possibilities > 2**30!  If you really desire this,'/ &
             1x,'   change the DO loop 7 statement and recompile.'// &
             1x,'   You may also need to alter the code to work with'/ &
             1x,'   REAL numbers rather than INTEGER numbers; Fortran'/ &
             1x,'   does not like to compute 2**j when j>30.')
 2100 format(1x,'WARNING: for some cases, a considerable performance'/ &
             1x,'   reduction has been observed when running a non-'/ &
             1x,'   optimal number of bits with the micro-GA.'/ &
             1x,'   If possible, use values for nposibl of 2**n,'/ &
             1x,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
!
      return
      end
!
!#######################################################################
      subroutine evalout(iskip,iend,ibest,fbar,best)
!
!  This subroutine evaluates the population, assigns fitness,
!  establishes the best individual, and outputs information.
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
      dimension paramsm(nparmax),paramav(nparmax),ibest(nchrmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
!
      fitsum=0.0d0
      best=-1.0d10
      do 29 n=1,nparam
         paramsm(n)=0.0d0
 29   continue
      jstart=1
      jend=npopsiz
      if(iskip.ne.0) jstart=iskip
      if(iend.ne.0) jend=iend
      do 30 j=jstart,jend
         call decode(j,parent,iparent)
         if(iskip.ne.0 .and. iend.ne.0 .and. iskip.eq.iend) &
         write(6,1075) j,(iparent(k,j),k=1,nchrome),(parent(kk,j),kk=1,nparam),0.0
!
!  Call function evaluator, write out individual and fitness, and add
!  to the summation for later averaging.
         call func(j,funcval)
         fitness(j)=funcval
         write(24,1075) j,(iparent(k,j),k=1,nchrome), &
                       (parent(kk,j),kk=1,nparam),fitness(j)
         fitsum=fitsum+fitness(j)
         do 22 n=1,nparam
            paramsm(n)=paramsm(n)+parent(n,j)
 22      continue
!
!  Check to see if fitness of individual j is the best fitness.
         if (fitness(j).gt.best) then
            best=fitness(j)
            jbest=j
            do 24 k=1,nchrome
               ibest(k)=iparent(k,j)
 24         continue
         endif
 30   continue
!
!  Compute parameter and fitness averages.
      fbar=fitsum/dble(npopsiz)
      do 23 n=1,nparam
         paramav(n)=paramsm(n)/dble(npopsiz)
 23   continue
!
!  Write output information
      if (npopsiz.eq.1) then
         write(24,1075) 1,(iparent(k,1),k=1,nchrome), &
                       (parent(k,1),k=1,nparam),fitness(1)
         write(24,*) ' Average Values:'
         write(24,1275) (parent(k,1),k=1,nparam),fbar
      else
         write(24,1275) (paramav(k),k=1,nparam),fbar
         write(24,1285) (parent(k,jbest),k=1,nparam),best
         rewind(49)
         write(49,*) (parent(k,jbest),k=1,nparam)
      endif
      write(6,1100) fbar
      write(24,1100) fbar
      write(6,1200) best
      write(24,1200) best
!
 !1075 format(i3,1x,64i1,4(1x,f7.4),1x,f8.5)
 !1075 format(i3,1x,32i1,4(1x,f7.4),1x,f8.5)
 1075 format(i3,1x,48i1,4(1x,f7.4),1x,f8.5)
 !1075 format(i3,1x,16i1,4(1x,f7.4),1x,f8.5)
 1100 format(1x,'Average Function Value of Generation=',f8.5)
 1200 format(1x,'Maximum Function Value              =',f8.5/)
 1275 format(/' Average Values:',10x,4(1x,f7.4),1x,f8.5)
 1285 format(' Best Values:',10x,4(1x,f7.4),1x,f8.5/)
      return
      end
!
!#######################################################################
      subroutine niche
!
!  Implement "niching" through Goldberg's multidimensional phenotypic
!  sharing scheme with a triangular sharing function.  To find the
!  multidimensional distance from the best individual, normalize all
!  parameter differences.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
!
!   Variable definitions:
!
!  alpha   = power law exponent for sharing function; typically = 1.0
!  del     = normalized multidimensional distance between ii and all
!            other members of the population
!            (equals the square root of del2)
!  del2    = sum of the squares of the normalized multidimensional
!            distance between member ii and all other members of
!            the population
!  nniche  = number of niched parameters
!  sigshar = normalized distance to be compared with del; in some sense,
!            1/sigshar can be viewed as the number of regions over which
!            the sharing function should focus, e.g. with sigshar=0.1,
!            the sharing function will try to clump in ten distinct
!            regions of the phase space.  A value of sigshar on the
!            order of 0.1 seems to work best.
!  share   = sharing function between individual ii and j
!  sumshar = sum of the sharing functions for individual ii
!
!      alpha=1.0
      sigshar=0.1d0
      nniche=0
      do 33 jj=1,nparam
         nniche=nniche+nichflg(jj)
 33   continue
      if (nniche.eq.0) then
         write(6,1900)
         write(24,1900)
         close(24)
         stop
      endif
      do 34 ii=1,npopsiz
         sumshar=0.0d0
         do 35 j=1,npopsiz
            del2=0.0d0
            do 36 k=1,nparam
               if (nichflg(k).ne.0) then
                  del2=del2+((parent(k,j)-parent(k,ii))/pardel(k))**2
               endif
 36         continue
            del=(dsqrt(del2))/dble(nniche)
            if (del.lt.sigshar) then
!               share=1.0-((del/sigshar)**alpha)
               share=1.0d0-(del/sigshar)
            else
               share=0.0d0
            endif
            sumshar=sumshar+share/dble(npopsiz)
 35      continue
         if (sumshar.ne.0.0d0) fitness(ii)=fitness(ii)/sumshar
 34   continue
!
 1900 format(1x,'ERROR: iniche=1 and all values in nichflg array = 0'/ &
            1x,'       Do you want to niche or not?')
!
      return
      end
!
!#######################################################################
      subroutine selectn(ipick,j,mate1,mate2)
!
!  Subroutine for selection operator.  Presently, tournament selection
!  is the only option available.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension fitness(indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt, &
                       itourny,ielite,icreep,itwopoint,iunifrm,iniche, &
                       iskip,iend,nchild,microga,kountmx
!
!  If tournament selection is chosen (i.e. itourny=1), then
!  implement "tournament" selection for selection of new population.
      if(itourny.eq.1) then
         call select(mate1,ipick)
         call select(mate2,ipick)
!        write(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
         do 46 n=1,nchrome
            ichild(n,j)=iparent(n,mate1)
            if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate2)
 46      continue
      endif
!
      return
      end
!
!#######################################################################
      subroutine crosovr(ncross,j,mate1,mate2)
!
!  Subroutine for crossover between the randomly selected pair.
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
!
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt, &
                       itourny,ielite,icreep,itwopoint,iunifrm,iniche, &
                        iskip,iend,nchild,microga,kountmx
!
       if (iunifrm.eq.0 .and. itwopoint.eq.0)  then
 !  Single-point crossover at a random chromosome point.
          call ran3(1,rand)
          if(rand.gt.pcross) goto 69
          ncross=ncross+1
          call ran3(1,rand)
          icross=2+dint(dble(nchrome-1)*rand)
          do 50 n=icross,nchrome
             ichild(n,j)=iparent(n,mate2)
             if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
  50      continue
       elseif (iunifrm.eq.0 .and. itwopoint.ne.0) then
 !  Two-point crossover at a random chromosome point.
          call ran3(1,rand)
          if(rand.gt.pcross) goto 69
          ncross=ncross+1
          !1st point
          call ran3(1,rand)
          icross=2+dint(dble(nchrome-1)*rand)
          !2nd point
          call ran3(1,rand)
          icross2=2+dint(dble(nchrome-1)*rand)
          if (icross .ne. icross2) then
          ! ordering icross positions
               if (icross2 .lt. icross) then
                 iaux=icross2
                 icross2=icross
                 icross=iaux
               endif
          do 70 n=1,icross
             ichild(n,j)=iparent(n,mate2)
             if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
  70      continue
          do 72 n=icross2,nchrome
             ichild(n,j)=iparent(n,mate1)
             if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate2)
  72      continue
           else
           !perform single point crossover if we get the same position for the two points
           !(particular case)
              do 80 n=icross,nchrome
             ichild(n,j)=iparent(n,mate2)
             if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
  80      continue
          endif
       else
 ! Perform uniform crossover between the randomly selected pair.
          do 60 n=1,nchrome
             call ran3(1,rand)
             if(rand.le.pcross) then
                ncross=ncross+1
                ichild(n,j)=iparent(n,mate2)
                if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
             endif
  60      continue
       endif
  69   continue
!
      return
      end
!
!#######################################################################
      subroutine mutate
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension nposibl(nparmax)
      dimension child(nparmax,indmax),ichild(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt, &
                       itourny,ielite,icreep,itwopoint,iunifrm,iniche, &
                       iskip,iend,nchild,microga,kountmx
!
!  This subroutine performs mutations on the children generation.
!  Perform random jump mutation if a random number is less than pmutate.
!  Perform random creep mutation if a different random number is less
!  than pcreep.
      nmutate=0
      ncreep=0
      do 70 j=1,npopsiz
         do 75 k=1,nchrome
!  Jump mutation (flip bit)
            call ran3(1,rand)
            if (rand.le.pmutate) then
               nmutate=nmutate+1
               if(ichild(k,j).eq.0) then
                  ichild(k,j)=1
               else
                  ichild(k,j)=0
               endif
               if (nowrite.eq.0) write(6,1300) j,k
               if (nowrite.eq.0) write(24,1300) j,k
            endif
 75      continue
!  Creep mutation (one discrete position away).
         if (icreep.ne.0) then
            do 76 k=1,nparam
               call ran3(1,rand)
               if(rand.le.pcreep) then
                  call decode(j,child,ichild)
                  ncreep=ncreep+1
                  creep=1.0d0
                  call ran3(1,rand)
                  if (rand.lt.0.5d0) creep=-1.0d0
                  child(k,j)=child(k,j)+g1(k)*creep
                  if (child(k,j).gt.parmax(k)) then
                     child(k,j)=parmax(k)-1.0d0*g1(k)
                  elseif (child(k,j).lt.parmin(k)) then
                     child(k,j)=parmin(k)+1.0d0*g1(k)
                  endif
                  call code(j,k,child,ichild)
                  if (nowrite.eq.0) write(6,1350) j,k
                  if (nowrite.eq.0) write(24,1350) j,k
               endif
 76         continue
         endif
 70   continue
      write(6,1250) nmutate,ncreep
      write(24,1250) nmutate,ncreep
!
 1250 format(/'  Number of Jump Mutations  =',i5/ &
              '  Number of Creep Mutations =',i5)
 1300 format('*** Jump mutation performed on individual  ',i4, &
             ', chromosome ',i3,' ***')
 1350 format('*** Creep mutation performed on individual ',i4, &
             ', parameter  ',i3,' ***')
!
      return
      end
!
!#######################################################################
      subroutine newgen(ielite,npossum,ig2sum,ibest)
!
!  Write child array back into parent array for new generation.  Check
!  to see if the best parent was replicated; if not, and if ielite=1,
!  then reproduce the best parent into a random slot.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
      dimension ibest(nchrmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
!
      if (npossum.lt.ig2sum) call possibl(child,ichild)
      kelite=0
      do 94 j=1,npopsiz
         jelite=0
         do 95 n=1,nchrome
            iparent(n,j)=ichild(n,j)
            if (iparent(n,j).eq.ibest(n)) jelite=jelite+1
            if (jelite.eq.nchrome) kelite=1
 95      continue
 94   continue
      if (ielite.ne.0 .and. kelite.eq.0) then
         call ran3(1,rand)
         irand=1d0+dint(dble(npopsiz)*rand)
         do 96 n=1,nchrome
            iparent(n,irand)=ibest(n)
 96      continue
         write(24,1260) irand
      endif
!
 1260 format('  Elitist Reproduction on Individual ',i4)
!
      return
      end
!
!#######################################################################
      subroutine gamicro(i,npossum,ig2sum,ibest)
!
!  Micro-GA implementation subroutine
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension ibest(nchrmax)
!
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
!
!  First, check for convergence of micro population.
!  If converged, start a new generation with best individual and fill
!  the remainder of the population with new randomly generated parents.
!
!  Count number of different bits from best member in micro-population
      icount=0
      do 81 j=1,npopsiz
         do 82 n=1,nchrome
            if(iparent(n,j).ne.ibest(n)) icount=icount+1
 82      continue
 81   continue
!
!  If icount less than 5% of number of bits, then consider population
!  to be converged.  Restart with best individual and random others.
      diffrac=dble(icount)/dble((npopsiz-1)*nchrome)
      if (diffrac.lt.0.05d0) then
      do 87 n=1,nchrome
         iparent(n,1)=ibest(n)
 87   continue
      do 88 j=2,npopsiz
         do 89 n=1,nchrome
            call ran3(1,rand)
            iparent(n,j)=1
            if(rand.lt.0.5d0) iparent(n,j)=0
 89      continue
 88   continue
      if (npossum.lt.ig2sum) call possibl(parent,iparent)
      write(6,1375) i
      write(24,1375) i
      endif
!
 1375 format(//'%%%%%%%  Restart micro-population at generation', &
             i5,'  %%%%%%%')
!
      return
      end
!
!#######################################################################
      subroutine select(mate,ipick)
!
!  This routine selects the better of two possible parents for mating.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
!
      if(ipick+1.gt.npopsiz) call shuffle(ipick)
      ifirst=ipick
      isecond=ipick+1
      ipick=ipick+2
      if(fitness(ifirst).gt.fitness(isecond)) then
         mate=ifirst
      else
         mate=isecond
      endif
!     write(3,*)'select',ifirst,isecond,fitness(ifirst),fitness(isecond)
!
      return
      end
!
!#######################################################################
      subroutine shuffle(ipick)
!
!  This routine shuffles the parent array and its corresponding fitness
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
!
      ipick=1
      do 10 j=1,npopsiz-1
         call ran3(1,rand)
         iother=j+1+dint(dble(npopsiz-j)*rand)
         do 20 n=1,nchrome
            itemp=iparent(n,iother)
            iparent(n,iother)=iparent(n,j)
            iparent(n,j)=itemp
 20      continue
         temp=fitness(iother)
         fitness(iother)=fitness(j)
         fitness(j)=temp
 10   continue
!
      return
      end
!
!#######################################################################
      subroutine decode(i,array,iarray)
!
!  This routine decodes a binary string to a real number.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
!
      l=1
      do 10 k=1,nparam
         iparam=0
         m=l
         do 20 j=m,m+ig2(k)-1
            l=l+1
            iparam=iparam+iarray(j,i)*(2**(m+ig2(k)-1-j))
 20      continue
         array(k,i)=g0(k)+g1(k)*dble(iparam)
 10   continue
!
      return
      end
!
!#######################################################################
      subroutine code(j,k,array,iarray)
!
!  This routine codes a parameter into a binary string.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
!
!  First, establish the beginning location of the parameter string of
!  interest.
      istart=1
      do 10 i=1,k-1
         istart=istart+ig2(i)
 10   continue
!
!  Find the equivalent coded parameter value, and back out the binary
!  string by factors of two.
      m=ig2(k)-1
      if (g1(k).eq.0.0d0) return
      iparam=nint((array(k,j)-g0(k))/g1(k))
      do 20 i=istart,istart+ig2(k)-1
         iarray(i,j)=0
         if ((iparam+1).gt.(2**m)) then
            iarray(i,j)=1
            iparam=iparam-2**m
         endif
         m=m-1
 20   continue
!     write(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
!
      return
      end
!
!#######################################################################
!
      subroutine possibl(array,iarray)
!
!  This subroutine determines whether or not all parameters are within
!  the specified range of possibility.  If not, the parameter is
!  randomly reassigned within the range.  This subroutine is only
!  necessary when the number of possibilities per parameter is not
!  optimized to be 2**n, i.e. if npossum < ig2sum.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax),nposibl(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
!
      do 10 i=1,npopsiz
         call decode(i,array,iarray)
         do 20 j=1,nparam
            n2ig2j=2**ig2(j)
            if(nposibl(j).ne.n2ig2j .and. array(j,i).gt.parmax(j)) then
               call ran3(1,rand)
               irand=dint(dble(nposibl(j))*rand)
               array(j,i)=g0(j)+dble(irand)*g1(j)
               call code(i,j,array,iarray)
               if (nowrite.eq.0) write(6,1000) i,j
               if (nowrite.eq.0) write(24,1000) i,j
            endif
 20      continue
 10   continue
!
 1000 format('*** Parameter adjustment to individual     ',i4, &
            ', parameter  ',i3,' ***')
!
      return
      end
!
!#######################################################################
      subroutine restart(i,istart,kount)
!
!  This subroutine writes restart information to the ga.restart file.
!
      implicit real*8 (a-h,o-z)
      save
!
      include 'params.inc'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt, &
                       itourny,ielite,icreep,itwopoint,iunifrm,iniche, &
                       iskip,iend,nchild,microga,kountmx

      kount=kount+1
      if(i.eq.maxgen+istart-1 .or. kount.eq.kountmx) then
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         write(25,*) i+1,npopsiz
         do 80 j=1,npopsiz
            write(25,1500) j,(iparent(l,j),l=1,nchrome)
 80      continue
         CLOSE (25)
         kount=0
      endif
!
 1500 format(i5,3x,30i2)
!
      return
      end
!
!#######################################################################
      subroutine ran3(idum,rand)
!
!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', "Numerical Recipes" p. 199.
!
      implicit real*8 (a-h,m,o-z)
      save
!      implicit real*4(m)
      parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
!     parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
!
!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.
      dimension ma(55)
      data iff /0/
      if (idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=mseed-dble(iabs(idum))
         mj=dmod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      rand=mj*fac
      return
      end
!
!#######################################################################
!
      subroutine func(j,funcval)
 !
       implicit real*8 (a-h,o-z)
       save
 !
       include 'params.inc'
       INCLUDE 'PARAM_UMAT.INC'
       dimension parent(nparmax,indmax)
       dimension iparent(nchrmax,indmax)
 !      dimension parent2(indmax,nparmax),iparent2(indmax,nchrmax)
 !
       common / ga2   / nparam,nchrome
       common / ga3   / parent,iparent
 !
       COMMON /KDATA/PTS
       DOUBLE PRECISION PTS(NPTS,2)
       PARAMETER(NTENS = 6, NSTATEV = 3, NPROPS = 14, NDI=3, NSHR=3)
       DIMENSION STRESS(NTENS),STATEV(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
 DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
       DOUBLE PRECISION gamma,aux,aux2
 !
 !  Benchmark function used in Goldberg and Richardson (1987) (maximum)
 !
 !      nvalley=6
 !      pi=4.0d0*datan(1.d0)
 !      funcval=1.0d0
 !      do 10 i=1,nparam
 !         f1=(sin(5.1d0*pi*parent(i,j) + 0.5d0))**nvalley
 !         f2=exp(-4.0d0*log(2.0d0)*((parent(i,j)-0.0667d0)**2)/0.64d0)
 !         funcval=funcval*f1*f2
 ! 10   continue
 !
 ! sphere function (minimum)
         !
         DO i1=1,NTENS
             DO j1=1,NTENS
         DDSDDE(i1,j1)=0.0D0
             ENDDO
             STRESS(i1)=0.0D0
         ENDDO
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
         !       do 11 i=1,nparam
         ! FREE MATERIAL PROPERTIES
         PROPS(2)=parent(1,j)
         PROPS(4)=parent(2,j)
         PROPS(5)=parent(3,j)
         PROPS(6)=parent(4,j) 
!
          funcval=0.d0
        do 10 j1=1,npts
           gamma=pts(j1,1)
           sexp=pts(j1,2)
           aux=aux+sexp
           DFGRD1(1,1)=gamma
           DFGRD1(2,2)=ONE/SQRT(gamma)
           DFGRD1(3,3)=ONE/SQRT(gamma)
         CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
         DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
         NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
         ! statev(3) is the 1st PK at direction(1,1)
         aux=(sexp-statev(3))!*(sexp**(-1.d0))
         aux2=aux*aux
         funcval=funcval+aux2
    10   continue
 !
 !        write(*,*) funcval
       funcval=1.d0/(1.d0+funcval)
  !     funcval=-funcaval
 !      write(*,*) funcval
 !      write(*,*)
 !
 ! rosenbrock function (minimum)
 !
 !funcval=0.d0
 !      do 10 i=1,nparam-1
 !         f1=parent(i+1,j)-parent(i,j)*parent(i,j)
 !         f2=parent(i,j)-1
 !         funcval=funcval+100.d0*f1**2.d0+f2**2.d0
 ! 10   continue
 !
 ! funcval=funcval*funcval
 ! funcval=1.d0/(1.d0+funcval)
       return
       end
!#######################################################################
