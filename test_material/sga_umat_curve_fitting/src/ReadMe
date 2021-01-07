This code initializes a random sample of individuals with different
parameters to be optimized using the genetic algorithm approach, i.e.
evolution via survival of the fittest.  The selection scheme used is
tournament selection with a shuffling technique for choosing random
pairs for mating.  The routine includes binary coding for the
individuals, jump mutation, creep mutation, and the option for
single-point or uniform crossover.  Niching (sharing) and an option
for the number of children per pair of parents has been added.  More
recently, an option for the use of a micro-GA has been added.


I have provided a sample subroutine "func", but ultimately
the user must supply this subroutine "func" which should be your
cost function.  You should be able to run the code with the
sample subroutine "func" and the provided ga.inp file and obtain
the optimal function value of 1.0000 at generation 187 with the uniform
crossover micro-GA enabled (this is 935 function evaluations).  Note that
because different computers may treat precision and truncation
differently, I have seen cases where two computers using the same
input produce different evolution histories (but still converge to the
optimal).

I still recommend using the micro-GA technique (microga=1)
with uniform crossover (iunifrm=1). However, if possible, I strongly
suggest that you use values of nposibl of 2**n (2, 4, 8, 16, 32, 64,
etc.). While my test function works fine for other values of nposibl,
I have encountered problems where the uniform crossover micro-GA has
difficulty with parameters having long bit strings and a non-2**n value
of nposibl, e.g. nposibl=1000, will have 10 bits assigned (for this case
I would suggest running nposibl=1024 rather than 1000);


_______________________________________________________________________

The code is presently set for a maximum population size of 200,
30 chromosomes (binary bits) and 2 parameters.  These values can be
changed in params.f as appropriate for your problem.  Correspondingly
you will have to change a few 'write' and 'format' statements if you
change nchrmax and/or nparmax.  In particular, if you change nchrome
and/or nparam, then you should change the 'format' statement numbers
1050, 1075, 1275, and 1500.  For example, if you have a problem with
4 parameters and 16 chromosomes (bits), then you should change these
format statements to be:
 1050 format(1x,' #      Binary Code',8x,'Param1 Param2 Param3',
     +          ' Param4 Fitness')
 1075 format(i3,1x,16i1,4(1x,f6.2),1x,f6.2)
 1275 format(/' Average Values:',10x,4(1x,f6.2),1x,f6.2/)
 1500 format(i5,3x,16i2)
 _______________________________________________________________________
###############################################################################

micro-GA Tip:

My favorite GA technique is still the micro-GA.  At this point, I recommend
using the micro-GA with uniform crossover and a small population size.  The
following inputs gave me excellent performance:

       microga = 1
       npopsiz = 5
       maxgen  = 100
       iunifrm = 1

I have also gotten good performance with the single-point crossover (iunifrm=0),
micro-GA.

If you decide to use the micro-GA, you will not need to worry about the
population sizing or creep mutation tips below.

See the Krishnakumar reference below for more information about micro-GA's.

###############################################################################

Population Sizing Tip:

I've had a lot of people ask me about population sizing, especially
people who are attempting large problems where 100 individuals is probably
not enough.  The true authority on the subject is David Goldberg, but here is
a crude population scaling law in my paper (based on Goldberg & Deb, 1992):

      npopsiz = order[(l/k)(2**k)] for binary coding

where l = nchrome and k is the average size of the schema of interest
(effectively the average number of bits per parameter, i.e. approximately
equal to nchrome/nparam, rounded to the nearest integer).  I find that when
I have uniform crossover and niching turned on (which I recommend doing),
that this scaling law is usually overkill, i.e. you can most likely get by
with populations at least twice as small.

Remember to make the parameter 'indmax' (in 'params.f') greater than or equal
to 'npopsiz'.

###############################################################################

Creep Mutation Probability Tip:

I generally like to have approximately the same number of creep mutations and
jump mutations per generation.  Using basic probabilistic arguments, it can be
shown that you will get approximately the same number of creep and jump
mutations when
                pcreep = (nchrome/nparam) * pmutate

where pmutate (the jump mutation probability) is 1/npopsiz.

###############################################################################

Suggested reading that I have found to be of use:

Goldberg, D. E., and Richardson, J., "Genetic Algorithms with
Sharing for Multimodal Function Optimization," Genetic Algorithms and their
Applications: Proceedings of the Second International Conference on Genetic
Algorithms, 1987, pp. 41-49.

Goldberg, D. E., "Genetic Algorithms in Search, Optimization and
Machine Learning," Addison-Wesley, 1989.

Goldberg, D. E., "A Note on Boltzmann Tournament Selection for
Genetic Algorithms and Population-Oriented Simulated Annealing," in:
Complex Systems, Vol. 4, Complex Systems Publications, Inc., 1990, pp.
445-460.

Goldberg, D. E., "Real-coded Genetic Algorithms, Virtual Alphabets,
and Blocking," in: Complex Systems, Vol. 5, Complex Systems Publications,
Inc., 1991, pp. 139-167.

Goldberg, D. E., and Deb, K., "A Comparitive Analysis of Selection
Schemes Used in Genetic Algorithms," in: Foundations of Genetic Algorithms,
ed. by Rawlins, G.J.E., Morgan Kaufmann Publishers, San Mateo, CA, pp.
69-93, 1991.

Goldberg, D. E., Deb, K., and Clark, J. H., "Genetic Algorithms,
Noise, and the Sizing of Populations," in: Complex Systems, Vol. 6, Complex
Systems Pub., Inc., 1992, pp. 333-362.

Krishnakumar, K., "Micro-Genetic Algorithms for Stationary and
Non-Stationary Function Optimization," SPIE: Intelligent Control and
Adaptive Systems, Vol. 1196, Philadelphia, PA, 1989.

Syswerda, G., "Uniform Crossover in Genetic Algorithms," in:
Proceedings of the Third International Conference on Genetic Algorithms,
Schaffer, J. (Ed.), Morgan Kaufmann Publishers, Los Altos, CA, pp. 2-9,
1989.

###############################################################################

If you are interested in my work (which may give some insights into how
and why I coded some aspects of my GA), I can mail copies of three papers of
mine.

G. Yang, L.E. Reinstein, S. Pai, Z. Xu, and D.L. Carroll, "A new genetic
algorithm technique in optimization of permanent 125-I prostate implants,"
Medical Physics, Vol. 25, No. 12, 1998, pp. 2308-2315.

Carroll, D. L., "Chemical Laser Modeling with Genetic Algorithms,"
AIAA J., Vol. 34,  2, 1996, pp.338-346.
     (A preprint version of this paper can now be downloaded in PDF format
      via my website:
      <http://cuaerospace.com/carroll/gatips.html> look for AIAA1996.pdf)

Carroll, D. L., "Genetic Algorithms and Optimizing Chemical Oxygen-Iodine
Lasers," Developments in Theoretical and Applied Mechanics, Vol. XVIII,
eds. H.B. Wilson, R.C. Batra, C.W. Bert, A.M.J. Davis, R.A. Schapery, D.S.
Stewart, and F.F. Swinson, School of Engineering, The University of Alabama,
1996, pp.411-424.
     (This paper can now be downloaded in PDF format via my website:
      <http://cuaerospace.com/carroll/gatips.html> look for SECTAM18.pdf)

###############################################################################

Disclaimer:  this program is not guaranteed to be free of error
(although it is believed to be free of error), therefore it should
not be relied on for solving problems where an error could result in
injury or loss.  If this code is used for such solutions, it is
entirely at the user's risk and the author disclaims all liability.
