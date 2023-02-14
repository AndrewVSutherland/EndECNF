# EndECNF
Software associated to the paper *[Computing the endomorphism ring of an elliptic curve over a number field](http://arxiv.org/abs/2301.11169)*, by John E. Cremona and Andrew V. Sutherland

This repository includes sample implementations of Algorithm 1 and Algorithm 2 described in the paper in the following computer algebra systems

- [Magma](http://magma.maths.usyd.edu.au/magma/): see the file [cmdisc.m](cmdisc.m), tested on Magma version 2.27-7
- [Pari/GP](https://pari.math.u-bordeaux.fr/): see the file [cmdisc.gp](cmdisc.gp), tested on Pari/GP version 2.15.2
- [SageMath](https://www.sagemath.org/): see the file [cmdisc.py](cmdisc.py), tested on SageMath version 9.7

These are unoptimized implementations that incorporate some but not all of the implementation suggestions noted in Remarks in the paper (see the individual source files for details).

All three implementations include a function `CMProfile` that can be used to compare the performance of the implementations.  This function was used to generate the timings for Algorithms 2 that are listed in Table 1 of the paper using median discriminants taken from the file [cmdiscs1000.txt](cmdiscs1000.txt), which contains a list of discriminants $D<0$ with $h(D) \le 1000$ that is known to be complete under GRH.

Also included are files used to profile the distributed implementations in Magma and SageMath:

 - [Magma](http://magma.maths.usyd.edu.au/magma/): see the file [magma_cm.m](magma_cm.m), tested on Magma version 2.27-7
 - [SageMath](https://www.sagemath.org/): see the file [sage_cm.py](sage_cm.py), tested on SageMath version 9.7

Each includes a function `TestCMDiscriminant` that can be used to generate the timings listed in the final two columns of Table 1 of the paper.
