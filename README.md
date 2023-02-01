# EndECNF
Software associated to the paper "Computing the endomorphism ring of an elliptic curve over a number field", by John E. Cremona and Andrew V. Sutherland

This repository includes sample implementations of Algorithm 1 and Algorithm 2 described in the paper in the following computer algebra systems

- [Magma](http://magma.maths.usyd.edu.au/magma/): see the file [cmdisc.m](cmdisc.m), tested on Magma version 2.27-7
- [Pari/GP](https://pari.math.u-bordeaux.fr/): see the file [cmdisc.gp](cmdisc.gp), tested on Pari/GP version 2.15.2
- [SageMath](https://www.sagemath.org/): see the file [cmdisc.py](cmdisc.py), tested on SageMath version 9.7

These are unoptimized implementations that incorporate some but not all of the implementation suggestions noted in Remarks in the paper (see the individual source files for details).

All three implementations include a function `CMProfile` that can be used to generate timings for sample discriminants, including all the timings listed in Table 1 of the paper.
