# mutation-rate-intervals

Support for "The statistics of k-mers from a sequence undergoing a simple
mutation process without spurious matches," Blanca, Harris, Koslicki and
Medvedev. 2020.

### Quick start for computing an r1 confidence interval for the sketching mutation model

```bash 
$ r1-from-minhash-jaccard.py L=4.5M k=21 S=5K C=0.95 m=100 J=0.20
# (may take about a minute)
L       k  s    conf  jHat     r1Low    r1High
4500000 21 5000 0.950 0.200000 0.048886 0.053160
```

### Quick start for computing an r1 confidence interval from an observation of N<sub>mut</sub>

```bash 
$ r1-from-nmut.py L=4.5M k=21 C=0.95 N=3M
L       k  conf  nMut    r1Low    r1High
4500000 21 0.950 3000000 0.050725 0.051215
```

### Prerequisites

* python3
* scipy
* numpy

For computing confidence intervals, only scipy is required. An optional
module, mpmath, will be used if present (as described below).

Numpy is only used by the simulation programs.

Two addition packages are used if present: mpmath and mmh3.

mpmath is a multi-precision package, used here to avoid numerical problems that
can occur for very low mutation probabilities (e.g. 1e-8). If the module is not
present standard python floating-point is used instead.

mmh3 is a wrapper for MurmurHash3, used here for hashing kmers for bottom
sketches. If the module is not present the hashing options in
simulate_nucleotide_errors are not available.

### Usage Overview

The package has six parts:
* Two front-send modules to compute theoretical confidence intervals:
r1-from-minhash-jaccard.py and r1-from-nmut.py.
* A module to compute the theoretical confidence intervals for the sketching
mutation model described as theorem 6 in the manuscript:
hypergeometric_slicer.py.
* A module to compute the theoretical confidence intervals for r1 from an
observation of Nmut described as theorem 5 in the manuscript:
kmer_mutation_formulas_v1.py.
* Two programs to generate simulations: simulate_unit_errors.py and
simulate_nucleotide_errors.py.
* A program to evaluate how well the simulations conform to the theoretical
confidence intervals of theorem 6: evaluate_hypergeometric_slicer.py.
* Earlier versions of the simulation programs, simulate_unit_errors_v1.py and
simulate_nucleotide_errors_v1.py. These perform both simulation and
evaluation for theorem 5. 

Above we describe only the confidence interval modules. The simulation and
evaluation programs are described in the reproducibility folder.

