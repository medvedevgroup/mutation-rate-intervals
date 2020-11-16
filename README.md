# mutation-rate-intervals

Support for "The statistics of k-mers from a sequence undergoing a simple
mutation process without spurious matches," Blanca, Harris, Koslicki and
Medvedev. 2020.

### Quick start for computing a confidence interval for the sketching
mutation model

In a python3 interactive shell:
```bash 
import hypergeometric_slicer as hgslicer

L = 4500000        # number of kmers in sequence (4.5M)
k = 21             # kmer size
sketchSize = 5000  # size of bottom sketches (5K)
jHat = 0.20        # jaccard index observed in an experiment
confidence = 0.95
numSlices = 100    # number of slices (see manuscript)

hgslicer.r1_confidence_interval(L,k,sketchSize,1-confidence,numSlices,jHat)
# (may take about a minute)
# result is ≈ (0.048886,0.053160)
```

### Quick start for computing a confidence interval for r1 from an
observation of Nmut

In a python3 interactive shell:
```bash 
import kmer_mutation_formulas_v1 as v1

L = 4500000        # number of kmers in sequence (4.5M)
k = 21             # kmer size
nMut = 3000000     # number of mutated kmers observed in an experiment
confidence = 0.95

alpha = 1-confidence
z = v1.probit(1-alpha/2)
q1 = v1.q_for_n_affected_high(L,k,nMut,z)
q2 = v1.q_for_n_affected_low (L,k,nMut,z)
r1Low  = v1.q_to_r1(k,q1)
r1High = v1.q_to_r1(k,q2)
(r1Low,r1High)
# result is ≈ (0.050725,0.051215)
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

The package has four parts:
* A module to compute the theoretical confidence intervals for the sketching
mutation model described as theorem 6 in the manuscript:
hypergeometric_slicer.py.
* A module to compute the theoretical confidence intervals for r1 from an
observation of Nmut described as theorem 5 in the manuscript:
kmer_mutation_formulas_v1.py.
* Two programs to generate simulations: simulate_unit_errors.py and
simulate_nucleotide_errors.py.
* A program to evaluate how well the simulations conform to the theoretical
confidence intervals: evaluate_hypergeometric_slicer.py.

Above we describe only the confidence interval module. The simulation and
evaluation programs are described in the reproducibility folder.

