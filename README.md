# mutation-rate-intervals

Support for "The statistics of k-mers from a sequence undergoing a simple
mutation process without spurious matches," Blanca, Harris, Koslicki and
Medvedev. 2020.

### Quick start for computing a confidence interval

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
# result is (0.048886239445729, 0.05316000812867323)
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

The package has three parts:
* A module to compute the theoretical confidence intervals described in the
manuscript: hypergeometric_slicer.py.
* Two programs to generate simulations: simulate_unit_errors.py and
simulate_nucleotide_errors.py.
* A program to evaluate how well the simulations conform to the theoretical
confidence intervals: evaluate_hypergeometric_slicer.py.

Here we describe only the confidence interval module. The simulation and
evaluation programs are described in the reproducibility folder.

(more to come)
