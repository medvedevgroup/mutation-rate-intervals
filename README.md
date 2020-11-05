# mutation-rate-intervals

Support for "Statistics of k-mers from a sequence undergoing a simple mutation
process without spurious matches," Blanca, Harris, Koslicki and Medvedev. 2020.

### Prerequisites

* python3
* scipy
* numpy

Two addition packages are used if present: mpmath and mmh3.

mpmath is a multi-precision package, used here to avoid numerical problems that
can occur for very low mutation probabilities (e.g. 1e-8). If the module is not
present standard python floating-point is used instead.

mmh3 is a wrapper for MurmurHash3, used here for hashing kmers for bottom
sketches. If the module is not present the hashing options in
simulate_nucleotide_errors are not available.

### Usage Overview

Simulate 1,000 trials for sequences of length 10,000 with k=21 and a mutation
rate of 10%:

```bash 
simulate_unit_errors.py --linear T=1K L=10K K=21 --poisson=10%
```

(more to come)

