# mutation-rate-intervals

This software calculates confidence intervals and hypothesis tests for various
variables under a simple nucleotide mutation process. In this model, a sequence
B evolves from a sequence A by independently mutating every nucleotide with
probability r1. What is observed may be
* the number of mutated k-mers
* the Jaccard between A and B
* the minsketch Jaccard estimate between A and B

In particular, the software can be used to compute the following:
* given a mutation rate r1, a hypothesis test for the observed number of
  mutated k-mers, the observed Jaccard, or the observed minhash Jaccard estimate
* Given the observed number of mutated k-mers, the observed Jaccard, or the
  observed minhash Jaccard estimate, a confidence interval for the mutation
  rate r1.

The two most common cases where this is applicable are
* There are two sequences that have evolved from a common ancestor and we
  observe their Jaccard similarity. What was the mutation rate r1?
* There is a read r generated with a known error rate r1. We would like to know
  if r was generated from a sequence s based on the observed Jaccard between
  them.

### Quick start

To compute an r1 confidence interval from an observed number of mutated k-mers:

```bash 
$ r1-from-nmut.py L=4.5M k=21 C=0.95 N=2997034
L       k  conf  nMut    r1Low    r1High
4499980 21 0.950 2997034 0.050636 0.051126
```

To compute an r1 confidence interval from an observed Jaccard index:

```bash 
$ r1-from-jaccard.py L=4.5M k=21 C=0.95 J=0.20
L       k  conf  jaccard  r1Low    r1High
4499980 21 0.950 0.200000 0.050725 0.051215
```

To compute an r1 confidence interval from an observed minhash Jaccard estimate
in the sketching mutation model:

```bash 
$ r1-from-minhash-jaccard.py L=4.5M k=21 S=5K C=0.95 m=100 J=0.20
# (may take about a minute)
L       k  s    conf  jHat     slices r1Low    r1High
4499980 21 5000 0.950 0.200000 100    0.048886 0.053160
```

### How to choose parameters

The model assumes that k is large enough so that the effect of repetitive k-mers is negligible and that there are only point mutations (see paper for details). In this situation, L is the sequence length and N is the number of k-mers of A that are not in B (in this model, this is the same as the number of k-mers of B that are not in A).
In reality, these assumptions may be violated. 
How to handle it depends on the use cases, as follows:
TODO

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
* Two command-line programs to compute theoretical confidence intervals:
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

Above we describe only the confidence interval command-line programs. The
simulation and evaluation programs are described in the reproducibility folder.

### Citation
* Antonio Blanca, Robert S. Harris, David Koslicki and Paul Medvedev, "The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches", submitted 

