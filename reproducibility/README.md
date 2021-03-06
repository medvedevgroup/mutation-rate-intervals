# reproducibility

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

Simulate 1,000 trials, with sketch sizes 100 and 500, then evaluate how well
they conformed to the theoretical 100-slice confidence interval:

```bash 
simulate_unit_errors.py --linear \
  T=1K L=10K K=21 --poisson=10% \
  --sketch=100,500 \
  > simulation.dat

cat simulation.dat \
  | evaluate_hypergeometric_slicer.py --confidence=95% --slices=100
```
Simulate 1,000 trials for a real genomic sequence, using a 128 bit
pseudo-random hash function, then evaluate how well they conformed to the
theoretical 100-slice confidence interval:

```bash 
cat ecoli.fa \
  | simulate_nucleotide_errors.py --linear \
      --hashbits=128  \
      T=1K K=21 --poisson=10% \
      --sketch=100,500 \
  > simulation.dat

cat simulation.dat \
  | evaluate_hypergeometric_slicer.py --confidence=95% --slices=100 --useL.A,B
```

