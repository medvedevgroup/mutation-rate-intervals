# mutation-rate-intervals

This software calculates confidence intervals and hypothesis tests for various
variables under a simple nucleotide mutation process. In this model, a sequence
B evolves from a sequence A by independently mutating every nucleotide with
probability r<sub>1</sub>. We then observe a variable X which may be one of the
following:
* The number of mutated k-mers
* The Jaccard between A and B
* The minsketch Jaccard estimate between A and B

The software can be used to compute the following:
* Given a mutation rate r<sub>1</sub> and a desired significance level C (e.g. 95%), an interval such that the observed X lies in it with probability C 
(i.e. a hypothesis test for X)
* Given the observed value of X and a desired confidence C (e.g. 95%), an interval such that the mutation rate r<sub>1</sub> lies in it with probability C 
(i.e. a confidence interval for r<sub>1</sub>).

The two most common cases where this is applicable are
* There are two sequences that have evolved from a common ancestor and we
  observe X. What was the mutation rate r<sub>1</sub>?
* There is a read R generated with a known error rate r<sub>1</sub>. We would like to know
  if R was generated from a sequence s based on the observed X between
  them.

### Quick start

To compute a hypothesis test for the observed number of mutated k-mers:
```bash 
$ r1-to-nmut-hypothesis.py L=4.5M k=21 C=0.95 r1=.05
L       k  sig   r1       nMutLow nMutHigh
4500000 21 0.950 0.050000 2959275 2975671
```

To compute a hypothesis test for the observed Jaccard index:
```bash 
$ r1-to-jaccard-hypothesis.py L=4.5M k=21 C=0.95 r1=.05
L       k  sig   r1       jLow     jHigh
4500000 21 0.950 0.050000 0.203905 0.206552
```

To compute a hypothesis test for the observed minhash Jaccard estimate:
```bash 
$ r1-to-minhash-jaccard-hypothesis.py L=4.5M k=21 S=5K C=0.95 r1=.05
L       k  s    sig   slices r1       jHatLow  jHatHigh
4500000 21 5000 0.950 100    0.050000 0.194000 0.216800
```

To compute an r<sub>1</sub> confidence interval from an observed number of mutated k-mers:
```bash 
$ r1-from-nmut.py L=4.5M k=21 C=0.95 N=2997034
L       k  conf  nMut    r1Low    r1High
4500000 21 0.950 2997034 0.050636 0.051126
```

To compute an r<sub>1</sub> confidence interval from an observed Jaccard index:
```bash 
$ r1-from-jaccard.py L=4.5M k=21 C=0.95 J=0.20
L       k  conf  jaccard  r1Low    r1High
4500000 21 0.950 0.200000 0.050725 0.051215
```

To compute an r<sub>1</sub> confidence interval from an observed minhash Jaccard estimate
in the sketching mutation model:
```bash 
$ r1-from-minhash-jaccard.py L=4.5M k=21 S=5K C=0.95 J=0.20
# (may take about a minute)
L       k  s    conf  slices jHat     r1Low    r1High
4500000 21 5000 0.950 100    0.200000 0.048886 0.053160
```

### How to choose L and other parameters
In reality, you may not know L. In such cases, we recommend that you estimate
it from what you know. For example, if what you know is that the number of
distinct (i.e. counting duplicates only once) k-mers in A is nA and in B is nB,
then you can set L = (nA + nB) / 2. You can also try setting L = min(nA, nB) or
L = max(nA, nB).   

You may also want to get a confidence interval on r<sub>1</sub> from the number
of mutated k-mers N, but you might only known the number of shared k-mers, i.e.
the number of k-mers in both A and B. If this number is n, then you can set
N = L - n.

Note that the programs consider L (uppercase) as the number of kmers in the
sequence, and l (lowercase) as the number of nucleotides, with l = L + k-1.

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

### Repo organization 

The package has six parts:
* Six command-line programs to compute theoretical hypothesis tests or confidence intervals:
r1-to-nmut-hypothesis.py,
r1-to-jaccard-hypothesis.py,
r1-to-minhash-jaccard-hypothesis.py,
r1-from-nmut.py,
r1-from-jaccard.py,
and r1-from-minhash-jaccard.py.
* A module to compute the theoretical confidence intervals for the sketching
mutation model described as theorem 3 in the manuscript:
hypergeometric_slicer.py.
* A module to compute the theoretical confidence intervals for r<sub>1</sub> from an
observation of Nmut described as corollary 1 in the manuscript:
kmer_mutation_formulas_v1.py.
* Two programs to generate simulations: simulate_unit_errors.py and
simulate_nucleotide_errors.py.
* A program to evaluate how well the simulations conform to the theoretical
confidence intervals of theorem 3: evaluate_hypergeometric_slicer.py.
* Earlier versions of the simulation programs, simulate_unit_errors_v1.py and
simulate_nucleotide_errors_v1.py. These perform both simulation and
evaluation for corollary 1. 

Above we describe only the confidence interval command-line programs. The
simulation and evaluation programs are described in the reproducibility folder.

### Usage Details

r1-to-nmut-hypothesis.py

```bash  
Compute hypothesis test for the number of mutated k-mers, given the mutation
rate r1.

usage: r1-to-nmut-hypothesis.py [options]
  --r1=<list>                   (R1=) (cumulative) mutation rate; <list> is a
                                comma-separated list of probabilities
  --length=<N>                  (l=) sequence length (number of NUCLEOTIDES in
                                the sequence)
                                (default is 1000 plus kmer size minus 1)
  L=<N>                         (L=) sequence length (number of KMERS in
                                the sequence)
                                (default is 1000)
  --k=<N>                       (K=) kmer size
                                (default is 21)
  --significance=<probability>  (C=) significance level
                                (default is 95%)
  --validate[=<N>]              (V=) run simulations to validate the interval;
                                N is the number of simulation trials; if N is
                                not provided, 10,000 simulations are run
                                (by default, no simulation is performed)
  --seed=<string>               random seed for simulations
  --progress=<number>           periodically report how many simulations we've
                                performed
```

r1-to-jaccard-hypothesis.py

```bash  
Compute hypothesis test for the jaccard index, given the mutation rate r1.

usage: r1-to-jaccard-hypothesis.py [options]
  --r1=<list>                   (R1=) (cumulative) mutation rate; <list> is a
                                comma-separated list of probabilities
  --length=<N>                  (l=) sequence length (number of NUCLEOTIDES in
                                the sequence)
                                (default is 1000 plus kmer size minus 1)
  L=<N>                         (L=) sequence length (number of KMERS in
                                the sequence)
                                (default is 1000)
  --k=<N>                       (K=) kmer size
                                (default is 21)
  --significance=<probability>  (C=) significance level
                                (default is 95%)
  --validate[=<N>]              (V=) run simulations to validate the interval;
                                N is the number of simulation trials; if N is
                                not provided, 10,000 simulations are run
                                (by default, no simulation is performed)
  --seed=<string>               random seed for simulations
  --progress=<number>           periodically report how many simulations we've
                                performed
```

r1-to-minhash-jaccard-hypothesis.py

```bash  
Compute hypothesis test for the minhash Jaccard estimate, given the mutation
rate r1.

usage: r1-to-minhash-jaccard-hypothesis.py [options]
  --r1=<list>                   (R1=) (cumulative) mutation rate; <list> is a
                                comma-separated list of probabilities
  --length=<N>                  (l=) sequence length (number of NUCLEOTIDES in
                                the sequence)
                                (default is 1000 plus kmer size minus 1)
  L=<N>                         (L=) sequence length (number of KMERS in
                                the sequence)
                                (default is 1000)
  --k=<N>                       (K=) kmer size
                                (default is 21)
  --sketch=<N>                  (S=) sketch size
                                (there is no default)
  --significance=<probability>  (C=) significance level
                                (default is 95%)
  --slices=<N>                  (m=) number of slices
                                (default is 100)
  --validate[=<N>]              (V=) run simulations to validate the interval;
                                N is the number of simulation trials; if N is
                                not provided, 10,000 simulations are run
                                (by default, no simulation is performed)
  --seed=<string>               random seed for simulations
  --progress=<number>           periodically report how many simulations we've
                                performed
```

r1-from-nmut.py

```bash  
Compute confidence interval for the mutation rate r1, given the observed number
of mutated k-mers.

usage: r1-from-nmut.py [options]
  --nmut=<list>               (N=) (cumulative) observed number of mutated
                              k-mers; <list> is a comma-separated list of
                              numbers
  --length=<N>                (l=) sequence length (number of NUCLEOTIDES in
                              the sequence)
                              (default is 1000 plus kmer size minus 1)
  L=<N>                       (L=) sequence length (number of KMERS in
                              the sequence)
                              (default is 1000)
  --k=<N>                     (K=) kmer size
                              (default is 21)
  --confidence=<probability>  (C=) size of confidence interval
                              (default is 95%)
  --validate[=<N>]            (V=) run simulations to validate the interval;
                              N is the number of simulation trials; if N is not
                              provided, 10,000 simulations are run
                              (by default, no simulation is performed)
  --seed=<string>             random seed for simulations
  --progress=<number>         periodically report how many simulations we've
                              performed
```

r1-from-jaccard.py

```bash  
Compute confidence interval for the mutation rate r1, given the observed
Jaccard.

usage: r1-from-jaccard.py [options]
  --jhat=<list>               (J=) (cumulative) observations of jaccard
                              index; <list> is a comma-separated list of
                              numbers between 0 and 1
  --length=<N>                (l=) sequence length (number of NUCLEOTIDES in
                              the sequence)
                              (default is 1000 plus kmer size minus 1)
  L=<N>                       (L=) sequence length (number of KMERS in
                              the sequence)
                              (default is 1000)
  --k=<N>                     (K=) kmer size
                              (default is 21)
  --confidence=<probability>  (C=) size of confidence interval
                              (default is 95%)
  --validate[=<N>]            (V=) run simulations to validate the interval;
                              N is the number of simulation trials; if N is not
                              provided, 10,000 simulations are run
                              (by default, no simulation is performed)
  --seed=<string>             random seed for simulations
  --progress=<number>         periodically report how many simulations we've
                              performed
```

r1-from-minhash-jaccard.py

```bash  
Compute confidence interval for the mutation rate r1, given the observed
minhash Jaccard estimate.

usage: r1-from-minhash-jaccard.py [options]
  --jhat=<list>               (J=) (cumulative) observed estimates of jaccard
                              index; <list> is a comma-separated list of
                              numbers between 0 and 1
  --length=<N>                (l=) sequence length (number of NUCLEOTIDES in
                              the sequence)
                              (default is 1000 plus kmer size minus 1)
  L=<N>                       (L=) sequence length (number of KMERS in
                              the sequence)
                              (default is 1000)
  --k=<N>                     (K=) kmer size
                              (default is 21)
  --sketch=<N>                (S=) sketch size
                              (there is no default)
  --confidence=<probability>  (C=) size of confidence interval
                              (default is 95%)
  --slices=<N>                (m=) number of slices
                              (default is 100)
  --validate[=<N>]            (V=) run simulations to validate the interval;
                              N is the number of simulation trials; if N is not
                              provided, 10,000 simulations are run
                              (by default, no simulation is performed)
  --seed=<string>             random seed for simulations
  --progress=<number>         periodically report how many simulations we've
                              performed
```

simulate_unit_errors.py

```bash  
usage: simulate_unit_errors [options]
  --k=<N>                   (K=) kmer size
                            (default is 21)
  --n=<N>                   (N= or L=) sequence length (number of KMERS in the
                            sequence)
                            (default is 100)
  --sketch=<N>              (S=) (cumulative) sketch size
                            (default is "no sketch")
  --sequences=<N>           (T=) number of sequence pairs to generate
                            (default is 1)
  --poisson=<probability>   (P=) (required) inject random sequencing errors
                            (substitutions); each base suffers a substitution
                            error with the given probability (poisson-like
                            noise)
  --bernoulli=<probability> (B= or E=) (required) inject random sequencing
                            errors (substitutions); exactly
                            round(L*<probability>) errors will occur in the
                            sequence (bernoulli-like noise)
  --linear                  L kmers from linear sequences of length L+k-1
                            (this is the default)
  --circular                L kmers from circular sequences of length L
                            (this is not currently supported)
  --nosort                  don't sort output
                            (by default output is sorted by nMutated)
  --stats=<filename>        write stats to a file
                            (by default stats are written to stderr)
  --seed=<string>           set random seed
  --progress=<number>       periodically report how many sequence pairs we've
                            tested

Conceptually, generate pairs of sequences of values in the unit interval and
report the distribution of the number of mutated kmers as well as other related
stats.

A 'sequence pair' is a random linear sequence and a mutated version of it. The
mutated version is consistent with a model where the sequence represents hash
values of kmers in a sequence of nucleotides, and the individual nucleotides
are subject to error, with the caveat that no duplications are possible.

This program doesn't actually generate sequence pairs. Instead, it generates
the positions of the mutations and derives the relevant statistics from those.
```

simulate_nucleotide_errors.py

```bash  
usage: cat fasta | simulate_nucleotide_errors [options]
  --k=<N>                   (K=) kmer size
                            (default is 28)
  --sketch=<N>              (S=) (cumulative) sketch size
                            (default is "no sketch")
  --sequences=<N>           (T=) number of mutated sequences to generate
                            (default is 1)
  --poisson=<probability>   (P=) (required) inject random sequencing errors
                            (substitutions); each base suffers a substitution
                            error with the given probability (poisson-like
                            noise)
  --bernoulli=<probability> (B= or E=) (required) inject random sequencing
                            errors (substitutions); exactly
                            round(L*<probability>) errors will occur in the
                            sequence (bernoulli-like noise)
  --linear                  L kmers from linear sequences of length L+k-1
                            (this is the default)
  --circular                L kmers from circular sequences of length L
  --nosort                  don't sort output
                            (by default output is sorted by nMutated)
  --stats=<filename>        write stats to a file
                            (by default stats are written to stderr)
  --mutated=<filename>      file to write the mutated sequences to
  --mutatedonly             just write out the mutated sequences and quit
  --seed=<string>           set random seed
  --hash=<int>              set seed for hash function (only used for sketches);
                            it is highly recommended that users specify the
                            hash seed
                            (default is a 'randomly' chosen hash seed)
  --hashbits=<int>          number of bits for hash function output; this can
                            be 16, 32, 64, 128, or "none"; "none" indicates
                            kmers should not be hashed
                            (default is no hashing)
  --progress=<number>       periodically report how many sequence pairs we've
                            tested

Repeatedly apply the specified mutation model to a single input sequence and
report the distribution of the number of mutated kmers as well as other
related stats.
```

evaluate_hypergeometric_slicer.py

```bash  
usage: cat <simulation_table> | evaluate_hypergeometric_slicer [options]
  --confidence=<p>          (C=) size of confidence interval; if present, this
                            overrides the value in the input file
                            (default is to get this from the input file)
  --slices=<N>              (m=) number of slices
                            (default is 100)
  --maxsketch=<N>           maximum sketch size; note that sketch sizes are
                            defined in the input file; this option causes us
                            to ignore large sketches
  --useL.A,B                use the column named "L.A,B" instead of L; input is
                            still required to contain an "L" column; the L.A,B
                            value is used as L in computing the confidence
                            interval but the original L is reported as L in the
                            output; L.A,B is reported as L in each record of
                            the details file
  --details=<filename>      write record-by-record details to a file
                            (by default we do not report these)
  --progress=<number>       periodically report how many input records we've
                            processed

typical input:
  #L      K  r     confidence q           nIntersection(s=100) nIntersection(s=500) ...
  4500000 21 0.100 0.95       0.890581011 7                    31                   ...
  4500000 21 0.100 0.95       0.890581011 4                    29                   ...
  4500000 21 0.100 0.95       0.890581011 7                    35                   ...
  4500000 21 0.100 0.95       0.890581011 5                    27                   ...
   ...

Columns L, k, r1, confidence, and q are required. At least one nIntersection
column is required. Sketch sizes are inferred fron the nIntersection column
headers.

nIntersection is the number of kmers (or kspans) that are in the intersection
of BS(A), BS(B), and BS(A union B) for a sketch of the given size. The Jaccard
estimate is nIntersection/s.
```


### Citation
If using this software, please cite
* Antonio Blanca, Robert S. Harris, David Koslicki and Paul Medvedev, "The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches", submitted 

