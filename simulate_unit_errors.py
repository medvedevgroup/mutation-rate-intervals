#!/usr/bin/env python3
"""
In sequence pairs chosen from the unit interval under a mutation model, compute
the number of mutated kmers and islands and associated stats.
"""

from sys          import argv,stdin,stdout,stderr,exit
from random       import Random,seed as random_seed,random as unit_random
from math         import sqrt,floor,ceil
from gzip         import open as gzip_open
from numpy.random import RandomState
from scipy.stats  import hypergeom
from hypergeometric_slicer \
                  import p_mutated,exp_n_mutated,var_n_mutated,estimate_r1_from_n_mutated, \
                         exp_n_island,var_n_island


def usage(s=None):
	message = """
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
the positions of the mutations and derives the relevant statistics from those."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	kmerSize           = 21
	kmerSequenceLength = 100
	sketchSizes        = None
	numSequences       = None
	noiseKind          = None
	pSubstitution      = None
	sequenceType       = "linear"
	sortBy             = "nMutated"
	statsFilename      = None
	prngSeed           = None
	reportProgress     = None
	debug              = []

	statsOfInterest = ["r1","k","L","trials","q",
		               "E[nMut].theory","StDev[nMut].theory",
		               "Mean[nMut].obs","StDev[nMut].obs","RMSE(StDev[nMut])",
		               "E[nIsl].theory","StDev[nIsl].theory",
		               "Mean[nIsl].obs","StDev[nIsl].obs","RMSE(StDev[nIsl])"]

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--kmer=")) or (arg.startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--set=")) or (arg.startswith("N=")) or (arg.startswith("L=")):
			kmerSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--sketch=")) or (arg.startswith("S=")):
			if (sketchSizes == None): sketchSizes = []
			sketchSizes += map(int_with_unit,argVal.split(","))
		elif (arg.startswith("--sequences=")) or (arg.startswith("T=")):
			numSequences = int_with_unit(argVal)
		elif (arg.startswith("--poisson=")) or (arg.startswith("--noise=")) or (arg.startswith("P=")):
			noiseKind = "poisson"
			pSubstitution = parse_probability(argVal)
		elif (arg.startswith("--bernoulli=")) or (arg.startswith("--error=")) or (arg.startswith("B=")) or (arg.startswith("E=")):
			noiseKind = "bernoulli"
			pSubstitution = parse_probability(argVal)
		elif (arg == "--linear"):
			sequenceType = "linear"
		elif (arg == "--circular"):
			sequenceType = "circular"
		elif (arg == "--nosort"):
			sortBy = None
		elif (arg.startswith("--stats=")):
			statsFilename = argVal
		elif (arg.startswith("--seed=")):
			prngSeed = argVal
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (pSubstitution == None):
		usage("you have to tell me the mutation probability")

	if (numSequences == None):
		numSequences = 1

	if (noiseKind == None):
		usage("you must specify either poisson or bernoulli error model")

	if (sequenceType == "circular"):
		# all the estimator code assumes linear sequences
		usage("circular sequences are not currently supported")

	if (sketchSizes != None):
		sketchSizes = list(set(sketchSizes))   # (remove duplicates)
		sketchSizes.sort()

	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			statsOfInterest += ["Mean[nIntersection(S=%d)].obs" % sketchSize,
			                    "Mean[Jaccard(S=%d)].obs" % sketchSize,
			                    "StDev[Jaccard(S=%d)].obs" % sketchSize]

	# set up randomness
	#
	# we use separate PRNGs for errors and sketches so that we'll get the same
	# simulation results for errors whether we do sketches or not (as long as
	# the user specifies a seed), and for a given sketch size, we'll get the
	# same simulation results regardless of what other sketch sizes we do

	sketchPrngs = None
	if (prngSeed != None):
		random_seed(prngSeed.encode("utf-8"))
		if ("prng" in debug):
			print("prng = %s" % prngSeed,file=stderr)

	if (sketchSizes != None):
		sketchPrngs = {}
		for sketchSize in sketchSizes:
			sketchPrngSeed = "%s.%s" % (prngSeed,sketchSize)
			sketchPrng = RandomState(seed=list(sketchPrngSeed.encode("utf-8")))
			sketchPrngs[sketchSize] = sketchPrng
			if ("prng" in debug):
				print("prng[S=%d] = %s" % (sketchSize,sketchPrngSeed),file=stderr)

	# set up model/generator

	if (noiseKind == "poisson") and (sequenceType == "linear"):
		mutationModel = PoissonModel \
		                  (kmerSequenceLength+kmerSize-1,kmerSize,pSubstitution,
		                   count_mutated_kmers_linear_naive if ("naive" in debug) else count_mutated_kmers_linear,
		                   count_islands_linear)
	elif (noiseKind == "bernoulli") and (sequenceType == "linear"):
		mutationModel = BernoulliModel \
		                  (kmerSequenceLength+kmerSize-1,kmerSize,pSubstitution,
		                   count_mutated_kmers_linear_naive if ("naive" in debug) else count_mutated_kmers_linear,
		                   count_islands_linear)
	elif (noiseKind == "poisson") and (sequenceType == "circular"):
		mutationModel = PoissonModel \
		                  (kmerSequenceLength,kmerSize,pSubstitution,
		                   count_mutated_kmers_circular_naive if ("naive" in debug) else count_mutated_kmers_circular,
		                   count_islands_circular)
	elif (noiseKind == "bernoulli") and (sequenceType == "circular"):
		mutationModel = BernoulliModel \
		                  (kmerSequenceLength,kmerSize,pSubstitution,
		                   count_mutated_kmers_circular_naive if ("naive" in debug) else count_mutated_kmers_circular,
		                   count_islands_circular)
	else:
		assert (False), "internal error"

	# generate sequences and collect stats

	q = p_mutated(kmerSize,pSubstitution)

	nErrorsObserved       = []
	nMutatedObserved      = []
	nIslandObserved       = []
	r1EstNMutatedObserved = []

	if (sketchSizes != None):
		nIntersectionObserved = {}
		jaccardObserved       = {}
		for sketchSize in sketchSizes:
			nIntersectionObserved[sketchSize] = []
			jaccardObserved[sketchSize]       = []
			mutationModel.set_sketch_prng(sketchSize,sketchPrngs[sketchSize])

	for seqNum in range(numSequences):
		if (reportProgress != None):
			if (1+seqNum <= 2) or ((1+seqNum) % reportProgress == 0):
				print("testing sequence %d" % (1+seqNum),file=stderr)

		# generate a (conceptual) sequence pair and collect stats

		mutationModel.generate()
		(nErrors,nMutated,nIsland) = mutationModel.count()
		nErrorsObserved  += [nErrors]
		nMutatedObserved += [nMutated]
		nIslandObserved  += [nIsland]

		r1EstNMutated = estimate_r1_from_n_mutated(kmerSequenceLength,kmerSize,nMutated)
		r1EstNMutatedObserved += [r1EstNMutated]

		# generate (conceptual) sketches and collect basic stats

		if (sketchSizes != None):
			for sketchSize in sketchSizes:
				nIntersection = mutationModel.simulate_sketch(kmerSequenceLength,nMutated,sketchSize)
				nIntersectionObserved[sketchSize] += [nIntersection]
				jaccardObserved[sketchSize]       += [float(nIntersection)/sketchSize]

	# report per-trial results

	if (sortBy == "nMutated"):
		order = [(nMutatedObserved[ix],ix) for ix in range(numSequences)]
		order.sort()
		order = [ix for (_,ix) in order]
	else: # if (sortBy == None):
		order = list(range(numSequences))

	header = ["L","K","r","trial","q","nErr","nMut","nIsl","r1est.nMut"]
	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			header += ["nIntersection(s=%d)" % sketchSize]
			header += ["j.est(nMut,s=%d)" % sketchSize]
	print("#%s" % "\t".join(header))

	for ix in range(numSequences):
		line = "\t".join(["%d","%d","%0.3f","%d","%0.9f","%d","%d","%d","%0.9f"]) \
		     % (kmerSequenceLength,
		        kmerSize,
		        pSubstitution,
		        1+order[ix],
		        q,
		        nErrorsObserved[order[ix]],
		        nMutatedObserved[order[ix]],
		        nIslandObserved[order[ix]],
		        r1EstNMutatedObserved[order[ix]])
		if (sketchSizes != None):
			for sketchSize in sketchSizes:
				line += "\t%d"    % nIntersectionObserved[sketchSize][order[ix]]
				line += "\t%0.9f" % jaccardObserved[sketchSize][order[ix]]
		print(line)

	# compute stats

	nMutatedMean      = sample_mean(nMutatedObserved)
	nMutatedStDev     = sqrt(sample_variance(nMutatedObserved))
	predNMutatedMean  = exp_n_mutated(kmerSequenceLength,kmerSize,pSubstitution)
	predNMutatedStDev = sqrt(var_n_mutated(kmerSequenceLength,kmerSize,pSubstitution))
	rmseNMutatedStDev = abs(nMutatedStDev-predNMutatedStDev)

	nIslandMean       = sample_mean(nIslandObserved)
	nIslandStDev      = sqrt(sample_variance(nIslandObserved))
	predNIslandMean   = exp_n_island(kmerSequenceLength,kmerSize,pSubstitution)
	predNIslandStDev  = sqrt(var_n_island(kmerSequenceLength,kmerSize,pSubstitution))
	rmseNIslandStDev  = abs(nIslandStDev-predNIslandStDev)

	if (sketchSizes != None):
		nIntersectionMean = {}
		jaccardEstMean = {}
		jaccardEstStDev = {}
		for sketchSize in sketchSizes:
			nIntersectionMean[sketchSize] = sample_mean(nIntersectionObserved[sketchSize])
			jaccardEstMean   [sketchSize] = sample_mean(jaccardObserved[sketchSize])
			jaccardEstStDev  [sketchSize] = sqrt(sample_variance(jaccardObserved[sketchSize]))

	# report stats

	statToText = {}
	statToText["r1"]                        = "%0.3f" % pSubstitution
	statToText["k"]                         = "%d"    % kmerSize
	statToText["L"]                         = "%d"    % kmerSequenceLength
	statToText["trials"]                    = "%d"    % numSequences
	statToText["q"]                         = "%0.9f" % q
	statToText["E[nMut].theory"]            = "%0.9f" % predNMutatedMean
	statToText["StDev[nMut].theory"]        = "%0.9f" % predNMutatedStDev
	statToText["Mean[nMut].obs"]            = "%0.9f" % nMutatedMean
	statToText["StDev[nMut].obs"]           = "%0.9f" % nMutatedStDev
	statToText["RMSE(StDev[nMut])"]         = "%0.9f" % rmseNMutatedStDev
	statToText["E[nIsl].theory"]            = "%0.9f" % predNIslandMean
	statToText["StDev[nIsl].theory"]        = "%0.9f" % predNIslandStDev
	statToText["Mean[nIsl].obs"]            = "%0.9f" % nIslandMean
	statToText["StDev[nIsl].obs"]           = "%0.9f" % nIslandStDev
	statToText["RMSE(StDev[nIsl])"]         = "%0.9f" % rmseNIslandStDev

	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			statToText["Mean[nIntersection(S=%d)].obs" % sketchSize] = "%0.9f" % nIntersectionMean[sketchSize]
			statToText["Mean[Jaccard(S=%d)].obs"       % sketchSize] = "%0.9f" % jaccardEstMean[sketchSize]
			statToText["StDev[Jaccard(S=%d)].obs"      % sketchSize] = "%0.9f" % jaccardEstStDev[sketchSize]

	if (statsFilename != None):
		if (statsFilename.endswith(".gz")) or (statsFilename.endswith(".gzip")):
			statsF = gzip_open(statsFilename,"wt")
		else:
			statsF = open(statsFilename,"wt")

		print("#%s" % "\t".join(statsOfInterest),file=statsF)
		statsLine = [statToText[stat] for stat in statsOfInterest]
		print("\t".join(statsLine),file=statsF)
		statsF.close()
	else:
		statW = max(len(stat) for stat in statsOfInterest)
		for stat in statsOfInterest:
			print("%*s = %s" % (statW,stat,statToText[stat]),file=stderr)


# PoissonModel--
#	Generate a sequence of poisson-type errors, and report the number of
#	errors, mutated kmers, and islands.

class PoissonModel(object):

	def __init__(self,ntSequenceLength,kmerSize,pSubstitution,mutatedKmerCounter,islandCounter):
		self.ntSequenceLength   = ntSequenceLength
		self.kmerSize           = kmerSize
		self.pSubstitution      = pSubstitution
		self.mutatedKmerCounter = mutatedKmerCounter
		self.islandCounter      = islandCounter
		self.sketchPrngs        = {}

	def set_sketch_prng(self,sketchSize,prng=None):
		self.sketchPrngs[sketchSize] = prng

	def count(self,regenerate=False):
		if (regenerate): self.generate()
		return (sum(self.errorSeq),
		        self.mutatedKmerCounter(self.kmerSize,self.errorSeq),
		        self.islandCounter(self.kmerSize,self.errorSeq))

	def generate(self):
		self.errorSeq = list(map(lambda _:1 if (unit_random()<self.pSubstitution) else 0,range(self.ntSequenceLength)))
		return self.errorSeq

	def simulate_sketch(self,kmerSequenceLength,nMutated,sketchSize):
		if not (0 < sketchSize <= kmerSequenceLength): raise ValueError
		prng = self.sketchPrngs[sketchSize] if (sketchSize in self.sketchPrngs) else None
		# Given sequence length L and N mutated kmers, we consider the kmers
		# in A union B to be numbered from 0 to L+N-1, and we consider the
		# *un*mutated kmers to be the first L-N of these
		#    <---unmutated--> <---mutated, in A--> <---mutated, in B-->
		#   +----------------+--------------------+--------------------+
		#   | 0        L-1-N | L-N            L-1 | L            L+N-1 |
		#   +----------------+--------------------+--------------------+
		# The L-N *un*mutated kmers are A intersection B. The hash function
		# would effectively choose a random set of s of all L+N kmers as bottom
		# sketch BS(A union B), where s is the sketch size. So conceptually, we
		# have an urn with L+N balls, s of which are 'red'. We draw L-N balls
		# and want to know how many are red. This is the size of the
		# intersection of BS(A), BS(B), and BS(A union B).
		L = kmerSequenceLength
		N = nMutated
		s = sketchSize
		if (N == L):    # hypergeom.rvs doesn't handle this case, a
			return 0    # .. case that seems perfectly legitimate
		nIntersection = hypergeom.rvs(L+N,s,L-N,random_state=prng)
		return nIntersection


# BernoulliModel--
#	Generate a sequence of bernoulli-type errors, and report the number of
#	errors, mutated kmers, and islands.

class BernoulliModel(object):

	def __init__(self,ntSequenceLength,kmerSize,pSubstitution,mutatedKmerCounter,islandCounter):
		self.ntSequenceLength   = ntSequenceLength
		self.kmerSize           = kmerSize
		self.pSubstitution      = pSubstitution
		self.mutatedKmerCounter = mutatedKmerCounter
		self.islandCounter      = islandCounter

	def generate(self):
		self.errorSeq = []
		nErrors = remainingErrors = round(pSubstitution * ntSequenceLength)
		for remainingSeqLen in range(ntSequenceLength,0,-1):
			if (remainingSeqLen*unit_random() >= remainingErrors):
				self.errorSeq += [0]
			else:
				self.errorSeq += [1]
				remainingErrors -= 1
		return self.errorSeq


# count_mutated_kmers_linear--
#	Given a sequence of errors, report the number of errors and mutated
#	kmers, treating the sequence as linear. The error sequence is a list of
#	1 (error) and 0 (non-error).
#
#	The algorithm maintains a sum of errors over a sliding window of length K.
#	wherever this sum is non-zero, the corresponding kmer is 'mutated'.

def count_mutated_kmers_linear(kmerSize,errorSeq):
	nMutated = 0
	errorsInKmer = 0
	for pos in range(len(errorSeq)+1):
		# invariant: errorsInKmer is sum of errorSeq pos-kmerSize through pos-1
		if (pos >= kmerSize):
			if (errorsInKmer > 0): nMutated += 1   # pos-kmerSize is 'mutated'
			errorsInKmer -= errorSeq[pos-kmerSize]
		if (pos == len(errorSeq)): break
		errorsInKmer += errorSeq[pos]
	if ("mutated" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nMutated),file=stderr)
	return nMutated


def count_mutated_kmers_linear_naive(kmerSize,errorSeq):
	seqLen = len(errorSeq) - (kmerSize-1)
	nMutated = 0
	for pos in range(seqLen):
		errorsInKmer = sum(errorSeq[pos:pos+kmerSize])
		if (errorsInKmer > 0): nMutated += 1       # pos is 'mutated'
	if ("mutated" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nMutated),file=stderr)
	return nMutated


# count_mutated_kmers_circular--
#	Given a sequence of errors, report the number of errors and mutated
#	kmers, treating the sequence as circular. The error sequence is a list of
#	1 (error) and 0 (non-error).
#
#	The algorithm maintains a sum of errors over a sliding window of length K.
#	wherever this sum is non-zero, the corresponding kmer is 'mutated'.

def count_mutated_kmers_circular(kmerSize,errorSeq):
	seqLen = len(errorSeq)
	nMutated = 0
	errorsInKmer = 0
	for pos in range(seqLen+kmerSize):
		# invariant: errorsInKmer is sum of errorSeq pos-kmerSize through pos-1
		if (pos >= kmerSize):
			if (errorsInKmer > 0): nMutated += 1   # pos-kmerSize is 'mutated'
			errorsInKmer -= errorSeq[pos-kmerSize]
		errorsInKmer += errorSeq[pos] if (pos < seqLen) else errorSeq[pos-seqLen]
	if ("mutated" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nMutated),file=stderr)
	return nMutated


def count_mutated_kmers_circular_naive(kmerSize,errorSeq):
	seqLen = len(errorSeq)
	extendedSeq = errorSeq + errorSeq[:kmerSize-1]

	nMutated = 0
	for pos in range(seqLen):
		errorsInKmer = sum(extendedSeq[pos:pos+kmerSize])
		if (errorsInKmer > 0): nMutated += 1       # pos is 'mutated'
	if ("mutated" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nMutated),file=stderr)
	return nMutated


# count_islands_linear--
#	Given a sequence of errors, report the number of islands, treating the
#	sequence as linear. The error sequence is a list of 1 (error) and 0
#	(non-error).
#
#	An ocean is a maximal run of k or more non-errors. An island is a maximal
#	segment containing at least 1 error and no ocean. If we assume a ghost
#	ocean to the left of the sequence start, every island has a first error
#	that immediately follows an ocean.

def count_islands_linear(kmerSize,errorSeq):
	assert (len(errorSeq) >= kmerSize)
	nIsland = 0
	errorFreeRunLength = kmerSize  # assume ghost ocean before start of sequence
	for pos in range(len(errorSeq)):
		if (errorSeq[pos] == 0):
			errorFreeRunLength += 1
		else: # if (errorSeq[pos] == 1):
			if (errorFreeRunLength >= kmerSize): nIsland += 1
			errorFreeRunLength = 0
	if ("islands" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nIsland),file=stderr)
	return nIsland


# count_islands_circular--
#	Given a sequence of errors, report the number of islands, treating the
#	sequence as circular. The error sequence is a list of 1 (error) and 0
#	(non-error).
#
#	See count_islands_linear for definition of oceans and islands. Usually
#	every island has a first error that immediately follows an ocean, but we
#	must carefully consider an ocean that wraps from the end of the sequence
#	to the start. And we also have to account for a special case of a sequence
#	that is entirely an island.

def count_islands_circular(kmerSize,errorSeq):
	assert (len(errorSeq) >= kmerSize)
	# count the length of any error free run at the end of the sequence; we
	# needn't count beyond the length that would make the suffix an ocean
	errorFreeRunLength = 0
	for pos in range(len(errorSeq)-1,-1,-1):
		if (errorSeq[pos] == 1):
			break
		else: # if (errorSeq[pos] == 0):
			errorFreeRunLength += 1
			if (errorFreeRunLength >= kmerSize):
				break

	# count islands by looking for a first error following any ocean;  we
	# include the count above in the count of any error-free run prior to the
	# first error we encounter, thus if an ocean wraps from the end to the
	# start, we properly detect it; if the sequence ends with an ocean and
	# starts with an island, we detect that because the run length counter will
	# start with k; if an island wraps from the end to the start, we will
	# eventually count it when we encounter the first error following an ocean;
	# but if there is no ocean the island won't get counted (see special case
	# below)
	nIsland = 0
	hasAnyErrors = False
	for pos in range(len(errorSeq)):
		if (errorSeq[pos] == 0):
			errorFreeRunLength += 1
		else: # if (errorSeq[pos] == 1):
			if (errorFreeRunLength >= kmerSize): nIsland += 1
			hasAnyErrors = True
			errorFreeRunLength = 0
	# check for the all-island special case
	if (hasAnyErrors) and (nIsland == 0):
		nIsland = 1
	if ("islands" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nIsland),file=stderr)
	return nIsland


# mean, variance, mean_squared_error--

def sample_mean(observed):
	return float(sum(observed)) / len(observed)

def sample_variance(observed):
	if (len(observed) <= 1): return 0.0
	m = sample_mean(observed)
	return float(sum([(n-m)**2 for n in observed])) / (len(observed)-1)

def mean_squared_error(observed,predicted):
	return float(sum([(n-predicted)**2 for n in observed])) / len(observed)


# parse_probability--
#	Parse a string as a probability

def parse_probability(s,strict=True):
	scale = 1.0
	if (s.endswith("%")):
		scale = 0.01
		s = s[:-1]

	try:
		p = float(s)
	except:
		try:
			(numer,denom) = s.split("/",1)
			p = float(numer)/float(denom)
		except:
			raise ValueError

	p *= scale

	if (strict) and (not 0.0 <= p <= 1.0):
		raise ValueError

	return p


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
