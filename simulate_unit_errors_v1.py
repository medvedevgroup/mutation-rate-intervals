#!/usr/bin/env python3
"""
In sequence pairs chosen from the unit interval under a mutation model, compute
the number of affected kmers and islands.
"""

from sys          import argv,stdin,stdout,stderr,exit
from random       import Random,seed as random_seed,random as unit_random
from math         import sqrt,floor,ceil
from gzip         import open as gzip_open
from numpy.random import RandomState
from scipy.stats  import hypergeom
from kmer_mutation_formulas_v1 \
                  import p_affected,exp_n_affected,var_n_affected,estimate_r1_from_n_affected, \
                         confidence_interval_r1_from_n_affected,in_confidence_interval_q_from_n_affected, \
                         exp_n_island,var_n_island,estimate_r1_from_n_island,impossible_n_island, \
                         confidence_interval_r1_from_n_island,in_confidence_interval_q_from_n_island, \
                         confidence_interval_jaccard_from_sketch_n_affected, \
                         in_confidence_interval_jaccard_from_sketch_n_affected


def usage(s=None):
	message = """
usage: simulate_unit_errors [options]
  --k=<N>                   (K=) kmer size
                            (default is 28)
  --n=<N>                   (N= or L=) sequence length
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
  --confidence=<p>          (C=) size of confidence interval
                            (default is 99%)
  --alpha1=<p>              (cumulative) alpha', which must be s.t.
                            1-confidence<alpha'<1
  --alpha1=alpha2           choose alpha' s.t. alpha''=alpha'
                            (this is the default)
  --noinverse               for confidence interval tests, do NOT use inverse
                            functions
                            (by default inverse functions are used)
  --nosort                  don't sort output
                            (by default output is sorted by nAffected)
  --stats=<filename>        write stats to a file
                            (by default stats are written to stderr)
  --seed=<string>           set random seed
  --progress=<number>       periodically report how many sequence pairs we've
                            tested

Conceptually, generate pairs of sequences in the unit interval and report the
distribution of the number of affected kmers as well as other related stats.

A 'sequence pair' is a random circular sequence and a mutated version of it.
The mutated version is consistent with a model where the sequence represents
hash values of kmers in a sequence of nucleotides, and the individual
nucleotides are subject to error, with the caveat that no duplications are
possible.

This program doesn't actually generate sequence pairs. Instead, it generates
the positions of the mutations."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	kmerSize           = 28
	kmerSequenceLength = 100
	sketchSizes        = None
	numSequences       = None
	noiseKind          = None
	pSubstitution      = None
	sequenceType       = "linear"
	confidence         = 0.99
	alpha1s            = None
	ciUseInverse       = True
	sortBy             = "nAffected"
	statsFilename      = None
	prngSeed           = None
	reportProgress     = None
	debug              = []

	statsOfInterest = ["r1","k","L","confidence","trials","q",
		               "E[nAff].theory","StDev[nAff].theory",
		               "inCI(r1est.nAff).obs",
		               "Mean[nAff].obs","StDev[nAff].obs","RMSE(StDev[nAff])","RMSE(r1est.nAff)",
		               "E[nIsl].theory","StDev[nIsl].theory",
		               "inCI(r1est.nIsl).obs",
		               "Mean[nIsl].obs","StDev[nIsl].obs","RMSE(StDev[nIsl])","RMSE(r1est.nIsl)",
		               "r1est.nIsl.impossible"]

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
		elif (arg.startswith("--confidence=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--alpha1=alpha2")):
			if (alpha1s == None): alpha1s = []
			alpha1s += [None]
		elif (arg.startswith("--alpha1=")):
			if (alpha1s == None): alpha1s = []
			for alpha1 in argVal.split(","):
				if (alpha1 == "alpha2"):
					alpha1s += ["alpha2"]
				else:
					alpha1s += [parse_probability(alpha1)]
		elif (arg == "--noinverse"):
			ciUseInverse = False
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

	if (sequenceType == "circular"):
		# all the estimator code assumes linear sequences
		usage("circular sequences are not currently supported")

	if (alpha1s != None):
		if (sketchSizes == None):
			print("WARNING: alpha1 is irrelevant without sketch sizes; ignoring it")
		else:
			for alpha1 in alpha1s:
				if (alpha1 == "alpha2"): continue
				if (0 < alpha1 < 1-confidence ): continue
				usage("alpha1=%s is not between 0 and 1-confidence (%s))" % (alpha1,1-confidence))
	if (alpha1s == None):
		alpha1s = ["alpha2"]
	if ("alpha2") in alpha1s:
		ix = alpha1s.index("alpha2")
		alpha1s[ix] = 1-sqrt(confidence)
	alpha1s.sort()

	if (sketchSizes != None):
		sketchSizes = list(set(sketchSizes))   # (remove duplicates)
		sketchSizes.sort()

	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			statsOfInterest += ["Mean[nIntersection(S=%d)].obs" % sketchSize,
			                    "Mean[Jaccard(S=%d)].obs" % sketchSize,
			                    "StDev[Jaccard(S=%d)].obs" % sketchSize]
			for alpha1 in alpha1s:
				statsOfInterest += ["inCI(jest.nAff,S=%d,alpha'=%.6f).obs" % (sketchSize,alpha1)]

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
		                   count_affected_kmers_linear_naive if ("naive" in debug) else count_affected_kmers_linear,
		                   count_islands_linear)
	elif (noiseKind == "bernoulli") and (sequenceType == "linear"):
		mutationModel = BernoulliModel \
		                  (kmerSequenceLength+kmerSize-1,kmerSize,pSubstitution,
		                   count_affected_kmers_linear_naive if ("naive" in debug) else count_affected_kmers_linear,
		                   count_islands_linear)
	elif (noiseKind == "poisson") and (sequenceType == "circular"):
		mutationModel = PoissonModel \
		                  (kmerSequenceLength,kmerSize,pSubstitution,
		                   count_affected_kmers_circular_naive if ("naive" in debug) else count_affected_kmers_circular,
		                   count_islands_circular)
	elif (noiseKind == "bernoulli") and (sequenceType == "circular"):
		mutationModel = BernoulliModel \
		                  (kmerSequenceLength,kmerSize,pSubstitution,
		                   count_affected_kmers_circular_naive if ("naive" in debug) else count_affected_kmers_circular,
		                   count_islands_circular)
	else:
		assert (False), "internal error"

	# generate sequences and collect stats

	nErrorsObserved        = []
	nAffectedObserved      = []
	nIslandObserved        = []
	r1EstNAffectedObserved = []
	r1EstNIslandObserved   = []
	numImpossibleNislands  = 0   # counts when nIsland can't be achieved with any r1

	if (sketchSizes != None):
		nIntersectionObserved = {}
		jaccardObserved       = {}
		for sketchSize in sketchSizes:
			nIntersectionObserved[sketchSize] = []
			jaccardObserved[sketchSize]       = []
			mutationModel.set_sketch_prng(sketchSize,sketchPrngs[sketchSize])

	for seqNum in range(numSequences):
		if (reportProgress != None):
			if (1+seqNum == 1) or ((1+seqNum) % reportProgress == 0):
				print("testing sequence %d" % (1+seqNum),file=stderr)

		# generate a (conceptual) sequence pair and collect stats

		mutationModel.generate()
		(nErrors,nAffected,nIsland) = mutationModel.count()
		nErrorsObserved   += [nErrors]
		nAffectedObserved += [nAffected]
		nIslandObserved   += [nIsland]

		r1EstNAffected = estimate_r1_from_n_affected(kmerSequenceLength,kmerSize,nAffected)
		r1EstNAffectedObserved += [r1EstNAffected]

		# note: when r1 is estimated from nIsland, there can be more than one
		# solution; we take the solution that is closest to r1 estimated from
		# nAffected

		r1EstNIslandList = estimate_r1_from_n_island(kmerSequenceLength,kmerSize,nIsland)
		if (len(r1EstNIslandList) == 0):
			r1EstNIsland = float("nan")
		elif (len(r1EstNIslandList) == 1):
			r1EstNIsland = r1EstNIslandList[0]
		else:
			r1EstNIsland = None
			for r1Est in r1EstNIslandList:
				diff = abs(r1Est-r1EstNAffected)
				if (r1EstNIsland == None) or (diff < bestDiff):
					r1EstNIsland = r1Est
					bestDiff = diff
		r1EstNIslandObserved += [r1EstNIsland]

		if (impossible_n_island(kmerSequenceLength,kmerSize,nIsland)):
			numImpossibleNislands += 1

		# generate (conceptual) sketches and collect basic stats

		if (sketchSizes != None):
			for sketchSize in sketchSizes:
				nIntersection = mutationModel.simulate_sketch(kmerSequenceLength,nAffected,sketchSize)
				nIntersectionObserved[sketchSize] += [nIntersection]
				jaccardObserved[sketchSize]       += [float(nIntersection)/sketchSize]

	# report per-trial results

	if (sortBy == "nAffected"):
		order = [(nAffectedObserved[ix],ix) for ix in range(numSequences)]
		order.sort()
		order = [ix for (_,ix) in order]
	else: # if (sortBy == None):
		order = list(range(numSequences))

	header = ["L","K","r","trial","nErr","nAff","nIsl","r1est.nAff","r1.est.nIsl"]
	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			header += ["nIntersection(s=%d)" % sketchSize]
			header += ["j.est(nAff,s=%d)" % sketchSize]
	print("#%s" % "\t".join(header))

	for ix in range(numSequences):
		line = "\t".join(["%d","%d","%0.3f","%d","%d","%d","%d","%0.9f","%0.9f"]) \
		     % (kmerSequenceLength,
		        kmerSize,
		        pSubstitution,
		        1+order[ix],
		        nErrorsObserved[order[ix]],
		        nAffectedObserved[order[ix]],
		        nIslandObserved[order[ix]],
		        r1EstNAffectedObserved[order[ix]],
		        r1EstNIslandObserved[order[ix]])
		if (sketchSizes != None):
			for sketchSize in sketchSizes:
				line += "\t%d"    % nIntersectionObserved[sketchSize][order[ix]]
				line += "\t%0.9f" % jaccardObserved[sketchSize][order[ix]]
		print(line)

	# compute stats

	alpha = 1 - confidence
	q = p_affected(kmerSize,pSubstitution)

	nAffectedMean      = sample_mean(nAffectedObserved)
	nAffectedStDev     = sqrt(sample_variance(nAffectedObserved))
	predNAffectedMean  = exp_n_affected(kmerSequenceLength,kmerSize,pSubstitution)
	predNAffectedStDev = sqrt(var_n_affected(kmerSequenceLength,kmerSize,pSubstitution))
	rmseNAffectedStDev = abs(nAffectedStDev-predNAffectedStDev)

	nIslandMean        = sample_mean(nIslandObserved)
	nIslandStDev       = sqrt(sample_variance(nIslandObserved))
	predNIslandMean    = exp_n_island(kmerSequenceLength,kmerSize,pSubstitution)
	predNIslandStDev   = sqrt(var_n_island(kmerSequenceLength,kmerSize,pSubstitution))
	rmseNIslandStDev   = abs(nIslandStDev-predNIslandStDev)

	rmseR1EstNAffected = sqrt(mean_squared_error(r1EstNAffectedObserved,pSubstitution))
	rmseR1EstNIsland   = sqrt(mean_squared_error(r1EstNIslandObserved,pSubstitution))

	(predR1EstNAffectedLow,predR1EstNAffectedHigh) \
	                   = confidence_interval_r1_from_n_affected(kmerSequenceLength,kmerSize,pSubstitution,alpha)
	inConfR1EstNAffected \
	                   = in_confidence_interval_q_from_n_affected(kmerSequenceLength,kmerSize,pSubstitution,alpha,
	                                                              nAffectedObserved,useInverse=ciUseInverse)

	(predR1EstNIslandLow,predR1EstNIslandHigh) \
	                   = confidence_interval_r1_from_n_island(kmerSequenceLength,kmerSize,pSubstitution,alpha)
	inConfR1EstNIsland = in_confidence_interval_q_from_n_island(kmerSequenceLength,kmerSize,pSubstitution,alpha,
	                                                            nIslandObserved,nAffectedObserved,useInverse=ciUseInverse)

	if (sketchSizes != None):
		nIntersectionMean = {}
		jaccardEstMean = {}
		jaccardEstStDev = {}
		predJaccardEstNAffectedLow = {}
		predJaccardEstNAffectedHigh = {}
		inConfJaccardEstNAffected = {}
		for sketchSize in sketchSizes:
			nIntersectionMean          [sketchSize] = sample_mean(nIntersectionObserved[sketchSize])
			jaccardEstMean             [sketchSize] = sample_mean(jaccardObserved[sketchSize])
			jaccardEstStDev            [sketchSize] = sqrt(sample_variance(jaccardObserved[sketchSize]))
			(predJaccardEstNAffectedLow[sketchSize],predJaccardEstNAffectedHigh[sketchSize])\
			                                        = confidence_interval_jaccard_from_sketch_n_affected \
			                                            (kmerSequenceLength,kmerSize,sketchSize,pSubstitution,
			                                             alpha,alpha1)
			for (alpha1Ix,alpha1) in enumerate(alpha1s):
				inConfJaccardEstNAffected[(sketchSize,alpha1Ix)] \
				                                    = in_confidence_interval_jaccard_from_sketch_n_affected \
				                                        (kmerSequenceLength,kmerSize,sketchSize,pSubstitution,
				                                         alpha,alpha1,
			                                             jaccardObserved[sketchSize],useInverse=ciUseInverse)

	# report stats

	statToText = {}
	statToText["r1"]                        = "%0.3f" % pSubstitution
	statToText["k"]                         = "%d"    % kmerSize
	statToText["L"]                         = "%d"    % kmerSequenceLength
	statToText["confidence"]                = "%0.3f" % confidence
	statToText["trials"]                    = "%d"    % numSequences
	statToText["q"]                         = "%0.9f" % q
	statToText["E[nAff].theory"]            = "%0.9f" % predNAffectedMean
	statToText["StDev[nAff].theory"]        = "%0.9f" % predNAffectedStDev
	statToText["CIlow(r1est.nAff).theory"]  = "%0.9f" % predR1EstNAffectedLow
	statToText["CIhigh(r1est.nAff).theory"] = "%0.9f" % predR1EstNAffectedHigh
	statToText["inCI(r1est.nAff).obs"]      = "%0.9f" % (float(inConfR1EstNAffected) / numSequences)
	statToText["Mean[nAff].obs"]            = "%0.9f" % nAffectedMean
	statToText["StDev[nAff].obs"]           = "%0.9f" % nAffectedStDev
	statToText["RMSE(StDev[nAff])"]         = "%0.9f" % rmseNAffectedStDev
	statToText["RMSE(r1est.nAff)"]          = "%0.9f" % rmseR1EstNAffected
	statToText["E[nIsl].theory"]            = "%0.9f" % predNIslandMean
	statToText["StDev[nIsl].theory"]        = "%0.9f" % predNIslandStDev
	statToText["CIlow(r1est.nIsl).theory"]  = "%0.9f" % predR1EstNIslandLow
	statToText["CIhigh(r1est.nIsl).theory"] = "%0.9f" % predR1EstNIslandHigh
	statToText["inCI(r1est.nIsl).obs"]      = "%0.9f" % (float(inConfR1EstNIsland) / numSequences)
	statToText["Mean[nIsl].obs"]            = "%0.9f" % nIslandMean
	statToText["StDev[nIsl].obs"]           = "%0.9f" % nIslandStDev
	statToText["RMSE(StDev[nIsl])"]         = "%0.9f" % rmseNIslandStDev
	statToText["RMSE(r1est.nIsl)"]          = "%0.9f" % rmseR1EstNIsland
	statToText["r1est.nIsl.impossible"]     = "%0.9f" % (float(numImpossibleNislands)/numSequences)

	if (sketchSizes != None):
		for sketchSize in sketchSizes:
			statToText["Mean[nIntersection(S=%d)].obs" % sketchSize] = "%0.9f" % nIntersectionMean[sketchSize]
			statToText["Mean[Jaccard(S=%d)].obs"       % sketchSize] = "%0.9f" % jaccardEstMean[sketchSize]
			statToText["StDev[Jaccard(S=%d)].obs"      % sketchSize] = "%0.9f" % jaccardEstStDev[sketchSize]
			statToText["CIlow(jest.nAffS=%d).theory"   % sketchSize] = "%0.9f" % predJaccardEstNAffectedLow[sketchSize]
			statToText["CIhigh(jest.nAffS=%d).theory"  % sketchSize] = "%0.9f" % predJaccardEstNAffectedHigh[sketchSize]
			for (alpha1Ix,alpha1) in enumerate(alpha1s):
				statToText["inCI(jest.nAff,S=%d,alpha'=%.6f).obs" % (sketchSize,alpha1)] \
				                                                     = "%0.9f" % (float(inConfJaccardEstNAffected[(sketchSize,alpha1Ix)]) / numSequences)

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
			print("%*s = %s" % (stat,statW,statToText[stat]),file=stderr)


# PoissonModel--
#	Generate a sequence of poisson-type errors, and report the number of
#	errors, affected kmers, and islands.

class PoissonModel(object):

	def __init__(self,ntSequenceLength,kmerSize,pSubstitution,affectedKmerCounter,islandCounter):
		self.ntSequenceLength    = ntSequenceLength
		self.kmerSize            = kmerSize
		self.pSubstitution       = pSubstitution
		self.affectedKmerCounter = affectedKmerCounter
		self.islandCounter       = islandCounter
		self.sketchPrngs         = {}

	def set_sketch_prng(self,sketchSize,prng=None):
		self.sketchPrngs[sketchSize] = prng

	def count(self,regenerate=False):
		if (regenerate): self.generate()
		return (sum(self.errorSeq),
		        self.affectedKmerCounter(self.kmerSize,self.errorSeq),
		        self.islandCounter(self.kmerSize,self.errorSeq))

	def generate(self):
		self.errorSeq = list(map(lambda _:1 if (unit_random()<self.pSubstitution) else 0,range(self.ntSequenceLength)))
		return self.errorSeq

	def simulate_sketch(self,kmerSequenceLength,nAffected,sketchSize):
		if not (0 < 2*sketchSize < kmerSequenceLength): raise ValueError
		prng = self.sketchPrngs[sketchSize] if (sketchSize in self.sketchPrngs) else None
		# Given sequence length L and N affected kmers, we consider the kmers
		# in A union B to be numbered from 0 to L+N-1, and we consider the
		# *un*affected kmers to be the first L-N of these
		#    <--unaffected--> <--affected, in A--> <--affected, in B-->
		#   +----------------+--------------------+--------------------+
		#   | 0        L-1-N | L-N            L-1 | L            L+N-1 |
		#   +----------------+--------------------+--------------------+
		# The L-N *un*affected kmers are A intersection B. The hash function
		# would effectively choose a random set of s of all L+N kmers as bottom
		# sketch BS(A union B), where s is the sketch size. So conceptually, we
		# have an urn with L+N balls, s of which are 'red'. We draw L-N balls
		# and want to know how many are red. This is the size of the
		# intersection of BS(A), BS(B), and BS(A union B).
		L = kmerSequenceLength
		N = nAffected
		s = sketchSize
		if (N == L):    # hypergeom.rvs doesn't handle this case, a
			return 0    # .. case that seems perfectly legitimate
		nIntersection = hypergeom.rvs(L+N,s,L-N,random_state=prng)
		return nIntersection


# BernoulliModel--
#	Generate a sequence of bernoulli-type errors, and report the number of
#	errors, affected kmers, and islands.

class BernoulliModel(object):

	def __init__(self,ntSequenceLength,kmerSize,pSubstitution,affectedKmerCounter,islandCounter):
		self.ntSequenceLength    = ntSequenceLength
		self.kmerSize            = kmerSize
		self.pSubstitution       = pSubstitution
		self.affectedKmerCounter = affectedKmerCounter
		self.islandCounter       = islandCounter

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


# count_affected_kmers_linear--
#	Given a sequence of errors, report the number of errors and affected
#	kmers, treating the sequence as linear. The error sequence is a list of
#	1 (error) and 0 (non-error).
#
#	The algorithm maintains a sum of errors over a sliding window of length K.
#	wherever this sum is non-zero, the corresponding kmer is 'affected'.

def count_affected_kmers_linear(kmerSize,errorSeq):
	nAffected = 0
	errorsInKmer = 0
	for pos in range(len(errorSeq)+1):
		# invariant: errorsInKmer is sum of errorSeq pos-kmerSize through pos-1
		if (pos >= kmerSize):
			if (errorsInKmer > 0): nAffected += 1   # pos-kmerSize is 'affected'
			errorsInKmer -= errorSeq[pos-kmerSize]
		if (pos == len(errorSeq)): break
		errorsInKmer += errorSeq[pos]
	if ("affected" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nAffected),file=stderr)
	return nAffected


def count_affected_kmers_linear_naive(kmerSize,errorSeq):
	seqLen = len(errorSeq) - (kmerSize-1)
	nAffected = 0
	for pos in range(seqLen):
		errorsInKmer = sum(errorSeq[pos:pos+kmerSize])
		if (errorsInKmer > 0): nAffected += 1       # pos is 'affected'
	if ("affected" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nAffected),file=stderr)
	return nAffected


# count_affected_kmers_circular--
#	Given a sequence of errors, report the number of errors and affected
#	kmers, treating the sequence as circular. The error sequence is a list of
#	1 (error) and 0 (non-error).
#
#	The algorithm maintains a sum of errors over a sliding window of length K.
#	wherever this sum is non-zero, the corresponding kmer is 'affected'.

def count_affected_kmers_circular(kmerSize,errorSeq):
	seqLen = len(errorSeq)
	nAffected = 0
	errorsInKmer = 0
	for pos in range(seqLen+kmerSize):
		# invariant: errorsInKmer is sum of errorSeq pos-kmerSize through pos-1
		if (pos >= kmerSize):
			if (errorsInKmer > 0): nAffected += 1   # pos-kmerSize is 'affected'
			errorsInKmer -= errorSeq[pos-kmerSize]
		errorsInKmer += errorSeq[pos] if (pos < seqLen) else errorSeq[pos-seqLen]
	if ("affected" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nAffected),file=stderr)
	return nAffected


def count_affected_kmers_circular_naive(kmerSize,errorSeq):
	seqLen = len(errorSeq)
	extendedSeq = errorSeq + errorSeq[:kmerSize-1]

	nAffected = 0
	for pos in range(seqLen):
		errorsInKmer = sum(extendedSeq[pos:pos+kmerSize])
		if (errorsInKmer > 0): nAffected += 1       # pos is 'affected'
	if ("affected" in debug):
		print("%s (%d)" % (" ".join(map(str,errorSeq)),nAffected),file=stderr)
	return nAffected


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
	return (float(sum([n**2 for n in observed])) / (len(observed)-1)) - m**2

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
