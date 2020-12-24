#!/usr/bin/env python3
"""
Induce errors in a fasta sequence, under a mutation model, and compute the
number of mutated kmers and related stats."""

from sys    import argv,stdin,stdout,stderr,exit
from random import Random,seed as random_seed, \
                   random as unit_random, \
                   choice as random_choice, \
                   sample as random_sample
from math   import sqrt,floor,ceil
from gzip   import open as gzip_open
from kmer_mutation_formulas_v1 \
            import p_mutated,exp_n_mutated,var_n_mutated,estimate_r1_from_n_mutated, \
                   confidence_interval_r1_from_n_mutated,in_confidence_interval_q_from_n_mutated


def usage(s=None):
	message = """
usage: cat fasta | simulate_nucleotide_errors [options]
  --k=<N>                   (K=) kmer size
                            (default is 28)
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
  --confidence=<p>          (C=) size of confidence interval
                            (default is 99%)
  --noinverse               for confidence interval tests, do NOT use inverse
                            functions
                            (by default inverse functions are used)
  --nosort                  don't sort output
                            (by default output is sorted by nMutated)
  --stats=<filename>        write stats to a file
                            (by default stats are written to stderr)
  --mutated=<filename>      file to write the mutated sequences to
  --mutateeonly             just write out the mutated sequences and quit
  --seed=<string>           set random seed
  --progress=<number>       periodically report how many sequence pairs we've
                            tested

Repeatedly apply the specified mutation model to a single input sequence and
report the distribution of the number of mutated kmers as well as other
related stats."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	kmerSize        = 28
	numSequences    = None
	noiseKind       = None
	pSubstitution   = None
	sequenceType    = "linear"
	confidence      = 0.99
	ciUseInverse    = True
	sortBy          = "nMutated"
	statsFilename   = None
	mutatedFilename = None
	mutateOnly      = False
	prngSeed        = None
	reportProgress  = None
	debug           = []

	statsOfInterest = ["name",
	                   "r1","k","L","confidence","trials","q",
		               "Mean[|A|].obs","Mean[|B|].obs","Mean[|A^B|].obs","Mean[|AuB|].obs","Mean[nMut.A,B].obs","Mean[L.A,B].obs",
		               "Mean[r1est.A,B].obs","inCI(r1est.A,B).obs"]

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--kmer=")) or (arg.startswith("K=")):
			kmerSize = int(argVal)
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
		elif (arg == "--noinverse"):
			ciUseInverse = False
		elif (arg == "--nosort"):
			sortBy = None
		elif (arg.startswith("--stats=")):
			statsFilename = argVal
		elif (arg.startswith("--mutated=")):
			mutatedFilename = argVal
		elif (arg == "--mutateonly"):
			mutateOnly = True
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

	if (noiseKind == "bernoulli"):
		# the presence of non-ACGT nucleotides isn't considered
		usage("the bernoulli noise model is not currently supported")

	if (sequenceType == "circular"):
		# all the estimator code assumes linear sequences
		usage("circular sequences are not currently supported")

	# set up randomness

	if (prngSeed != None):
		random_seed(prngSeed.encode("utf-8"))

	# open a file to receive the mutated sequences

	mutatedF = None
	if (mutateOnly) and (mutatedFilename == None):
		mutatedF = stdout
	else:
		if (mutatedFilename != None):
			if (mutatedFilename.endswith(".gz")) or (mutatedFilename.endswith(".gzip")):
				mutatedF = gzip_open(mutatedFilename,"wt")
			else:
				mutatedF = open(mutatedFilename,"wt")

	# fetch the *single* input sequence

	numSequencesSeen = 0
	for (seqName,seq) in fasta_sequences(stdin):
		numSequencesSeen += 1
		assert (numSequencesSeen < 2), "there was more than one sequence in the input"
		seqLen = len(seq)

	assert (numSequencesSeen == 1), "there were no sequences in the input"

	ntSequenceLength = len(seq)
	assert (ntSequenceLength >= kmerSize), "input sequence length (%d) is shorter than the kmer size (%d)" % (ntSequenceLength,kmerSize)

	distinctKmersA = kmer_set(seq,kmerSize)
	numDistinctKmersA = len(distinctKmersA)

	# set up model/generator

	if (noiseKind == "poisson") and (sequenceType == "linear"):
		kmerSequenceLength = ntSequenceLength - (kmerSize-1)
		mutationModel = PoissonModel \
		                  (seq,kmerSize,pSubstitution,
		                   count_mutated_kmers_linear)
	elif (noiseKind == "bernoulli") and (sequenceType == "linear"):
		kmerSequenceLength = ntSequenceLength - (kmerSize-1)
		mutationModel = BernoulliModel \
		                  (seq,kmerSize,pSubstitution,
		                   count_mutated_kmers_linear)
	elif (noiseKind == "poisson") and (sequenceType == "circular"):
		kmerSequenceLength = ntSequenceLength
		mutationModel = PoissonModel \
		                  (seq,kmerSize,pSubstitution,
		                   count_mutated_kmers_circular)
	elif (noiseKind == "bernoulli") and (sequenceType == "circular"):
		kmerSequenceLength = ntSequenceLength
		mutationModel = BernoulliModel \
		                  (seq,kmerSize,pSubstitution,
		                   count_mutated_kmers_circular)
	else:
		assert (False), "internal error"

	# generate sequences and collect stats

	alpha = 1 - confidence

	nErrorsObserved               = []
	nMutatedObserved              = []
	r1EstNMutatedObserved         = []
	nDistinctAObserved            = []
	nDistinctBObserved            = []
	nDistinctIntersectionObserved = []
	nDistinctUnionObserved        = []
	nMutatedABObserved            = []
	kmerSequenceLengthABObserved  = []
	r1EstABObserved               = []
	inConfR1EstABObserved         = []

	for seqNum in range(numSequences):
		if (reportProgress != None):
			if (1+seqNum == 1) or ((1+seqNum) % reportProgress == 0):
				print("testing sequence %d" % (1+seqNum),file=stderr)

		# generate a mutated sequence and collect stats

		mutatedSeq = mutationModel.generate()
		if (mutatedF != None):
			write_fasta(mutatedF,seqName+"_mutation_%d)"%(1+seqNum),mutatedSeq)
		(nErrors,nMutated) = mutationModel.count()
		nErrorsObserved  += [nErrors]
		nMutatedObserved += [nMutated]

		r1EstNMutated = estimate_r1_from_n_mutated(kmerSequenceLength,kmerSize,nMutated)
		r1EstNMutatedObserved += [r1EstNMutated]

		distinctKmersB = kmer_set(mutatedSeq,kmerSize)
		numDistinctKmersB = len(distinctKmersB)
		nDistinctKmersIntersection = len(distinctKmersA.intersection(distinctKmersB))
		nDistinctKmersUnion        = len(distinctKmersA.union(distinctKmersB))
		nDistinctAObserved            += [numDistinctKmersA]
		nDistinctBObserved            += [numDistinctKmersB]
		nDistinctIntersectionObserved += [nDistinctKmersIntersection]
		nDistinctUnionObserved        += [nDistinctKmersUnion]

		kmerSequenceLengthAB = (numDistinctKmersA+numDistinctKmersB)/2.0
		nMutatedAB          = kmerSequenceLengthAB - nDistinctKmersIntersection
		r1EstAB = estimate_r1_from_n_mutated(kmerSequenceLengthAB,kmerSize,nMutatedAB)
		nMutatedABObserved           += [nMutatedAB]
		kmerSequenceLengthABObserved += [kmerSequenceLengthAB]
		r1EstABObserved              += [r1EstAB]
		inConfR1EstAB = in_confidence_interval_q_from_n_mutated(kmerSequenceLengthAB,kmerSize,pSubstitution,alpha,
		                                                        nMutatedAB,useInverse=ciUseInverse)
		inConfR1EstABObserved += [inConfR1EstAB]

	# report per-trial results

	if (sortBy == "nMutated"):
		order = [(nDistinctIntersectionObserved[ix],ix) for ix in range(numSequences)]
		order.sort()
		order.reverse()
		order = [ix for (_,ix) in order]
	else: # if (sortBy == None):
		order = list(range(numSequences))

	header = ["L","K","r","trial","nErr","nMut","r1est.nMut","|A|","|B|","|A^B|","|AuB|","nMut.A,B","L.A,B","r1est.A,B","inCI(r1est.A,B)"]
	print("#%s" % "\t".join(header))

	for ix in range(numSequences):
		line = "\t".join(["%d","%d","%0.3f","%d","%d","%d","%0.9f","%d","%d","%d","%d","%0.1f","%0.1f","%0.9f","%d"]) \
		     % (kmerSequenceLength,                       # L
		        kmerSize,                                 # K
		        pSubstitution,                            # r
		        1+order[ix],                              # trial
		        nErrorsObserved[order[ix]],               # nErr
		        nMutatedObserved[order[ix]],              # nMut
		        r1EstNMutatedObserved[order[ix]],         # r1est.nMut
		        nDistinctAObserved[order[ix]],            # |A|
		        nDistinctBObserved[order[ix]],            # |B|
		        nDistinctIntersectionObserved[order[ix]], # |A^B|
		        nDistinctUnionObserved[order[ix]],        # |AuB|
		        nMutatedABObserved[order[ix]],            # nMut.A,B
		        kmerSequenceLengthABObserved[order[ix]],  # L.A,B
		        r1EstABObserved[order[ix]],               # r1est.A,B
		        inConfR1EstABObserved[order[ix]])         # inCI(r1est.A,B)"]
		print(line)

	if (mutatedF != None) and (mutatedF != stdout):
		mutatedF.close()

	if (mutateOnly):
		exit()

	# compute stats

	q = p_mutated(kmerSize,pSubstitution)

	nMutatedMean        = sample_mean(nMutatedObserved)
	nMutatedStDev       = sqrt(sample_variance(nMutatedObserved))
	predNMutatedMean    = exp_n_mutated(kmerSequenceLength,kmerSize,pSubstitution)
	predNMutatedStDev   = sqrt(var_n_mutated(kmerSequenceLength,kmerSize,pSubstitution))
	rmseNMutatedStDev   = abs(nMutatedStDev-predNMutatedStDev)
	rmseR1EstNMutated   = sqrt(mean_squared_error(r1EstNMutatedObserved,pSubstitution))

	(predR1EstNMutatedLow,predR1EstNMutatedHigh) \
	                    = confidence_interval_r1_from_n_mutated(kmerSequenceLength,kmerSize,pSubstitution,alpha)
	inConfR1EstNMutated = in_confidence_interval_q_from_n_mutated(kmerSequenceLength,kmerSize,pSubstitution,alpha,
	                                                              nMutatedObserved,useInverse=ciUseInverse)

	nDistinctAMean       = sample_mean(nDistinctAObserved)
	nDistinctBMean       = sample_mean(nDistinctBObserved)
	nDistinctIntersectionMean \
	                     = sample_mean(nDistinctIntersectionObserved)
	nDistinctUnionMean   = sample_mean(nDistinctUnionObserved)
	nMutatedABMean       = sample_mean(nMutatedABObserved)
	kmerSequenceLengthABMean \
	                     = sample_mean(kmerSequenceLengthABObserved)
	r1EstABMean          = sample_mean(r1EstABObserved)

	# report stats

	statToText = {}
	statToText["name"]                          = seqName
	statToText["r1"]                            = "%0.3f" % pSubstitution
	statToText["k"]                             = "%d"    % kmerSize
	statToText["L"]                             = "%d"    % kmerSequenceLength
	statToText["confidence"]                    = "%0.3f" % confidence
	statToText["trials"]                        = "%d"    % numSequences
	statToText["q"]                             = "%0.9f" % q
	statToText["E[nMut].theory"]                = "%0.9f" % predNMutatedMean
	statToText["StDev[nMut].theory"]            = "%0.9f" % predNMutatedStDev
	statToText["CIlow(r1est.nMut).theory"]      = "%0.9f" % predR1EstNMutatedLow
	statToText["CIhigh(r1est.nMut).theory"]     = "%0.9f" % predR1EstNMutatedHigh
	statToText["inCI(r1est.nMut).obs"]          = "%0.9f" % (float(inConfR1EstNMutated) / numSequences)
	statToText["Mean[nMut].obs"]                = "%0.9f" % nMutatedMean
	statToText["StDev[nMut].obs"]               = "%0.9f" % nMutatedStDev
	statToText["RMSE(StDev[nMut])"]             = "%0.9f" % rmseNMutatedStDev
	statToText["RMSE(r1est.nMut)"]              = "%0.9f" % rmseR1EstNMutated
	statToText["Mean[|A|].obs"]                 = "%d"    % nDistinctAMean
	statToText["Mean[|B|].obs"]                 = "%d"    % nDistinctBMean
	statToText["Mean[|A^B|].obs"]               = "%d"    % nDistinctIntersectionMean
	statToText["Mean[|AuB|].obs"]               = "%d"    % nDistinctUnionMean
	statToText["Mean[nMut.A,B].obs"]            = "%d"    % nMutatedABMean
	statToText["Mean[L.A,B].obs"]               = "%d"    % kmerSequenceLengthABMean
	statToText["Mean[r1est.A,B].obs"]           = "%0.9f" % r1EstABMean
	statToText["inCI(r1est.A,B).obs"]           = "%0.9f" % (float(sum(inConfR1EstABObserved)) / numSequences)

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
#	errors and mutated kmers.

ntToMutations = {"A":"CGT","C":"AGT","G":"ACT","T":"ACG",
		         "a":"cgt","c":"agt","g":"act","t":"acg"}

class PoissonModel(object):

	def __init__(self,seq,kmerSize,pSubstitution,mutatedKmerCounter):
		self.seq                = seq
		self.mutatedSeq         = None
		self.kmerSize           = kmerSize
		self.pSubstitution      = pSubstitution
		self.mutatedKmerCounter = mutatedKmerCounter

	def count(self,regenerate=False):
		if (regenerate): self.generate()
		return (sum([(self.seq[ix]!=self.mutatedSeq[ix]) for ix in range(len(self.seq))]),
		        self.mutatedKmerCounter(self.seq,self.mutatedSeq,self.kmerSize))

	def generate(self):
		errorSeq = list(map(lambda _:1 if (unit_random()<self.pSubstitution) else 0,range(len(self.seq))))
		errorPositions = [pos for (pos,err) in enumerate(errorSeq)
		                  if (err==1) and (self.seq[pos] in ntToMutations)]
		self.mutatedSeq = self.apply_errors(errorPositions)
		return self.mutatedSeq

	def apply_errors(self,errorPositions):
		# Create a mutated copy of a sequence.
		#
		# Change the nucleotide in each position of a given list to one of
		# the other three nucleotides; if the position originally contains
		# something other than a valid nucleotide, change it to one of the four
		# nucleotides.
		mutatedSeq = list(self.seq)
		for pos in errorPositions:
			nuc = self.seq[pos]
			mutations = ntToMutations[nuc] if (nuc in ntToMutations) else "ACGT"
			mutatedSeq[pos] = random_choice(mutations)
		return "".join(mutatedSeq)


# BernoulliModel--
#	Generate a sequence of bernoulli-type errors, and report the number of
#	errors and mutated kmers.

class BernoulliModel(object):

	def __init__(self,seq,kmerSize,pSubstitution,mutatedKmerCounter):
		self.seq                = seq
		self.mutatedSeq         = None
		self.kmerSize           = kmerSize
		self.pSubstitution      = pSubstitution
		self.mutatedKmerCounter = mutatedKmerCounter

	def generate(self):
		# $$$ this needs to consider that some positions don't have valid ACGT,
		#     and exclude them from the nErrors computation and from eligibility
		#     as an error position
		nErrors = round(pSubstitution * ntSequenceLength)
		if (nErrors == 0):
			self.mutatedSeq = self.seq
		else:
			errorPositions = random_sample(range(len(self.seq)),nErrors)
			self.mutatedSeq = self.apply_errors(errorPositions)
		return self.mutatedSeq


# count_mutated_kmers_linear--
#	Given a sequence pair, report the number of errors and mutated kmers,
#	treating the sequence as linear.

# $$$ needs to exclude invalid kmers
def count_mutated_kmers_linear(seq,mutatedSeq,kmerSize):
	assert (len(seq) == len(mutatedSeq))
	assert (len(seq) >= kmerSize)
	nKmers = len(seq) - (kmerSize-1)
	nMutated = 0
	for pos in range(nKmers):
		if (seq[pos:pos+kmerSize] != mutatedSeq[pos:pos+kmerSize]):
			nMutated += 1       # pos is mutated
	return nMutated


# count_mutated_kmers_circular--
#	Given a sequence pair, report the number of errors and mutated kmers,
#	treating the sequence as circular.

# $$$ needs to exclude invalid kmers
def count_mutated_kmers_circular(seq,mutatedSeq,kmerSize):
	assert (len(seq) == len(mutatedSeq))
	assert (len(seq) >= kmerSize)
	sExtended = seq + seq[:kmerSize-1]
	mExtended = mutatedSeq + mutatedSeq[:kmerSize-1]
	nKmers = len(seq)
	nMutated = 0
	for pos in range(nKmers):
		if (sExtended[pos:pos+kmerSize] != mExtended[pos:pos+kmerSize]):
			nMutated += 1       # pos is mutated
	return nMutated


# mean, variance, mean_squared_error--

def sample_mean(observed):
	return float(sum(observed)) / len(observed)

def sample_variance(observed):
	if (len(observed) <= 1): return 0.0
	m = sample_mean(observed)
	assert (False), "this formula is wrong"
	return (float(sum([n**2 for n in observed])) / (len(observed)-1)) - m**2

def mean_squared_error(observed,predicted):
	return float(sum([(n-predicted)**2 for n in observed])) / len(observed)


# kmer_set--
#	Report the set of distinct kmers in a sequence; kmers containing non-ACGT
#	are excluded

def kmer_set(seq,kmerSize):
	# note: this assumes linear sequences, not circular
	nKmers = len(seq) - (kmerSize-1)
	return set([seq[pos:pos+kmerSize] for pos in range(nKmers)
	           if (is_valid_kmer(seq[pos:pos+kmerSize]))])

def is_valid_kmer(seq):
	validCount = sum([(nt in "ACGTacgt") for nt in seq])
	return (validCount == len(seq))


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = line[1:].strip().split()[0]
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


# write_fasta--
#	Write a fasta sequence to a file

def write_fasta(f,name,seq,wrapLength=100):
	print(">%s" % name,file=f)
	if (wrapLength == None):
		print(seq,file=f)
	else:
		for i in range(0,len(seq),wrapLength):
			print(seq[i:i+wrapLength],file=f)


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
