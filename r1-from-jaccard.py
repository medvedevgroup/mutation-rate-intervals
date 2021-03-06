#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil,log10
import kmer_mutation_formulas_thm5 as thm5
import mutation_model_simulator as mms
import hypergeometric_slicer as hgslicer


def usage(s=None):
	message = """
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
                              performed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	jaccardObserved    = []
	ntSequenceLength   = None
	kmerSequenceLength = None
	kmerSize           = 21
	confidence         = 0.95
	numSimulations     = None
	prngSeed           = None
	reportProgress     = None
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--help","-help","--h","-h"]):
			usage()
		elif (arg.lower().startswith("--jhat=")) or (arg.upper().startswith("J=")):
			jaccardObserved += list(map(parse_probability,argVal.split(",")))
		elif (arg.startswith("--length=")) or (arg.startswith("l=")):
			ntSequenceLength = int_with_unit(argVal)
		elif(arg.startswith("L=")):
			kmerSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--kmer=")) or (arg.upper().startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--confidence=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--validate=")) or (arg.startswith("V=")):
			numSimulations = int_with_unit(argVal)
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

	if (jaccardObserved == []):
		usage("you have to give me at least one jaccard estimate")

	if (prngSeed != None) and (numSimulations == None):
		print("WARNING, seed is ignored since --validate was not enabled",file=stderr)

	if (ntSequenceLength != None) and (kmerSequenceLength != None):
		if (kmerSequenceLength != ntSequenceLength + kmerSize-1):
			usage("nucleotide and kmer sequence lengths are inconsistent\nyou only need to specify one of them")
	elif (kmerSequenceLength != None):
		ntSequenceLength = kmerSequenceLength + (kmerSize-1)
	elif (ntSequenceLength == None):
		ntSequenceLength = 1000 + (kmerSize-1)

	if ("nocache" in debug):
		hgslicer.useCache = False

	if ("nojmonotonicity" in debug):
		hgslicer.doJMonotonicityCheck = False
	else:
		hgslicer.doJMonotonicityCheck = True

	if ("nsanity" in debug):
		hgslicer.doNLowSanityCheck  = True
		hgslicer.doNHighSanityCheck = True

	# compute the confidence interval(s)

	L = ntSequenceLength - (kmerSize-1)
	k = kmerSize
	alpha = 1 - confidence
	z = thm5.probit(1-alpha/2)

	header = ["L","k","conf","jaccard","r1Low","r1High"]
	if (numSimulations != None):
		header += ["r1Hat","validate"]
	print("\t".join(header))

	for (jaccardIx,jaccard) in enumerate(jaccardObserved):
		nMut = L * (1-jaccard)/(1+jaccard)
		q1 = thm5.q_for_n_mutated_high(L,k,nMut,z)
		q2 = thm5.q_for_n_mutated_low (L,k,nMut,z)
		r1Low  = thm5.q_to_r1(k,q1)
		r1High = thm5.q_to_r1(k,q2)
		line = ["%d\t%d\t%.3f\t%.6f\t%.6f\t%.6f" % (L,k,confidence,jaccard,r1Low,r1High)]
		if (numSimulations != None):
			numDigits = max(2,int(ceil(log10(numSimulations))))
			r1Hat = hgslicer.jaccard_to_r1(k,jaccard)
			prngSeedForJaccard = ("%s_%.9f_%d" % (prngSeed,jaccard,jaccardIx)) if (prngSeed != None) else None
			successRate = mms.r1_simulations(numSimulations,L,k,r1Hat,None,r1Low,r1High,
			                                 prngSeed=prngSeedForJaccard,
			                                 reportProgress=reportProgress)
			line += ["%.6f" % r1Hat]
			line += ["%.*f" % (numDigits,successRate["no sketch"])]
		print("\t".join(line))


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
