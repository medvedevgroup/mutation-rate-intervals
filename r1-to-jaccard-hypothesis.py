#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import floor,ceil,log10
import kmer_mutation_formulas_v1 as v1
import mutation_model_simulator as mms
import hypergeometric_slicer as hgslicer


def usage(s=None):
	message = """
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
                                performed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	r1Values           = []
	ntSequenceLength   = None
	kmerSequenceLength = None
	kmerSize           = 21
	confidence         = 0.95
	numSimulations     = None
	prngSeed           = None
	reportProgress     = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.lower().startswith("--r1=")) or (arg.upper().startswith("R1=")):
			r1Values += list(map(parse_probability,argVal.split(",")))
		elif (arg.startswith("--length=")) or (arg.startswith("l=")):
			ntSequenceLength = int_with_unit(argVal)
		elif(arg.startswith("L=")):
			kmerSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--kmer=")) or (arg.upper().startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--significance=")) or (arg.startswith("C=")):
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

	if (r1Values == []):
		usage("you have to give me at least one r1 probability")

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

	if ("jmonotonicity" in debug):
		hgslicer.doJMonotonicityCheck = True

	if ("nsanity" in debug):
		hgslicer.doNLowSanityCheck  = True
		hgslicer.doNHighSanityCheck = True

	# compute the interval(s)

	L = ntSequenceLength - (kmerSize-1)
	k = kmerSize
	alpha = 1 - confidence
	z = v1.probit(1-alpha/2)

	header = ["L","k","sig","r1","jLow","jHigh"]
	if (numSimulations != None):
		header += ["jaccard","validate"]
	print("\t".join(header))

	for (r1Ix,r1) in enumerate(r1Values):
		q = v1.r1_to_q(k,r1)
		qLow  = v1.n_low (L,k,q,z) / L
		qHigh = v1.n_high(L,k,q,z) / L
		jLow  = (1-qHigh)/(1+qHigh)
		jHigh = (1-qLow )/(1+qLow )
		line = ["%d\t%d\t%.3f\t%.6f\t%.6f\t%.6f" % (L,k,confidence,r1,jLow,jHigh)]
		if (numSimulations != None):
			numDigits = max(2,int(ceil(log10(numSimulations))))
			jaccard = hgslicer.r1_to_jaccard(k,r1)
			prngSeedForR1 = ("%s_%.9f_%d" % (prngSeed,r1,r1Ix)) if (prngSeed != None) else None
			successRate = mms.jaccard_simulations(numSimulations,L,k,jaccard,None,jLow,jHigh,
			                                      prngSeed=prngSeedForR1,
			                                      reportProgress=reportProgress)
			line += ["%.6f" % jaccard]
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
