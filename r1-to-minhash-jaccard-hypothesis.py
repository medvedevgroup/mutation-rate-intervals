#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil,log10
import mutation_model_simulator as mms


def usage(s=None):
	message = """
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
                                performed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	r1Values           = []
	ntSequenceLength   = None
	kmerSequenceLength = None
	kmerSize           = 21
	sketchSize         = None
	confidence         = 0.95
	numSlices          = 100
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
		elif (arg.startswith("--sketch=")) or (arg.startswith("S=")):
			sketchSize = int_with_unit(argVal)
		elif (arg.startswith("--significance=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--slices=")) or (arg.lower().startswith("m=")): \
			numSlices = int(argVal)
		elif (arg.startswith("--validate=")) or (arg.startswith("V=")):
			numSimulations = int_with_unit(argVal)
		elif (arg.startswith("--seed=")):
			prngSeed = argVal
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (r1Values == []):
		usage("you have to give me at least one jaccard estimate")

	if (sketchSize == None):
		usage("you have to tell me the sketch size")

	if (prngSeed != None) and (numSimulations == None):
		print("WARNING, seed is ignored since --validate was not enabled",file=stderr)

	if (ntSequenceLength != None) and (kmerSequenceLength != None):
		if (kmerSequenceLength != ntSequenceLength + kmerSize-1):
			usage("nucleotide and kmer sequence lengths are inconsistent\nyou only need to specify one of them")
	elif (kmerSequenceLength != None):
		ntSequenceLength = kmerSequenceLength + (kmerSize-1)
	elif (ntSequenceLength == None):
		ntSequenceLength = 1000 + (kmerSize-1)

	# compute the interval(s)

	L = ntSequenceLength - (kmerSize-1)
	k = kmerSize
	s = sketchSize
	alpha = 1 - confidence
	m = numSlices

	header = ["L","k","s","sig","slices","r1","jHatLow","jHatHigh"]
	if (numSimulations != None):
		header += ["jHat","validate"]
	print("\t".join(header))

	for (r1Ix,r1) in enumerate(r1Values):
		(jLow,jHigh) = hgslicer.jaccard_bounds(L,k,r1,s,alpha,m)
		line = ["%d\t%d\t%d\t%.3f\t%d\t%.6f\t%.6f\t%.6f" % (L,k,s,confidence,m,r1,jLow,jHigh)]
		if (numSimulations != None):
			numDigits = max(2,int(ceil(log10(numSimulations))))
			jHat = hgslicer.r1_to_jaccard(k,r1)
			prngSeedForR1 = ("%s_%.9f_%d" % (prngSeed,r1,r1Ix)) if (prngSeed != None) else None
			successRate = mms.jaccard_simulations(numSimulations,L,k,jHat,[s],jLow,jHigh,
			                                      prngSeed=prngSeedForR1,
			                                      reportProgress=reportProgress)
			line += ["%.6f" % jHat]
			line += ["%.*f" % (numDigits,successRate[s])]
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
