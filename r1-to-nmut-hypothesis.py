#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import floor,ceil,log10
import kmer_mutation_formulas_thm5 as thm5
import mutation_model_simulator as mms


def usage(s=None):
	message = """
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
                                performed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

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

		if (arg in ["--help","-help","--h","-h"]):
			usage()
		elif (arg.lower().startswith("--r1=")) or (arg.upper().startswith("R1=")):
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

	# compute the interval(s)

	L = ntSequenceLength - (kmerSize-1)
	k = kmerSize
	alpha = 1 - confidence
	z = thm5.probit(1-alpha/2)

	header = ["L","k","sig","r1","nMutLow","nMutHigh"]
	if (numSimulations != None):
		header += ["nMut","validate"]
	print("\t".join(header))

	for (r1Ix,r1) in enumerate(r1Values):
		q = thm5.r1_to_q(k,r1)
		nMutLow  = max(0,floor(thm5.n_low (L,k,q,z)))
		nMutHigh = min(L,ceil (thm5.n_high(L,k,q,z)))
		line = ["%d\t%d\t%.3f\t%.6f\t%d\t%d" % (L,k,confidence,r1,nMutLow,nMutHigh)]
		if (numSimulations != None):
			numDigits = max(2,int(ceil(log10(numSimulations))))
			nMutHat = int(round(L*q))
			prngSeedForR1 = ("%s_%.9f_%d" % (prngSeed,r1,r1Ix)) if (prngSeed != None) else None
			successRate = mms.nmut_simulations(numSimulations,L,k,nMutHat,None,nMutLow,nMutHigh,
			                                   prngSeed=prngSeedForR1,
			                                   reportProgress=reportProgress)
			line += ["%d" % nMutHat]
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
