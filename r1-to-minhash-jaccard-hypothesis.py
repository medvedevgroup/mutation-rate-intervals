#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil
import hypergeometric_slicer as hgslicer


def usage(s=None):
	message = """
Compute hypothesis test for the minhash Jaccard estimate, given the mutation
rate r1.

usage: r1-to-minhash-jaccard-hypothesis.py [options]
  --r1=<list>                   (R1=) (cumulative) mutation rate; <list> is a
                                comma-separated list of probabilities
  --length=<N>                  (L=) sequence length (number of NUCLEOTIDES in
                                the sequence)
                                (default is 1K)
  --k=<N>                       (K=) kmer size
                                (default is 21)
  --sketch=<N>                  (S=) sketch size
                                (there is no default)
  --significance=<probability>  (C=) significance level
                                (default is 95%)
  --slices=<N>                  (m=) number of slices
                                (default is 100)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	r1Values         = []
	ntSequenceLength = 1*1000
	kmerSize         = 21
	sketchSize       = None
	confidence       = 0.95
	numSlices        = 100

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.lower().startswith("--r1=")) or (arg.upper().startswith("R1=")):
			r1Values += list(map(parse_probability,argVal.split(",")))
		elif (arg.startswith("--set=")) or (arg.startswith("L=")):
			ntSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--kmer=")) or (arg.upper().startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--sketch=")) or (arg.startswith("S=")):
			sketchSize = int_with_unit(argVal)
		elif (arg.startswith("--significance=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--slices=")) or (arg.lower().startswith("m=")): \
			numSlices = int(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (r1Values == []):
		usage("you have to give me at least one jaccard estimate")

	if (sketchSize == None):
		usage("you have to tell me the sketch size")

	# compute the interval(s)

	L = ntSequenceLength - (kmerSize-1)
	k = kmerSize
	s = sketchSize
	alpha = 1 - confidence
	m = numSlices

	print("\t".join(["L","k","s","sig","slices","r1","jHatLow","jHatHigh"]))
	for r1 in r1Values:
		(jLow,jHigh) = hgslicer.jaccard_bounds(L,k,r1,s,alpha,m)
		print("%d\t%d\t%d\t%.3f\t%d\t%.6f\t%.6f\t%.6f" % (L,k,s,confidence,m,r1,jLow,jHigh))


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
