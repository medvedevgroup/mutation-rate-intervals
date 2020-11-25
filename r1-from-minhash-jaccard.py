#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil
import hypergeometric_slicer as hgslicer


def usage(s=None):
	message = """
usage: r1-from-minhash-jaccard.py [options]
  --jhat=<list>               (J=) (cumulative) observed estimates of jaccard
                              index; <list> is a comma-separated list of
                              numbers between 0 and 1
  --length=<N>                (N= or L=) number of kmers in the sequence
                              (default is 1K)
  --k=<N>                     (K=) kmer size
                              (default is 21)
  --sketch=<N>                (S=) sketch size
                              (there is no default)
  --confidence=<probability>  (C=) size of confidence interval
                              (default is 95%)
  --slices=<N>                (m=) number of slices
                              (default is 100)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	jaccardObserved    = []
	kmerSequenceLength = 1*1000
	kmerSize           = 21
	sketchSize         = None
	confidence         = 0.95
	numSlices          = 100

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.lower().startswith("--jhat=")) or (arg.upper().startswith("J=")):
			jaccardObserved += list(map(parse_probability,argVal.split(",")))
		elif (arg.startswith("--set=")) or (arg.startswith("N=")) or (arg.startswith("L=")):
			kmerSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--kmer=")) or (arg.upper().startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--sketch=")) or (arg.startswith("S=")):
			sketchSize = int_with_unit(argVal)
		elif (arg.startswith("--confidence=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--slices=")) \
		 or (arg.startswith("m=")) or (arg.startswith("--m=")) \
		 or (arg.startswith("M=")) or (arg.startswith("--M=")):
			numSlices = int(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (jaccardObserved == []):
		usage("you have to give me at least one jaccard estimate")

	if (sketchSize == None):
		usage("you have to tell me the sketch size")

	# compute the confidence interval(s)

	L = kmerSequenceLength
	k = kmerSize
	s = sketchSize
	alpha = 1 - confidence
	m = numSlices

	print("\t".join(["L","k","s","alpha","jHat","r1Low","r1High"]))
	for jHat in jaccardObserved:
		(r1Low,r1High) = hgslicer.r1_confidence_interval(L,k,s,alpha,m,jHat)
		print("%d\t%d\t%d\t%.3f\t%.6f\t%.6f\t%.6f" % (L,k,s,alpha,jHat,r1Low,r1High))


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
