#!/usr/bin/env python3

from sys  import argv,stdin,stdout,stderr,exit
from math import ceil
import kmer_mutation_formulas_v1 as v1


def usage(s=None):
	message = """
usage: r1-from-nmut.py [options]
  --nmut=<list>               (N=) (cumulative) observed estimates of jaccard
                              index; <list> is a comma-separated list of
                              numbers
  --length=<N>                (L=) number of kmers in the sequence
                              (default is 1K)
  --k=<N>                     (K=) kmer size
                              (default is 21)
  --confidence=<probability>  (C=) size of confidence interval
                              (default is 95%)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global reportProgress,debug

	# parse the command line

	nMutationObserved    = []
	kmerSequenceLength = 1*1000
	kmerSize           = 21
	confidence         = 0.95

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.lower().startswith("--nmut=")) or (arg.upper().startswith("N=")):
			nMutationObserved += list(map(int_with_unit,argVal.split(",")))
		elif (arg.startswith("--set=")) or (arg.startswith("L=")):
			kmerSequenceLength = int_with_unit(argVal)
		elif (arg.startswith("--kmer=")) or (arg.upper().startswith("K=")):
			kmerSize = int(argVal)
		elif (arg.startswith("--confidence=")) or (arg.startswith("C=")):
			confidence = parse_probability(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (nMutationObserved == []):
		usage("you have to give me at least one nMut observation")

	# compute the confidence interval(s)

	L = kmerSequenceLength
	k = kmerSize
	alpha = 1 - confidence
	z = v1.probit(1-alpha/2)

	print("\t".join(["L","k","conf","nMut","r1Low","r1High"]))
	for nMut in nMutationObserved:
		q1 = v1.q_for_n_affected_high(L,k,nMut,z)
		q2 = v1.q_for_n_affected_low (L,k,nMut,z)
		r1Low  = v1.q_to_r1(k,q1)
		r1High = v1.q_to_r1(k,q2)
		print("%d\t%d\t%.3f\t%d\t%.6f\t%.6f" % (L,k,confidence,nMut,r1Low,r1High))


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
