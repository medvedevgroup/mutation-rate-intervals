#!/usr/bin/env python3
"""
Formulas relating to kmers in strings subjected to independent/uniform
nucleotide substitutions.

k:  Kmer length.
L:  Sequence length; specifically, the number of complete KMERS in the sequence.
    The corresponding nucleotide sequence length would be L+k-1.
r1: Nucleotide substitution rate.
q:  1-(1-r1)^k, the probability that a kmer is mutated, i.e. that a kmer
    contains a least one substitution.

For info on the brentq solver, which is used herein, see
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html"""

from sys            import stderr,exit
from math           import sqrt
from scipy.optimize import brentq
from scipy.stats    import norm as scipy_norm

try:
	from mpmath import mp as mpmath,mpf
	mpmath.dps = 50
except ModuleNotFoundError:
	mpf = lambda v:float(v)

#==========
# formulas for Nmutated
#==========

def p_mutated(k,r1):
	return r1_to_q(k,r1)

def r1_to_q(k,r1):
	#return 1-(1-r1)**k
	r1 = mpf(r1)
	q = 1-(1-r1)**k
	return float(q)


def p_mutated_inverse(k,q):
	return q_to_r1(k,r1)

def q_to_r1(k,q):
	if not (0 <= q <= 1): return float("nan")
	#return 1-(1-q)**(1.0/k)
	q = mpf(q)
	r1 = 1-(1-q)**(1.0/k)
	return float(r1)


def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q


def var_n_mutated(L,k,r1,q=None):
	# there are computational issues in the variance formula that we solve here
	# by the use of higher-precision arithmetic; the problem occurs when r is
	# very small; for example, with L=10,k=2,r1=1e-6 standard precision
	# gives varN<0 which is nonsense; by using the mpf type, we get the correct
	# answer which is about 0.000038
	if (r1 == 0): return 0.0
	r1 = mpf(r1)
	if (q == None): # we assume that if q is provided, it is correct for r1
		q = r1_to_q(k,r1)
	varN = L*(1-q)*(q*(2*k+(2/r1)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
	     + (2*(1-q)/(r1**2))*((1+(k-1)*(1-q))*r1-q)
	assert (varN>=0.0)
	return float(varN)


# estimate_r1_from_n_mutated:
#   q = 1-(1-r1)^k  and e[nMutated] = qL
# so
#   qHat = nMutated/L
#   r1Est = 1-kth_root(1-qHat) = 1-kth_root(1-nMutated/L) 

def estimate_q_from_n_mutated(L,nMutated):
	return float(nMutated)/L


def estimate_r1_from_n_mutated(L,k,nMutated):
	return 1 - (1-float(nMutated)/L)**(1.0/k)


def confidence_interval_r1_from_n_mutated(L,k,r1,alpha):
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)
	varN = var_n_mutated(L,k,r1,q=q)
	(nLow,nHigh) = confidence_interval(L,q,varN,z)
	r1Low  = q_to_r1(k,nLow/L)
	r1High = q_to_r1(k,nHigh/L)
	return (r1Low,r1High)


def in_confidence_interval_q_from_n_mutated(L,k,r1,alpha,nMutatedObserved,useInverse=True):
	# nMutatedObserved argument can be a single value or a list
	if (not isinstance(nMutatedObserved,list)):
		nMutatedObserved = [nMutatedObserved]
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)

	if (useInverse):
		numInCI = 0
		for nMutated in nMutatedObserved:
			q1 = q_for_n_mutated_high(L,k,nMutated,z)   # nHigh(q1) == nMut
			q2 = q_for_n_mutated_low (L,k,nMutated,z)   # nLow (q2) == nMut
			if (q1 < q < q2):
				numInCI += 1
		return numInCI
	else:
		# we expect this to give exactly the same results as the useInverse case
		qLow  = n_low (L,k,q,z) / L
		qHigh = n_high(L,k,q,z) / L
		numInCI = 0
		for nMutated in nMutatedObserved:
			qHat = float(nMutated) / L
			if (qLow < qHat < qHigh):
				numInCI += 1
		return numInCI


# q_for_n_mutated_high--
#	find q s.t. nHigh(q) == nMut
#
# Note: nMut==0 is a special case. When q=0 the formula for variance has a zero
# in the denominator and thus fails to compute. Hoever, the limit of that
# formula as q goes to zero is zero (and in fact, it is easy to see that the
# variance is truly zero when q=0). This means the formulas for nLow and nHigh,
# e.g. L*q-z*sqrt(varN), are zero when q=0. Thus if nMut=0, 0 is the q for
# which nHigh(q) == nMut.
#
# nMut==L is another special case. There are two solutions in this case, one
# of which is q=1. We are interested in the other solutions

def q_for_n_mutated_high(L,k,nMut,z):
	if (nMut == 0): return 0.0   # special case, see note above
	qRight = 1 if (nMut<L) else 1-1e-5
	qLeft = 1e-5
	attemptsLeft = 10
	while (n_high(L,k,qLeft,z) >= nMut):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (n_high(L,k,qLeft,z) >= nMut):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_n_mutated_high(L=%s,k=%d,nMut=%s)" % (L,k,nMut),file=stderr)
		print("n_high(...,qLeft=%s)=%s" % (qLeft,n_low(L,k,qLeft,z)),file=stderr)
		raise ValueError

	# at this point,
	#	n_high(L,k,qLeft,z)  - nMut < 0
	#	n_high(L,k,qRight,z) - nMut > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: n_high(L,k,q,z)-nMut
	return brentq(func,qLeft,qRight)


# q_for_n_mutated_low--
#	find q s.t. nLow(q) == nMut

def q_for_n_mutated_low(L,k,nMut,z):
	qRight = 1
	qLeft = 1e-5
	attemptsLeft = 10
	while (n_low(L,k,qLeft,z) >= nMut):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (n_low(L,k,qLeft,z) >= nMut):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_n_mutated_low(L=%s,k=%d,nMut=%s)" % (L,k,nMut),file=stderr)
		print("n_low(...,qLeft=%s)=%s" % (qLeft,n_low(L,k,qLeft,z)),file=stderr)
		raise ValueError

	# at this point,
	#	n_low(L,k,qLeft,z)  - nMut < 0
	#	n_low(L,k,qRight,z) - nMut > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: n_low(L,k,q,z)-nMut
	return brentq(func,qLeft,qRight)


# confidence interval for Nmutated

def n_low(L,k,q,z):
	r1 = q_to_r1(k,q)
	varN = var_n_mutated(L,k,r1)
	return L*q - z*sqrt(varN)

def n_high(L,k,q,z):
	r1 = q_to_r1(k,q)
	varN = var_n_mutated(L,k,r1)
	return L*q + z*sqrt(varN)

#==========
# formulas for Nisland
#==========

def exp_n_island(L,k,r1):
	q = r1_to_q(k,r1)
	return L*r1*(1-q) + q - r1*(1-q)


def exp_n_island_max(L,k):
	# maximum value of E[Nisland]
	return 1 + float(L-2)/(k+1) * ((float(L-2)*k)/((L-1)*(k+1)))**k


def exp_n_island_argmax_r1(L,k):
	# value of r1 which maximizes E[Nisland]
	return float(L+k-1)/((L-1)*(k+1))


def var_n_island(L,k,r1):
	q = r1_to_q(k,r1)
	return L*r1*(1-q)*(1-r1*(1-q)*(2*k+1)) \
	     + (k**2)*(r1**2)*(1-q)**2 \
	     + k*r1*(3*r1+2)*(1-q)**2 \
	     + (1-q)*((1-q)*(r1**2)-q-r1)


# estimate_r1_from_n_island:
#	We want to find r1 s.t. exp_n_island(r1) == nIsland. Note that there are
#	usually two solutions.
#
# We look for solutions to f(r) = E[Nisland(L,k,r)] - nIsland == 0. Note that
# E[Nisland(L,k,0)] = 0 and E[Nisland(L,k,1)] = 1, and that the derivative of
# E[Nisland(L,k,r)] wrt r is positive at r=0, crosses zero at some 0<r'<1, and
# returns to zero at r=1. Thus E[Nisland] peaks at r' between 0 and 1. So long
# as the observed nIsland is more than 1 and less than this peak, we'll have
# one solution between 0 and r' and another between r' and 1.

def estimate_r1_from_n_island(L,k,nIsland):
	if (nIsland < 0): return ()
	if (nIsland == 0): return (0.0,)
	assert (nIsland >= 1)
	rPeak = exp_n_island_argmax_r1(L,k)
	if (nIsland >= exp_n_island_max(L,k)): return (rPeak,)
	# at this point,
	#	E[Nisland(L,k,0)     - nIsland < 0
	#	E[Nisland(L,k,rPeak) - nIsland > 0
	#	E[Nisland(L,k,1)     - nIsland < 0
	# so we can use the Brent's method to find solutions in those two bracketed
	# intervals
	func = lambda r1: exp_n_island(L,k,r1)-nIsland
	soln1 = brentq(func,0.0,rPeak)
	soln2 = brentq(func,rPeak,1.0)
	return (soln1,soln2)


def impossible_n_island(L,k,nIsland):
	if (nIsland < 0): return True
	if (nIsland == 0): return False
	return (nIsland >= exp_n_island_max(L,k))


def estimate_q_from_n_island(L,k,nIsland):
	r1Solutions = estimate_r1_from_n_island(L,k,nIsland)
	return map(lambda r1:r1_to_q(k,r1),r1Solutions)


def confidence_interval_r1_from_n_island(L,k,r1,alpha):
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)
	varN = var_n_island(L,k,r1)
	(nLow,nHigh) = confidence_interval(L,q,varN,z)
	r1Low  = q_to_r1(k,nLow/L)
	r1High = q_to_r1(k,nHigh/L)
	return (r1Low,r1High)


def in_confidence_interval_q_from_n_island(L,k,r1,alpha,nIslandObserved,nMutatedObserved,useInverse=True):
	return float("nan") # not implemented

#==========
# formulas relating to confidence intervals
#==========

# confidence_interval--

def confidence_interval(L,q,varN,z):
	ciMiddle = L*q
	ciHalfWidth  = z*sqrt(varN)
	return (ciMiddle-ciHalfWidth,ciMiddle+ciHalfWidth)


# probit--
#
# see https://stackoverflow.com/questions/20626994/how-to-calculate-the-inverse-of-the-normal-cumulative-distribution-function-in-p

def probit(p):
	return scipy_norm.ppf(p)

