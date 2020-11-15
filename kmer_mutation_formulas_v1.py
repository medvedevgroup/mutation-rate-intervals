#!/usr/bin/env python3
"""
Formulas relating to kmers in strings subjected to independent/uniform
nucleotide substitutions.

k:  Kmer length.
L:  Sequence length; specifically, the number of complete kmers in the sequence.
    The corresponding nucleotide sequence length would be L+k-1.
s:  Sketch size.
r1: Nucleotide substitution rate.
q:  1-(1-r1)^k, the probability that a kmer is 'affected'. A kmer is affected
    if it contains a least one substitution.

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
# formulas for Naffected
#==========

def p_affected(k,r1):
	return r1_to_q(k,r1)

def r1_to_q(k,r1):
	#return 1-(1-r1)**k
	r1 = mpf(r1)
	q = 1-(1-r1)**k
	return float(q)


def p_affected_inverse(k,q):
	return q_to_r1(k,r1)

def q_to_r1(k,q):
	if not (0 <= q <= 1): return float("nan")
	#return 1-(1-q)**(1.0/k)
	q = mpf(q)
	r1 = 1-(1-q)**(1.0/k)
	return float(r1)


def exp_n_affected(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q


def var_n_affected(L,k,r1,q=None):
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


# estimate_r1_from_n_affected:
#   q = 1-(1-r1)^k  and e[nAffected] = qL
# so
#   qHat = nAffected/L
#   r1Est = 1-kth_root(1-qHat) = 1-kth_root(1-nAffected/L) 

def estimate_q_from_n_affected(L,nAffected):
	return float(nAffected)/L


def estimate_r1_from_n_affected(L,k,nAffected):
	return 1 - (1-float(nAffected)/L)**(1.0/k)


def confidence_interval_r1_from_n_affected(L,k,r1,alpha):
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)
	varN = var_n_affected(L,k,r1,q=q)
	(nLow,nHigh) = confidence_interval(L,q,varN,z)
	r1Low  = q_to_r1(k,nLow/L)
	r1High = q_to_r1(k,nHigh/L)
	return (r1Low,r1High)


def in_confidence_interval_q_from_n_affected(L,k,r1,alpha,nAffectedObserved,useInverse=True):
	# nAffectedObserved argument can be a single value or a list
	if (not isinstance(nAffectedObserved,list)):
		nAffectedObserved = [nAffectedObserved]
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)

	if (useInverse):
		numInCI = 0
		for nAffected in nAffectedObserved:
			q1 = q_for_n_affected_high(L,k,nAffected,z)   # nHigh(q1) == nAff
			q2 = q_for_n_affected_low (L,k,nAffected,z)   # nLow (q2) == nAff
			if (q1 < q < q2):
				numInCI += 1
		return numInCI
	else:
		# we expect this to give exactly the same results as the useInverse case
		qLow  = n_low (L,k,q,z) / L
		qHigh = n_high(L,k,q,z) / L
		numInCI = 0
		for nAffected in nAffectedObserved:
			qHat = float(nAffected) / L
			if (qLow < qHat < qHigh):
				numInCI += 1
		return numInCI


# q_for_n_affected_high--
#	find q s.t. nHigh(q) == nAff
#
# Note: nAff==0 is a special case. When q=0 the formula for variance has a zero
# in the denominator and thus fails to compute. Hoever, the limit of that
# formula as q goes to zero is zero (and in fact, it is easy to see that the
# variance is truly zero when q=0). This means the formulas for nLow and nHigh,
# e.g. L*q-z*sqrt(varN), are zero when q=0. Thus if nAff=0, 0 is the q for
# which nHigh(q) == nAff.
#
# nAff==L is another special case. There are two solutions in this case, one
# of which is q=1. We are interested in the other solutions

def q_for_n_affected_high(L,k,nAff,z):
	if (nAff == 0): return 0.0   # special case, see note above
	qRight = 1 if (nAff<L) else 1-1e-5
	qLeft = 1e-5
	attemptsLeft = 10
	while (n_high(L,k,qLeft,z) >= nAff):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (n_high(L,k,qLeft,z) >= nAff):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_n_affected_high(L=%s,k=%d,nAff=%s)" % (L,k,nAff),file=stderr)
		print("n_high(...,qLeft=%s)=%s" % (qLeft,n_low(L,k,qLeft,z)),file=stderr)
		raise ValueError

	# at this point,
	#	n_high(L,k,qLeft,z)  - nAff < 0
	#	n_high(L,k,qRight,z) - nAff > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: n_high(L,k,q,z)-nAff
	return brentq(func,qLeft,qRight)


# q_for_n_affected_low--
#	find q s.t. nLow(q) == nAff

def q_for_n_affected_low(L,k,nAff,z):
	qRight = 1
	qLeft = 1e-5
	attemptsLeft = 10
	while (n_low(L,k,qLeft,z) >= nAff):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (n_low(L,k,qLeft,z) >= nAff):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_n_affected_low(L=%s,k=%d,nAff=%s)" % (L,k,nAff),file=stderr)
		print("n_low(...,qLeft=%s)=%s" % (qLeft,n_low(L,k,qLeft,z)),file=stderr)
		raise ValueError

	# at this point,
	#	n_low(L,k,qLeft,z)  - nAff < 0
	#	n_low(L,k,qRight,z) - nAff > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: n_low(L,k,q,z)-nAff
	return brentq(func,qLeft,qRight)


# confidence interval for Naffected

def n_low(L,k,q,z):
	r1 = q_to_r1(k,q)
	varN = var_n_affected(L,k,r1)
	return L*q - z*sqrt(varN)

def n_high(L,k,q,z):
	r1 = q_to_r1(k,q)
	varN = var_n_affected(L,k,r1)
	return L*q + z*sqrt(varN)


#==========
# formulas for sketch jaccard (from Naffected)
#==========

def confidence_interval_jaccard_from_sketch_n_affected(L,k,s,r1,alpha,alpha1=None):
	z = probit(1-alpha/2)
	q = r1_to_q(k,r1)
	if (alpha1 == None): alpha1 = 1 - sqrt(1-alpha)   # alpha1 s.t. alpha2=alpha1
	jLow  = j_low (L,k,q,s,alpha,alpha1)
	jHigh = j_high(L,k,q,s,alpha,alpha1)
	return (jLow,jHigh)


def in_confidence_interval_jaccard_from_sketch_n_affected(L,k,s,r1,alpha,alpha1,jaccardObserved,useInverse=True):
	# note: if (alpha1 == None) we use alpha1 s.t. alpha2=alpha1
	assert(useInverse), "non-inverse method has not been implemented"
	# jaccardObserved argument can be a single value or a list
	if (not isinstance(jaccardObserved,list)):
		jaccardObserved = [jaccardObserved]
	q = r1_to_q(k,r1)
	numInCI = 0
	for jaccard in jaccardObserved:
		q1 = q_for_j_low (L,k,s,jaccard,alpha,alpha1)   # jLow (q1) == jaccard
		q2 = q_for_j_high(L,k,s,jaccard,alpha,alpha1)   # jHigh(q2) == jaccard
		if (q1 < q < q2):
			numInCI += 1
	return numInCI


# q_for_j_low--
#	find q s.t. jLow(q) == jaccard
#
# jLow(q) decreases from 1 at q=0, crosses below zero at some q1<1, then
# increases back to zero at some q2<1. For q>q2 jLow(q) is poorly defined
# because its definition violates the valid range of psi.
#
# So for 0<j<=1 jLow(q)=j has a unique solution. For j=0 there are either two
# solutions (q1 and q2) or infinitely many if one decides to treat the interval
# of poor definition as zero.
#
# Special cases:
#	jLow == 1:
#		jLow(0) = 1 regardless of the other parameters.
#	jLow == 0:
#		We are only interested in the smallest q s.t. jLow(q) = 0, which is
#		called q1 above. Looking at the definition of jLow(q) we define
#		u = psi(qHigh(q)) so that jLow(q) is u-z2*sqrt((u(1-u))/s). Solving
#		jLow(q)=0 gives us (after we ignore the u=0 solution)
#			psi(qHigh(q)) = u = z2^2/(s+z2^2)
#		or equivalently
#			qHigh(q) = s/(s+2*z2^2)
#		We then solve for q using the inverse of qHigh(); actually using the
#		inverse of nHigh. 

def q_for_j_low(L,k,s,jaccard,alpha,alpha1=None):
	assert (0<alpha<1)
	if (alpha1 == None): alpha1 = 1 - sqrt(1-alpha)   # alpha1 s.t. alpha2=alpha1

	if (jaccard == 1):
		return 0.0

	if (jaccard == 0):
		alpha2 = (alpha-alpha1) / (1-alpha1)
		z1 = probit(1-alpha1/2)
		z2 = probit(1-alpha2/2)
		qHigh = s/(s + 2*z2**2)
		return q_for_n_affected_high(L,k,L*qHigh,z1)

	qLeft  = 1e-5
	qRight = 1-1e-5

	attemptsLeft = 10
	while (j_low(L,k,qLeft,s,alpha,alpha1) <= jaccard):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (j_low(L,k,qLeft,s,alpha,alpha1) <= jaccard):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_j_low(L=%s,k=%d,s=%d,jaccard=%s)" % (L,k,s,jaccard),file=stderr)
		print("j_low(...,qLeft=%s)=%s" % (qLeft,j_low(L,k,qLeft,s,alpha,alpha1)),file=stderr)
		raise ValueError

	attemptsLeft = 10
	while (j_low(L,k,qRight,s,alpha,alpha1) >= jaccard):
		qRight = (1+qRight)/2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (j_low(L,k,qRight,s,alpha,alpha1) >= jaccard):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_j_low(L=%s,k=%d,s=%d,jaccard=%s)" % (L,k,s,jaccard),file=stderr)
		print("j_low(...,qRight=%s)=%s" % (qRight,j_low(L,k,qRight,s,alpha,alpha1)),file=stderr)
		raise ValueError

	# at this point,
	#	j_low(L,k,qLeft,s,alpha,alpha1)  - jaccard < 0
	#	j_low(L,k,qRight,s,alpha,alpha1) - jaccard > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: j_low(L,k,q,s,alpha,alpha1)-jaccard
	return brentq(func,qLeft,qRight)


# q_for_j_high--
#	find q s.t. jHigh(q) == jaccard
#
# jHigh(q) is poorly defined at q=0 because its definition violates the valid
# range of psi. As q increases, jHigh(q) continues to be poorly defined until
# some q1<1 where it is 1. It then increases above 1, then decreases back to
# cross from 1 to less than 1 at some q2<1. It then continues to decrease until
# it reaches zero at q=1.
#
# So for 0<=j<1 jHigh(q)=j has a unique solution. For j=1 there are either two
# solutions (q1 and q2) or infinitely many if one decides to treat the interval
# of poor definition as 1.
#
# Special cases:
#	jHigh == 0:
#		jHigh(1) = 0 regardless of the other parameters.
#	jHigh == 1:
#		We are only interested in the largest q s.t. jHigh(q) = 1, which is
#		called q2 above. Looking at the definition of jHigh(q) we define
#		u = psi(qLow(q)) so that jHigh(q) is u+z2*sqrt((u(1-u))/s). Solving
#		jHigh(q)=1 gives us (after we ignore the u=1 solution)
#			psi(qLow(q)) = u = s/(s+z2^2)
#		or equivalently
#			qLow(q) = z2^2/(2*s+z2^2)
#		We then solve for q using the inverse of qLow(); actually using the
#		inverse of nLow. 

def q_for_j_high(L,k,s,jaccard,alpha,alpha1=None):
	assert (0<alpha<1)
	if (alpha1 == None): alpha1 = 1 - sqrt(1-alpha)   # alpha1 s.t. alpha2=alpha1

	if (jaccard == 0):
		return 1.0

	if (jaccard == 1):
		alpha2 = (alpha-alpha1) / (1-alpha1)
		z1 = probit(1-alpha1/2)
		z2 = probit(1-alpha2/2)
		qLow = (z2**2)/(2*s + z2**2)
		return q_for_n_affected_low(L,k,L*qLow,z1)

	qLeft  = 1e-5
	qRight = 1-1e-5

	attemptsLeft = 10
	while (j_high(L,k,qLeft,s,alpha,alpha1) <= jaccard):
		qLeft /= 2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (j_high(L,k,qLeft,s,alpha,alpha1) <= jaccard):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_j_high(L=%s,k=%d,s=%d,jaccard=%s)" % (L,k,s,jaccard),file=stderr)
		print("j_high(...,qLeft=%s)=%s" % (qLeft,j_high(L,k,qLeft,s,alpha,alpha1)),file=stderr)
		raise ValueError

	attemptsLeft = 10
	while (j_high(L,k,qRight,s,alpha,alpha1) >= jaccard):
		qRight = (1+qRight)/2
		attemptsLeft -= 1
		if (attemptsLeft < 0): break
	if (j_high(L,k,qRight,s,alpha,alpha1) >= jaccard):
		# this is just laziness, it really means our assumptions about the
		# solution space are wrong
		print("q_for_j_high(L=%s,k=%d,s=%d,jaccard=%s)" % (L,k,s,jaccard),file=stderr)
		print("j_high(...,qRight=%s)=%s" % (qRight,j_high(L,k,qRight,s,alpha,alpha1)),file=stderr)
		raise ValueError

	# at this point,
	#	j_high(L,k,qLeft,s,alpha,alpha1)  - jaccard < 0
	#	j_high(L,k,qRight,s,alpha,alpha1) - jaccard > 0
	# so we can use the Brent's method to find the solution in the bracketed
	# interval
	func = lambda q: j_high(L,k,q,s,alpha,alpha1)-jaccard
	return brentq(func,qLeft,qRight)


# confidence interval for Jaccard from Naffected

def psi(x):                 # only valid for 0<=x<=1
	assert (0<=x<=1)
	return (1-x)/(1+x)

def j_low(L,k,q,s,alpha,alpha1):
	assert (0<alpha<1)
	assert (0<alpha1<alpha)
	alpha2 = (alpha-alpha1) / (1-alpha1)
	z1 = probit(1-alpha1/2)
	z2 = probit(1-alpha2/2)
	qHigh = n_high(L,k,q,z1)/L
	psiQHigh = psi(min(1,qHigh))
	return psiQHigh - z2*sqrt(psiQHigh*(1-psiQHigh)/s)

def j_high(L,k,q,s,alpha,alpha1):
	assert (0<alpha<1)
	assert (0<alpha1<alpha)
	alpha2 = (alpha-alpha1) / (1-alpha1)
	z1 = probit(1-alpha1/2)
	z2 = probit(1-alpha2/2)
	qLow = n_low(L,k,q,z1)/L
	psiQLow = psi(max(0,qLow))
	return psiQLow + z2*sqrt(psiQLow*(1-psiQLow)/s)

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


def in_confidence_interval_q_from_n_island(L,k,r1,alpha,nIslandObserved,nAffectedObserved,useInverse=True):
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

