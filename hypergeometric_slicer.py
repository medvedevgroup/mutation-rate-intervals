#!/usr/bin/env python3
"""
This implements the "hypergeometric slicer" for the confidence interval of the
Jaccard index estimated from the sketching nucleotide mutation model, as
described in "The statistics of kmers from a sequence undergoing a simple
mutation process without spurious matches," Blanca, Harris, Koslicki and
Medvedev.

k:  Kmer length.
L:  Sequence length; specifically, the number of complete KMERS in the sequence.
    The corresponding nucleotide sequence length would be L+k-1.
s:  Sketch size.
m:  Number of slices.
r1: Nucleotide substitution rate.
q:  1-(1-r1)^k, the probability that a kmer is mutated, i.e. that a kmer
    contains a least one substitution."""

from sys         import argv,stdin,stdout,stderr,exit
from math        import sqrt,log2,floor,ceil
from scipy.stats import hypergeom,norm as scipy_norm

try:
	from mpmath import mp as mpmath,mpf
	mpmath.dps = 50
except ModuleNotFoundError:
	# it can be installed with pip: "python3 -m pip install --user mpmath"
	mpf = lambda v:float(v)

# module name
#   v1 clipped n_low() and n_high() at 0 and L, respectively; but caller was
#	.. unaware of whether clipping had occurred.
#   v2 passes an indication of whether clipping occurred back to beta_low() and
#	.. beta_high() so that only the first clipped result contributes to the
#	.. corresponding sum

moduleName = "slicer.v2"

#==========
# 'hypergeometric slicer' formulas for sketch jaccard (from Nmutated)
#
# nota bene: Because these can be expensive to compute, we back each function
#            with a cache. This implementation would probably be cleaner if it
#            used the python memoization paradigm.
#==========

useCache = True
useNLowSanityCheck = True
useNHighSanityCheck = True
showQLeftSearch = False
showQRightSearch = False
#showZetaCalls = False


zeta_cache = {}
def zeta(L,s,Nmutated,a):
	# zeta(n,a)=Pr[Z>=a|N=Nmutated]
	# i.e. the tail of a hypergeometric distribution
	# in the manuscript, zeta is F
	#
	# see
	#   https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
	#
	# hypergeom.cdf(k,M,n,N)
	#   k = number of red balls drawn = a-1 (not a)
	#   M = total number of balls     = L+Nmutated
	#   n = number of red balls       = L-Nmutated
	#   N = number of draws           = s 

	#if (showZetaCalls):
	#	callStr = "zeta(%s,%s,%s,%s) = 1-hypergeom.cdf(%s,%s,%s,%s) = %.12f" \
	#	        % (L,s,Nmutated,a,
	#	           a-1,L+Nmutated,L-Nmutated,s,
	#               1 - hypergeom.cdf(a-1,L+Nmutated,L-Nmutated,s))
	#	cacheKey = (L,s,Nmutated,a)
	#	if (cacheKey in zeta_cache): callStr += " (from cache)"
	#	print(callStr,file=stderr)

	if (useCache):
		cacheKey = (L,s,Nmutated,a)
		if (cacheKey in zeta_cache):
			return zeta_cache[cacheKey]

	p = 1 - hypergeom.cdf(a-1,L+Nmutated,L-Nmutated,s)

	if (useCache):
		zeta_cache[cacheKey] = p
	return p


n_low_cache = {}
def n_low(L,k,q,m,i):
	assert (m>=2)
	assert (0<=i<=m)
	if (useCache):
		cacheKey = (L,k,q,m,i)
		if (cacheKey in n_low_cache):
			return n_low_cache[cacheKey]

	clipped = False
	if (i == 0):
		nLow = 0.0                        # (practical limit)
		clipped = True
	else:
		alphai = float(i) / m
		zi = probit(1-alphai/2)
		r1 = q_to_r1(k,q)
		varN = var_n_mutated(L,k,r1,q=q)
		sigma = sqrt(varN)
		nLow = L*q - zi*sigma
		if (nLow < 0):
			nLow = 0.0                    # (practical limit)
			clipped = True

	if (useCache):
		n_low_cache[cacheKey] = (nLow,clipped)
	return (nLow,clipped)


n_high_cache = {}
def n_high(L,k,q,m,i):
	assert (m>=2)
	assert (0<=i<=m)
	if (useCache):
		cacheKey = (L,k,q,m,i)
		if (cacheKey in n_high_cache):
			return n_high_cache[cacheKey]

	clipped = False
	if (i == 0):
		nHigh = float(L)                  # (practical limit)
		clipped = True
	else:
		alphai = float(i) / m
		zi = probit(1-alphai/2)
		r1 = q_to_r1(k,q)
		varN = var_n_mutated(L,k,r1,q=q)
		sigma = sqrt(varN)
		nHigh = L*q + zi*sigma
		if (nHigh > L):
			nHigh = float(L)              # (practical limit)
			clipped = True

	if (useCache):
		n_high_cache[cacheKey] = (nHigh,clipped)
	return (nHigh,clipped)


def precompute_n_high_low(L,k,q,m):
	# no need to cache this since n_low and n_high are cached

	nLow  = {}
	nHigh = {}

	for i in range(m+1):
		(nLow [i],_) = n_low (L,k,q,m,i)
		(nHigh[i],_) = n_high(L,k,q,m,i)

	# (sanity check)

	if (useNLowSanityCheck):
		assert (nLow[0] == 0), \
			   "for L=%d,k=%d,q=%.9f,m=%d n_low(0) = %s (expected 0)" \
			 % (L,k,q,m,nLow[0])
		for i in range(1,m+1):
			if (nLow[i-1] == 0) and (nLow[i] == 0):
				continue   # (nLow went out of bounds, avoid sanity check )
			assert (nLow[i-1] < nLow[i]), \
				   "for L=%d,k=%d,q=%.9f,m=%d n_low(%d) >= n_low(%d) (%s >= %s)" \
				 % (L,k,q,m,i-1,i,nLow[i-1],nLow[i])
		assert (abs(nLow[m] - L*q) < 1e-10)

	if (useNHighSanityCheck):
		assert (abs(nHigh[m] - L*q) < 1e-10)
		for i in range(m-1,-1,-1):
			if (nHigh[i+1] == L) and (nHigh[i] == L):
				continue   # (nHigh went out of bounds, avoid sanity check )
			assert (nHigh[i+1] < nHigh[i]), \
				   "for L=%d,k=%d,q=%.9f,m=%d n_high(%d) >= n_high(%d) (%s >= %s)" \
				 % (L,k,q,m,i+1,i,nHigh[i+1],nHigh[i])
		assert (nHigh[0] == L), \
			   "for L=%d,k=%d,q=%.9f,m=%d n_high(0) = %s (expected %d)" \
			 % (L,k,q,m,nHigh[0],L)

	return (nLow,nHigh)


beta_low_cache = {}
def beta_low(L,k,q,s,m,a):
	# in the manuscript, beta is B
	assert (m>=2)
	if (useCache):
		cacheKey = (L,k,q,s,m,a)
		if (cacheKey in beta_low_cache):
			return beta_low_cache[cacheKey]

	episilon = 0.0  # placeholder for an error term we're ignoring

	hadAClippedLow = hadAClippedHigh = False
	betaLow = 0.0
	for i in range(1,m+1):
		#betaLow += zeta(L,s,ceil(n_low(L,k,q,m,i)),a) \
		#        +  zeta(L,s,ceil(n_high(L,k,q,m,i-1)),a)

		(nLow, lowIsClipped)  = n_low (L,k,q,m,i)
		(nHigh,highIsClipped) = n_high(L,k,q,m,i-1)

		if (not lowIsClipped) or (not hadAClippedLow):
			betaLow += zeta(L,s,ceil(nLow),a)
		if (lowIsClipped): 
			hadAClippedLow = True

		if (not highIsClipped) or (not hadAClippedHigh):
			betaLow += zeta(L,s,ceil(nHigh),a)
		if (highIsClipped): 
			hadAClippedHigh = True

	betaLow /= (2*m + 2*episilon)

	if (useCache):
		beta_low_cache[cacheKey] = betaLow
	return betaLow


beta_high_cache = {}
def beta_high(L,k,q,s,m,a):
	# in the manuscript, beta is B
	assert (m>=2)
	if (useCache):
		cacheKey = (L,k,q,s,m,a)
		if (cacheKey in beta_high_cache):
			return beta_high_cache[cacheKey]

	episilon = 0.0  # placeholder for an error term we're ignoring

	hadAClippedLow = hadAClippedHigh = False
	betaHigh = 0.0
	for i in range(1,m+1):
		#betaHigh += zeta(L,s,floor(n_low(L,k,q,m,i-1)),a) \
		#         +  zeta(L,s,floor(n_high(L,k,q,m,i)),a)

		(nLow, lowIsClipped)  = n_low (L,k,q,m,i-1)
		(nHigh,highIsClipped) = n_high(L,k,q,m,i)

		if (not lowIsClipped) or (not hadAClippedLow):
			betaHigh += zeta(L,s,floor(nLow),a)
		if (lowIsClipped): 
			hadAClippedLow = True

		if (not highIsClipped) or (not hadAClippedHigh):
			betaHigh += zeta(L,s,floor(nHigh),a)
		if (highIsClipped): 
			hadAClippedHigh = True

	betaHigh /= (2*m - 2*episilon)

	if (useCache):
		beta_high_cache[cacheKey] = betaHigh
	return betaHigh


a_min_cache = {}
def a_min(L,k,q,s,alpha,m):
	# aMin = max{a : a/2 > 1-betaHigh(a)}
	assert (0<alpha<1)
	assert (m>=2)
	assert (s>=1)

	if (useCache):
		cacheKey = (L,k,q,s,alpha,m)
		if (cacheKey in a_min_cache):
			return a_min_cache[cacheKey]

	aMin = a_min_search(L,k,q,s,alpha,m)
	if (aMin == None): aMin = 0  # (no suitable a was found)

	if (useCache):
		a_min_cache[cacheKey] = aMin
	return aMin


def a_min_search(L,k,q,s,alpha,m):
	# binary search to find  max{a : alpha/2 > 1-betaHigh(a)}
	# we assume betaHigh(a) is non-increasing, that it decreases (or does not
	# .. increase) as a increases
	maxIterations = ceil(log2(s))
	halfAlpha = alpha/2

	betaHigh = beta_high(L,k,q,s,m,s)
	if (halfAlpha > 1-betaHigh):
		return s
	betaHigh = beta_high(L,k,q,s,m,0)
	if (not halfAlpha > 1-betaHigh):
		return None              # (no suitable a exists)

	# invariant:
	#   aLo < aHi  and  1-beta_high(aLo) < alpha/2 <= 1-beta_high(aHi)

	aLo = 0
	aHi = s
	iterationNum = 0
	while (aLo < aHi-1):
		iterationNum += 1
		assert (iterationNum <= maxIterations), "internal error"
		aMid = (aLo + aHi) // 2  # (truncated division)
		betaHigh = beta_high(L,k,q,s,m,aMid)
		if (halfAlpha > 1-betaHigh):
			aLo = aMid           # 1-beta_high(new aLo) < alpha/2
		else:
			aHi = aMid           # alpha/2 <= 1-beta_high(new aHi)

	return aLo


a_max_cache = {}
def a_max(L,k,q,s,alpha,m):
	# aMax = min{a : a/2 > betaLow(a)}
	assert (0<alpha<1)
	assert (m>=2)
	assert (s>=1)

	if (useCache):
		cacheKey = (L,k,q,s,alpha,m)
		if (cacheKey in a_max_cache):
			return a_max_cache[cacheKey]

	aMax = a_max_search(L,k,q,s,alpha,m)
	if (aMax == None): aMax = s  # (no suitable a was found)

	if (useCache):
		a_max_cache[cacheKey] = aMax
	return aMax


def a_max_search(L,k,q,s,alpha,m):
	# binary search to find min{a : alpha/2 > betaLow(a)}
	# we assume betaLow(a) is non-increasing, that it decreases (or does not
	# .. increase) as a increases
	maxIterations = ceil(log2(s))
	halfAlpha = alpha/2

	betaLow = beta_low(L,k,q,s,m,0)
	if (halfAlpha > betaLow):
		return 0
	betaLow = beta_low(L,k,q,s,m,s)
	if (not halfAlpha > betaLow):
		return None              # (no suitable a exists)

	# invariant:
	#   aLo < aHi  and  beta_low(aLo) >= alpha/2 > beta_low(aHi)

	aLo = 0
	aHi = s
	iterationNum = 0
	while (aLo < aHi-1):
		iterationNum += 1
		assert (iterationNum <= maxIterations), "internal error"
		aMid = (aLo + aHi) // 2  # (truncated division)
		betaLow = beta_low(L,k,q,s,m,aMid)
		if (halfAlpha > betaLow):
			aHi = aMid           # alpha/2 > beta_low(new aHi)
		else:
			aLo = aMid           # beta_low(new aLo) >= alpha/2

	return aHi


def j_low(L,k,q,s,alpha,m):
	# no need to cache this since a_min is cached
	aMin = a_min(L,k,q,s,alpha,m)
	return aMin / float(s)


def j_high(L,k,q,s,alpha,m):
	# no need to cache this since a_max is cached
	aMax = a_max(L,k,q,s,alpha,m)
	return aMax / float(s)


jaccard_bounds_cache = {}
def jaccard_bounds(L,k,r1,s,alpha,m):
	# confidence interval of jaccard estimated from a sketch of Nmutated
	if (useCache):
		cacheKey = (L,k,r1,s,alpha,m)
		if (cacheKey in jaccard_bounds_cache):
			return jaccard_bounds_cache[cacheKey]

	q = r1_to_q(k,r1)
	jLow  = j_low (L,k,q,s,alpha,m)
	jHigh = j_high(L,k,q,s,alpha,m)

	if (useCache):
		jaccard_bounds_cache[cacheKey] = (jLow,jHigh)
	return (jLow,jHigh)


def truth_in_jaccard_bounds(L,k,r1,s,alpha,m,jaccardObserved):
	# number of times the true q falls in the confidence interval of observed
	# jHat(s); this shortcuts the computation by counting the number of times
	# the observed jHat falls in the confidence interval of the true q; the
	# count should be equivalent

	# (the jaccardObserved argument can be a single value or a list)
	if (not isinstance(jaccardObserved,list)):
		jaccardObserved = [jaccardObserved]

	(jLow,jHigh) = jaccard_bounds(L,k,r1,s,alpha,m)
	numInCI = 0
	for jHat in jaccardObserved:
		if (jLow <= jHat <= jHigh):
			numInCI += 1
	return numInCI


def r1_confidence_interval(L,k,s,alpha,m,jaccardObserved):
	# (the jaccardObserved argument can be a single value or a list) find
	# r1Left and r1Right s.t. j_low'(r1Left) = jHat = j_high'(r1Right), within
	# epsilon, where j_low'(r1) = j_low(r1_to_q(r1)) and j_high'(r1) =
	# j_high(r1_to_q(r1))
	returnAsList = True
	if (not isinstance(jaccardObserved,list)):
		jaccardObserved = [jaccardObserved]
		returnAsList = False

	intervals = []
	for jHat in jaccardObserved:
		(qLeft,qRight) = q_confidence_interval(L,k,s,alpha,m,jHat)
		intervals += [(q_to_r1(k,qLeft),q_to_r1(k,qRight))]

	if (not returnAsList):
		return intervals[0]
	else:
		return intervals


q_confidence_interval_cache = {}
def q_confidence_interval(L,k,s,alpha,m,jHat,epsilon=1e-6):
	# the preferred user interface is to call this via r1_confidence_interval()
	#
	# find minimum qLeft and maximum qRight s.t. j_low(qLeft) = jHat =
	# j_high(qRight), within epsilon
	if (useCache):
		cacheKey = (L,k,s,alpha,m,jHat)
		if (cacheKey in q_confidence_interval_cache):
			return q_confidence_interval_cache[cacheKey]

	qLeft  = q_left_search (L,k,s,alpha,m,jHat,epsilon=epsilon/2)
	qRight = q_right_search(L,k,s,alpha,m,jHat,epsilon=epsilon/2)

	if (useCache):
		q_confidence_interval_cache[cacheKey] = (qLeft,qRight)
	return (qLeft,qRight)


def q_left_search(L,k,s,alpha,m,jHat,epsilon=0.5e-6):
	# find minimum qLeft s.t. j_low(qLeft) = jHat, within epsilon
	#
	# we assume that jLow(q) is decreasing (actually, non- increasing) as q
	# increases from 0 to 1
	maxIterations = 1 + ceil(-log2(epsilon))

	jLow = j_low(L,k,0.0,s,alpha,m)     # (this corresponds to qLo = 0.0)
	if   (jLow == jHat): return 0.0
	elif (jLow <  jHat): return 0.0     # (no suitable q exists)

	jLow = j_low(L,k,1.0,s,alpha,m)     # (this corresponds to qHi = 1.0)
	if (jLow > jHat): return 1.0        # (no suitable q exists)

	# invariant:
	#   qLo < qHi  and  j_low(qLo) > jHat >= j_low(qHi)

	qLo = 0.0
	qHi = 1.0
	iterationNum = 0
	while (qLo < qHi-epsilon):
		iterationNum += 1
		assert (iterationNum <= maxIterations), "internal error"
		qMid = (qLo + qHi) / 2
		jLow = j_low(L,k,qMid,s,alpha,m)
		if (showQLeftSearch):
			print("iter %d: qLo=%.12f qHi=%.12f j_low(%.12f)=%.12f" % (iterationNum,qLo,qHi,qMid,jLow),file=stderr)
		elif (jLow <= jHat):
			qHi = qMid           # jHat >= j_low(new qHi)
		else: # if (jLow > jHat):
			qLo = qMid           # j_low(new qLo) > jHat

	return qHi


def q_right_search(L,k,s,alpha,m,jHat,epsilon=0.5e-6):
	# find maximum qRight s.t. j_high(qRight) = jHat, within epsilon
	#
	# we assume that jHigh(q) is decreasing (actually, non- increasing) as q
	# increases from 0 to 1
	maxIterations = 1 + ceil(-log2(epsilon))

	jHigh = j_high(L,k,0.0,s,alpha,m)   # (this corresponds to qLo = 0.0)
	if (jHigh < jHat): return 0.0       # (no suitable q exists)

	jHigh = j_high(L,k,1.0,s,alpha,m)   # (this corresponds to qHi = 1.0)
	if   (jHigh == jHat): return 1.0
	elif (jHigh >  jHat): return 1.0    # (no suitable q exists)

	# invariant:
	#   qLo < qHi  and  j_high(qLo) >= jHat > j_high(qHi)

	qLo = 0.0
	qHi = 1.0
	iterationNum = 0
	while (qLo < qHi-epsilon):
		iterationNum += 1
		assert (iterationNum <= maxIterations), "internal error"
		qMid = (qLo + qHi) / 2
		jHigh = j_high(L,k,qMid,s,alpha,m)
		if (showQRightSearch):
			print("iter %d: qLo=%.12f qHi=%.12f j_high(%.12f)=%.12f" % (iterationNum,qLo,qHi,qMid,jHigh),file=stderr)
		elif (jHigh < jHat):
			qHi = qMid           # jHat > j_high(new qHi)
		else: # if (jHigh >= jHat):
			qLo = qMid           # j_high(new qLo) >= jHat

	return qLo

#==========
# formulas for Nmutated
#==========

def p_mutated(k,r1):
	return r1_to_q(k,r1)


def p_mutated_inverse(k,q):
	return q_to_r1(k,r1)


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
	q = mpf(q)
	varN = L*(1-q)*(q*(2*k+(2/r1)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
	     + (2*(1-q)/(r1**2))*((1+(k-1)*(1-q))*r1-q)
	assert (varN>=0.0), \
	       "for L=%d,k=%d,r1=%.9f,q=%.9f var_n_mutated evaluated as %s" \
	     % (L,k,r1,q,varN)

	return float(varN)


def estimate_r1_from_n_mutated(L,k,Nmutated):
	q = estimate_q_from_n_mutated(L,Nmutated)
	return q_to_r1(k,q)


def estimate_q_from_n_mutated(L,Nmutated):
	# e[Nmutated] = qL  ==>  q = e[Nmutated]/L
	return Nmutated / float(L)

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


def var_n_island(L,k,r1,q=None):
	# there are computational issues in the variance formula; see the note in
	# var_n_mutated()
	if (r1 == 0): return 0.0
	r1 = mpf(r1)
	if (q == None): # we assume that if q is provided, it is correct for r1
		q = r1_to_q(k,r1)
	q = mpf(q)
	varN = L*r1*(1-q)*(1-r1*(1-q)*(2*k+1)) \
	     + (k**2)*(r1**2)*((1-q)**2) \
	     + k*r1*(3*r1+2)*((1-q)**2) \
	     + (1-q)*((1-q)*(r1**2)-q-r1)
	assert (varN>=0.0), \
	       "for L=%d,k=%d,r1=%.9f,q=%.9f var_n_island evaluated as %s" \
	     % (L,k,r1,q,varN)

	return float(varN)

#==========
# base formulas
#==========

def r1_to_q(k,r1):
	#return 1-(1-r1)**k
	r1 = mpf(r1)
	q = 1-(1-r1)**k
	return float(q)


def r1_to_jaccard(k,r1):
	if not (0 <= r1 <= 1): return float("nan")
	return q_to_jaccard(r1_to_q(k,r1))


def q_to_r1(k,q):
	if not (0 <= q <= 1): return float("nan")
	#return 1-(1-q)**(1.0/k)
	q = mpf(q)
	r1 = 1-(1-q)**(1.0/k)
	return float(r1)


def q_to_jaccard(q):
	if not (0 <= q <= 1): return float("nan")
	#return (1-q)/(1+q)
	q = mpf(q)
	jaccard = (1-q)/(1+q)
	return float(jaccard)


def jaccard_to_r1(k,jaccard):
	if not (0 <= jaccard <= 1): return float("nan")
	return q_to_r1(k,jaccard_to_q(jaccard))


def jaccard_to_q(jaccard):
	if not (0 <= jaccard <= 1): return float("nan")
	#return (1-jaccard)/(1+jaccard)
	jaccard = mpf(jaccard)
	q = (1-jaccard)/(1+jaccard)
	return float(q)


# probit--
#
# see https://stackoverflow.com/questions/20626994/how-to-calculate-the-inverse-of-the-normal-cumulative-distribution-function-in-p
#
# nota bene: Because this might be expensive to compute, we back it with a
#            cache. This implementation would probably be cleaner if it used
#            the python memoization paradigm.

probit_cache = {}
def probit(p):
	cacheKey = p
	if (cacheKey in probit_cache):
		return probit_cache[cacheKey]

	z = scipy_norm.ppf(p)

	probit_cache[cacheKey] = z
	return z
