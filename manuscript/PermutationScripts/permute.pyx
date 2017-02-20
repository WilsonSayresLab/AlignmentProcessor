'''This is module for performing permutations tests on two data sets.
It will return the counts for mean and medians before they are divided
by 10,000.'''

import random
from statistics import mean
from statistics import median

def twoSided(a, x, lena):
	cdef long meancount = 0
	cdef long mediancount = 0
	cdef long n = 0
	# Determine observed values for each pair
	cdef double obsmean = mean(x)/mean(a)
	cdef double obsmedian = median(x)/median(a)
	# Define types for permuted mean and median
	cdef double pmean
	cdef double pmedian
	# Join data sets
	join = a + x
	while n < 10000:
		n += 1
		# Shuffle and determine permuted means and medians
		random.shuffle(join, random.random)
		pmean = (mean(join[lena + 1:]))/(mean(join[:lena + 1]))
		pmedian = (median(join[lena + 1:]))/(median(join[:lena + 1]))
		# Add 1 to counts if permuted ratio is higher than absolute value
		# of observed
		if (pmean) >= abs(obsmean):
			meancount += 1
		if (pmedian) >= abs(obsmedian):
			mediancount += 1
	return meancount, mediancount

def oneSided(a, x, lena):
	cdef long lmeancount = 0
	cdef long lmediancount = 0
	cdef long rmeancount = 0
	cdef long rmediancount = 0
	cdef long n = 0
	# Determine observed values for each pair
	cdef double obsmean = mean(x)/mean(a)
	cdef double obsmedian = median(x)/median(a)
	# Define types for permuted mean and median
	cdef double pmean
	cdef double pmedian
	# Join data sets
	join = a + x
	while n < 10000:
		n += 1
		# Shuffle and determine permuted means and medians
		random.shuffle(join, random.random)
		pmean = (mean(join[lena + 1:]))/(mean(join[:lena + 1]))
		pmedian = (median(join[lena + 1:]))/(median(join[:lena + 1]))
		# Add 1 to left tail counts if permuted ratio is less than observed
		if (pmean) < (obsmean):
			lmeancount += 1
		if (pmedian) < (obsmedian):
			lmediancount += 1
		# Add 1 to right tail counts if permuted ratio is higher than observed
		if (pmean) > (obsmean):
			rmeancount += 1
		if (pmedian) > (obsmedian):
			rmediancount += 1
	return lmeancount, lmediancount, rmeancount, rmediancount
