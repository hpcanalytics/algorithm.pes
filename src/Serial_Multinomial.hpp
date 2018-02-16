/* 
 * File    : Serial_Multinomial.hpp
 * Purpose : Serial/Sequential implementation of binomial and multinomial distribution
 * Author  : Md Hasanuzzaman Bhuiyan
 * Email   : mhb@vbi.vt.edu
 * Date	   : August 31, 2012
 */
 
#ifndef _SERIAL_MULTINOMIAL_HPP
#define	_SERIAL_MULTINOMIAL_HPP

#include <cstdlib>
#include <cmath>
#include <limits>

using namespace std;

typedef long int LI;
typedef double	 LD;


// ========================= UNIFORM ================================

// Return uniformly distributed LD in range [0.0, 1.0].
inline LD Uniform1()
{
	return ((LD)rand()) / (LD) RAND_MAX;
}

// ====================== END : UNIFORM =============================




// ================ Binomial Inverse Transform Method ===============

// Return a binomial random number with parameter n and p
// using inverse transform method (exact binomial distrbuition calculation, NOT approximate calculation)
// for large values of <n>, it splits <n> into smaller <sn>, calculate binomial for all small <sn> and add their results

LI BinomialInv(LI n, LD p)
{
	LD Pn, Pi, F, U, fp, c, threshold;		
	LI i, sn, sum_sn, BinomialValue = 0;	// BinomialValue will contain the final value of Binomial distribution
											// sn -> small n [n is splitted into small values (sn) so that (1-p)^sn does not become 0]
	c  = - log(numeric_limits<LD>::min());
	threshold = c / (2.0 * p);			
	sn = (LI)floor(threshold);
	Pn = pow((LD)(1.0-p), sn);
	fp = p / (1.0-p);
	
	for(sum_sn = sn; sum_sn <= n; sum_sn += sn)
	{
		F = Pi = Pn;	
		U = Uniform1();
		for (i=0; F < U; i++) 
		{
			Pi *= (fp * (sn-i)) / (i+1);
			F  += Pi;
		}
		BinomialValue += i;
	}
	
	sn = n - (sum_sn - sn);		// last remaining part of n
	if(sn > 0)
	{
		F = Pi = pow((LD)(1.0-p), sn);
		U = Uniform1();
		for (i=0; F < U; i++) 
		{
			Pi *= (fp * (sn-i)) / (i+1);
			F  += Pi;
		}
		BinomialValue += i;
	}
	return BinomialValue;
}

// ============= END : Binomial Inverse Transform Method ============



// ========================== Binomial =============================

// Return a binomial random number with parameter n and p
// using Inverse transform (Inv)

LI Binomial(LI n, LD p)
{
    if (p == 0.0) return 0;
    if (p == 1.0) return n;
    if (p < 0.5)
		return BinomialInv(n, p);
	return n - BinomialInv(n, (LD)(1-p));
}

// ======================== END : Binomial ==========================




// ========================= Multinomial ============================

// Generate k multinomial random variables N (a vector of k long integers)
// using parameters n and parobability vector P of size k.

LI *Multinomial(LI n, LD *P, LI k, LI *N)
{
	LI i;
	LD cp = 0.0;	// cumulative prob. - set in the next line

	for (i=0; i<k; i++)		// normalize vector P if the 
		cp += P[i];			// probabilties do not sum to 1.
	for (i=0; i<k-1; i++) 
	{
		if (cp > 0) 
		{
			N[i] = Binomial(n, P[i]/cp);
			n	 = n  - N[i];
			cp 	 = cp - P[i];
		}
		else N[i] = 0;
	}
	N[k-1] = n;				// the last event of multinomial
	return N;
}

// ===================== END : Multinomial ==========================


#endif	/* _SERIAL_MULTINOMIAL_HPP */