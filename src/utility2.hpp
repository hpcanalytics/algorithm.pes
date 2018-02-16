/* 
 * File  : utility2.hpp
 * Author: Md Hasanuzzaman Bhuiyan
 * Email : mhb@vbi.vt.edu
 * Created on February 13, 2012
 */

#ifndef _UTILITY2_HPP
#define	_UTILITY2_HPP

#include <climits>
#include <cfloat>

// Reference: http://www.cplusplus.com/reference/clibrary/climits/
//			: http://www.cplusplus.com/reference/clibrary/cfloat/
//			: http://www.cplusplus.com/reference/std/limits/numeric_limits/
//			: http://msdn.microsoft.com/en-us/library/296az74e%28v=vs.80%29.aspx


// ================= returns max value of character (127) ===============================

inline char maxValue (char x)
{
	return CHAR_MAX; 
}

// ================= returns max value of unsigned character (255) ======================

inline unsigned char maxValue (unsigned char x)
{
	return UCHAR_MAX;
}

// ================= returns max value of integer (2147483647) ==========================

inline int maxValue (int x)
{
	return INT_MAX;
}

// ================= returns max value of unsigned integer (4294967295) =================

inline unsigned int maxValue (unsigned int x)
{
	return UINT_MAX;
}

// ================= returns max value of long (2147483647) =============================

inline long maxValue (long x)
{
	return LONG_MAX;
}

// ================= returns max value of unsigned long (4294967295) ====================

inline unsigned long maxValue (unsigned long x)
{
	return ULONG_MAX;
}

// ================= returns max value of double (approx. 1E+37) ========================

inline double maxValue (double x)
{
	return DBL_MAX;
}

#endif	/* _UTILITY2_HPP */
