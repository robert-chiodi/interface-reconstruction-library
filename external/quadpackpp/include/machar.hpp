/* machar.hpp
 *
 * Class that provides key floating-point parameters based on the
 * original MACHAR routine.
 *
 * Copyright (C) 2010 Jerry Gagelman
 * Copyright (C) 2006 John Burkardt
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef QUADPACKPP_MACHAR_HPP
#define QUADPACKPP_MACHAR_HPP

/** \brief Key prarameters for the \a Real floating-point type.

 The class constructor adapts source code from
 <http://people.sc.fsu.edu/~jburkardt/c_src/machar/machar.html>.
 The original algorithm is described in
 - William Cody, <em>Algorithm 665: MACHAR, a subroutine to dynamically
 determine machine parameters,</em> ACM Transactions on Mathematical Software,
 Volume 14, (1988) pp. 303-311.
 */
template <class Real>
class Machar
{
protected:
	Real eps_;	///< Machine epsilon
	Real xmin_;	///< Underflow threshold
	Real xmax_;	///< Overflow threshold

	Real inline abs(Real x) { return (x < Real(0)) ? -(x) : x; }
	Real inline max(Real a, Real b) { return a > b ? a : b; }
	Real inline min(Real a, Real b) { return a < b ? a : b; }

public:
	Machar();
	~Machar() {}

	Real inline get_eps() { return eps_; }
	Real inline get_xmin() { return xmin_; }
	Real inline get_xmax() { return xmax_; }
};

/* Naming conventions from the original FORTRAN program have been maintained.
 * A brief description is as follows:
 *
 C  IBETA   - the radix for the floating-point representation.
 C
 C  IT      - the number of base IBETA digits in the floating-point
 C            significand.
 C
 C  IRND    - 0 if floating-point addition chops,
 C            1 if floating-point addition rounds, but not in the
 C              IEEE style,
 C            2 if floating-point addition rounds in the IEEE style,
 C            3 if floating-point addition chops, and there is
 C              partial underflow,
 C            4 if floating-point addition rounds, but not in the
 C              IEEE style, and there is partial underflow,
 C            5 if floating-point addition rounds in the IEEE style,
 C              and there is partial underflow.
 C
 C  NGRD    - the number of guard digits for multiplication with
 C            truncating arithmetic.  It is
 C            0 if floating-point arithmetic rounds, or if it
 C              truncates and only  IT  base  IBETA digits
 C              participate in the post-normalization shift of the
 C              floating-point significand in multiplication;
 C            1 if floating-point arithmetic truncates and more
 C              than  IT  base  IBETA  digits participate in the
 C              post-normalization shift of the floating-point
 C              significand in multiplication.
 C
 C  MACHEP  - the largest negative integer such that
 C            1.0 + FLOAT(IBETA)**MACHEP != 1.0, except that
 C            MACHEP is bounded below by  -(IT+3).
 C
 C  NEGEP   - the largest negative integer such that
 C            1.0 - FLOAT(IBETA)**NEGEP != 1.0, except that
 C            NEGEP is bounded below by  -(IT+3).
 C
 C  IEXP    - the number of bits (decimal places if IBETA = 10)
 C            reserved for the representation of the exponent
 C            (including the bias or sign) of a floating-point
 C            number.
 C
 C  MINEXP  - the largest in magnitude negative integer such that
 C            FLOAT(IBETA)**MINEXP is positive and normalized.
 C
 C  MAXEXP  - the smallest positive power of  BETA  that overflows.
 C
 C  EPS     - the smallest positive floating-point number such
 C            that  1.0+EPS != 1.0. In particular, if either
 C            IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
 C            Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2.
 C
 C  EPSNEG  - A small positive floating-point number such that
 C            1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
 C            or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEP.
 C            Otherwise,  EPSNEG = (IBETA**NEGEP)/2.  Because
 C            NEGEP is bounded below by -(IT+3), EPSNEG may not
 C            be the smallest number that can alter 1.0 by
 C            subtraction.
 C
 C  XMIN    - the smallest non-vanishing normalized floating-point
 C            power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP.
 C
 C  XMAX    - the largest finite floating-point number.  In
 C            particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
 C            Note - on some machines  XMAX  will be only the
 C            second, or perhaps third, largest number, being
 C            too small by 1 or 2 units in the last digit of
 C            the significand.
 */
template <class Real>
Machar<Real>::Machar()
{
	long int
	ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	Real epsneg;

	Real one = (Real)1.0;
	Real two = one + one;
	Real zero = (Real)0.0;

	/* Determine IBETA and BETA. */
	Real a = one;
	Real temp, temp1;
	do {
		a = a + a;
		temp = a + one;
		temp1 = temp - a;
	} while ( temp1 - one == zero );

	Real b = one;
	do {
		b = b + b;
		temp = a + b;
	} while ( temp - a == zero );

	Real beta = temp - a;
	ibeta = (long int)beta;
	Real betah = beta / two;
	Real betain = one / beta;

	/* Determine IRND, IT. */
	it = 0;
	b = one;
	do {
		it = it + 1;
		b = b * beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);

	irnd = 0;
	temp = a + betah;
	temp1 = temp - a;
	if ( temp1 != zero )
		irnd = 1;

	Real tempa = a + beta;
	temp = tempa + betah;

	if (irnd == 0 && (temp - tempa != zero ))
		irnd = 2;

	/* Determine NEGEP, EPSNEG. */
	negep = it + 3;
	a = one;
	for (long int i = 1; i <= negep; ++i )
		a = a * betain;

	tempa = a;
	temp = one - a;
	while (temp - one == zero ) {
		a = a * beta;
		negep = negep - 1;
		temp = one - a;
	}

	negep = -negep;
	epsneg = a;

	if ( ibeta != 2 && irnd != 0 ) {
		a = ( a * (one + a) ) / two;
		temp = one - a;
		if ( temp - one != zero )
			epsneg = a;
	}

	/* Determine MACHEP, EPS. */
	machep = -it - 3;
	a = tempa;
	temp = one + a;
	while ( temp - one == zero ) {
		a = a * beta;
		machep = machep + 1;
		temp = one + a;
	}

	eps_ = a;

	/* Determine NGRD. */
	ngrd = 0;
	temp = one + eps_;
	if ( irnd == 0 && ( temp * one - one ) != zero )
		ngrd = 1;

	/* Determine IEXP, MINEXP and XMIN.
	 * Loop to determine largest I such that (1/BETA) ** (2**(I))
	 * does not underflow.  Exit from loop is signaled by an underflow.
	 */
	long int i = 0, k = 1, nxres = 0, iz = 0, mx = 0;
	Real y, z = betain, t = one + eps_;

	for ( ; ; ) {
		y = z;
		z = y * y;
		/* Check for underflow */
		a = z * one;
		temp = z * t;
		if ( ( a + a == zero ) || z > y ) {
			break;
		}
		temp1 = temp * betain;
		if ( temp1 * beta == z ) {
			iexp = i + 1;
			mx = k + k;
			break;
		}
		i = i + 1;
		k = k + k;
	}

	if ( ibeta == 10 ) {
		/* For decimal machines only */
		iexp = 2;
		iz = ibeta;
		while ( k >= iz ) {
			iz = iz * ibeta;
			iexp = iexp + 1;
		}
		mx = iz + iz - 1;
	}

	/* Loop to determine MINEXP, XMIN.
	 * Exit from loop is signaled by an underflow.
	 */
	for ( ; ; )
	{
		xmin_ = y;
		y = y * betain;
		a = y * one;
		temp = y * t;
		if ( ( a + a == zero ) || y >=  xmin_ )
			break;
		k = k + 1;
		temp1 = temp * betain;
		if ( temp1 * beta  == y ) {
			nxres = 3;
			xmin_ = y;
			break;
		}
	}

	minexp = -k;

	/* Determine MAXEXP, XMAX. */
	if ( mx <= ( k + k - 3 ) && ibeta != 10 ) {
		mx = mx + mx;
		iexp = iexp + 1;
	}

	maxexp = mx + minexp;

	/* Adjust IRND to reflect partial underflow. */
	irnd = irnd + nxres;

	/* Adjust for IEEE-style machines. */
	if ( ( irnd >= 2 ) || ( irnd == 5 ) )
		maxexp = maxexp - 2;

	/* Adjust for non-IEEE machines with partial underflow. */
	if ( ( irnd == 3 ) || ( irnd == 4 ) ) {
		maxexp = maxexp - it;
	}

	/* Adjust for machines with implicit leading bit in binary
	 * significand and machines with radix point at extreme
	 * right of significand.
	 */
	if ( ( ibeta == 2 ) && ( maxexp + minexp == 0 ) )
		maxexp = maxexp - 1;

	if ( maxexp + minexp > 20 )
		maxexp = maxexp - 1;

	if ( a != y )
		maxexp = maxexp - 2;

	xmax_ = one - epsneg;
	if ( xmax_ * one != xmax_ )
		xmax_ = one - beta * epsneg;

	xmax_ /= ( beta * beta * beta * xmin_ );
	i = maxexp + minexp + 3;

	if ( i > 0 ) {
		for (long int j = 1; j <= i; j++ ) {
			if ( ibeta == 2 )
				xmax_ = xmax_ + xmax_;
			else
				xmax_ = xmax_ * beta;
		}
	}
}

#endif	// QUADPACKPP_MACHAR_HPP
