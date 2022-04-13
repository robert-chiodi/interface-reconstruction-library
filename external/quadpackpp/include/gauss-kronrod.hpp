/* gauss-kronrod.hpp
 *
 * Gauss-Kronrod quadrature using templated floating-point type;
 * based on "qk*.c" codes from <http://www.gnu.org/software/gsl/>.
 *
 * Copyright (C) 2010, 2011 Jerry Gagelman
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

#ifndef QUADPACKPP_GAUSS_KRONROD_HPP
#define QUADPACKPP_GAUSS_KRONROD_HPP

#include "machar.hpp"
#include "function.h"

/** \brief Gauss-Kronrod quadrature rule(s) plus error-estimate data for
 adaptive quadrature routines.

 The Gauss-Kronrod abscissae consist of 2m+1 points \f$x_1 < \cdots < x_{2m+1}\f$
 in the interval (-1, 1) used for a low-order and a high-order quadrature rule:
 \f[ Q_m^G f = \sum_{k=1}^m a_k f(x_{2k}), \qquad
 Q_m^{GK} f = \sum_{k=1}^{2m+1} b_k f(x_k).
 \f]

 The weights and abscissae are available through member functions, however
 they are stored according to compact QUADPACK convention. Due to symmetry,
 the positive abscissae \f$x_{2m+1}\f$, \f$x_{2m}\f$,..., \f$x_m\f$ are returned
 as values <em>xgk(0)</em>, <em>xgk(1)</em>,..., <em>xgk(m+1)</em> respectively
 of the member function \ref xgk(). Note the reverse order. The corresponding
 weights \f$b_{2m+1}\f$, \f$b_{2m}\f$,..., \f$b_m\f$ are returned by the
 respective values of \ref wgk(). The weights \f$a_m\f$, \f$a_{m-1}\f$,...,
 corresponding to the even-indexed \f$x_{2m}\f$, \f$x_{2m-2}\f$, ...., are
 returned by the values of \ref wg() in their reverse order.

 <h3>Computational details</h3>

 The even-indexed abscissae \f$x_2\f$, ..., \f$x_{2m}\f$ are the zeros of the m-th
 Legendre polynomial \f$P_m\f$. The odd indexed points are zeros of a polynomial
 that is represented as a Chebyshev sum,
 \f[ E_{m+1} = T_{m+1} + c_{m-1}T_{m-1} + c_{m-3} T_{m-3} + \cdots,
 \f]
 whose coefficients are defined by explicit formulae in
 - Giovanni Monegato, <em>Some remarks on the construction of extended Gaussian
 quadrature rules,</em> Math. Comp., Vol. 32 (1978) pp. 247-252.
 <a href="http://www.jstor.org/stable/2006272">[jstor]</a>.

 The zeros of both of these polynomials are computed by Newton's method. Upper
 bounds for their round-off errors, as functions of machine epsilon, are
 incorporated in the stopping criteria for for the root finders.

 The weights \f$a_1\f$, ..., \f$a_m\f$ are Gauss-Legendre weights.
 The \f$b_k\f$ are given by the formulae
 \f[ b_{2k} = a_k + \frac{2 p_m}{(2m+1)t_{m+1} P_m'(x_{2k}) E_{m+1}(x_{2k})},
 \qquad k = 1,\ldots,m, \f]
 and
 \f[ b_{2k+1} = \frac{2 p_m}{(2m+1) t_{m+1} P_m(x_{2k+1}) E_{m+1}'(x_{2k+1})},
 \qquad k = 0, \ldots, m, \f]
 where \f$p_m\f$ and \f$t_{m+1}\f$ are the leading coefficients of the
 polynomials \f$P_m\f$ and \f$T_{m+1}\f$ respectively. These are from
 - Giovanni Monegato, <em>A note on extended Gaussian quadrature rules,</em>
 Math. Comp., Vol. 30 (1976) pp. 812-817.
 <a href="http://www.jstor.org/stable/2005400">[jstor]</a>.
 */
template <class Real>
class GaussKronrod : public Machar<Real>
{
private:
	size_t m_;  // Gauss-Legendre degree
	size_t n_;  // size of Gauss-Kronrod arrays
	Real*  xgk_;  // Gauss-Kronrod abscissae
	Real*  wg_;   // Gauss-Legendre weights
	Real*  wgk_;  // Gauss-Kronrod weights

	Real *coefs;  // Chebyshev coefficients
	Real *zeros;  // zeros of Legendre polynomial
	Real *fv1, *fv2;  // scratch space for error estimator

	Real rescale_error(Real err, const Real result_abs,
							 const Real result_asc);

	void legendre_zeros();
	void chebyshev_coefs();
	void gauss_kronrod_abscissae();
	void gauss_kronrod_weights();

	Real legendre_err(int deg, Real x, Real& err);
	Real legendre_deriv(int deg, Real x);
	Real chebyshev_series(Real x, Real& err);
	Real chebyshev_series_deriv(Real x);

public:
	//! Initializes class for (2m+1)-point Gauss-Kronrod quadrature.
	GaussKronrod(size_t m = 10);
	~GaussKronrod();

	//! Approximates \f$\int_a^b f\,dx\f$ using the Gauss-Kronrod rule.
	void qk(FtnBase<Real>& f, Real a, Real b,
			  Real& result, Real& abserr, Real& resabs, Real& resasc);

	//! Size of arrays of Gauss-Kronrod abscissae and weights.
	size_t size() { return n_; };

	//! Array of Gauss-Kronrod abscissae in (0, 1); QUADPACK convention.
	Real xgk(int k)
	{
		return (0 <= k && k < n_) ? xgk_[k] : Real(0);
	}

	//! Array of corresponding Gauss-Kronrod weights; QUADPACK convention.
	Real wgk(int k)
	{
		return (0 <= k && k < n_) ? wgk_[k] : Real(0);
	}

	//! Gauss-Legendre weights for odd indexed abscissae; QUADPACK convention.
	Real wg(int k)
	{
		return (0 <= k && k < n_/2) ? wg_[k] : Real(0);
	}
};

template <class Real>
GaussKronrod<Real>::GaussKronrod(size_t m) : Machar<Real>()
{
	m_ = m;
	n_ = m_ + 1;
	xgk_ = new Real[n_];
	wg_  = new Real[n_ / 2];
	wgk_ = new Real[n_];
	coefs = new Real[n_ + 1];
	zeros = new Real[m_ + 2];
	fv1 = new Real[n_];
	fv2 = new Real[n_];

	legendre_zeros();
	chebyshev_coefs();
	gauss_kronrod_abscissae();
	gauss_kronrod_weights();
}

template <class Real>
GaussKronrod<Real>::~GaussKronrod()
{
	delete[] xgk_;
	delete[] wgk_;
	delete[] wg_;
	delete[] coefs;
	delete[] zeros;
	delete[] fv1;
	delete[] fv2;
}

/**
 Computes the zeros of the Legendre polynomial \f$P_m\f$. Upon exit, these
 are stored consecutively as elements of the array \c zeros[] <b>indexed by
 1,..., m.</b>
 */
template <class Real>
void GaussKronrod<Real>::legendre_zeros()
{
	Real* temp = new Real[m_+1];
	zeros[0] = Real(-1);
	zeros[1] = Real(1);
	Real delta, epsilon;

	for (int k = 1; k <= m_; ++k) {
		// loop to locate zeros of P_k interlacing z_0,...,z_k
		for (int j = 0; j < k; ++j) {
			// Newton's method for P_k :
			// initialize solver at midpoint of (z_j, z_{j+1})
			delta = 1;
			Real x_j = (zeros[j] + zeros[j+1]) / 2;
			Real P_k = legendre_err(k, x_j, epsilon);
			while (this->abs(P_k) > epsilon &&
					 this->abs(delta) > this->eps_)
			{
				delta = P_k / legendre_deriv(k, x_j);
				x_j -= delta;
				P_k = legendre_err(k, x_j, epsilon);
			}
			temp[j] = x_j;
		}

		// copy roots tmp_0,...,tmp_{k-1} to z_1,...z_k:
		zeros[k+1] = zeros[k];
		for (int j = 0; j < k; ++j)
			zeros[j+1] = temp[j];

	}
	delete[] temp;
}

/**
 Computes coefficients of polynomial \f$E_{m+1}\f$ in the array \c coefs[].
 */
template <class Real>
void GaussKronrod<Real>::chebyshev_coefs()
{
	size_t ell = (m_ + 1)/2;
	Real* alpha = new Real[ell+1];
	Real* f = new Real[ell+1];

	/* Care must be exercised in initalizing the constants in the definitions.
	 * Compilers interpret expressions like "(2*k + 1.0)/(k + 2.0)" as floating
	 * point precision, before casting to Real.
	 */
	f[1] = Real(m_+1)/Real(2*m_ + 3);
	alpha[0] = Real(1); // coefficient of T_{m+1}
	alpha[1] = -f[1];

	for (int k = 2; k <= ell; ++k) {
		f[k] = f[k-1] * (2*k - 1) * (m_ + k) / (k * (2*m_ + 2*k + 1));
		alpha[k] = -f[k];
		for (int i = 1; i < k; ++i)
			alpha[k] -= f[i] * alpha[k-i];
	}

	for (int k = 0; k <= ell; ++k) {
		coefs[m_ + 1 - 2*k] = alpha[k];
		if (m_  >= 2*k)
			coefs[m_ - 2*k] = Real(0);
	}

	delete[] alpha;
	delete[] f;
}

/**
 Computes Gauss-Legendre weights \c wg_[] and Gauss-Kronrod weights \c wgk_[].
 */
template <class Real>
void GaussKronrod<Real>::gauss_kronrod_weights()
{
	Real err;
	/* Gauss-Legendre weights:
	 */
	for (int k = 0; k < n_ / 2; ++k)
	{
		Real x = xgk_[2*k + 1];
		wg_[k] = (Real(-2) /
					 ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)));
	}

	/* The ratio of leading coefficients of P_n and T_{n+1} is computed
	 * from the recursive formulae for the respective polynomials.
	 */
	Real F_m = Real(2) / Real(2*m_ + 1);
	for (int k = 1; k <= m_; ++k)
		F_m *= (Real(2*k) / Real(2*k - 1));

	/* Gauss-Kronrod weights:
	 */
	for (size_t k = 0; k < n_; ++k)
	{
		Real x = xgk_[k];
		if (k % 2 == 0)
		{
			wgk_[k] = F_m / (legendre_err(m_, x, err) * chebyshev_series_deriv(x));
		}
		else
		{
			wgk_[k] = (wg_[k/2] +
						  F_m / (legendre_deriv(m_, x) * chebyshev_series(x, err)));
		}
	}
}

/**
 Computes the zeros of the polynomial \f$E_{m+1}\f$, using the fact that these
 interlace the zeros of the Legendre polynomial \f$P_m\f$, which are stored in
 the array \c zeros[]. Appropriate elements of \c zeros[] are then copied
 into \c xgk_[].
 */
template <class Real>
void GaussKronrod<Real>::gauss_kronrod_abscissae()
{
	Real delta, epsilon;

	for (int k = 0; 2*k < n_; ++k)
	{
		delta = 1;
		// Newton's method for E_{n+1} :
		Real x_k = (zeros[m_-k] + zeros[m_+1-k])/Real(2);
		Real E = chebyshev_series(x_k, epsilon);
		while (this->abs(E) > epsilon &&
				 this->abs(delta) > this->eps_)
		{
			delta = E / chebyshev_series_deriv(x_k);
			x_k -= delta;
			E = chebyshev_series(x_k, epsilon);
		}
		xgk_[2*k] = x_k;
		// copy adjacent Legendre-zero into the array:
		if (2*k+1 < n_)
			xgk_[2*k+1] = zeros[m_-k];
	}
}

/**
 Recursive definition of the Legendre polynomials
 \f[ (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x),
 \f]
 and estimate of the rounding error,
 \f[ E_{k+1} = \frac{(2k+1)|x|E_k + kE_{k-1}}{2(k+1)},
 \f]
 are from the routine <tt>gsl_sf_legendre_Pl_e</tt> distributed with GSL.
 */
template <class Real>
Real GaussKronrod<Real>::legendre_err(int n, Real x, Real& err)
{
	if (n == 0) {
		err = Real(0);
		return Real(1);
	}
	else if (n == 1) {
		err = Real(0);
		return x;
	}

	Real P0 = Real(1), P1 = x, P2;
	Real E0 = this->eps_;
	Real E1 = this->abs(x) * this->eps_;
	for (int k = 1; k < n; ++k)
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
		err = ((2*k + 1) * this->abs(x) * E1 + k * E0) / (2*(k + 1));
		P0 = P1; P1 = P2;
		E0 = E1; E1 = err;
	}
	return P2;
}

/**
 Three-term recursion identity for the Legendre derivatives:
 \f[ P_{k+1}'(x) = (2k+1) P_k(x) + P_{k-1}'(x).
 \f]
 */
template <class Real>
Real GaussKronrod<Real>::legendre_deriv(int n, Real x)
{
	if (n == 0)
		return Real(0);
	else if (n == 1)
		return Real(1);

	Real P0 = Real(1), P1 = x, P2;
	Real dP0 = Real(0), dP1 = Real(1), dP2;
	for (int k = 1; k < n; ++k)
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / (k + Real(1));
		dP2 = (2*k + 1) * P1 + dP0;
		P0 = P1; P1 = P2;
		dP0 = dP1; dP1 = dP2;
	}
	return dP2;
}

/**
 Evaluation of the polynomial \f$E_{m+1}\f$ is using the Clenshaw method and
 (truncation) error estimate is taken from the routine
 <tt>gsl_cheb_eval_err</tt> distributed with GSL.
*/
template <class Real>
Real GaussKronrod<Real>::chebyshev_series(Real x, Real& err)
{
	Real d1(0), d2(0);
	Real absc = this->abs(coefs[0]); // final term for truncation error
	Real y2 = 2 * x; // linear term for Clenshaw recursion

	for (int k = n_; k >= 1; --k) {
		Real temp = d1;
		d1 = y2 * d1 - d2 + coefs[k];
      d2 = temp;
		absc += this->abs(coefs[k]);
	}

	err = absc * this->eps_;
	return x * d1 - d2 + coefs[0]/2;
}

/**
 Derivatives of Chebyshev polynomials satisfy the identity \f$T_n' = nU_{n-1}\f$,
 where the \f$U_k\f$ are Chebyshev polynomials of the second kind. The derivative
 of the polynomial \f$E_{m+1}\f$ is implemented using the Clenshaw algorithm for
 the latter polynomials.
 */
template <class Real>
Real
GaussKronrod<Real>::chebyshev_series_deriv(Real x)
{
	Real d1(0), d2(0);
	Real y2 = 2 * x; // linear term for Clenshaw recursion

	for (int k = n_; k >= 2; --k) {
		Real temp = d1;
		d1 = y2 * d1 - d2 + k * coefs[k];
      d2 = temp;
	}

	return y2 * d1 - d2 + coefs[1];
}

/**
 QUADPACK's nonlinear formula for the absolute error.
 */
template <class Real>
Real GaussKronrod<Real>::rescale_error (Real err, const Real result_abs,
													 const Real result_asc)
{
	err = this->abs(err);

	if (result_asc != Real(0) && err != Real(0))
	{
		// cast 1.5 as Real number
		Real exponent = Real(3)/Real(2);
		Real scale = pow((200 * err / result_asc), exponent);

		if (scale < Real(1))
		{
			err = result_asc * scale ;
		}
		else
		{
			err = result_asc ;
		}
	}

	if (result_abs > this->xmin_ / (50 * this->eps_))
	{
      Real min_err = 50 * this->eps_ * result_abs ;

      if (min_err > err)
		{
			err = min_err ;
		}
	}

	return err ;
}

template <class Real>
void GaussKronrod<Real>::qk(FtnBase<Real>& f, Real a, Real b,
									 Real& result, Real& abserr,
									 Real& resabs, Real& resasc)
{
	const Real center = (a + b) / 2;
	const Real half_length = (b - a) / 2;
	const Real abs_half_length = this->abs(half_length);
	// const Real f_center = f.function(center, f.params);
	const Real f_center = f(center);

	Real result_gauss = Real(0);
	Real result_kronrod = f_center * wgk_[n_ - 1];
	Real result_abs = this->abs(result_kronrod);
	Real result_asc = Real(0);
	Real mean = Real(0), err = Real(0);

	int j;

	 if (n_ % 2 == 0)
	 {
		 result_gauss = f_center * wg_[n_/2 - 1];
	 }

	for (j = 0; j < (n_ - 1) / 2; j++)
	{
      int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
      Real abscissa = half_length * xgk_[jtw];
//      Real fval1 = f.function( center - abscissa , f.params);
//      Real fval2 = f.function( center + abscissa , f.params);
		Real fval1 = f(center - abscissa);
		Real fval2 = f(center + abscissa);
      Real fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg_[j] * fsum;
      result_kronrod += wgk_[jtw] * fsum;
      result_abs += wgk_[jtw] * (this->abs(fval1) + this->abs(fval2));
	}

	for (j = 0; j < n_ / 2; j++)
	{
      int jtwm1 = j * 2;
      Real abscissa = half_length * xgk_[jtwm1];
//      Real fval1 = f.function( center - abscissa , f.params);
//      Real fval2 = f.function( center + abscissa , f.params);
 		Real fval1 = f(center - abscissa);
		Real fval2 = f(center + abscissa);
		fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk_[jtwm1] * (fval1 + fval2);
      result_abs += wgk_[jtwm1] * (this->abs(fval1) + this->abs(fval2));
	};

	mean = result_kronrod / 2;

	result_asc = wgk_[n_ - 1] * this->abs(f_center - mean);

	for (j = 0; j < n_ - 1; j++)
	{
      result_asc += wgk_[j] * (this->abs(fv1[j] - mean) +
										 this->abs(fv2[j] - mean));
	}

	/* scale by the width of the integration region */

	err = (result_kronrod - result_gauss) * half_length;

	result_kronrod *= half_length;
	result_abs *= abs_half_length;
	result_asc *= abs_half_length;

	result = result_kronrod;
	resabs = result_abs;
	resasc = result_asc;
	abserr = rescale_error (err, result_abs, result_asc);
}

#endif	// QUADPACKPP_GAUSS_KRONROD_HPP
