/* workspace.hpp
 *
 * Adaptive quadrature for a templated floating-point type based on the
 * routine QAG implemented in <http://www.gnu.org/software/gsl/>.
 *
 * Copyright (C) 2010 Jerry Gagelman
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#ifndef QUADPACKPP_WORKSPACE_HPP
#define QUADPACKPP_WORKSPACE_HPP

#include "errcodes.h"
#include "gauss-kronrod.hpp"

/** \brief The <tt>gsl_integration_workspace</tt> structure with member
 function \ref qag() for adaptive quadrature.
 */
template <class Real>
class Workspace : public GaussKronrod<Real>
{
public:
	Workspace();
	//! Initialize workspace for \a limit refinement intervals.
	Workspace(size_t limit);
	//! Initialize for \a limit refinement intervals and (2m+1)-point quadrature.
	Workspace(size_t limit, size_t m);
	~Workspace();

	/** \brief Implementation of <tt>gsl_integration_qag</tt> using class data.

	 The arguments correspond to the <a href=
	 "http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html">
	 original GSL</a> function except that <em>limit</em> and <em>key</em> are
	 use the workspace size and Gauss-Kronrod rule initializing the class.
	 The stopping criterion is
	 <center>
	 |\a result \f$-\int f\,dx\f$| \f$\leq\f$
	 max{\a epsabs, \a epsrel\f$|\int f\,dx|\f$},
	 </center>
	 and \a abserr is the internal estimate of the right hand side.
	 */
	int qag(FtnBase<Real>& f, Real a, Real b,
			  Real epsabs, Real epsrel, Real& result, Real& abserr);

private:
	void	allocate(size_t limit);

	/* data from gsl_integration_workspace */
	size_t limit;
	size_t size;
	size_t nrmax;
	size_t i_work;
	size_t maximum_level;
	Real   *alist;
	Real   *blist;
	Real   *rlist;
	Real   *elist;
	size_t *order;
	size_t *level;

	/* auxillary functions for adaptive quadrature */
	void  append_interval(Real a1, Real b1, Real area1, Real error1);
	void  initialise(Real a, Real b);
	void  set_initial_result(Real result, Real error);
	void  qpsrt();
	void  sort_results();
	void  retrieve(Real& a, Real& b, Real& r, Real& e);
	int   subinterval_too_small (Real a1, Real a2, Real b2);
	Real  sum_results();
	void  update(Real a1, Real b1, Real area1, Real error1,
					 Real a2, Real b2, Real area2, Real error2);

};

template <class Real>
void
Workspace<Real>::allocate(size_t n)
{
	alist = new Real[n];
	blist = new Real[n];
	rlist = new Real[n];
	elist = new Real[n];
	order = new size_t[n];
	level = new size_t[n];

	size = 0;
	limit = n;
	maximum_level = 0;
}

template <class Real>
Workspace<Real>::Workspace() : GaussKronrod<Real>()
{
	allocate(1);
}

template <class Real>
Workspace<Real>::Workspace(size_t n) : GaussKronrod<Real>()
{
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

template <class Real>
Workspace<Real>::Workspace(size_t n, size_t m) : GaussKronrod<Real>(m)
{
	if (n == 0) n = 1;	// ensure that workspace has positive size
	allocate(n);
}

template <class Real>
Workspace<Real>::~Workspace()
{
	delete[] alist;
	delete[] blist;
	delete[] rlist;
	delete[] elist;
	delete[] order;
	delete[] level;
}

template <class Real>
void
Workspace<Real>::append_interval (Real a1, Real b1, Real area1, Real error1)
{
	alist[size] = a1;
	blist[size] = b1;
	rlist[size] = area1;
	elist[size] = error1;
	order[size] = size;
	level[size] = 0;

	size++;
}

template <class Real>
void
Workspace<Real>::initialise (Real a, Real b)
{
	size = 0;
	nrmax = 0;
	i_work = 0;
	alist[0] = a;
	blist[0] = b;
	rlist[0] = Real(0);
	elist[0] = Real(0);
	order[0] = 0;
	level[0] = 0;

	maximum_level = 0;
}

template <class Real>
void
Workspace<Real>::sort_results ()
{
	size_t i;

	for (i = 0; i < size; i++)
	{
      size_t i1 = order[i];
      Real e1 = elist[i1];
      size_t i_max = i1;
      size_t j;

      for (j = i + 1; j < size; j++)
		{
			size_t i2 = order[j];
			Real e2 = elist[i2];

			if (e2 >= e1)
			{
				i_max = i2;
				e1 = e2;
			}
		}

      if (i_max != i1)
		{
			order[i] = order[i_max];
			order[i_max] = i1;
		}
	}

	i_work = order[0] ;
}

template <class Real>
void
Workspace<Real>::qpsrt ()
{
	const size_t last = size - 1;

	Real errmax ;
	Real errmin ;
	int i, k, top;

	size_t i_nrmax = nrmax;
	size_t i_maxerr = order[i_nrmax] ;

	/* Check whether the list contains more than two error estimates */

	if (last < 2)
	{
      order[0] = 0 ;
      order[1] = 1 ;
      i_work = i_maxerr ;
      return ;
	}

	errmax = elist[i_maxerr] ;

	/* This part of the routine is only executed if, due to a difficult
	 integrand, subdivision increased the error estimate. In the normal
	 case the insert procedure should start after the nrmax-th largest
	 error estimate. */

	while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]])
	{
      order[i_nrmax] = order[i_nrmax - 1] ;
      i_nrmax-- ;
	}

	/* Compute the number of elements in the list to be maintained in
	 descending order. This number depends on the number of
	 subdivisions still allowed. */

	if(last < (limit/2 + 2))
	{
      top = last ;
	}
	else
	{
      top = limit - last + 1;
	}

	/* Insert errmax by traversing the list top-down, starting
	 comparison from the element elist(order(i_nrmax+1)). */

	i = i_nrmax + 1 ;

	/* The order of the tests in the following line is important to
	 prevent a segmentation fault */

	while (i < top && errmax < elist[order[i]])
	{
      order[i-1] = order[i] ;
      i++ ;
	}

	order[i-1] = i_maxerr ;

	/* Insert errmin by traversing the list bottom-up */

	errmin = elist[last] ;

	k = top - 1 ;

	while (k > i - 2 && errmin >= elist[order[k]])
	{
      order[k+1] = order[k] ;
      k-- ;
	}

	order[k+1] = last ;

	/* Set i_max and e_max */

	i_maxerr = order[i_nrmax] ;

	i_work = i_maxerr ;
	nrmax = i_nrmax ;
}

template <class Real>
void
Workspace<Real>::set_initial_result (Real result, Real error)
{
	size = 1;
	rlist[0] = result;
	elist[0] = error;
}

template <class Real>
void
Workspace<Real>::retrieve (Real& a, Real& b, Real& r, Real& e)
{
	a = alist[i_work] ;
	b = blist[i_work] ;
	r = rlist[i_work] ;
	e = elist[i_work] ;
}

template <class Real>
Real
Workspace<Real>::sum_results ()
{
	size_t k;
	Real result_sum = Real(0);

	for (k = 0; k < size; k++)
	{
      result_sum += rlist[k];
	}
	return result_sum;
}

template <class Real>
int
Workspace<Real>::subinterval_too_small (Real a1, Real a2, Real b2)
{
	const Real e = this->eps_;
	const Real u = this->xmin_;

	Real tmp = (1 + 100 * e) * (this->abs(a2) + 1000 * u);
	int status = (this->abs(a1) <= tmp && this->abs(b2) <= tmp);
	return status;
}

template <class Real>
void
Workspace<Real>::update (Real a1, Real b1, Real area1, Real error1,
								 Real a2, Real b2, Real area2, Real error2)
{
	const size_t i_max = i_work ;
	const size_t i_new = size ;

	const size_t new_level = level[i_max] + 1;

	/* append the newly-created intervals to the list */

	if (error2 > error1)
	{
      alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
      rlist[i_max] = area2;
      elist[i_max] = error2;
      level[i_max] = new_level;

      alist[i_new] = a1;
      blist[i_new] = b1;
      rlist[i_new] = area1;
      elist[i_new] = error1;
      level[i_new] = new_level;
	}
	else
	{
      blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
      rlist[i_max] = area1;
      elist[i_max] = error1;
      level[i_max] = new_level;

      alist[i_new] = a2;
      blist[i_new] = b2;
      rlist[i_new] = area2;
      elist[i_new] = error2;
      level[i_new] = new_level;
	}

	size++;

	if (new_level > maximum_level)
	{
      maximum_level = new_level;
	}

	qpsrt () ;
}

/* -------------------------------------------------------------------------- */
template <class Real>
int Workspace<Real>::qag(FtnBase<Real>& f, Real a, Real b,
							Real epsabs, Real epsrel, Real& result, Real& abserr)
{
	Real area, errsum;
	Real result0, abserr0, resabs0, resasc0;
	Real tolerance;
	size_t iteration = 0;
	int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

	Real round_off;

	/* Initialize results */

	initialise (a, b);

	result = Real(0);
	abserr = Real(0);

	if (epsabs <= Real(0) && epsrel < 50 * this->eps_)
	{
      static const char*
		message = "tolerance cannot be acheived with given epsabs and epsrel";
		throw message;
		return GSL_EBADTOL;
	}

	/* perform the first integration */

	this->qk(f, a, b, result0, abserr0, resabs0, resasc0);

	set_initial_result (result0, abserr0);

	/* Test on accuracy */

	tolerance = this->max (epsabs, epsrel * this->abs (result0));

	/* need IEEE rounding here to match original quadpack behavior
	 *
	 * round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);
	 */
	round_off = 50 * this->eps_ * resabs0;

	if (abserr0 <= round_off && abserr0 > tolerance)
	{
      result = result0;
      abserr = abserr0;

      static const char*
		message = "cannot reach tolerance because of roundoff error on first attempt";
		throw message;
		return GSL_EROUND;
	}
	else if ((abserr0 <= tolerance && abserr0 != resasc0)
				|| abserr0 == Real(0))
	{
      result = result0;
      abserr = abserr0;

      return GSL_SUCCESS;
	}
	else if (limit == 1)
	{
      result = result0;
      abserr = abserr0;

      static const char*
		message = "a maximum of one iteration was insufficient";
		throw message;
		return GSL_EMAXITER;
	}

	area = result0;
	errsum = abserr0;

	iteration = 1;

	do
	{
      Real a1, b1, a2, b2;
      Real a_i, b_i, r_i, e_i;
      Real area1 = Real(0), area2 = Real(0), area12 = Real(0);
      Real error1 = Real(0), error2 = Real(0), error12 = Real(0);
      Real resasc1, resasc2;
      Real resabs1, resabs2;

      /* Bisect the subinterval with the largest error estimate */

      retrieve(a_i, b_i, r_i, e_i);

      a1 = a_i;
      b1 = (a_i + b_i) / Real(2);
      a2 = b1;
      b2 = b_i;

      this->qk(f, a1, b1, area1, error1, resabs1, resasc1);
      this->qk(f, a2, b2, area2, error2, resabs2, resasc2);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_i);
      area += area12 - r_i;

		/*
		 * ISSUE: Should checking for roundoff errors using QUADPACK's
		 * tolerances for double precision be disabled for multiple
		 * precision arithmetic.
		 *
		 */

      if (resasc1 != error1 && resasc2 != error2)
		{
			Real delta = r_i - area12;

			if (this->abs(delta) <= 1.0e-5 * this->abs(area12)
				 && error12 >= 0.99 * e_i)
			{
				roundoff_type1++;
			}
			if (iteration >= 10 && error12 > e_i)
			{
				roundoff_type2++;
			}
		}

		tolerance = this->max (epsabs, epsrel * this->abs(area));

      if (errsum > tolerance)
		{
			if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
			{
				error_type = 2;   /* round off error */
			}

			/* set error flag in the case of bad integrand behaviour at
			 a point of the integration range */

			if (subinterval_too_small (a1, a2, b2))
			{
				error_type = 3;
			}
		}

      update (a1, b1, area1, error1, a2, b2, area2, error2);

      retrieve(a_i, b_i, r_i, e_i);

      iteration++;

	}
	while (iteration < limit && !error_type && errsum > tolerance);

	result = sum_results();
	abserr = errsum;

	if (errsum <= tolerance)
	{
      return GSL_SUCCESS;
	}
	else if (error_type == 2)
	{
      static const char*
		message = "roundoff error prevents tolerance from being achieved";
		throw message;
		return GSL_EROUND;
	}
	else if (error_type == 3)
	{
      static const char*
		message = "bad integrand behavior found in the integration interval";
		throw message;
		return GSL_ESING;
	}
	else if (iteration == limit)
	{
      static const char*
		message = "maximum number of subdivisions reached";
		throw message;
		return GSL_EMAXITER;
	}
	else
	{
      static const char*
		message = "could not integrate function";
		throw message;
		return GSL_EFAILED;
	}
}

#endif	// QUADPACKPP_WORKSPACE_HPP
