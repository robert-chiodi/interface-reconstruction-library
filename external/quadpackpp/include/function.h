/* function.h
 * Part of the quadpack++.
 */

#ifndef QUADPACKPP_FUNCTION_H
#define QUADPACKPP_FUNCTION_H

//! C++ version of <tt>gsl_function</tt> passed to quadrature routines.
template<typename Real>
class FtnBase {
public:
	virtual Real operator() (Real x) =0;
};

/** \brief C++ version of <tt>gsl_function</tt> with constructors.

 User-defined functions for numerical routines may require additional parameters
 depending on the application. However function-pointers, when passed as
 arguments to C-routines, depend on their "signature" type. To accomodate a
 universal signature, the GSL uses the struct <tt>gsl_function</tt> with the
 member <tt>void* params</tt> for all possible parameter/structure types.

 The approach in quadpack++ is analgous except that the parameter type is
 templated. The universality is accomodated in that only the base type,
 \ref FtnBase, is passed as an argument to routines.

 <b>Note:</b> the parameter \ref params_ is must be a pointer instead of a
 refererence since this object is declared outside of the scope.
 */
template<typename Real, class param_t>
class Function : public FtnBase<Real> {
public:
	//! Signature defining function \ref FtnBase::operator()
	typedef Real defn_t(Real, param_t*);

	//! Reference to definition
	defn_t& function_;
	//! Pointer to current parameters.
	param_t* params_;

	//! Overload the \ref FtnBase::operator() with \ref function_.
	virtual Real operator() (Real x)
	{
		return function_(x, params_);
	}

	//! Constructor without parameters.
	Function(defn_t& function) : function_(function), params_(0) {}
	//! General constructor.
	Function(defn_t& function, param_t* params) : function_(function), params_(params) {}
	~Function() {}
};

#endif // QUADPACKPP_FUNCTION_H
