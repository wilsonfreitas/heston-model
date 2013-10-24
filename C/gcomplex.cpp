/** 
 * @file gcomplex.cpp
 * @brief Complex number support
 * @author Anshul Nigham
 * @date 2006-11-02
 */

#include <gsl/gsl_complex.h>
#include "gcomplex.h"

namespace gslcpp {
	gcomplex::gcomplex()
	{
		data_ = gsl_complex_rect(0.0,0.0);
	}
	
	gcomplex::gcomplex(gsl_complex z) : data_(z) { }
	
	gcomplex gcomplex::rect(double x, double y)
	{
		gsl_complex z = gsl_complex_rect(x,y);
		return gcomplex(z);
	}
	
	gcomplex gcomplex::polar(double r, double theta)
	{
		gsl_complex z = gsl_complex_polar(r,theta);
		return gcomplex(z);
	}
	
	double gcomplex::real() const
	{
		return GSL_REAL(data_);
	}
	
	double gcomplex::imag() const
	{
		return GSL_IMAG(data_);
	}
	
	void gcomplex::set(double x, double y)
	{
		GSL_SET_COMPLEX(&data_, x, y);
	}
	
	void gcomplex::set_real(double a)
	{
		GSL_SET_REAL(&data_, a);
	}
	
	void gcomplex::set_imag(double a)
	{
		GSL_SET_IMAG(&data_, a);
	}
	
	// elementary functions
	double gcomplex::arg() const
	{
		return gsl_complex_arg(data_);
	}
	
	double gcomplex::abs() const
	{
		return gsl_complex_abs(data_);
	}
	
	double gcomplex::abs2() const
	{
		return gsl_complex_abs2(data_);
	}
	
	double gcomplex::logabs() const
	{
		return gsl_complex_logabs(data_);
	}
	
	gcomplex gcomplex::exp() const
	{
		return gcomplex(gsl_complex_exp(data_));
	}
	
	gcomplex gcomplex::exp(const gcomplex& c)
	{
		return gcomplex(gsl_complex_exp(c.data_));
	}
	
	gcomplex gcomplex::log() const
	{
		return gcomplex(gsl_complex_log(data_));
	}
	
	gcomplex gcomplex::log(const gcomplex& c)
	{
		return gcomplex(gsl_complex_log(c.data_));
	}
	
	gcomplex gcomplex::pow(double r) const {
		return gcomplex(gsl_complex_pow_real(data_, r));
	}
	
	gcomplex gcomplex::pow(const gcomplex& c, double r) {
		return gcomplex(gsl_complex_pow_real(c.data_, r));
	}
	
	gcomplex gcomplex::sqrt(const gcomplex& c) {
		return gcomplex(gsl_complex_sqrt(c.data_));
	}
	
	
	// operators
//	ostream& operator<<(ostream& out, const gcomplex& c)
//	{
//		out << "(" << c.real() << ", " << c.imag() << ")";
//		return out;
//	}
//	istream& operator>>(istream& in, gcomplex& c)
//	{
//		return in;
//	}
	gcomplex gcomplex::operator+(const gcomplex& c) const
	{
		return gcomplex(gsl_complex_add(data_,c.data_));
	}
	gcomplex gcomplex::operator-(const gcomplex& c) const
	{
		return gcomplex(gsl_complex_sub(data_,c.data_));
	}
	gcomplex gcomplex::operator*(const gcomplex& c) const
	{
		return gcomplex(gsl_complex_mul(data_,c.data_));
	}
	gcomplex gcomplex::operator/(const gcomplex& c) const
	{
		return gcomplex(gsl_complex_div(data_,c.data_));
	}
	gcomplex gcomplex::operator+(double x) const
	{
		return gcomplex(gsl_complex_add_real(data_,x));
	}
	gcomplex gcomplex::operator-(double x) const
	{
		return gcomplex(gsl_complex_sub_real(data_,x));
	}
	gcomplex gcomplex::operator*(double x) const
	{
		return gcomplex(gsl_complex_mul_real(data_,x));
	}
	gcomplex gcomplex::operator/(double x) const
	{
		return gcomplex(gsl_complex_div_real(data_,x));
	}
	
	gcomplex gcomplex::add_imag(double y) const
	{
		return gcomplex(gsl_complex_add_imag(data_,y));
	}
	gcomplex gcomplex::sub_imag(double y) const
	{
		return gcomplex(gsl_complex_sub_imag(data_,y));
	}
	gcomplex gcomplex::mul_imag(double y) const
	{
		return gcomplex(gsl_complex_mul_imag(data_,y));
	}
	gcomplex gcomplex::div_imag(double y) const
	{
		return gcomplex(gsl_complex_div_imag(data_,y));
	}
	gcomplex gcomplex::conjugate() const
	{
		return gcomplex(gsl_complex_conjugate(data_));
	}
	gcomplex gcomplex::inverse() const
	{
		return gcomplex(gsl_complex_inverse(data_));
	}
	gcomplex gcomplex::operator-(void) const
	{
		return gcomplex(gsl_complex_negative(data_));
	}
	
	// trigonometric functions
	gcomplex gcomplex::sin() const
	{
		return gcomplex(gsl_complex_sin(data_));
	}
	gcomplex gcomplex::cos() const
	{
		return gcomplex(gsl_complex_cos(data_));
	}
	gcomplex gcomplex::tan() const
	{
		return gcomplex(gsl_complex_tan(data_));
	}
	gcomplex gcomplex::sec() const
	{
		return gcomplex(gsl_complex_sec(data_));
	}
	gcomplex gcomplex::csc() const
	{
		return gcomplex(gsl_complex_csc(data_));
	}
	gcomplex gcomplex::cot() const
	{
		return gcomplex(gsl_complex_cot(data_));
	}
	
	// inverse trignometric functions
	gcomplex gcomplex::arcsin() const
	{
		return gcomplex(gsl_complex_arcsin(data_));
	}
	gcomplex gcomplex::arcsin(double x)
	{
		return gcomplex(gsl_complex_arcsin_real(x));
	}
	gcomplex gcomplex::arccos() const
	{
		return gcomplex(gsl_complex_arccos(data_));
	}
	gcomplex gcomplex::arccos(double x)
	{
		return gcomplex(gsl_complex_arccos_real(x));
	}
	gcomplex gcomplex::arctan() const
	{
		return gcomplex(gsl_complex_arctan(data_));
	}
	gcomplex gcomplex::arcsec() const
	{
		return gcomplex(gsl_complex_arcsec(data_));
	}
	gcomplex gcomplex::arcsec(double x)
	{
		return gcomplex(gsl_complex_arcsec_real(x));
	}
	gcomplex gcomplex::arccsc() const
	{
		return gcomplex(gsl_complex_arccsc(data_));
	}
	gcomplex gcomplex::arccsc(double x)
	{
		return gcomplex(gsl_complex_arccsc_real(x));
	}
	gcomplex gcomplex::arccot() const
	{
		return gcomplex(gsl_complex_arccot(data_));
	}
	
	// hyperbolic functions
	gcomplex gcomplex::sinh() const
	{
		return gcomplex(gsl_complex_sinh(data_));
	}
	gcomplex gcomplex::cosh() const
	{
		return gcomplex(gsl_complex_cosh(data_));
	}
	gcomplex gcomplex::tanh() const
	{
		return gcomplex(gsl_complex_tanh(data_));
	}
	gcomplex gcomplex::sech() const
	{
		return gcomplex(gsl_complex_sech(data_));
	}
	gcomplex gcomplex::csch() const
	{
		return gcomplex(gsl_complex_csch(data_));
	}
	gcomplex gcomplex::coth() const
	{
		return gcomplex(gsl_complex_coth(data_));
	}
	
	// inverse hyperbolic functions
	gcomplex gcomplex::arcsinh() const
	{
		return gcomplex(gsl_complex_arcsinh(data_));
	}
	gcomplex gcomplex::arccosh() const
	{
		return gcomplex(gsl_complex_arccosh(data_));
	}
	gcomplex gcomplex::arccosh(double x)
	{
		return gcomplex(gsl_complex_arccosh_real(x));
	}
	gcomplex gcomplex::arctanh() const
	{
		return gcomplex(gsl_complex_arctanh(data_));
	}
	gcomplex gcomplex::arctanh(double x)
	{
		return gcomplex(gsl_complex_arctanh_real(x));
	}
	gcomplex gcomplex::arcsech() const
	{
		return gcomplex(gsl_complex_arcsech(data_));
	}
	gcomplex gcomplex::arccsch() const
	{
		return gcomplex(gsl_complex_arccsch(data_));
	}
	gcomplex gcomplex::arccoth() const
	{
		return gcomplex(gsl_complex_arccoth(data_));
	}
}