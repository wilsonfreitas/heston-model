#ifndef  GSLCPP_GCOMPLEX_H
#define  GSLCPP_GCOMPLEX_H

/** 
 * @file gcomplex.h
 * @brief Complex number support
 * @author Anshul Nigham
 * @date 2006-11-02
 */

#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

namespace gslcpp {
	
	/** 
	 * @brief Class to represent complex numbers
	 *
	 * Note: In this documentation, an object of this class is referred to as \f$z\f$.
	 *
	 * \todo Check to ensure that gsl_complex struct does not contain pointers and is safe to copy
	 */
	class gcomplex {
	private:
		gsl_complex data_;
		gcomplex(gsl_complex z);
		
	public:
		/** 
		 * @brief Default constructor. Sets the complex number to \f$z=0 + 0i\f$
		 */
		gcomplex();
		/** 
		 * @brief Returns a complex number initialized from rectangular coordinates
		 * 
		 * @param x Real part
		 * @param y Imaginary part
		 * 
		 * @return New complex number \f$z=x + yi\f$
		 */
		static gcomplex rect(double x, double y);
		/** 
		 * @brief Returns a complex number initialized from polar coordinates
		 * 
		 * @param r The radius
		 * @param theta The angle
		 * 
		 * @return Returns the complex number \f$z=re^{i\theta} = r(\cos(\theta)+i\sin(\theta))\f$
		 */
		static gcomplex polar(double r, double theta);
		/** 
		 * @brief Returns the real part.
		 */
		double real() const;
		/** 
		 * @brief Returns the imaginary part.
		 */
		double imag() const;
		/** 
		 * @brief Set the complex number value based on rectangular coordinates
		 * 
		 * @param x Real part
		 * @param y Imaginary part
		 */
		void set(double x, double y);
		/** 
		 * @brief Set the real part
		 */
		void set_real(double a);
		/** 
		 * @brief Set the imaginary part
		 */
		void set_imag(double a);
		/** 
		 * @brief Returns the argument of the complex number.
		 */
		double arg() const;
		/** 
		 * @brief Returns the absolute value of the complex number.
		 */
		double abs() const;
		/** 
		 * @brief Returns the square of the absolute value of the complex number.
		 */
		double abs2() const;
		/** 
		 * @brief Returns the natural logarithm of the magnitude of the complex number.
		 */
		double logabs() const;
		/** 
		 * @brief Returns the natural logarithm of the magnitude of the complex number.
		 */
		gcomplex exp() const;
		static gcomplex exp(const gcomplex& c);
		gcomplex log() const;
		static gcomplex log(const gcomplex& c);
		gcomplex pow(double r) const;
		static gcomplex pow(const gcomplex& c, double r);
		static gcomplex sqrt(const gcomplex& c);
		
//		ostream& operator<<(ostream& o, const gcomplex& c); // output
//		istream& operator>>(istream& i, gcomplex& c);
		/** 
		 * @brief For complex c, returns a new complex number \f$z_n = z+c\f$
		 */
		gcomplex operator+(const gcomplex& c) const;
		/** 
		 * @brief For complex c, returns a new complex number \f$z_n = z-c\f$
		 */
		gcomplex operator-(const gcomplex& c) const;
		/** 
		 * @brief For complex c, returns a new complex number \f$z_n = z*c\f$
		 */
		gcomplex operator*(const gcomplex& c) const;
		/** 
		 * @brief For complex c, returns a new complex number \f$z_n = z/c\f$
		 */
		gcomplex operator/(const gcomplex& c) const;
		/** 
		 * @brief For real x, returns a new complex number \f$z_n = z+x\f$
		 */
		gcomplex operator+(double x) const;
		/** 
		 * @brief For real x, returns a new complex number \f$z_n = z-x\f$
		 */
		gcomplex operator-(double x) const;
		/** 
		 * @brief For real x, returns a new complex number \f$z_n = z*x\f$
		 */
		gcomplex operator*(double x) const;
		/** 
		 * @brief For real x, returns a new complex number \f$z_n = z+x\f$
		 */
		gcomplex operator/(double x) const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = z+iy\f$
		 */
		gcomplex add_imag(double y) const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = z-iy\f$
		 */
		gcomplex sub_imag(double y) const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = z*iy\f$
		 */
		gcomplex mul_imag(double y) const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = z/iy\f$
		 */
		gcomplex div_imag(double y) const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = z^* = x - iy\f$
		 */
		gcomplex conjugate() const;
		/** 
		 * @brief For real y, returns a new complex number \f$z_n = 1/z\f$
		 */
		gcomplex inverse() const;
		/** 
		 * @brief Returns a new complex number \f$z_n = -z\f$
		 */
		gcomplex operator-(void) const;
		
		// trig functions
		
		/** 
		 * @brief Returns the complex sine of the complex number \f$z, \sin(z) = (e^{iz} - e^{-iz})/(2i)\f$.
		 */
		gcomplex sin() const;
		/** 
		 * @brief Returns the complex cosine of the complex number \f$z, \cos(z) = (e^{iz} + e^{-iz})/2\f$.
		 */
		gcomplex cos() const;
		/** 
		 * @brief Returns the complex tangent of the complex number \f$z, \tan(z) = \sin(z)/\cos(z)\f$.
		 */
		gcomplex tan() const;
		/** 
		 * @brief Returns the complex secant of the complex number \f$z, \sec(z) = 1/\cos(z)\f$.
		 */
		gcomplex sec() const;
		/** 
		 * @brief Returns the complex cosecant of the complex number \f$z, \csc(z) = 1/\sin(z)\f$.
		 */
		gcomplex csc() const;
		/** 
		 * @brief Returns the complex cotangent of the complex number \f$z, \cot(z) = 1/\tan(z)\f$.
		 */
		gcomplex cot() const;
		
		// inverse trig functions
		
		/** 
		 * @brief Returns the complex arcsine of \f$z, \arcsin(z)\f$. The branch cuts are on the real axis, less than -1 and greater than 1.
		 */
		gcomplex arcsin() const;
		/** 
		 * @brief Returns the complex arcsine of the real number \f$x, \arcsin(x)\f$. For x between -1 and 1, the function returns a real value in the range \f$[-\pi/2,\pi/2]\f$. For x less than -1 the result has a real part of \f$-\pi/2\f$ and a positive imaginary part. For x greater than 1 the result has a real part of \f$\pi/2\f$ and a negative imaginary part.
		 */
		static gcomplex arcsin(double x);
		/** 
		 * @brief Returns the complex arccosine of \f$z, \arccos(z)\f$. The branch cuts are on the real axis, less than -1 and greater than 1.
		 */
		gcomplex arccos() const;
		/** 
		 * @brief Returns the complex arccosine of the real number \f$x, \arccos(x)\f$. For x between -1 and 1, the function returns a real value in the range \f$[0,\pi/2]\f$. For x less than -1 the result has a real part of \f$\pi\f$ and a negative imaginary part. For x greater than 1 the result is purely imaginary and positive.
		 */
		static gcomplex arccos(double x);
		/** 
		 * @brief Returns the complex arctangent of \f$z, \arctan(z)\f$. The branch cuts are on the imaginary axis, below \f$-i\f$ and above \f$-i\f$.
		 */
		gcomplex arctan() const;
		/** 
		 * @brief Returns the complex arcsecant of \f$z, mathrm{arcsec}(z) = \arccos(1/z)\f$.
		 */
		gcomplex arcsec() const;
		/** 
		 * @brief Returns the complex arcsecant of real number \f$x, \mathrm{arcsec}(x) = \arccos(1/x)\f$. 
		 */
		static gcomplex arcsec(double x);
		/** 
		 * @brief Returns the complex arccosecant of \f$z, \mathrm{arccsc}(z) = \arcsin(1/z)\f$.
		 */
		gcomplex arccsc() const;
		/** 
		 * @brief Returns the complex arccosecant of real number \f$x, \mathrm{arccsc}(x) = \arcsin(1/x)\f$. 
		 */
		static gcomplex arccsc(double x);
		/** 
		 * @brief Returns the complex arccotangent of \f$z, \mathrm{arccot}(z) = \arctan(1/z)\f$. 
		 */
		gcomplex arccot() const;
		
		// hyperbolic functions
		
		/** 
		 * @brief Returns the complex hyperbolic sine of the complex hyperbolic number \f$z, \sin(z) = (e^{z} - e^{-z})/2\f$.
		 */
		gcomplex sinh() const;
		/** 
		 * @brief Returns the complex hyperbolic cosine of the complex hyperbolic number \f$z, \cos(z) = (e^{z} + e^{-z})/2\f$.
		 */
		gcomplex cosh() const;
		/** 
		 * @brief Returns the complex hyperbolic tangent of the complex hyperbolic number \f$z, \tanh(z) = \sinh(z)/\cosh(z)\f$.
		 */
		gcomplex tanh() const;
		/** 
		 * @brief Returns the complex hyperbolic secant of the complex hyperbolic number \f$z, \mathrm{sech}(z) = 1/\cosh(z)\f$.
		 */
		gcomplex sech() const;
		/** 
		 * @brief Returns the complex hyperbolic cosecant of the complex hyperbolic number \f$z, \mathrm{csch}(z) = 1/\sinh(z)\f$.
		 */
		gcomplex csch() const;
		/** 
		 * @brief Returns the complex hyperbolic cotangent of the complex hyperbolic number \f$z, \mathrm{coth}(z) = 1/\tanh(z)\f$.
		 */
		gcomplex coth() const;
		
		// inverse hypoerbolic functions
		
		/** 
		 * @brief Returns the complex hyperbolic arcsine of \f$z, \mathrm{arcsinh}(z)\f$. The branch cuts are on the imaginary axis, below \f$-i\f$ and above \f$-i\f$.
		 */
		gcomplex arcsinh() const;
		/** 
		 * @brief Returns the complex hyperbolic arccosine of \f$z, \mathrm{arccosh}(z)\f$. The branch cuts are on the real axis, less than 1.
		 */
		gcomplex arccosh() const;
		/** 
		 * @brief Returns the complex hyperbolic arccosine of the real number \f$x, \mathrm{arccosh}(x)\f$. 
		 */
		static gcomplex arccosh(double x);
		/** 
		 * @brief Returns the complex hyperbolic arctangent of \f$z, \mathrm{arctanh}(z)\f$. The branch cuts are on the reak axis, less than -1 and greater than 1.
		 */
		gcomplex arctanh() const;
		/** 
		 * @brief Returns the complex hypoerbolic arctangent of real number \f$x, \mathrm{arctanh}(x)\f$. 
		 */
		static gcomplex arctanh(double x);
		/** 
		 * @brief Returns the complex hyperbolic arcsecant of \f$z, \mathrm{arcsech}(z) = \mathrm{arccosh}(1/z)\f$.
		 */
		gcomplex arcsech() const;
		/** 
		 * @brief Returns the complex hyperbolic arccosecant of \f$z, \mathrm{arccsch}(z) = \mathrm{arcsinh}(1/z)\f$.
		 */
		gcomplex arccsch() const;
		/** 
		 * @brief Returns the complex hyperbolic arccotangent of \f$z, \mathrm{arccoth}(z) = \mathrm{arctanh}(1/z)\f$. 
		 */
		gcomplex arccoth() const;
	};
}

#endif   /* ----- #ifndef GSLCPP_GCOMPLEX_H  ----- */