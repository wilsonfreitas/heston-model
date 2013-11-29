

#include <Python.h>
#include "heston.h"

static PyObject *
ucall(PyObject *self, PyObject *args)
{
	double F;
	double K;
	double tau;
	double v;
	double vbar;
	double lambda;
	double eta;
	double rho;
	
	if (!PyArg_ParseTuple(args, "dddddddd", &F, &K, &tau, &v, &vbar, &lambda, &eta, &rho))
		return NULL;  
	
	double h = heston_ucall(F, K, tau, v, vbar, lambda, eta, rho);
	return PyFloat_FromDouble(h);
}

static PyObject *
call(PyObject *self, PyObject *args)
{
	double S;
	double K;
	double tau;
	double r;
	double q;
	double v;
	double vbar;
	double lambda;
	double eta;
	double rho;
	
	if (!PyArg_ParseTuple(args, "dddddddddd", &S, &K, &tau, &r, &q, &v, &vbar, &lambda, &eta, &rho))
		return NULL;  
	double h = heston_call(S, K, tau, r, q, v, vbar, lambda, eta, rho);
	return PyFloat_FromDouble(h);
}

PyMethodDef methods[] = {
    {"ucall", (PyCFunction)ucall, METH_VARARGS, "Returns undiscounted value of an European call option."},
    {"call", (PyCFunction)call, METH_VARARGS, "Returns value of an European call option."},
    {NULL, NULL, 0, NULL}
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC 
init_heston(void)
{
    (void) Py_InitModule("_heston", methods);
}


