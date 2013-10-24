
from distutils.core import setup, Extension

setup(name="_heston", version="0.01",
      ext_modules=[ Extension("_heston", ["_heston.c", "heston.c"], 
		libraries=["m", "gsl", "gslcblas"],
		include_dirs=["/usr/local/include"], 
		library_dirs=["/usr/local/lib"]) ])

