// fjutil.hh
#ifndef HEPPYY_FJUTIL_HH
#define HEPPYY_FJUTIL_HH

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <vector>
#include <fastjet/PseudoJet.hh>
#include <numpy/arrayobject.h>
#include <stdexcept>

namespace alian
{
	void process_numpy_array(PyObject *array);

	// Static initializer class
	class NumpyInitializer
	{
	public:
		NumpyInitializer()
		{
			if (_do_init() != 0)
			{
				throw std::runtime_error("NumPy initialization failed");
			}
		}

	private:
		int _do_init()
		{
			if (_import_array() < 0)
			{
				PyErr_Print();
				PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
				return -1;
			}
			return 0;
		}
	};

	// Static instance to trigger initialization
	static NumpyInitializer numpy_initializer;

	std::vector<fastjet::PseudoJet> numpy_pxpypz_to_pseudojets(PyObject *px, PyObject *py, PyObject *pz, double m);
	std::vector<fastjet::PseudoJet> numpy_ptetaphi_to_pseudojets(PyObject *pt, PyObject *eta, PyObject *phi, double m);
}

#endif