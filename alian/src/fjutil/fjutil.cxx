// fjutil.cxx
#include "fjutil.hh"
#include <iostream>

namespace alian
{
	void process_numpy_array(PyObject *array)
	{
		PyArrayObject *np_array = reinterpret_cast<PyArrayObject *>(array);
		if (PyArray_NDIM(np_array) != 1)
		{
			PyErr_SetString(PyExc_TypeError, "Array must be one-dimensional");
			return;
		}

		double *data = static_cast<double *>(PyArray_DATA(np_array));
		npy_intp size = PyArray_SIZE(np_array);

		for (npy_intp i = 0; i < size; ++i)
		{
			std::cout << "Element " << i << ": " << data[i] << std::endl;
		}
	}

	// function to transform three numpy arrays into a vector of PseudoJets
	std::vector<fastjet::PseudoJet> numpy_pxpypz_to_pseudojets(PyObject *px, PyObject *py, PyObject *pz, double m)
	{
		PyArrayObject *np_px = reinterpret_cast<PyArrayObject *>(px);
		PyArrayObject *np_py = reinterpret_cast<PyArrayObject *>(py);
		PyArrayObject *np_pz = reinterpret_cast<PyArrayObject *>(pz);

		if (PyArray_NDIM(np_px) != 1 || PyArray_NDIM(np_py) != 1 || PyArray_NDIM(np_pz) != 1)
		{
			PyErr_SetString(PyExc_TypeError, "Arrays must be one-dimensional");
			return std::vector<fastjet::PseudoJet>();
		}

		npy_intp size = PyArray_SIZE(np_px);
		if (size != PyArray_SIZE(np_py) || size != PyArray_SIZE(np_pz))
		{
			PyErr_SetString(PyExc_ValueError, "Arrays must have the same size");
			return std::vector<fastjet::PseudoJet>();
		}

		float *px_data = static_cast<float *>(PyArray_DATA(np_px));
		float *py_data = static_cast<float *>(PyArray_DATA(np_py));
		float *pz_data = static_cast<float *>(PyArray_DATA(np_pz));

		std::vector<fastjet::PseudoJet> particles;
		particles.reserve(size);

		for (npy_intp i = 0; i < size; ++i)
		{
			double E = sqrt(px_data[i] * px_data[i] + py_data[i] * py_data[i] + pz_data[i] * pz_data[i] + m * m);
			particles.emplace_back(px_data[i], py_data[i], pz_data[i], E);
		}

		return particles;
	}

	// function to transform three numpy arrays pt,eta,phi into a vector of PseudoJets
	std::vector<fastjet::PseudoJet> numpy_ptetaphi_to_pseudojets(PyObject *pt, PyObject *eta, PyObject *phi, double m)
	{
		PyArrayObject *np_pt = reinterpret_cast<PyArrayObject *>(pt);
		PyArrayObject *np_eta = reinterpret_cast<PyArrayObject *>(eta);
		PyArrayObject *np_phi = reinterpret_cast<PyArrayObject *>(phi);

		if (PyArray_NDIM(np_pt) != 1 || PyArray_NDIM(np_eta) != 1 || PyArray_NDIM(np_phi) != 1)
		{
			PyErr_SetString(PyExc_TypeError, "Arrays must be one-dimensional");
			return std::vector<fastjet::PseudoJet>();
		}

		npy_intp size = PyArray_SIZE(np_pt);
		if (size != PyArray_SIZE(np_eta) || size != PyArray_SIZE(np_phi))
		{
			PyErr_SetString(PyExc_ValueError, "Arrays must have the same size");
			return std::vector<fastjet::PseudoJet>();
		}

		float *pt_data  = static_cast<float *>(PyArray_DATA(np_pt));
		float *eta_data = static_cast<float *>(PyArray_DATA(np_eta));
		float *phi_data = static_cast<float *>(PyArray_DATA(np_phi));

		std::vector<fastjet::PseudoJet> particles;
		particles.reserve(size);

		for (npy_intp i = 0; i < size; ++i)
		{
			if (pt_data[i] <= 0.001) // this is 1 MeV (!)
			{
				// PyErr_SetString(PyExc_ValueError, "pt must be positive");
				// return std::vector<fastjet::PseudoJet>();
				continue;
			}
			double px = pt_data[i] * cos(phi_data[i]);
			double py = pt_data[i] * sin(phi_data[i]);
			double pz = pt_data[i] * sinh(eta_data[i]);
			double E = sqrt(px * px + py * py + pz * pz + m * m);

			particles.emplace_back(px, py, pz, E);
			particles.back().set_user_index(i);
		}

		return particles;
	}
}