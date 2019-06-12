#ifndef _EIGEN_NUMPY_H_
#define _EIGEN_NUMPY_H_

#include <boost/python.hpp>

#if PY_VERSION_HEX >= 0x03000000
void* SetupEigenConverters();
#else
void SetupEigenConverters();
#endif

#endif
