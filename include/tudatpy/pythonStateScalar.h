#ifndef TUDATPY_PYTHON_STATE_SCALAR_H
#define TUDATPY_PYTHON_STATE_SCALAR_H

#include <string>

namespace tudatpy
{

enum class PythonStateScalarType { double_precision, quad_precision };

void setPythonStateScalarType( const std::string& requestedType );

PythonStateScalarType getPythonStateScalarType( );

std::string getPythonStateScalarTypeName( );

bool quadPrecisionPythonBindingsAvailable( );

}  // namespace tudatpy

#endif  // TUDATPY_PYTHON_STATE_SCALAR_H
