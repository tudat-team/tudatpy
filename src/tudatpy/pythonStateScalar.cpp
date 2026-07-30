#include "pythonStateScalar.h"

#include <stdexcept>

namespace tudatpy
{
namespace
{

PythonStateScalarType activePythonStateScalarType = PythonStateScalarType::double_precision;

}  // namespace

void setPythonStateScalarType( const std::string& requestedType )
{
    if( requestedType == "double" )
    {
        activePythonStateScalarType = PythonStateScalarType::double_precision;
    }
    else if( requestedType == "quad" )
    {
        if( !quadPrecisionPythonBindingsAvailable( ) )
        {
            throw std::runtime_error( "This TudatPy kernel was built without quad-precision Python bindings." );
        }
        activePythonStateScalarType = PythonStateScalarType::quad_precision;
    }
    else
    {
        throw std::runtime_error( "Unknown TudatPy state scalar type '" + requestedType + "'." );
    }
}

PythonStateScalarType getPythonStateScalarType( )
{
    return activePythonStateScalarType;
}

std::string getPythonStateScalarTypeName( )
{
    return activePythonStateScalarType == PythonStateScalarType::quad_precision ? "quad" : "double";
}

bool quadPrecisionPythonBindingsAvailable( )
{
#if TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE
    return true;
#else
    return false;
#endif
}

}  // namespace tudatpy
