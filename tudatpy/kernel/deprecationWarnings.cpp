#include "deprecationWarnings.h"

namespace tudatpy {

void printDeprecationWarning(
        const std::string& oldName,
        const std::string& newName,
        const std::string& description )
{
    std::cerr<<"Warning, the function "<<oldName<<
               " is deprecated, and will be removed in a future version. Please use the functionally equivalent "<<
               newName<<" instead"<<std::endl;
    if( description != "" )
    {
        std::cerr<<description<<std::endl;
    }
}


} // namespace tudatpy

