#include "tudat/basics/deprecationWarnings.h"

namespace tudat
{

namespace utilities
{

void printDeprecationError(
        const std::string& name,
        const std::string& descriptionPage )
{
    std::string errorMessage = "Error, the " + name + " function/class is no longer available. To understand how to modify your code to the new version, go to the page " +
            descriptionPage + ", or open an issue on Github (https://github.com/tudat-team/tudatpy) if this page is not sufficient.";
    throw std::runtime_error( errorMessage );
}

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

}

}

