#ifndef TUDAT_APPLICATIONOUTPUT_H
#define TUDAT_APPLICATIONOUTPUT_H

namespace tudat_applications
{

//! Get path for output directory.
static inline std::string getOutputPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    return filePath_.substr( 0, filePath_.length( ) -
                                std::string( "applicationOutput.h" ).length( ) );
}

}

#endif // TUDAT_APPLICATIONOUTPUT_H
