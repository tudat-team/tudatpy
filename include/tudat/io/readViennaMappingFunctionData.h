/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READVMFDATA_H
#define TUDAT_READVMFDATA_H

#include <vector>
#include <string>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace input_output
{

struct VMFData {
public:
    std::map< double, std::vector< double > > delayData;
    std::map< double, std::vector< double > > meteoData;
    std::map< double, std::vector< double > > gradientData;

    void validate( );

    void getFullDataSet( std::map< double, Eigen::VectorXd >& processedTroposphereData,
                         std::map< double, Eigen::VectorXd >& processedMeteoData );
};

void readVMFFile( const std::string& fileName,
                  std::map< std::string, VMFData >& vmfData,
                  const bool fileHasMeteo,
                  const bool fileHasGradient );

void readVMFFiles( const std::vector< std::string >& fileName,
                   std::map< std::string, VMFData >& vmfData,
                   const bool fileHasMeteo,
                   const bool fileHasGradient );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_READVMFDATA_H
