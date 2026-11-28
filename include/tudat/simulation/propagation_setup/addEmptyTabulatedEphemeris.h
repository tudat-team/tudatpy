/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ADDEMPTYTABULATEDEPHEMERIS_H
#define TUDAT_ADDEMPTYTABULATEDEPHEMERIS_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
void addEmptyTabulatedEphemeris( const simulation_setup::SystemOfBodies& bodies,
                                 const std::string& bodyName,
                                 const std::string& ephemerisOrigin = "",
                                 const bool isPartOfMultiArc = false )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting empty tabulated ephemeris for body " + bodyName + ", no such body found" );
    }
    std::string ephemerisOriginToUse = ( ephemerisOrigin == "" ) ? bodies.getFrameOrigin( ) : ephemerisOrigin;
    std::shared_ptr< ephemerides::Ephemeris > tabulatedEphemeris =
            std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >( ),
                    ephemerisOriginToUse,
                    bodies.getFrameOrientation( ) );
    if( !isPartOfMultiArc )
    {
        bodies.at( bodyName )->setEphemeris( tabulatedEphemeris );
    }
    else
    {
        bodies.at( bodyName )
                ->setEphemeris( std::make_shared< ephemerides::MultiArcEphemeris >(
                        std::map< double, std::shared_ptr< ephemerides::Ephemeris > >( { { 0.0, tabulatedEphemeris } } ),
                        ephemerisOriginToUse,
                        bodies.getFrameOrientation( ) ) );
    }

    bodies.processBodyFrameDefinitions< StateScalarType, TimeType >( );
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_ADDEMPTYTABULATEDEPHEMERIS_H
