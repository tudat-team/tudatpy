/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Basics/timeType.h"
#include "Tudat/SimulationSetup/PropagationSetup/setNumericallyIntegratedStates.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

//! Function checking feasibility of resetting the translational dynamics
void checkTranslationalStatesFeasibility(
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    // Check feasibility of epheme                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ris origins.
    for( simulation_setup::NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), bodyIterator->first ) ==
                bodiesToIntegrate.end( ) )
        {
            std::string ephemerisOrigin
                    = bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( );

            if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), ephemerisOrigin )
                != bodiesToIntegrate.end( ) )
            {
                std::cerr << "Warning, found non-integrated body with an integrated body as ephemeris origin" +
                             bodyIterator->second->getEphemeris( )->getReferenceFrameOrigin( ) + " " +
                             bodyIterator->first << std::endl;
            }
        }

    }

    // Check whether each integrated body exists, and whether it has a TabulatedEphemeris
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        std::string bodyToIntegrate = bodiesToIntegrate.at( i );

        if( bodyMap.count( bodyToIntegrate ) == 0 )
        {
            if( bodyMap.at( bodyToIntegrate )->getEphemeris( ) == NULL )
            {
                throw std::runtime_error( "Error when checking translational dynamics feasibility of body " +
                                          bodyToIntegrate + "no such body found" );
            }
        }
        else
        {
            if( bodyMap.at( bodyToIntegrate )->getEphemeris( ) == NULL )
            {
                throw std::runtime_error( "Error when checking translational dynamics feasibility of body " +
                                          bodyToIntegrate + " no ephemeris found" );
            }
        }
    }

}

//! Function to create an interpolator for the new translational state of a body.
template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return boost::make_shared<
        interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >( stateMap, 6 );
}

//! Function to create an interpolator for the new translational state of a body.
template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >
createStateInterpolator( const std::map< double, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return boost::make_shared<
        interpolators::LagrangeInterpolator< double,
                                             Eigen::Matrix< long double, 6, 1 > > >( stateMap, 6 );
}

//! Function to create an interpolator for the new translational state of a body.
template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > >
createStateInterpolator( const std::map< Time, Eigen::Matrix< long double, 6, 1 > >& stateMap )
{
    return boost::make_shared<
        interpolators::LagrangeInterpolator<
            Time, Eigen::Matrix< long double, 6, 1 >, long double > >( stateMap, 6 );
}


//! Function to create an interpolator for the new translational state of a body.
template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > >
createStateInterpolator( const std::map< Time, Eigen::Matrix< double, 6, 1 > >& stateMap )
{
    return boost::make_shared<
        interpolators::LagrangeInterpolator<
            Time, Eigen::Matrix< double, 6, 1 >, long double > >( stateMap, 6 );
}


template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( stateMap, 6 );
}

template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< double, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::Matrix< long double, 7, 1 > > >( stateMap, 6 );
}

template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< Time, Eigen::Matrix< double, 7, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< double, 7, 1 >, long double > >( stateMap, 6 );
}

template< >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 7, 1 > > >
createRotationalStateInterpolator( const std::map< Time, Eigen::Matrix< long double, 7, 1 > >& stateMap )
{
    return boost::make_shared< interpolators::LagrangeInterpolator< Time, Eigen::Matrix< long double, 7, 1 >, long double > >( stateMap, 6 );
}


} // namespace propagators

} // namespace tudat
