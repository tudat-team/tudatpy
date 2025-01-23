/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SETNUMERICALLYINTEGRATEDSTATES_H
#define TUDAT_SETNUMERICALLYINTEGRATEDSTATES_H

#include "tudat/basics/utilities.h"
#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

#include "tudat/math/interpolators/lagrangeInterpolator.h"


namespace tudat
{

namespace propagators
{


template< typename StateScalarType, typename TimeType >
void addEmptyTabulatedEphemeris(
    const simulation_setup::SystemOfBodies& bodies, const std::string& bodyName, const std::string& ephemerisOrigin = "",
    const bool isPartOfMultiArc = false )
{
    if( bodies.count( bodyName ) ==  0 )
    {
        throw std::runtime_error( "Error when setting empty tabulated ephemeris for body " + bodyName + ", no such body found" );
    }
    std::string ephemerisOriginToUse = ( ephemerisOrigin == "" ) ? bodies.getFrameOrigin( ) : ephemerisOrigin;
    std::shared_ptr< ephemerides::Ephemeris > tabulatedEphemeris =
        std::make_shared< ephemerides::TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
            std::shared_ptr< interpolators::OneDimensionalInterpolator
                < TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >( ), ephemerisOriginToUse, bodies.getFrameOrientation( ) );
    if( !isPartOfMultiArc )
    {
        bodies.at( bodyName )->setEphemeris( tabulatedEphemeris );
    }
    else
    {
        bodies.at( bodyName )->setEphemeris( std::make_shared< ephemerides::MultiArcEphemeris >(
            std::map< double, std::shared_ptr< ephemerides::Ephemeris > >(
                { { 0.0, tabulatedEphemeris } } ), ephemerisOriginToUse, bodies.getFrameOrientation( ) ) );
    }

    bodies.processBodyFrameDefinitions( );
}



template< typename StateScalarType, typename TimeType >
void addEmptyTabulatedRotationalEphemeris(
    const simulation_setup::SystemOfBodies& bodies, const std::string& bodyName, const std::string& bodyFixedFrameName = ""  )
{
    if( bodies.count( bodyName ) ==  0 )
    {
        throw std::runtime_error( "Error when setting empty tabulated rotational ephemeris for body " + bodyName + ", no such body found" );
    }
    std::string bodyFixedFrameNameToUse = ( bodyFixedFrameName == "" ) ? ( bodyName + "_fixed" ) : bodyFixedFrameName;

    bodies.at( bodyName )->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< StateScalarType, TimeType > >(
        std::shared_ptr< interpolators::OneDimensionalInterpolator
            < TimeType, Eigen::Matrix< StateScalarType, 7, 1 > > >( ), bodies.getFrameOrientation( ), bodyFixedFrameNameToUse ) );

}

//! Function to create an interpolator for the new translational state of a body.
/*!
 * Function to create an interpolator for the new translational state of a body.
 * \param stateMap New state history, w.r.t. the required ephemeris origin.
 * \return Lagrange interpolator (order 6) that produces the required continuous state.
 */
template< typename TimeType, typename StateScalarType >
std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
createStateInterpolator(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& stateMap );

//! Function to reset the tabulated ephemeris of a body
/*!
 * Function to reset the tabulated ephemeris of a body
 * \param ephemerisInput New state history that is to be set
 * \param tabulatedEphemeris Ephemeris in which the ephemerisInput is to be set.
 */
template< typename StateTimeType, typename StateScalarType, typename EphemerisTimeType, typename EphemerisScalarType  >
void resetIntegratedEphemerisOfBody(
        const std::map< StateTimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& ephemerisInput,
        const std::shared_ptr< ephemerides::TabulatedCartesianEphemeris< EphemerisScalarType, EphemerisTimeType > > tabulatedEphemeris )
{
    std::map< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 6, 1 > > castEphemerisInput;
    utilities::castMatrixMap< StateTimeType, StateScalarType, EphemerisTimeType, EphemerisScalarType, 6, 1 >(
                ephemerisInput, castEphemerisInput );
    
    std::shared_ptr< interpolators::OneDimensionalInterpolator< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 6, 1 > > >
            ephemerisInterpolator = createStateInterpolator( castEphemerisInput );
    tabulatedEphemeris->resetInterpolator( ephemerisInterpolator );
}

//! Function to reset the tabulated ephemeris of a body
/*!
 * Function to reset the tabulated ephemeris of a body, this requires the requested body to possess
 * an ephemeris of type TabulatedCartesianEphemeris< StateScalarType, TimeType >
 * \param bodies List of bodies used in simulations.
 * \param ephemerisInput New state history of the body
 * \param bodyToIntegrate Name of body for which the ephemeris is to be reset.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerisOfBody(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& ephemerisInput,
        const std::string& bodyToIntegrate )
{
    using namespace tudat::interpolators;
    using namespace tudat::ephemerides;
    
    // If body does not have ephemeris, give error message.
    if( bodies.at( bodyToIntegrate )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when resetting integrated ephemeris of body " +
                                  bodyToIntegrate + "no ephemeris found" );
    }
    
    // If current ephemeris is not already a tabulated ephemeris, give error message.
    else if( !isTabulatedEphemeris( bodies.at( bodyToIntegrate )->getEphemeris( ) ) )
    {
        throw std::runtime_error( "Error when resetting integrated ephemeris of body " +
                                  bodyToIntegrate + " no tabulated ephemeris found" );
        
    }
    // Else, update existing tabulated ephemeris
    else
    {
        if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                    bodies.at( bodyToIntegrate )->getEphemeris( ) ) != nullptr )
        {
            std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
                    ephemerisInterpolator = createStateInterpolator( ephemerisInput );
            std::shared_ptr< TabulatedCartesianEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
                    std::dynamic_pointer_cast< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                        bodies.at( bodyToIntegrate )->getEphemeris( ) );
            tabulatedEphemeris->resetInterpolator( ephemerisInterpolator );
        }
        else
        {
            if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
                        bodies.at( bodyToIntegrate )->getEphemeris( ) ) != nullptr )
            {
                resetIntegratedEphemerisOfBody(
                            ephemerisInput, std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
                                bodies.at( bodyToIntegrate )->getEphemeris( ) ) );
            }
            else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >(
                         bodies.at( bodyToIntegrate )->getEphemeris( ) ) != nullptr )
            {
                resetIntegratedEphemerisOfBody(
                            ephemerisInput, std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >(
                                bodies.at( bodyToIntegrate )->getEphemeris( ) ) );
            }
            else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >(
                         bodies.at( bodyToIntegrate )->getEphemeris( ) ) != nullptr )
            {
                resetIntegratedEphemerisOfBody(
                            ephemerisInput, std::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >(
                                bodies.at( bodyToIntegrate )->getEphemeris( ) ) );
            }
            else if( std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >(
                         bodies.at( bodyToIntegrate )->getEphemeris( ) ) != nullptr )
            {
                resetIntegratedEphemerisOfBody(
                            ephemerisInput, std::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >(
                                bodies.at( bodyToIntegrate )->getEphemeris( ) ) );
            }
            else
            {
                throw std::runtime_error( "Error, no tabulated ephemeris found when resetting ephemeris" );
            }
        }
    }
}

//! Function to convert output of translational motion to input for the ephemeris.
/*!
 * Function to convert output of translational motion from the numerical integrator to the required
 * input for the ephemeris.  It extracts the state history of a single body from the full list of
 * integrated states Additionally, it changes the origin of the reference frame in which the states
 * are given, by using the integrationToEphemerisFrameFunction input variable.
 * \param bodyIndex Index of integrated body for which the state is to be retrieved
 * \param startIndex Index in entries of equationsOfMotionNumericalSolution where the translational states start.
 * \param equationsOfMotionNumericalSolution Full numerical solution of numerical integrator,
 * already converted to Cartesian states (w.r.t. the integration origin of the body of bodyIndex)
 * \param ephemerisTable State history of body bodyIndex w.r.t. the origin with which its ephemeris is defined
 * (returned by reference).
 * \param integrationToEphemerisFrameFunction Function to provide the state of the ephemeris origin
 * of the current body w.r.t. its integration origin.
*/
template< typename TimeType, typename StateScalarType >
void convertNumericalSolutionToEphemerisInput(
        const int bodyIndex,
        const int startIndex,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
        equationsOfMotionNumericalSolution,
        std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& ephemerisTable,
        const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >
        integrationToEphemerisFrameFunction = nullptr )
{
    // If no integrationToEphemerisFrameFunction is provided, origin is already correct; only
    // extract required indices.
    if( integrationToEphemerisFrameFunction == 0 )
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
             bodyIterator = equationsOfMotionNumericalSolution.begin( );
             bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {
            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( startIndex + 6 * bodyIndex, 0, 6, 1 );
        }
    }
    // Else, extract indices and add required translation from integrationToEphemerisFrameFunction
    else
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
             bodyIterator = equationsOfMotionNumericalSolution.begin( );
             bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {
            
            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( startIndex + 6 * bodyIndex, 0, 6, 1 ) -
                    integrationToEphemerisFrameFunction( bodyIterator->first );
        }
    }
}

//! Function to extract the numerical solution for the translational dynamics of a single body from full propagation history.
/*!
 * Function to extract the numerical solution for the translational dynamics of a single body from full propagation history.
 * Function can perform frame translation if required.
 * \param bodiesToIntegrate List of names of bodies which are numericall integrated (in the order in
 * which they are in the equationsOfMotionNumericalSolution map.
 * \param translationalStateStartIndex Index in entries of equationsOfMotionNumericalSolution where the translational states start
 * \param bodyForWhichToRetrieveState Name of body for which the states are to be extracted
 * \param equationsOfMotionNumericalSolution Numerical solution of dynamics, with translational results in Cartesian elements
 * w.r.t. integratation origins.
 * \param ephemerisInput State history of requested body (returned by reference)
 * \param bodyIndex Index of bodyForWhichToRetrieveState in bodiesToIntegrate (returned by reference)
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris
 * origins of each body w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void getSingleBodyStateHistoryFromPropagationOutpiut(
        const std::vector< std::string >& bodiesToIntegrate,
        const int translationalStateStartIndex,
        const std::string& bodyForWhichToRetrieveState,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& ephemerisInput,
        int& bodyIndex,
        const std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
        integrationToEphemerisFrameFunctions =
        std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    // Get index of current body to be updated in bodiesToIntegrate.
    std::vector< std::string >::const_iterator bodyFindIterator = std::find(
                bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), bodyForWhichToRetrieveState );
    if( bodyFindIterator == bodiesToIntegrate.end( ) )
    {
        throw std::runtime_error( "Error when creating and setting ephemeris after integration, cannot find body " +
                                  bodyForWhichToRetrieveState );
    }
    bodyIndex = std::distance( bodiesToIntegrate.begin( ), bodyFindIterator );
    
    // Get frame origin function if applicable
    std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > integrationToEphemerisFrameFunction = nullptr;
    if( integrationToEphemerisFrameFunctions.count( bodiesToIntegrate.at( bodyIndex ) ) > 0 )
    {
        integrationToEphemerisFrameFunction =
                integrationToEphemerisFrameFunctions.at( bodiesToIntegrate.at( bodyIndex ) );
    }
    
    // Create and reset interpolator.
    convertNumericalSolutionToEphemerisInput(
                bodyIndex, translationalStateStartIndex, equationsOfMotionNumericalSolution, ephemerisInput, integrationToEphemerisFrameFunction );
}

//! Create and reset ephemerides interpolator
/*!
 * Creates and resets the interpolator for the ephemerides of the integrated bodies from the
 * numerical integration results.
 * \param bodies List of bodies used in simulations.
 * \param bodiesToIntegrate List of names of bodies which are numericall integrated (in the order in
 * which they are in the equationsOfMotionNumericalSolution map.
 * \param startIndex Index in entries of equationsOfMotionNumericalSolution where the translational states start.
 * \param ephemerisUpdateOrder Order in which to update the ephemeris objects.
 * \param equationsOfMotionNumericalSolution Numerical solution of translational equations of
 * motion, in Cartesian elements w.r.t. integratation origins.
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris
 * origins of each body w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void createAndSetInterpolatorsForEphemerides(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const int startIndex,
        const std::vector< std::string >& ephemerisUpdateOrder,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
        integrationToEphemerisFrameFunctions =
        std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    using namespace tudat::interpolators;
    
    std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > ephemerisInput;
    int bodyIndex;
    
    // Iterate over all bodies that are integrated numerically and create state interpolator.
    for( unsigned int i = 0; i < ephemerisUpdateOrder.size( ); i++ )
    {
        ephemerisInput.clear( );
        getSingleBodyStateHistoryFromPropagationOutpiut(
                    bodiesToIntegrate, startIndex, ephemerisUpdateOrder.at( i ), equationsOfMotionNumericalSolution,
                    ephemerisInput, bodyIndex, integrationToEphemerisFrameFunctions );
        resetIntegratedEphemerisOfBody(
                    bodies, ephemerisInput, bodiesToIntegrate.at( bodyIndex ) );
    }
}

//! Resets the ephemerides of the integrated bodies from the numerical integration results.
/*!
 * Resets the ephemerides of the integrated bodies from the numerical integration results, and
 * performs associated computation for ephemeris-dependent environment variables.
 * \param bodies List of bodies used in simulations.
 * \param bodiesToIntegrate List of names of bodies which are numerically integrated (in the order in
 * which they are in the equationsOfMotionNumericalSolution map.
 * \param startIndexAndSize Pair with start index and total (contiguous) size of integrated states in entries of
 * equationsOfMotionNumericalSolution
 * \param ephemerisUpdateOrder Order in which to update the ephemeris objects (empty if arbitrary).
 * \param equationsOfMotionNumericalSolution Numerical solution of translational equations of
 * motion, in Cartesian elements w.r.t. integratation origins.
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris
 * origins of each body w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerides(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::vector< std::string >& bodiesToIntegrate,
        const std::pair< unsigned int, unsigned int > startIndexAndSize,
        std::vector< std::string > ephemerisUpdateOrder = std::vector< std::string >( ),
        const std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
        integrationToEphemerisFrameFunctions =
        std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    // Set update order arbitrarily if no order is provided.
    if( ephemerisUpdateOrder.size( ) == 0 )
    {
        ephemerisUpdateOrder = bodiesToIntegrate;
    }
    // Check input consistency
    else if( ephemerisUpdateOrder.size( ) != bodiesToIntegrate.size( ) )
    {
        throw std::runtime_error( "Error when resetting ephemerides, input vectors have inconsistent size" );
    }
    
    if( static_cast< unsigned int >( equationsOfMotionNumericalSolution.begin( )->second.rows( ) )
            < startIndexAndSize.first + startIndexAndSize.second )
    {
        throw std::runtime_error( "Error when resetting ephemerides, input solution inconsistent with start index and size." );
    }
    
    if( startIndexAndSize.second != 6 * bodiesToIntegrate.size( ) )
    {
        throw std::runtime_error( "Error when resetting ephemerides, number of bodies inconsistent with input size." );
        
    }
    
    // Create interpolators from numerical integration results (states) at discrete times.
    createAndSetInterpolatorsForEphemerides(
                bodies, bodiesToIntegrate, startIndexAndSize.first, ephemerisUpdateOrder,
                equationsOfMotionNumericalSolution, integrationToEphemerisFrameFunctions );
}

//! Resets the ephemerides of the integrated bodies from the numerical multi-arc integration results.
/*!
 * Resets the ephemerides of the integrated bodies from the numerical multi-arc integration results, and
 * performs associated computation for ephemeris-dependent environment variables.
 * \param bodies List of bodies used in simulations.
 * \param equationsOfMotionNumericalSolution Numerical multi-arc solution of translational equations of
 * motion, in Cartesian elements w.r.t. integratation origins (one vector entry represents one arc).
 * \param arcStartTimes List of start times of the propagation arcs.
 * \param bodiesToIntegrate List of names of bodies which are numerically integrated (in the order in
 * which they are in the equationsOfMotionNumericalSolution map.
 * \param startIndexAndSize Pair with start index and total (contiguous) size of integrated states in entries of
 * equationsOfMotionNumericalSolution
 * \param ephemerisUpdateOrder Order in which to update the ephemeris objects (empty if arbitrary).
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris
 * origins of each body w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void resetMultiArcIntegratedEphemerides(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >&
        equationsOfMotionNumericalSolution,
        const std::vector< double > arcStartTimes,
        const std::vector< std::vector< std::string > >& bodiesToIntegrate,
        const std::pair< unsigned int, unsigned int > startIndexAndSize,
        std::vector< std::vector< std::string > > ephemerisUpdateOrder = std::vector< std::vector< std::string > >( ),
        const std::map< std::string, std::vector< std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > >&
        multiArcIntegrationToEphemerisFrameFunctions =
        std::map< std::string, std::vector< std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > >( ) )
{
    using namespace tudat::interpolators;
    using namespace tudat::ephemerides;

    // Check input consistency
    if( ephemerisUpdateOrder.size( ) == 0 )
    {
        ephemerisUpdateOrder = bodiesToIntegrate;
    }


    //    // Check that, for each arc, the bodies to integrate and ephemeris update order are identical (for now).
    //    for ( unsigned int i = 1 ; i < bodiesToIntegrate.size( ) ; i++ )
    //    {
    //        for ( unsigned int k = 0 ; k < bodiesToIntegrate.at( i ).size( ) ; k++ )
    //        {
    //            if ( ( bodiesToIntegrate.at( 0 ).at( k ) != bodiesToIntegrate.at( i ).at( k ) ) ||
    //                 ( ephemerisUpdateOrder.at( 0 ).at( k ) != ephemerisUpdateOrder.at( i ).at( k ) ) )
    //            {
    //                throw std::runtime_error( "Error when resetting multi-arc integrated ephemerides, the bodies to integrate and/or "
    //                                          " ephemeris update order are inconsistent in-between arcs." );
    //            }
    //        }
    //    }

    //    std::cout << "bodies to integrate size per arc!" << "\n\n";
    //    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    //    {
    //        std::cout << "arc " << i <<  " = " << bodiesToIntegrate.at( i ).size( ) << "\n\n";
    //    }

    // Update ephemerides
    std::map< std::string, std::vector< std::shared_ptr< Ephemeris > > > arcEphemerisListPerBody;
    std::map< std::string, std::vector< double > > arcStartingTimesPerBody;
    std::map< std::string, int > counterArcPerBody;
    for ( unsigned int arc = 0 ; arc < arcStartTimes.size( ) ; arc++ )
    {
        for ( unsigned int i = 0 ; i < ephemerisUpdateOrder.at( arc ).size( ) ; i++ )
        {
            // Find index of current body in bodiesToIntegrate.
            std::vector< std::string >::const_iterator bodyFindIterator = std::find(
                        bodiesToIntegrate.at( arc ).begin( ), bodiesToIntegrate.at( arc ).end( ), ephemerisUpdateOrder.at( arc ).at( i ) );
            int bodyIndex = std::distance( bodiesToIntegrate.at( arc ).begin( ), bodyFindIterator );

            // Update counter for nb of arc associated with given body.
            if ( counterArcPerBody.count( bodiesToIntegrate.at( arc ).at( bodyIndex ) ) == 0 )
            {
                counterArcPerBody[ bodiesToIntegrate.at( arc ).at( bodyIndex ) ] = 0;
            }
            else
            {
                counterArcPerBody[ bodiesToIntegrate.at( arc ).at( bodyIndex ) ] += 1;
            }
            int counterArcCurrentBody = counterArcPerBody[ bodiesToIntegrate.at( arc ).at( bodyIndex ) ];

            // Check environment consistency
            std::shared_ptr<MultiArcEphemeris> currentBodyEphemeris =
                    std::dynamic_pointer_cast< MultiArcEphemeris >(
                        bodies.at( bodiesToIntegrate.at( arc ).at( bodyIndex ) )->getEphemeris( ) );
            if ( currentBodyEphemeris == nullptr )
            {
                throw std::runtime_error("Error when resetting ephemeris of body " + bodiesToIntegrate.at( arc ).at( bodyIndex ) +
                                         ", original ephemeris is of incompatible type");
            }

            //            std::vector<std::shared_ptr<Ephemeris> > arcEphemerisList;
            //            for ( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
            //            {
            std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >
                    integrationToEphemerisFrameFunction = nullptr;

            // Create transformation function, if needed.
            if ( multiArcIntegrationToEphemerisFrameFunctions.count( bodiesToIntegrate.at( arc ).at( bodyIndex ) ) > 0 )
            {
                if ( multiArcIntegrationToEphemerisFrameFunctions.at( bodiesToIntegrate.at( arc ).at( bodyIndex ) ).at( counterArcCurrentBody /*arc*/ ) != nullptr )
                {
                    integrationToEphemerisFrameFunction =
                            multiArcIntegrationToEphemerisFrameFunctions.at( bodiesToIntegrate.at( arc ).at( bodyIndex ) ).at( counterArcCurrentBody /*arc*/ );
                }
            }

            std::map<TimeType, Eigen::Matrix<StateScalarType, 6, 1> > currentArcSolution;
            convertNumericalSolutionToEphemerisInput(
                        bodyIndex, startIndexAndSize.first,
                        equationsOfMotionNumericalSolution.at( arc ), currentArcSolution, integrationToEphemerisFrameFunction );

            // Create interpolator.
            std::shared_ptr<OneDimensionalInterpolator< TimeType, Eigen::Matrix<StateScalarType, 6, 1 > > >
                    ephemerisInterpolator = createStateInterpolator( currentArcSolution );

            arcEphemerisListPerBody[ bodiesToIntegrate.at( arc ).at( bodyIndex ) ].push_back(
                        std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                            ephemerisInterpolator, currentBodyEphemeris->getReferenceFrameOrigin( ),
                            currentBodyEphemeris->getReferenceFrameOrientation( ) ) );

            arcStartingTimesPerBody[ bodiesToIntegrate.at( arc ).at( bodyIndex ) ].push_back( arcStartTimes.at( arc ) );

            //                arcEphemerisList.push_back( std::make_shared< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
            //                        ephemerisInterpolator, currentBodyEphemeris->getReferenceFrameOrigin( ),
            //                        currentBodyEphemeris->getReferenceFrameOrientation( ) ) );
            //            }
            currentBodyEphemeris->resetSingleArcEphemerides( arcEphemerisListPerBody.at( bodiesToIntegrate.at( arc ).at( bodyIndex ) ) /*arcEphemerisList*/,
                                                             arcStartingTimesPerBody.at( bodiesToIntegrate.at( arc ).at( bodyIndex ) ) /*arcStartTimes*/ );

        }
    }

    // Having set new ephemerides, update body properties depending on ephemerides.
    for( auto bodyIterator : bodies.getMap( )  )
    {
        bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
    }
}

//! Function to reset the tabulated rotational ephemeris of a body
/*!
 * Function to reset the tabulatedrotational  ephemeris of a body, this requires the requested body to possess
 * a rotational ephemeris of type TabulatedRotationalEphemeris< StateScalarType, TimeType >
 * \param bodies List of bodies used in simulations.
 * \param rotationalEphemerisInterpolator New rotational state history of the body
 * \param bodyToIntegrate Name of body for which the rotational ephemeris is to be reset.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedRotationalEphemerisOfBody(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > > >
        rotationalEphemerisInterpolator,
        const std::string bodyToIntegrate )
{
    using namespace tudat::interpolators;
    using namespace tudat::ephemerides;
    
    if( bodies.at( bodyToIntegrate )->getRotationalEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error, no rotational ephemeris detected for body " +
                                  bodyToIntegrate + " when resetting ephemeris" );
    }
    // If current ephemeris is not already a tabulated ephemeris, create new ephemeris.
    else if( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< StateScalarType, TimeType > >(
                 bodies.at( bodyToIntegrate )->getRotationalEphemeris( ) ) == nullptr )
    {
        throw std::runtime_error( "Error when resetting integrated rotational ephemeris of body, rotation model type is incompatible " );
        
    }
    // Else, update existing tabulated ephemeris
    else
    {
        std::shared_ptr< TabulatedRotationalEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
                std::dynamic_pointer_cast< TabulatedRotationalEphemeris<  StateScalarType, TimeType > >(
                    bodies.at( bodyToIntegrate )->getRotationalEphemeris( ) );
        tabulatedEphemeris->reset( rotationalEphemerisInterpolator );
    }
}

//! Function to convert output of rotational motion to input for the rotational ephemeris.
/*!
 * Function to convert output of rotational motion from the numerical integrator to the required
 * input for the rotational ephemeris.  It extracts the state history of a single body from the full list of
 * integrated states.
 * \param startIndex Index in entries of equationsOfMotionNumericalSolution where the rotational states start.
 * \param bodyIndex Index of integrated body for which the state is to be retrieved
 * \param ephemerisTable State history of body bodyIndex w.r.t. the origin with which its ephemeris is defined
 * (returned by reference).
 * \param equationsOfMotionNumericalSolution Full numerical solution of numerical integrator,
 * already converted to Cartesian states (w.r.t. the integration origin of the body of bodyIndex)
*/
template< typename TimeType, typename StateScalarType >
void convertNumericalSolutionToRotationalEphemerisInput(
        const int startIndex,
        const int bodyIndex,
        std::map< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > >& ephemerisTable,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution )
{
    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator bodyIterator =
         equationsOfMotionNumericalSolution.begin( ); bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
    {
        ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( startIndex + 7 * bodyIndex, 0, 7, 1 );
    }
}

//! Function to create an interpolator for the new translational state of a body.
/*!
 * Function to create an interpolator for the new translational state of a body.
 * \param stateMap New state history, w.r.t. the required ephemeris origin.
 * \return Lagrange interpolator (order 6) that produces the required continuous state.
 */
template< typename TimeType, typename StateScalarType >
std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
createStateInterpolator(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& stateMap );

//! Function to create an interpolator for the new rotational state of a body.
/*!
 * Function to create an interpolator for the new rotational state of a body.
 * \param stateMap New rotational state history.
 * \return Lagrange interpolator (order 6) that produces the required continuous rotational state.
 */
template< typename TimeType, typename StateScalarType >
std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > > >
createRotationalStateInterpolator(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > >& stateMap );

//! Function to reset the tabulated rotational ephemeris of a body
/*!
 * Function to reset the tabulated rotational ephemeris of a body
 * \param bodies List of bodies used in simulations.
 * \param bodiesToIntegrate List of names of bodies for which rotational state is numerically integrated
 * \param startIndex Index in the state vector where the rotational state starts.
 * \param equationsOfMotionNumericalSolution New rotational state history that is to be set
 */
template< typename TimeType, typename StateScalarType >
void createAndSetInterpolatorsForRotationalEphemerides(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const int startIndex,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution )
{
    using namespace tudat::interpolators;
    
    // Iterate over all bodies that are integrated numerically and create state interpolator.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        std::map< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > > ephemerisInput;
        convertNumericalSolutionToRotationalEphemerisInput( startIndex, i, ephemerisInput, equationsOfMotionNumericalSolution );

        // Create interpolator.
        std::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 7, 1 > > >
                ephemerisInterpolator = createRotationalStateInterpolator(
                    ephemerisInput );
        
        resetIntegratedRotationalEphemerisOfBody( bodies, ephemerisInterpolator, bodiesToIntegrate.at( i ) );
    }
}


//! Resets the rotational ephemerides of a set of bodies from the numerical integration results.
/*!
 * Resets the rotational ephemerides of a set of bodies from the numerical integration results, and
 * performs associated computation for ephemeris-dependent environment variables.
 * \param bodies List of bodies used in simulations.
 * \param equationsOfMotionNumericalSolution Numerical solution of rotational equations of motion
 * \param bodiesToIntegrate List of names of bodies which are numerically integrated
 * \param startIndexAndSize Pair with start index and total (contiguous) size of integrated states in entries of
 * equationsOfMotionNumericalSolution
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedRotationalEphemerides(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::vector< std::string >& bodiesToIntegrate,
        const std::pair< unsigned int, unsigned int > startIndexAndSize )
{
    // Create interpolators from numerical integration results (states) at discrete times.
    createAndSetInterpolatorsForRotationalEphemerides(
                bodies, bodiesToIntegrate, startIndexAndSize.first, equationsOfMotionNumericalSolution );
    
    // Having set new ephemerides, update body properties depending on ephemerides.
    for( auto bodyIterator : bodies.getMap( )  )
    {
        //NOTE: Inefficient, should be done once following full integration.
        bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
    }
}

//! Resets the mass models of the integrated bodies from the numerical integration results.
/*!
 * Resets the mass models of the integrated bodies from the numerical integration results.
 * \param bodies List of bodies used in simulations.
 * \param equationsOfMotionNumericalSolution Numerical solution of the body masses.
 * \param bodiesToIntegrate List of names of bodies for which mass is numerically integrated (in the order in
 * which they are in the equationsOfMotionNumericalSolution map.
 * \param startIndexAndSize Pair with start index and total (contiguous) size of integrated states in entries of
 * equationsOfMotionNumericalSolution
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedBodyMass(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::vector< std::string >& bodiesToIntegrate ,
        const std::pair< unsigned int, unsigned int > startIndexAndSize )
{
    if( startIndexAndSize.second != bodiesToIntegrate.size( ) )
    {
        throw std::runtime_error( "Error when resetting body masses, number of bodies inconsistent with input size." );
    }
    
    // Iterate over all bodies for which mass is propagated.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        std::map< double, double > currentBodyMassMap;
        
        // Create mass map with double entries.
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
             stateIterator = equationsOfMotionNumericalSolution.begin( );
             stateIterator != equationsOfMotionNumericalSolution.end( ); stateIterator++ )
        {
            currentBodyMassMap[ static_cast< double >( stateIterator->first ) ] =
                    static_cast< double >( stateIterator->second( startIndexAndSize.first + i ) );
        }
        
        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
        
        // Create and set interpolator.
        bodies.at( bodiesToIntegrate.at( i ) )->setBodyMassFunction(
                    std::bind(
                        static_cast< double( LocalInterpolator::* )( const double ) >
                        ( &LocalInterpolator::interpolate ),
                        std::make_shared< interpolators::LagrangeInterpolatorDouble >( currentBodyMassMap, 6 ), std::placeholders::_1 ) );
    }
}

//! Base class for settings how numerically integrated states are processed
/*!
 *  Base class for defining settings on how numerically integrated states are to be processed in the
 *  environment. This class is pure virtual, and a derived class must be provided for each state
 *  type.
 */
template< typename TimeType, typename StateScalarType >
class IntegratedStateProcessor
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state that is to be set in environment.
     * \param startIndexAndSize Start index of current state type in integrated state vector, and
     * its size in the vector (first and second entry of pair).
     */
    IntegratedStateProcessor( const IntegratedStateType stateType,
                              const std::pair< int, int > startIndexAndSize ):
        stateType_( stateType ), startIndexAndSize_( startIndexAndSize ){ }

    //! Virtual destructor.
    virtual ~IntegratedStateProcessor( ){ }

    //! Type of state that is to be set in environment.
    IntegratedStateType stateType_;

    //! Start index of current state type in integrated state vector, and its size in the vector
    /*!
     *  Start index of current state type in integrated state vector, and its size in the vector
     * (first and second entry of pair).
     */
    std::pair< int, int > startIndexAndSize_;
};


//! Base class for settings how numerically integrated single-arc states are processed
/*!
 *  Base class for defining settings on how numerically integrated single-arc states are to be processed in the
 *  environment. This class is pure virtual, and a derived class must be provided for each state
 *  type.
 */
template< typename TimeType, typename StateScalarType >
class SingleArcIntegratedStateProcessor : public IntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state that is to be set in environment.
     * \param startIndexAndSize Start index of current state type in integrated state vector, and
     * its size in the vector (first and second entry of pair).
     */
    SingleArcIntegratedStateProcessor( const IntegratedStateType stateType,
                                       const std::pair< int, int > startIndexAndSize,
                                       const simulation_setup::SystemOfBodies& bodies,
                                       const std::vector< std::string >& bodiesToIntegrate ) :
        IntegratedStateProcessor< TimeType, StateScalarType >( stateType, startIndexAndSize ),
        bodies_( bodies ), bodiesToIntegrate_( bodiesToIntegrate ){ }

    //! Virtual destructor.
    virtual ~SingleArcIntegratedStateProcessor( ){ }

    //! Function that processes the entries of the stateType_ in the full numericalSolution
    /*!
     * Function that processes the entries of the stateType_ in the full numericalSolution
     * \param numericalSolution Full numerical solution, in global representation (see
     * convertToOutputSolution function in associated SingleStateTypeDerivative derived class.
     */
    virtual void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType,
            Eigen::Dynamic, 1 > >& numericalSolution ) = 0;

    //! List of bodies used in simulations.
    simulation_setup::SystemOfBodies bodies_;

    //! List of bodies for which the translational state is numerically integrated.
    /*!
     * List of bodies for which the translational state is numerically integrated. Order in this
     * vector is the same as the order in state vector.
     */
    std::vector< std::string > bodiesToIntegrate_;
};


//! Base class for settings how numerically integrated single-arc states are processed
/*!
 *  Base class for defining settings on how numerically integrated single-arc states are to be processed in the
 *  environment. This class is pure virtual, and a derived class must be provided for each state
 *  type.
 */
template< typename TimeType, typename StateScalarType >
class MultiArcIntegratedStateProcessor : public IntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state that is to be set in environment.
     * \param startIndexAndSize Start index of current state type in integrated state vector, and
     * its size in the vector (first and second entry of pair).
     */
    MultiArcIntegratedStateProcessor( const IntegratedStateType stateType,
                                      const std::pair< int, int > startIndexAndSize,
                                      const simulation_setup::SystemOfBodies& bodies,
                                      const std::vector< double > arcStartTimes ) :
        IntegratedStateProcessor< TimeType, StateScalarType >( stateType, startIndexAndSize ),
        bodies_( bodies ), arcStartTimes_( arcStartTimes ){ }

    //! Virtual destructor.
    virtual ~MultiArcIntegratedStateProcessor( ){ }

    virtual void processIntegratedMultiArcStates(
            const std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >& numericalSolution,
            const std::vector< double >& arcStartTimes ) = 0;


    //! List of bodies used in simulations.
    simulation_setup::SystemOfBodies bodies_;

    //! List of bodies for which the translational state is numerically integrated (per arc).
    /*!
     * List of bodies for which the translational state is numerically integrated (per arc). Order in this
     * vector is the same as the order in state vector.
     */
    std::vector< std::vector< std::string > > bodiesToIntegrate_;

    // List of arc-wise starting times.
    std::vector< double > arcStartTimes_;
};


//! Class used for processing numerically integrated translational states
/*!
 *  Class used for processing numerically integrated translational states, updates ephemeris object
 *  of each integrated body, performing frame translations if needed.
 */
template< typename TimeType, typename StateScalarType >
class TranslationalStateIntegratedStateProcessor: public SingleArcIntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param startIndex Index in the state vector where the translational state starts.
     * \param bodies List of bodies used in simulations.
     * \param bodiesToIntegrate List of bodies for which the translational state is numerically
     * integrated. Order in this vector is the same as the order in state vector.
     * \param centralBodies List of origing w.r.t. which the translational states are integrated
     * (with the same order as bodiesToIntegrate).
     * \param frameManager Object to get state of one body w.r.t. another body.
     */
    TranslationalStateIntegratedStateProcessor(
            const int startIndex,
            const simulation_setup::SystemOfBodies& bodies,
            const std::vector< std::string >& bodiesToIntegrate,
            const std::vector< std::string >& centralBodies,
            const std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager ):
        SingleArcIntegratedStateProcessor< TimeType, StateScalarType >(
            translational_state, std::make_pair( startIndex, 6 * bodiesToIntegrate.size( ) ), bodies, bodiesToIntegrate )
    {
        // Get update orders.
        ephemerisUpdateOrder_ = determineEphemerisUpdateorder(
                    this->bodiesToIntegrate_, centralBodies,
                    frameManager->getEphemerisOrigins( bodiesToIntegrate ) );

        // Get required frame origin translations.
        integrationToEphemerisFrameFunctions_ = ephemerides::getTranslationFunctionsFromIntegrationFrameToEphemerisFrame< StateScalarType, TimeType >(
                    centralBodies, this->bodiesToIntegrate_, frameManager );
    }

    ~TranslationalStateIntegratedStateProcessor( ){ }

    //! Function processing single-arc translational state, resetting bodies' ephemerides with new states
    /*!
     * Function processing single-arc translational state, resetting bodies' ephemerides with new states in numericalSolution
     * variable. It extracts and converts the states to the required frames, and updates the associated
     * ephemerides.
     * \param numericalSolution Full numerical solution, in global representation (see
     * convertToOutputSolution function in NBodyStateDerivative class.
     */
    void processIntegratedStates(
            const std::map< TimeType,
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution )
    {
        resetIntegratedEphemerides< TimeType, StateScalarType >(
                    this->bodies_, numericalSolution, this->bodiesToIntegrate_, this->startIndexAndSize_, ephemerisUpdateOrder_,
                    integrationToEphemerisFrameFunctions_ );
    }

    std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > getIntegrationToEphemerisFrameFunctions( )
    {
        return integrationToEphemerisFrameFunctions_;
    }

    std::vector< std::string > getBodiesToIntegrate( )
    {
        return this->bodiesToIntegrate_;
    }

    std::vector< std::string > getEphemerisUpdateOrder( )
    {
        return ephemerisUpdateOrder_;
    }

private:

    //! Order in which to update the ephemeris objects
    std::vector< std::string > ephemerisUpdateOrder_;

    //! Function to provide the states of the ephemeris origins of each body w.r.t. their respective
    //! integration origins.
    std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >
    integrationToEphemerisFrameFunctions_;

};


//! Class used for processing numerically integrated translational states
/*!
*  Class used for processing numerically integrated translational states, updates ephemeris object
*  of each integrated body, performing frame translations if needed.
*/
template< typename TimeType, typename StateScalarType >
class MultiArcTranslationalStateIntegratedStateProcessor: public MultiArcIntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param startIndex Index in the state vector where the translational state starts.
     * \param bodies List of bodies used in simulations.
     * \param bodiesToIntegrate List of bodies for which the translational state is numerically
     * integrated. Order in this vector is the same as the order in state vector.
     * \param centralBodies List of origing w.r.t. which the translational states are integrated
     * (with the same order as bodiesToIntegrate).
     * \param frameManager Object to get state of one body w.r.t. another body.
     */
    MultiArcTranslationalStateIntegratedStateProcessor(
            const simulation_setup::SystemOfBodies& bodies,
            std::vector< double >& arcStartTimes,
            const std::vector< std::shared_ptr<
            TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > > > singleArcTranslationalStateProcessors ):
        MultiArcIntegratedStateProcessor< TimeType, StateScalarType >(
            translational_state, singleArcTranslationalStateProcessors.at( 0 )->startIndexAndSize_, bodies, arcStartTimes ),
        singleArcTranslationalStateProcessors_( singleArcTranslationalStateProcessors )
    {

        // All bodies to integrate as keys, and the values are maps containing the general arc indices and arc indices specifically for the given body.
        std::map< std::string, std::map< int, int > > bodiesToIntegrateArcsAndIndices;

        unsigned int counterMultiArcStates = 0;
        for ( unsigned int i = 0 ; i < this->arcStartTimes_.size( ) ; i++ )
        {
            ephemerisUpdateOrder_.push_back(singleArcTranslationalStateProcessors_.at(i)->getEphemerisUpdateOrder());
            this->bodiesToIntegrate_.push_back(singleArcTranslationalStateProcessors_.at(i)->getBodiesToIntegrate());

            for ( unsigned int k = 0 ; k < singleArcTranslationalStateProcessors_.at( i )->getBodiesToIntegrate( ).size( ) ; k++ ) {
                std::string currentBodyToIntegrate = singleArcTranslationalStateProcessors_.at( i )->getBodiesToIntegrate( ).at( k );
                if ( bodiesToIntegrateArcsAndIndices.count( currentBodyToIntegrate ) == 0 )
                {
                    std::map< int, int > currentBodyArcAndIndex;
                    currentBodyArcAndIndex[ i ] = 0;
                    bodiesToIntegrateArcsAndIndices[ currentBodyToIntegrate ] = currentBodyArcAndIndex;
                }
                else
                {
                    int counterArcs = bodiesToIntegrateArcsAndIndices[ currentBodyToIntegrate ].rbegin( )->second;
                    bodiesToIntegrateArcsAndIndices.at( currentBodyToIntegrate )[ i ] = counterArcs + 1;
                }
            }

            counterMultiArcStates += singleArcTranslationalStateProcessors_.at( i )->getBodiesToIntegrate( ).size( );
        }
        this->startIndexAndSize_ = std::make_pair(  this->startIndexAndSize_.first, 6 * counterMultiArcStates );

        for ( unsigned int i = 0 ; i < this->arcStartTimes_.size( ) ; i++ )
        {
            std::shared_ptr<
                    TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > > currentTranslationalStateProcessor =
                    singleArcTranslationalStateProcessors_.at( i );

            std::map< std::string, std::function< Eigen::Matrix<StateScalarType, 6, 1 >( const TimeType ) > > singleArcIntegrationToEphemerisFrameFunctions
                    = currentTranslationalStateProcessor->getIntegrationToEphemerisFrameFunctions( );
            for ( auto itr: singleArcIntegrationToEphemerisFrameFunctions )
            {
                int currentArcIndexForBody = bodiesToIntegrateArcsAndIndices.at( itr.first ).at( i );
                if ( multiArcIntegrationToEphemerisFrameFunctions_.count( itr.first ) == 0 )
                {
                    std::vector< std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > ephemerisToFrameFunctions;
                    ephemerisToFrameFunctions.resize( bodiesToIntegrateArcsAndIndices.at( itr.first ).size( ) /*this->arcStartTimes_.size( )*/ );
                    ephemerisToFrameFunctions[ currentArcIndexForBody ] = itr.second;
                    multiArcIntegrationToEphemerisFrameFunctions_[ itr.first ] = ephemerisToFrameFunctions;
                }
                else
                {
                    multiArcIntegrationToEphemerisFrameFunctions_.at( itr.first ).at( currentArcIndexForBody ) = itr.second;
                }

            }
        }
    }

    ~MultiArcTranslationalStateIntegratedStateProcessor( ){ }


    //! Function processing multi-arc translational state, resetting bodies' ephemerides with new states
    /*!
     * Function processing multi-arc translational state, resetting bodies' ephemerides with new states in numericalSolution
     * variable. It extracts and converts the states to the required frames, and updates the associated
     * ephemerides.
     * \param numericalSolution Vector of numerical solutions (one vector per arc), in global representation (see
     * convertToOutputSolution function in NBodyStateDerivative class.
     * \param arcStartTimes List of start times of the propagation arcs.
     */
    void processIntegratedMultiArcStates(
            const std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >& numericalSolution,
            const std::vector< double >& arcStartTimes )
    {
        resetMultiArcIntegratedEphemerides< TimeType, StateScalarType >(
                    this->bodies_, numericalSolution, this->arcStartTimes_,
                    this->bodiesToIntegrate_, this->startIndexAndSize_, ephemerisUpdateOrder_, multiArcIntegrationToEphemerisFrameFunctions_ /*integrationToEphemerisFrameFunctions_*/ );
    }

private:
    //! Order in which to update the ephemeris objects
    std::vector< std::vector< std::string > > ephemerisUpdateOrder_;

    std::map< std::string, std::vector< std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > >
    multiArcIntegrationToEphemerisFrameFunctions_;

    std::vector< std::shared_ptr<
    TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > > > singleArcTranslationalStateProcessors_;

};


//! Class used for processing numerically integrated rotational states
/*!
 *  Class used for processing numerically integrated rotational states, updates rotational ephemeris object
 *  of each integrated body.
 */
template< typename TimeType, typename StateScalarType >
class RotationalStateIntegratedStateProcessor: public SingleArcIntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param startIndex Index in the state vector where the rotational state starts.
     * \param bodies List of bodies used in simulations.
     * \param bodiesToIntegrate List of bodies for which the rotational state is numerically
     * integrated. Order in this vector is the same as the order in state vector.
     */
    RotationalStateIntegratedStateProcessor(
            const int startIndex,
            const simulation_setup::SystemOfBodies& bodies,
            const std::vector< std::string >& bodiesToIntegrate ):
        SingleArcIntegratedStateProcessor< TimeType, StateScalarType >(
            rotational_state, std::make_pair( startIndex, 7 * bodiesToIntegrate.size( ) ),
            bodies, bodiesToIntegrate )
    { }

    //! Function processing rotational state in the full numericalSolution
    /*!
     * Function that processes the entries of the rotational state in the full numericalSolution,
     * extracts the states for each body, and updates the associated rotational ephemerides.
     * \param numericalSolution Full numerical solution, in global representation (see
     * convertToOutputSolution function in RotationalMotionStateDerivative class.
     */
    void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution )
    {
        resetIntegratedRotationalEphemerides< TimeType, StateScalarType >(
                    this->bodies_, numericalSolution, this->bodiesToIntegrate_, this->startIndexAndSize_ );
    }

private:

};

//! Class used for processing numerically integrated masses of bodies.
template< typename TimeType, typename StateScalarType >
class BodyMassIntegratedStateProcessor: public SingleArcIntegratedStateProcessor< TimeType, StateScalarType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param startIndex Index in the state vector where the translational state starts.
     * \param bodies List of bodies used in simulations.
     * \param bodiesToIntegrate List of bodies for which the mass is numerically
     * integrated. Order in this vector is the same as the order in state vector.
     */
    BodyMassIntegratedStateProcessor(
            const int startIndex,
            const simulation_setup::SystemOfBodies& bodies,
            const std::vector< std::string >& bodiesToIntegrate ):
        SingleArcIntegratedStateProcessor<  TimeType, StateScalarType >(
            body_mass_state, std::make_pair( startIndex, static_cast< int >( bodiesToIntegrate.size( ) ) ), bodies, bodiesToIntegrate )
    { }

    //! Destructor
    ~BodyMassIntegratedStateProcessor( ){ }

    //! Function processing mass state in the full numericalSolution
    /*!
     * Function that processes the entries of the propagated mass in the full numericalSolution, resetting bodies' mass models
     * \param numericalSolution Full numerical solution of state, in global representation (representation is constant
     * for mass).
     */
    void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution )
    {
        resetIntegratedBodyMass( this->bodies_, numericalSolution, this->bodiesToIntegrate_, this->startIndexAndSize_ );
    }

private:

};


template< typename StateScalarType, typename TimeType >
void checkRotationalStatesFeasibility(
    const std::vector< std::string >& bodiesToIntegrate,
    const simulation_setup::SystemOfBodies& bodies,
    const bool setIntegratedResult )
{
    // Check whether each integrated body exists, and whether it has a TabulatedEphemeris
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        std::string bodyToIntegrate = bodiesToIntegrate.at( i );

        if( bodies.count( bodyToIntegrate ) == 0 )
        {
            throw std::runtime_error( "Error when checking rotational dynamics feasibility of body " +
                                      bodyToIntegrate + " no such body found" );
        }
        else
        {
            if( setIntegratedResult )
            {
                if( bodies.at( bodyToIntegrate )->getRotationalEphemeris( ) == nullptr )
                {
                    addEmptyTabulatedRotationalEphemeris< StateScalarType, TimeType >( bodies, bodyToIntegrate );
                }
                else if( !ephemerides::isTabulatedRotationalEphemeris( bodies.at( bodyToIntegrate )->getRotationalEphemeris( )  ) )
                {
                    throw std::runtime_error( "Error when checking rotational dynamics feasibility of body " +
                                              bodyToIntegrate + " rotational ephemeris exists, but is not tabulated." );
                }
            }
        }
    }
}

template< typename StateScalarType, typename TimeType >
void checkTranslationalStatesFeasibility(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const bool setIntegratedResult = false,
        const bool isPartOfMultiArc = false )

{
    // Check feasibility of epheme                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ris origins.
    for( auto bodyIterator : bodies.getMap( )  )
    {
        if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), bodyIterator.first ) ==
            bodiesToIntegrate.end( ) )
        {
            if( bodyIterator.second->getEphemeris( ) != nullptr )
            {
                std::string ephemerisOrigin
                    = bodyIterator.second->getEphemeris( )->getReferenceFrameOrigin( );

                if( std::find( bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), ephemerisOrigin )
                    != bodiesToIntegrate.end( ) )
                {
                    std::cerr << "Warning, found non-integrated body with an integrated body as ephemeris origin " +
                                 bodyIterator.second->getEphemeris( )->getReferenceFrameOrigin( ) + " " +
                                 bodyIterator.first << std::endl;
                }
            }
        }

    }

    // Check whether each integrated body exists, and whether it has a TabulatedEphemeris
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        std::string bodyToIntegrate = bodiesToIntegrate.at( i );
        std::string centralBody = centralBodies.at( i );

        if( bodies.count( bodyToIntegrate ) == 0 )
        {
            throw std::runtime_error( "Error when checking translational dynamics feasibility of body " +
                                      bodyToIntegrate + " no such body found" );
        }
        else
        {
            if( setIntegratedResult )
            {
                if ( bodies.at( bodyToIntegrate )->getEphemeris( ) == nullptr )
                {
                    addEmptyTabulatedEphemeris<StateScalarType, TimeType>( bodies, bodyToIntegrate, centralBody, isPartOfMultiArc );
                }
                else if ( !ephemerides::isTabulatedEphemeris( bodies.at( bodyToIntegrate )->getEphemeris( ) ) && !isPartOfMultiArc )
                {
                    throw std::runtime_error( "Error when checking translational dynamics feasibility of body " +
                                              bodyToIntegrate + " ephemeris exists, but is not tabulated." );
                }
                else if ( ( std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( bodies.at( bodyToIntegrate )->getEphemeris( ) ) == nullptr )
                    && isPartOfMultiArc )
                {
                    throw std::runtime_error( "Error when checking translational dynamics feasibility of body " +
                                              bodyToIntegrate + " ephemeris exists, but is not multi-arc." );
                }
            }
        }

        if( bodiesToIntegrate.at( i ) ==  bodies.getFrameOrigin( ) )
        {
            throw std::runtime_error( "Error when propagating translational dynamics, cannnot propagate " + bodiesToIntegrate.at( i ) +
                                      ", as it is the global origin" );
        }
    }
}



template< typename StateScalarType, typename TimeType >
void checkPropagatedStatesFeasibility(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const bool isPartOfMultiArc )
{
    // Check dynamics type.
    switch( propagatorSettings->getStateType( ) )
    {
    case hybrid:
    {
        // Check input consistency
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType, TimeType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( multiTypePropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, multi-type propagator settings are inconsistent when make state processors" );
        }

        // Iterate over each propagated state type
        for( auto typeIterator :  multiTypePropagatorSettings->propagatorSettingsMap_ )
        {
            // Multi-type in multi-type not allowed (yet)
            if( typeIterator.first != hybrid )
            {
                for( unsigned int i = 0; i < typeIterator.second.size( ); i++ )
                {
                    if(  typeIterator.second.at( i ) == nullptr )
                    {
                        std::string errorMessage = "Error in when processing hybrid propagator settings, propagator entry " +
                                std::to_string( i ) + " is not defined.";
                        throw std::runtime_error( errorMessage );
                    }

                    //  Create state processor
                    checkPropagatedStatesFeasibility( typeIterator.second.at( i ), bodies, isPartOfMultiArc );
                }
            }
            else
            {
                throw std::runtime_error(
                            "Error when when checking dynamics feasibility, cannot handle multi-type propagator inside multi-type propagator" );
            }
        }
        break;
    }
    case translational_state:
    {
        // Check input feasibility
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >
                translationalPropagatorSettings = std::dynamic_pointer_cast
                < TranslationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( translationalPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type for translational dynamics is inconsistent when checking dynamics feasibility" );
        }
        checkTranslationalStatesFeasibility< StateScalarType, TimeType >(
                    translationalPropagatorSettings->bodiesToIntegrate_,
                    translationalPropagatorSettings->centralBodies_, bodies,
                    translationalPropagatorSettings->getOutputSettings( )->getSetIntegratedResult( ),
                    isPartOfMultiArc );
        break;
    }
    case rotational_state:
    {
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType, TimeType > > rotationalPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( rotationalPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type for rotational dynamics is inconsistent when checking dynamics feasibility" );
        }
        checkRotationalStatesFeasibility< StateScalarType, TimeType >(
                    rotationalPropagatorSettings->bodiesToIntegrate_,
                    bodies, rotationalPropagatorSettings->getOutputSettings( )->getSetIntegratedResult( ) );
        break;
    }
    case body_mass_state:
    {
        // Check input feasibility
        std::shared_ptr< MassPropagatorSettings< StateScalarType, TimeType > >
                massPropagatorSettings = std::dynamic_pointer_cast
                < MassPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if( massPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type for mass dynamics is inconsistent when checking dynamics feasibility" );
        }

        break;
    }
    case custom_state:
    {
        // Check input feasibility
        std::shared_ptr< CustomStatePropagatorSettings< StateScalarType > >
                customPropagatorSettings = std::dynamic_pointer_cast
                < CustomStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( customPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type for custom dynamics is inconsistent when checking dynamics feasibility" );
        }

        break;
    }
    default:
        throw std::runtime_error( "Error, integrated state type " +
                                  std::to_string( propagatorSettings->getStateType( ) ) +
                                  " not recognized when checking dynamics feasibility");
    }
}


//! Function to create list objects for processing numerically integrated results.
/*!
 * Function to create list objects for processing numerically integrated results, so that all
 * results are set in the environment models (i.e. resetting tabulated ephemeris for translational
 * dynamics).
 * \param propagatorSettings Settings for the propagation that is used
 * \param bodies List of body objects that represents the environemnt
 * \param frameManager Object for providinf conversion functions between different ephemeris origins.
 * \param startIndex Index of state vector where the state entries handled with propagatorSettings
 *        starts.
 * \return List objects for processing numerically integrated results.
 */
template< typename TimeType, typename StateScalarType >
std::map< IntegratedStateType,
std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
createIntegratedStateProcessors(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager,
        const int startIndex = 0 )
{
    std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > integratedStateProcessors;

    // Check dynamics type.
    switch( propagatorSettings->getStateType( ) )
    {
    case hybrid:
    {
        // Check input consistency
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType, TimeType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if ( multiTypePropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, multi-type propagator settings are inconsistent when make state processors" );
        }

        // Iterate over each propagated state type
        std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                singleTypeIntegratedStateProcessors;
        int currentStartIndex = 0;
        for( auto typeIterator :  multiTypePropagatorSettings->propagatorSettingsMap_ )
        {
            // Multi-type in multi-type not allowed (yet)
            if ( typeIterator.first != hybrid )
            {
                for ( unsigned int i = 0; i < typeIterator.second.size( ); i++ )
                {
                    if( typeIterator.second.at( i ) == nullptr )
                    {
                        std::string errorMessage = "Error in when processing hybrid propagator settings, propagator entry " +
                                std::to_string( i ) + " is not defined.";
                        throw std::runtime_error( errorMessage );
                    }

                    //  Create state processor
                    singleTypeIntegratedStateProcessors = createIntegratedStateProcessors< TimeType, StateScalarType >(
                                typeIterator.second.at( i ), bodies, frameManager, currentStartIndex );

                    if ( singleTypeIntegratedStateProcessors.size( ) > 1 )
                    {
                        throw std::runtime_error( "Error when making hybrid integrated result processors, multiple types found" );
                    }
                    else if ( ( singleTypeIntegratedStateProcessors.size( ) == 0 ) &&  ( typeIterator.first != custom_state ) )
                    {
                        throw std::runtime_error( "Error when making hybrid integrated result processors, no types found" );
                    }
                    else if ( ( singleTypeIntegratedStateProcessors.size( ) == 1 ) )
                    {
                        integratedStateProcessors[ singleTypeIntegratedStateProcessors.begin( )->first ] =
                                singleTypeIntegratedStateProcessors.begin( )->second;
                    }

                    currentStartIndex += typeIterator.second.at( i )->getConventionalStateSize( );
                }
            }
            else
            {
                throw std::runtime_error( "Error when making integrated state processors, cannot handle hybrid propagator inside hybrid propagator" );
            }
        }
        break;
    }

    case translational_state:
    {

        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >
                translationalPropagatorSettings =
                std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if ( translationalPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type is inconsistent in createIntegratedStateProcessors" );
        }
        // Create state processors
        integratedStateProcessors[ translational_state ] =
                std::make_shared< TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > >(
                    startIndex, bodies, translationalPropagatorSettings->bodiesToIntegrate_,
                    translationalPropagatorSettings->centralBodies_, frameManager );
        break;
    }

    case rotational_state:
    {
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType, TimeType > > rotationalPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );

        integratedStateProcessors[ rotational_state ] =
                std::make_shared< RotationalStateIntegratedStateProcessor< TimeType, StateScalarType > >(
                    startIndex, bodies, rotationalPropagatorSettings->bodiesToIntegrate_ );
        break;
    }

    case body_mass_state:
    {
        // Check input feasibility
        std::shared_ptr< MassPropagatorSettings< StateScalarType, TimeType > > massPropagatorSettings =
                std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
        if ( massPropagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error, input type is inconsistent in createIntegratedStateProcessors" );
        }

        // Create mass processors.
        integratedStateProcessors[ body_mass_state ] =
                std::make_shared< BodyMassIntegratedStateProcessor< TimeType, StateScalarType > >(
                    startIndex, bodies, massPropagatorSettings->bodiesWithMassToPropagate_ );
        break;
    }
    case custom_state:
    {
        break;
    }
    default:
        throw std::runtime_error( "Error, could not process integrated state type " +
                                  std::to_string( propagatorSettings->getStateType( ) ) );
    }

    return integratedStateProcessors;
}

//! Function to create list objects for processing numerically integrated multi-arc results.
/*!
 * Function to create list objects for processing numerically integrated multi-arc results, so that all
 * multi-arc results are set in the environment models (i.e. resetting tabulated ephemeris for translational
 * dynamics).
 * \param propagatorSettings Settings for the propagation that is used
 * \param bodies List of body objects that represents the environemnt
 * \param frameManager Object for providinf conversion functions between different ephemeris origins.
 * \param startIndex Index of state vector where the state entries handled with propagatorSettings
 *        starts.
 * \return List objects for processing numerically integrated results.
 */
template< typename TimeType, typename StateScalarType >
std::map< IntegratedStateType,
std::shared_ptr< MultiArcIntegratedStateProcessor< TimeType, StateScalarType > > >
createMultiArcIntegratedStateProcessors(
        const simulation_setup::SystemOfBodies& bodies,
        std::vector< double > arcStartTimes,
        const std::map< IntegratedStateType, std::vector< std::shared_ptr<
        SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > > singleArcIntegratedStatesProcessors )
{
    std::map< IntegratedStateType, std::shared_ptr< MultiArcIntegratedStateProcessor< TimeType, StateScalarType > > > integratedStateProcessors;

    for ( auto stateTypeItr : singleArcIntegratedStatesProcessors )
    {
        // Check input consistency
        if ( arcStartTimes.size( ) != stateTypeItr.second.size( ) )
        {
            throw std::runtime_error( "Error when creating multi-arc integrated state processors, input sizes are inconsistent." +
            std::to_string( arcStartTimes.size( ) ) + "_" + std::to_string( stateTypeItr.second.size( )) );
        }

        // Check dynamics type.
        switch( stateTypeItr.first )
        {
        case hybrid:
        {
            throw std::runtime_error( "Error when creating multi-arc integrated state processors, hybrid type should not be detected." );
            break;
        }
        case translational_state:
        {
            std::vector< std::shared_ptr<
                    TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > > > singleArcTranslationalStateProcessors;
            for ( unsigned int i = 0 ; i < stateTypeItr.second.size( ) ; i++ )
            {
                singleArcTranslationalStateProcessors.push_back(
                            std::dynamic_pointer_cast< TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > >( stateTypeItr.second.at( i ) ) );
            }
            integratedStateProcessors[ translational_state ] =
                    std::make_shared< MultiArcTranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > >(
                        bodies, arcStartTimes, singleArcTranslationalStateProcessors );
            break;
        }
        case rotational_state:
        {
            throw std::runtime_error( "Error, cannot yet reset multi-arc rotational model" );
            break;
        }
        case body_mass_state:
        {
            throw std::runtime_error( "Error, cannot yet reset multi-arc mass model" );
            break;
        }
        case custom_state:
        {
            break;
        }
        default:
            throw std::runtime_error( "Error, could not process integrated state type " + std::to_string( stateTypeItr.first ) );
        }

    }

    return integratedStateProcessors;
}


//! Function resetting dynamical properties of environment from numerical dynamics solution
/*!
 * Function to reset the dynamical properties of the environment from the numerically integrated
 * dynamics solution
 * \param equationsOfMotionNumericalSolution Solution produced by the numerical integration, in the
 * 'conventional form'
 * \sa SingleStateTypeDerivative::convertToOutputSolution
 * \param integratedStateProcessors List of objects (per dynamics type) used to process integrated
 * results into environment
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedStates(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< IntegratedStateType, std::shared_ptr<
        SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > integratedStateProcessors )
{
    for( typename std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor<
         TimeType, StateScalarType > > >::const_iterator updateIterator = integratedStateProcessors.begin( );
         updateIterator != integratedStateProcessors.end( ); updateIterator++ )
    {
        updateIterator->second->processIntegratedStates( equationsOfMotionNumericalSolution );
    }
}

} // namespace propagators

} // namespace tudat

#endif // TUDAT_SETNUMERICALLYINTEGRATEDSTATES_H
