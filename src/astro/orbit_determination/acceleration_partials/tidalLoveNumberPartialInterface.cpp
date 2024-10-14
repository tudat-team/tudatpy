/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

#include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"

namespace tudat
{

namespace orbit_determination
{

void getMaximumUsedDegreeAndOrder(
        const int maximumDegree, const int maximumOrder, const int evaluationDegree,
        int& maximumUsedDegree, int& maximumUsedOrder )
{
    if( maximumDegree >= evaluationDegree )
    {
        maximumUsedDegree = evaluationDegree;
    }
    else
    {
        maximumUsedDegree = 0;
    }

    if( maximumOrder >= evaluationDegree )
    {
        maximumUsedOrder = evaluationDegree;
    }
    else
    {
        maximumUsedOrder = maximumOrder;
    }
}

//! Function to obtain the indices of given list of body names in deformingBodies_ member vector
std::vector< int > TidalLoveNumberPartialInterface::getSelectedDeformingBodyIds(
        const std::vector< std::string >& selectedBodyNames )
{
    std::vector< int > selectedBodyIds;

    // Iterate over all input bodies
    for( unsigned int i = 0; i < selectedBodyNames.size( ); i++ )
    {
        // Check if current body is present in deformingBodies_, and return index (or throw exception if not found)
        std::vector< std::string >::iterator findIterator =
                std::find( deformingBodies_.begin( ), deformingBodies_.end( ), selectedBodyNames.at( i ) );
        if( findIterator != deformingBodies_.end( ) )
        {
            selectedBodyIds.push_back( std::distance( deformingBodies_.begin( ), findIterator ) );
        }
//        else
//        {
//            throw std::runtime_error( "Error when looking for deforming body " + selectedBodyNames.at( i ) +
//                                      " in TidalLoveNumberPartialInterface, body not found." );
//        }
    }
    return selectedBodyIds;
}

//! Function to calculate the partial of spherical harmonic acceleration wrt complex tidal love numbers.
std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > TidalLoveNumberPartialInterface::
calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    deformedBodyGravitationalParameter_ = deformedBodyGravitationalParameterFunction_( );

    // Set and calculate states and rotation needed for calculation of partial.
    updateCurrentTidalBodyStates( deformingBodyIndices );

    // Calculate partial of the coefficients wrt love number.
    std::vector< Eigen::Vector2d > realCoefficientPartials =
            calculateCoefficientPartialWrtRealTidalLoveNumber(
                degree, orders, deformingBodyIndices, maximumDegree, maximumOrder );

    return calculateSphericalHarmonicCoefficientPartialMatrix(
                realCoefficientPartials, complexLoveNumberScaler_ );
}

//! Function to calculate the partial of spherical harmonic acceleration wrt real tidal love numbers.
std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > TidalLoveNumberPartialInterface::
calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    deformedBodyGravitationalParameter_ = deformedBodyGravitationalParameterFunction_( );

    // Set and calculate states and rotation needed for calculation of partial.
    updateCurrentTidalBodyStates( deformingBodyIndices );

    // Calculate partial of the coefficients wrt love number.
    std::vector< Eigen::Vector2d > realCoefficientPartials = calculateCoefficientPartialWrtRealTidalLoveNumber(
                degree, orders, deformingBodyIndices, maximumDegree, maximumOrder );

    return calculateSphericalHarmonicCoefficientPartialMatrix(
                realCoefficientPartials, realLoveNumberScaler_ );
}

std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > TidalLoveNumberPartialInterface::
calculateSphericalHarmonicCoefficientsPartialWrtModeCoupledTidalLoveNumbers(
    const std::vector< std::pair< int, int > > parameterDegreeAndOrderIndices,
    const std::map< int, std::vector< int > >& ordersPerDegree,
    const std::vector< int >& deformingBodyIndices,
    const int maximumDegree,
    const int maximumOrder )
{
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > partials;
    partials.reserve( parameterDegreeAndOrderIndices.size( ) );
    deformedBodyGravitationalParameter_ = deformedBodyGravitationalParameterFunction_( );

    // Set and calculate states and rotation needed for calculation of partial.
    updateCurrentTidalBodyStates( deformingBodyIndices );

    // Compute list of partials for all orders at requested degree
    std::map< int, std::vector<Eigen::Vector2d> > partialsPerDegree;
    for( auto it : ordersPerDegree )
    {
        partialsPerDegree[ it.first ] = calculateCoefficientPartialWrtRealTidalLoveNumber(
            it.first, it.second, deformingBodyIndices, maximumDegree, maximumOrder );
    }

    // Add all partials into single list (ordered as per parameterDegreeAndOrderIndices)
    for( unsigned int i = 0; i < parameterDegreeAndOrderIndices.size( ); i++ )
    {
        partials.push_back( partialsPerDegree.at( parameterDegreeAndOrderIndices.at( i ).first ).at( parameterDegreeAndOrderIndices.at( i ).second ) );
    }
    return partials;
}

//! Function to set a dependency of this partial object w.r.t. a given double parameter.
std::pair< int, std::pair< int, int > > TidalLoveNumberPartialInterface::setParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        const int maximumDegree,
        const int maximumOrder )
{
    return std::make_pair( 0, std::make_pair( 0, 0 ) );
}

//! Function to set a dependency of this partial object w.r.t. a given vector parameter.
std::pair< int, std::pair< int, int > > TidalLoveNumberPartialInterface::setParameterPartialFunction(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        const int maximumDegree,
        const int maximumOrder )
{
    using namespace estimatable_parameters;

    // Define default values for return values
    int maximumUsedDegree = -1;
    int maximumUsedOrder = -1;
    int numberOfRows = 0;

    // Check if parameter is tidal property of deformed body
    if( parameter->getParameterName( ).second.first == deformedBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case full_degree_tidal_love_number:
        {
            // Cast parameter object to required type.
            std::shared_ptr< FullDegreeTidalLoveNumber > coefficientsParameter =
                    std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameter );
            if( coefficientsParameter == nullptr )
            {
                throw std::runtime_error(
                            "Error when setting partil function of full_degree_tidal_love_number in TidalLoveNumberPartialInterface, input is inconsistent" );
            }

            // Retrieve the maximum degree and order that are to be used for this model
            getMaximumUsedDegreeAndOrder(
                        maximumDegree, maximumOrder, coefficientsParameter->getDegree( ), maximumUsedDegree,  maximumUsedOrder );

            // Set parameter partial function if relevant terms are found
            if( maximumUsedDegree > 0 )
            {
                std::vector< int > selectedDeformingBodies = getSelectedDeformingBodyIds(
                            coefficientsParameter->getDeformingBodies( ) );

                // Check if deforming bodies correspond to bodies in model.
                if( selectedDeformingBodies.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) &&
                         deformingBodies_.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) )
                {
                    // Add partial function if it is not yet set.
                    if( parameterVectorPartialFunctions_.count(
                                std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ) == 0 )
                    {
                        if( coefficientsParameter->useComplexComponents( ) )
                        {
                            // Calculate partials for complex love number
                              parameterVectorPartialFunctions_[ std::make_pair(
                                        parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    std::bind( &TidalLoveNumberPartialInterface::
                                                 calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumber, this,
                                                 coefficientsParameter->getDegree( ), selectedDeformingBodies,
                                                 maximumUsedDegree, maximumUsedOrder );
                        }
                        else
                        {
                            // Calculate partial for real love number
                            parameterVectorPartialFunctions_[ std::make_pair(
                                        parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    std::bind( &TidalLoveNumberPartialInterface::
                                                 calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumber, this,
                                                 coefficientsParameter->getDegree( ), selectedDeformingBodies,
                                                 maximumUsedDegree, maximumUsedOrder );
                        }
                    }

                    numberOfRows = coefficientsParameter->getParameterSize( );
                }
            }

            break;
        }
        case single_degree_variable_tidal_love_number:
        {
            // Cast parameter object to required type.
            std::shared_ptr< SingleDegreeVariableTidalLoveNumber > coefficientsParameter =
                    std::dynamic_pointer_cast< SingleDegreeVariableTidalLoveNumber >( parameter );
            if( coefficientsParameter == nullptr )
            {
                throw std::runtime_error(
                            "Error when setting partil function of single_degree_variable_tidal_love_number in TidalLoveNumberPartialInterface, input is inconsistent" );
            }

            // Retrieve the maximum degree and order that are to be used for this model
            getMaximumUsedDegreeAndOrder(
                        maximumDegree, maximumOrder, coefficientsParameter->getDegree( ), maximumUsedDegree, maximumUsedOrder );

            // Set parameter partial function if relevant terms are found
            if( maximumUsedDegree > 0 )
            {
                // Check if deforming bodies correspond to bodies in model.
                std::vector< int > selectedDeformingBodies = getSelectedDeformingBodyIds(
                            coefficientsParameter->getDeformingBodies( ) );
                if( selectedDeformingBodies.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) &&
                        coefficientsParameter->getDeformingBodies( ).size( ) == deformingBodies_.size( ) )
                {
                    // Add partial function if it is not yet set.
                    if( parameterVectorPartialFunctions_.count(
                                std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ) == 0 )
                    {
                        if( coefficientsParameter->useComplexComponents( ) )
                        {
                            // Calculate partials for complex love number
                            parameterVectorPartialFunctions_[ std::make_pair(
                                        parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    std::bind( &TidalLoveNumberPartialInterface::
                                                 calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers, this,
                                                 coefficientsParameter->getDegree( ), coefficientsParameter->getOrders( ),
                                                 selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
                        }
                        else
                        {
                            // Calculate partial for real love number
                            parameterVectorPartialFunctions_[ std::make_pair(
                                        parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                    std::bind( &TidalLoveNumberPartialInterface::
                                                 calculateSphericalHarmonicCoefficientsPartialWrtRealTidalLoveNumbers, this,
                                                 coefficientsParameter->getDegree( ), coefficientsParameter->getOrders( ),
                                                 selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
                        }
                    }

                    numberOfRows = coefficientsParameter->getParameterSize( );
                }
            }

            break;
        }
        case mode_coupled_tidal_love_numbers:
        {
            // Cast parameter object to required type.
            std::shared_ptr< ModeCoupledTidalLoveNumber > coefficientsParameter =
                std::dynamic_pointer_cast< ModeCoupledTidalLoveNumber >( parameter );
            if( coefficientsParameter == nullptr )
            {
                throw std::runtime_error(
                    "Error when setting partil function of mode_coupled_tidal_love_numbers in TidalLoveNumberPartialInterface, input is inconsistent" );
            }

            // Retrieve the maximum degree and order that are to be used for this model
            getMaximumUsedDegreeAndOrder(
                maximumDegree, maximumOrder, coefficientsParameter->getMaximumForcingDegree( ), maximumUsedDegree, maximumUsedOrder );

            // Set parameter partial function if relevant terms are found
            if( maximumUsedDegree > 0 )
            {
                // Check if deforming bodies correspond to bodies in model.
                std::vector< int > selectedDeformingBodies = getSelectedDeformingBodyIds(
                    coefficientsParameter->getDeformingBodies( ) );
                if( selectedDeformingBodies.size( ) == coefficientsParameter->getDeformingBodies( ).size( ) &&
                    coefficientsParameter->getDeformingBodies( ).size( ) == deformingBodies_.size( ) )
                {
                    // Add partial function if it is not yet set.
                    if( parameterVectorPartialFunctions_.count(
                        std::make_pair( parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ) == 0 )
                    {
//                        if( coefficientsParameter->useComplexComponents( ) )
//                        {
//                            // Calculate partials for complex love number
//                            parameterVectorPartialFunctions_[ std::make_pair(
//                                parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
//                                std::bind( &TidalLoveNumberPartialInterface::
//                                           calculateSphericalHarmonicCoefficientsPartialWrtComplexTidalLoveNumbers, this,
//                                           coefficientsParameter->getDegree( ), coefficientsParameter->getOrders( ),
//                                           selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
//                        }
//                        else
//                        {
                            // Calculate partial for real love number
                            parameterVectorPartialFunctions_[ std::make_pair(
                                parameter, std::make_pair( maximumUsedDegree, maximumUsedOrder ) ) ] =
                                std::bind( &TidalLoveNumberPartialInterface::
                                           calculateSphericalHarmonicCoefficientsPartialWrtModeCoupledTidalLoveNumbers, this,
                                           coefficientsParameter->getParameterForcingDegreeAndOrderIndices( ),
                                           coefficientsParameter->getForcingOrdersPerDegree( ),
                                           selectedDeformingBodies, maximumUsedDegree, maximumUsedOrder );
//                        }
                    }



                    numberOfRows = coefficientsParameter->getParameterSize( );
                }
            }

            break;
        }
        default:
            break;
        }
    }
    return std::make_pair( numberOfRows, std::make_pair( maximumUsedDegree, maximumUsedOrder ) );

}


//! Function to pre-calculate all states of bodies involved.
void TidalLoveNumberPartialInterface::updateCurrentTidalBodyStates( const std::vector< int >& deformingBodiesToUpdate )
{
    positionOfTidallyDeformedBody_ = deformedBodyPositionFunction_( );

    for( unsigned int i = 0; i < deformingBodiesToUpdate.size( ); i++ )
    {
        positionsOfDeformingBodies_[ deformingBodiesToUpdate[ i ] ] =
                deformingBodyStateFunctions_[ deformingBodiesToUpdate[ i ] ]( );
    }
}

//! Function to set current state of single body and derived quantities
void TidalLoveNumberPartialInterface::setCurrentTidalBodyStates( const int degree, const int order, const int body )
{
    {
        // Calculate properties of currently considered body
        massRatio_ = deformingBodyGravitationalParameters_[ body ]( ) / deformedBodyGravitationalParameter_;
        relativeDeformingBodyPosition_ = rotationToTidallyDeformedBody_ * (
                    positionsOfDeformingBodies_[ body ] - positionOfTidallyDeformedBody_ );
        relativeDeformingBodySphericalPosition_  =coordinate_conversions::
                convertCartesianToSpherical( relativeDeformingBodyPosition_ );
        radiusRatio_ = deformedBodyReferenceRadius_ / relativeDeformingBodySphericalPosition_.x( );
        iLongitude_ = std::complex< double >( 0.0, relativeDeformingBodySphericalPosition_.z( ) );
        sineOfLatitude_ = std::sin( mathematical_constants::PI / 2.0 - relativeDeformingBodySphericalPosition_.y( ) );
    }
}

//! Function to compute Love number partials from pre-computed partials and provided scaling values
std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >  TidalLoveNumberPartialInterface::
calculateSphericalHarmonicCoefficientPartialMatrix(
        const std::vector< Eigen::Vector2d >& coefficientPartialsPerOrder,
        const std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, Eigen::Matrix< double, 2, Eigen::Dynamic > >&
        coefficientPartialScalers )
{
    std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > coefficientsPartial;

    // Check input consistency
    if( coefficientPartialScalers.first.rows( ) !=  coefficientPartialScalers.second.rows( ) )
    {
        throw std::runtime_error(
                    "Error when making spherical harmonic coefficient partial from coefficient partial in TidalLoveNumberPartialInterface, scale matrices not of same size." );
    }

    // COmpute scaled partials
    for( unsigned int m = 0; m < coefficientPartialsPerOrder.size( ); m++ )
    {
        coefficientsPartial.push_back(
                    ( coefficientPartialScalers.first * coefficientPartialsPerOrder[ m ].x( ) +
                      coefficientPartialScalers.second * coefficientPartialsPerOrder[ m ].y( ) ) );
    }

    return coefficientsPartial;
}

//! Function to calculate the partial of spherical harmonic coefficients w.r.t. real part of tidal love numbers.
std::vector< Eigen::Vector2d > TidalLoveNumberPartialInterface::calculateCoefficientPartialWrtRealTidalLoveNumber(
        const int degree,
        const std::vector< int >& orders,
        const std::vector< int >& deformingBodyIndices,
        const int maximumDegree,
        const int maximumOrder )
{
    // Set default values of return variable
    std::vector< Eigen::Vector2d > realCoefficientPartials;
    realCoefficientPartials.resize( orders.size( ) );
    for( unsigned int i = 0; i < orders.size( ); i ++ )
    {
        realCoefficientPartials[ i ] = Eigen::Vector2d::Zero( );
    }

    // Check if any dependency exists
    if( degree <= maximumDegree )
    {
        std::complex< double > unitLoveNumberCoefficientVariations;

        // Iterate over all bodies causing deformation and calculate and add associated corrections
        for( unsigned int i = 0; i < deformingBodyIndices.size( ); i++ )
        {
            for( unsigned int m = 0; m < orders.size( ); m++ )
            {
                if( orders.at( m ) <= maximumOrder )
                {
                    // Update member variables to current body, degree and order
                    setCurrentTidalBodyStates( degree, orders.at( m ), deformingBodyIndices.at( i ) );

                    // Calculate and add coefficients.
                    unitLoveNumberCoefficientVariations =
                            gravitation::calculateSolidBodyTideSingleCoefficientSetCorrectionFromAmplitude(
                                std::complex< double >( 1.0, 0.0 ), massRatio_,
                                basic_mathematics::raiseToIntegerPower( radiusRatio_, degree + 1 ) ,
                                basic_mathematics::computeLegendrePolynomialExplicit(
                                    degree, orders.at( m ), sineOfLatitude_ ),
                                static_cast< double >( orders.at( m ) ) * iLongitude_, degree, orders.at( m ) );
                    realCoefficientPartials[ m ].x( ) += unitLoveNumberCoefficientVariations.real( );
                    realCoefficientPartials[ m ].y( ) -= unitLoveNumberCoefficientVariations.imag( );
                }
            }
        }
    }
    return realCoefficientPartials;
}

}

}
