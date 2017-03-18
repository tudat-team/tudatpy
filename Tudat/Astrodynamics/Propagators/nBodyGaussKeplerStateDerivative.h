/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H
#define TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"

namespace tudat
{

namespace propagators
{

Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double semiLatusRectum,
        const double distance,
        const double meanMotion,
        const double orbitalAngularMomentum );

Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter );

Eigen::Vector6d computeGaussPlanetaryEquationsForKeplerElements(
        const Eigen::Vector6d& currentOsculatingKeplerElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter );

template< typename StateScalarType = double, typename TimeType = double >
class NBodyGaussKeplerStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyGaussKeplerStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                               const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                               const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, gauss_keplerian, bodiesToIntegrate )
    {
        currentTrueAnomalies_.resize( bodiesToIntegrate.size( ) );
        originalAccelerationModelsPerBody_ = this->accelerationModelsPerBody_ ;

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_ );
        this->createAccelerationModelList( );

    }

    //! Destructor
    ~NBodyGaussKeplerStateDerivative( ){ }


    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, false );


        Eigen::Vector3d currentAccelerationInRswFrame;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentAccelerationInRswFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                        currentCartesianLocalSoluton_.segment( i * 6, 6 ) ) *
                    stateDerivative.block( i * 6 + 3, 0, 3, 1 ).template cast< double >( );

                stateDerivative.block( i * 6, 0, 6, 1 ) = computeGaussPlanetaryEquationsForKeplerElements(
                            ( Eigen::Vector6d( )<<stateOfSystemToBeIntegrated.block( i * 6, 0, 5, 1 ).template cast< double >( ),
                              currentTrueAnomalies_.at( i ) ).finished( ), currentAccelerationInRswFrame,
                            centralBodyGravitationalParameters_.at( i )( ) ).template cast< StateScalarType >( );
        }

    }


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( cartesianSolution.rows( ) );

        // Subtract frame origin and Keplerian states from inertial state.
        Eigen::Matrix< StateScalarType, 6, 1 > currentCartesianState;
        Eigen::Matrix< StateScalarType, 6, 1 > currentKeplerianState;

        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianState = cartesianSolution.block( i * 6, 0, 6, 1 );
            currentKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements (
                        currentCartesianState, static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ) );
            StateScalarType eccentricAnomaly = orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                        currentKeplerianState( 5 ), currentKeplerianState( 1 ) );
            StateScalarType meanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomaly, currentKeplerianState( 1 ) );
            currentKeplerianState( 5 ) = meanAnomaly;
            currentState.segment( i * 6, 6 ) = currentKeplerianState;
        }

        return currentState;

    }


    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        // Add Keplerian state to perturbation from Encke algorithm to get Cartesian state in local frames.
        Eigen::Matrix< StateScalarType, 6, 1 > currentKeplerianState;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentKeplerianState = internalSolution.block( i * 6, 0, 6, 1 );
            StateScalarType currentEccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(
                        currentKeplerianState( 1 ), currentKeplerianState( 5 ) );
            StateScalarType currentTrueAnomaly = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                        currentEccentricAnomaly, currentKeplerianState( 1 ) );
            currentKeplerianState( 5 ) = currentTrueAnomaly;

            currentTrueAnomalies_[ i ] = currentTrueAnomaly;
            currentCartesianLocalSoluton.segment( i * 6, 6 ) = orbital_element_conversions::convertKeplerianToCartesianElements(
                        currentKeplerianState, static_cast< StateScalarType >( centralBodyGravitationalParameters_.at( i )( ) ) );

        }

        currentCartesianLocalSoluton_ = currentCartesianLocalSoluton.template cast< double >( );
    }

    basic_astrodynamics::AccelerationMap getFullAccelerationsMap( )
    {
        return originalAccelerationModelsPerBody_;
    }

private:

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters_;

    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_
    std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    basic_astrodynamics::AccelerationMap originalAccelerationModelsPerBody_;

    Eigen::VectorXd currentCartesianLocalSoluton_;

    std::vector< double > currentTrueAnomalies_;

};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_NGAUSSKEPLERSTATEDERIVATIVE_H
