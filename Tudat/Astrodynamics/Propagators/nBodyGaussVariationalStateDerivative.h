/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NGAUSSSTATEDERIVATIVE_H
#define TUDAT_NGAUSSSTATEDERIVATIVE_H

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

template< typename StateScalarType = double, typename TimeType = double >
class NBodyGaussStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
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
    NBodyGaussStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                               const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                               const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, gauss, bodiesToIntegrate )
    {
        originalAccelerationModelsPerBody_ = this->accelerationModelsPerBody_ ;

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_ );
    }

    //! Destructor
    ~NBodyGaussStateDerivative( ){ }


    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, false );

        Eigen::Vector6d bodyCartesianState;
        Eigen::Vector6d currentKeplerianState;

        Eigen::Vector3d currentAccelerationInRswFrame;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentKeplerianState = stateOfSystemToBeIntegrated.block( i * 6, 0, 6, 1 ).template cast< double >( );
            double currentEccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(
                        currentKeplerianState( 1 ), currentKeplerianState( 5 ) );
            double currentTrueAnomaly = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                        currentEccentricAnomaly, currentKeplerianState( 1 ) );
            currentKeplerianState( 5 ) = currentTrueAnomaly;


            bodyCartesianState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        currentKeplerianState,
                        centralBodyGravitationalParameters_.at( i )( ) );

            currentAccelerationInRswFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                        bodyCartesianState ) * stateDerivative.block( i * 6 + 3, 0, 3, 1 ).template cast< double >( );
            stateDerivative.block( i * 6, 0, 6, 1 ) = computeGaussPlanetaryEquationsForKeplerElements(
                        currentKeplerianState, currentAccelerationInRswFrame,
                        centralBodyGravitationalParameters_.at( i )( ) );
        }
    }


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( cartesianSolution.rows( ) );

        // Subtract frame origin and Keplerian states from inertial state.
        Eigen::Vector6d currentCartesianState;
        Eigen::Vector6d currentKeplerianState;

        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianState = cartesianSolution.block( i * 6, 0, 6, 1 ).template cast< double >( );
            currentKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements (
                        currentCartesianState, centralBodyGravitationalParameters_.at( i )( ) );
            double eccentricAnomaly = orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                        currentKeplerianState( 5 ), currentKeplerianState( 1 ) );
            double meanAnomaly = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
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
        Eigen::Vector6d currentKeplerianState;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentKeplerianState = internalSolution.block( i * 6, 0, 6, 1 ).template cast< double >( );
            double currentEccentricAnomaly = orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly(
                        currentKeplerianState( 1 ), currentKeplerianState( 5 ) );
            double currentTrueAnomaly = orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                        currentEccentricAnomaly, currentKeplerianState( 1 ) );
            currentKeplerianState( 5 ) = currentTrueAnomaly;

            currentCartesianLocalSoluton.segment( i * 6, 6 ) = orbital_element_conversions::convertKeplerianToCartesianElements(
                    currentKeplerianState, centralBodyGravitationalParameters_.at( i )( ) );
        }
        //std::cout<<"From Keplerian: "<<std::setprecision( 16 )<<time<<" "<<currentKeplerianState.transpose( )<<std::endl<<std::endl;
    }

    basic_astrodynamics::AccelerationMap getFullAccelerationsMap( )
    {
        return originalAccelerationModelsPerBody_;
    }

private:

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters_;


    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_/
    std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    basic_astrodynamics::AccelerationMap originalAccelerationModelsPerBody_;};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NGAUSSSTATEDERIVATIVE_H
