/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H
#define TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"

namespace tudat
{

namespace propagators
{

Eigen::Vector6d computeGaussPlanetaryEquationsForModifiedEquinoctialElements(
        const Eigen::Vector6d& osculatingModifiedEquinoctialElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter );


template< typename StateScalarType = double, typename TimeType = double >
class NBodyGaussModifiedEquinictialStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
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
     *  \param initialKeplerElements Kepler elements of bodiesToIntegrate, at initial propagation time
     */
    NBodyGaussModifiedEquinictialStateDerivative(
            const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
            const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const std::vector< std::string >& bodiesToIntegrate ,
            const std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& initialKeplerElements ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, gauss_modified_equinoctial, bodiesToIntegrate )
    {
        originalAccelerationModelsPerBody_ = this->accelerationModelsPerBody_ ;

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_ );
        this->createAccelerationModelList( );

        flipSingularities_.resize( bodiesToIntegrate.size( ) );
        for( unsigned int i = 0; i < initialKeplerElements.size( ); i++ )
        {
            if( initialKeplerElements.at( i )( orbital_element_conversions::inclinationIndex ) >  mathematical_constants::PI / 2.0 )
            {
                flipSingularities_[ i ] = true;
            }
            else
            {
                flipSingularities_[ i ] = false;
            }

            if( initialKeplerElements.at( i )( orbital_element_conversions::inclinationIndex ) >
                    mathematical_constants::PI  )
            {
                throw std::runtime_error( "Error when setting N Body Gauss MEE propagator, initial inclination of body is larger than pi." );
            }
        }

    }

    //! Destructor
    ~NBodyGaussModifiedEquinictialStateDerivative( ){ }


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
                        currentCartesianLocalSoluton_.segment( i * 6, 6 ).template cast< double >( ) ) *
                    stateDerivative.block( i * 6 + 3, 0, 3, 1 ).template cast< double >( );

            stateDerivative.block( i * 6, 0, 6, 1 ) = computeGaussPlanetaryEquationsForModifiedEquinoctialElements(
                        stateOfSystemToBeIntegrated.block( i * 6, 0, 6, 1 ).template cast< double >( ), currentAccelerationInRswFrame,
                        centralBodyGravitationalParameters_.at( i )( ) ).template cast< StateScalarType >( );
        }

    }


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( cartesianSolution.rows( ) );

        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentState.segment( i * 6, 6 ) = orbital_element_conversions::convertCartesianToModifiedEquinoctialElements< StateScalarType >(
                        cartesianSolution.block( i * 6, 0, 6, 1 ), static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ), flipSingularities_.at( i ) );
        }

        return currentState;

    }


    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianLocalSoluton.segment( i * 6, 6 ) =
                    orbital_element_conversions::convertModifiedEquinoctialToCartesianElements< StateScalarType >(
                        internalSolution.block( i * 6, 0, 6, 1 ), static_cast< StateScalarType >(
                            centralBodyGravitationalParameters_.at( i )( ) ), flipSingularities_.at( i ) );
        }

        currentCartesianLocalSoluton_ = currentCartesianLocalSoluton;
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

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentCartesianLocalSoluton_;

    std::vector< bool > flipSingularities_;

};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_NGAUSSMODIFIEDEQUINOCTIALSTATEDERIVATIVE_H
