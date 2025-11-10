/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_AERODYNAMIC_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ACCELERATION_H

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicForce.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/aerodynamics/flightConditions.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes primitive types as arguments to
 * perform the calculations. Therefore, these quantities (dynamic pressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param aerodynamicCoefficients Aerodynamic coefficients in right-handed reference frame.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given (assuming coefficients in positive direction).
 */
Eigen::Vector3d computeAerodynamicAcceleration( const double dynamicPressure,
                                                const double referenceArea,
                                                const Eigen::Vector3d& aerodynamicCoefficients,
                                                const double vehicleMass );

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and coefficients are retrieved.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given (assuming coefficients in positive direction).
 */
Eigen::Vector3d computeAerodynamicAcceleration( const double dynamicPressure,
                                                AerodynamicCoefficientInterfacePointer coefficientInterface,
                                                const double vehicleMass );

//! Class for calculation of aerodynamic accelerations.
/*!
 * Class for calculation of aerodynamic accelerations.
 * \sa AccelerationModel.
 */
class AerodynamicAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    AerodynamicAcceleration( const std::shared_ptr< AtmosphericFlightConditions > flightConditions,
                             const std::function< double( ) > currentMass ):
        flightConditions_( flightConditions ), currentMass_( currentMass ), dragComponentScaling_( 1.0 ), liftComponentScaling_( 1.0 ),
        sideComponentScaling_( 1.0 )
    {
        coefficientInterface_ = flightConditions_->getAerodynamicCoefficientInterface( );
        aerodynamicCoefficientFrame_ = coefficientInterface_->getForceCoefficientsFrame( );
        aerodynamicCompleteCoefficientFrame_ = getCompleteFrameForCoefficients( aerodynamicCoefficientFrame_ );
        coefficientMultiplier_ = areCoefficientsInNegativeDirection( aerodynamicCoefficientFrame_ ) == true ? -1.0 : 1.0;
    }

    //! Destructor
    virtual ~AerodynamicAcceleration( ) {}

    //! Update member variables used by the aerodynamic acceleration model.
    /*!
     * Updates member variables used by the aerodynamic acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentTime_ = currentTime;
            currentForceCoefficients_ = coefficientInterface_->getCurrentForceCoefficients( );
            currentForceCoefficients_ = coefficientMultiplier_ *
                    ( flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                              aerodynamicCompleteCoefficientFrame_, reference_frames::inertial_frame ) *
                      currentForceCoefficients_ );

            currentUnscaledAcceleration_ = computeAerodynamicAcceleration( flightConditions_->getCurrentDynamicPressure( ),
                                                                           coefficientInterface_->getReferenceArea( ),
                                                                           currentForceCoefficients_,
                                                                           currentMass_( ) );
            scaleAerodynamicAcceleration( );
        }
    }

    void scaleAerodynamicAcceleration( )
    {
        currentAcceleration_ = currentUnscaledAcceleration_;

        currentUnscaledAccelerationInAerodynamicFrame_ =
                flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                        reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) *
                currentAcceleration_;

        if( isScalingModelSet_ )
        {
            // if scalings are defined via functions (--> arc-wise!) update the scaling factor to current time
            if( (dragComponentScalingFunction_ != nullptr) )
            {
                dragComponentScaling_ = dragComponentScalingFunction_( currentTime_ );
            }
            if( (sideComponentScalingFunction_ != nullptr) )
            {
                sideComponentScaling_ = sideComponentScalingFunction_( currentTime_ );
            }
            if( (dragComponentScalingFunction_ != nullptr) )
            {
                liftComponentScaling_ = liftComponentScalingFunction_( currentTime_ );
            }

            currentAccelerationInAerodynamicFrame_( 0 ) = currentUnscaledAccelerationInAerodynamicFrame_( 0 ) * dragComponentScaling_;
            currentAccelerationInAerodynamicFrame_( 1 ) = currentUnscaledAccelerationInAerodynamicFrame_( 1 ) * sideComponentScaling_;
            currentAccelerationInAerodynamicFrame_( 2 ) = currentUnscaledAccelerationInAerodynamicFrame_( 2 ) * liftComponentScaling_;

            currentAcceleration_ = flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                                           reference_frames::aerodynamic_frame, reference_frames::inertial_frame ) *
                    currentAccelerationInAerodynamicFrame_;
        }
    }

    std::shared_ptr< AtmosphericFlightConditions > getFlightConditions( ) const
    {
        return flightConditions_;
    }

    std::shared_ptr< AerodynamicCoefficientInterface > getCoefficientInterface( ) const
    {
        return coefficientInterface_;
    }

    Eigen::Vector3d getCurrentForceCoefficientsInAerodynamicFrame( ) const
    {
        return flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                       reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) *
                ( -currentForceCoefficients_ );
    }

    double getCurrentMass( ) const
    {
        return currentMass_( );
    }

    void enableScaling( )
    {
        isScalingModelSet_ = true;
    }

    void setDragComponentScaling( double dragComponentScaling )
    {
        enableScaling( );
        dragComponentScaling_ = dragComponentScaling;
    }

    void setDragComponentScalingFunction( std::function<double(double)> dragComponentScalingFunction )
    {
        enableScaling( );
        dragComponentScalingFunction_ = dragComponentScalingFunction;
    }

    void setSideComponentScaling( double sideComponentScaling )
    {
        enableScaling( );
        sideComponentScaling_ = sideComponentScaling;
    }

    void setSideComponentScalingFunction( std::function<double(double)> sideComponentScalingFunction )
    {
        enableScaling( );
        sideComponentScalingFunction_ = sideComponentScalingFunction;
    }

    void setLiftComponentScaling( double liftComponentScaling )
    {
        enableScaling( );
        liftComponentScaling_ = liftComponentScaling;
    }

    void setLiftComponentScalingFunction( std::function<double(double)> liftComponentScalingFunction )
    {
        enableScaling( );
        liftComponentScalingFunction_ = liftComponentScalingFunction;
    }

    void setComponentScaling( const double scalingValue, const int index )
    {
        switch( index )
        {
            case 0:
                setDragComponentScaling( scalingValue );
                break;
            case 1:
                setSideComponentScaling( scalingValue );
                break;
            case 2:
                setLiftComponentScaling( scalingValue );
                break;
            default:
                throw std::runtime_error( "Error when setting aerodynamic component scaling factor, index not supported" +
                                          std::to_string( index ) );
        }
    }

    // setter interface for arc-wise parameter
    void setComponentScalingFunction( std::function<double(double)> scalingFunction, const int index )
    {
        switch( index )
        {
            case 0:
                setDragComponentScalingFunction( scalingFunction );
            break;
            case 1:
                setSideComponentScalingFunction( scalingFunction );
            break;
            case 2:
                setLiftComponentScalingFunction( scalingFunction );
            break;
            default:
                throw std::runtime_error( "Error when setting aerodynamic component scaling function, index not supported" +
                                          std::to_string( index ) );
        }
    }


    double getDragComponentScaling( )
    {
        return dragComponentScaling_;
    }

    std::function<double(double)> getDragComponentScalingFunction( )
    {
        return dragComponentScalingFunction_;
    }

    double getSideComponentScaling( )
    {
        return sideComponentScaling_;
    }

    std::function<double(double)> getSideComponentScalingFunction( )
    {
        return sideComponentScalingFunction_;
    }

    double getLiftComponentScaling( )
    {
        return liftComponentScaling_;
    }

    std::function<double(double)> getLiftComponentScalingFunction( )
    {
        return liftComponentScalingFunction_;
    }

    double getComponentScaling( const int index )
    {
        switch( index )
        {
            case 0:
                return getDragComponentScaling( );
                break;
            case 1:
                return getSideComponentScaling( );
                break;
            case 2:
                return getLiftComponentScaling( );
                break;
            default:
                throw std::runtime_error( "Error when retrieving aerodynamic component scaling factor, index not supported" +
                                          std::to_string( index ) );
        }
    }

    std::function<double(double)> getComponentScalingFunction( const int index )
    {
        switch( index )
        {
            case 0:
                return getDragComponentScalingFunction( );
            break;
            case 1:
                return getSideComponentScalingFunction( );
            break;
            case 2:
                return getLiftComponentScalingFunction( );
            break;
            default:
                throw std::runtime_error( "Error when retrieving aerodynamic component scaling factor, index not supported" +
                                          std::to_string( index ) );
        }
    }

    Eigen::Vector3d getCurrentUnscaledAcceleration( )
    {
        return currentUnscaledAcceleration_;
    }

    Eigen::Vector3d getCurrentUnscaledAccelerationInAerodynamicFrame( )
    {
        return currentUnscaledAccelerationInAerodynamicFrame_;
    }

private:
    std::shared_ptr< AtmosphericFlightConditions > flightConditions_;

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_;

    std::function< double( ) > currentMass_;

    Eigen::Vector3d currentForceCoefficients_;

    double coefficientMultiplier_;

    AerodynamicCoefficientFrames aerodynamicCoefficientFrame_;

    reference_frames::AerodynamicsReferenceFrames aerodynamicCompleteCoefficientFrame_;

    // new acceleration scaling
    Eigen::Vector3d currentUnscaledAcceleration_;

    Eigen::Vector3d currentUnscaledAccelerationInAerodynamicFrame_;

    Eigen::Vector3d currentAccelerationInAerodynamicFrame_;

    bool isScalingModelSet_;

    double dragComponentScaling_;

    std::function<double(double)> dragComponentScalingFunction_;

    double liftComponentScaling_;

    std::function<double(double)> liftComponentScalingFunction_;

    double sideComponentScaling_;

    std::function<double(double)> sideComponentScalingFunction_;

};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_AERODYNAMIC_ACCELERATION_H
