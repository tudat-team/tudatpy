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
                             flightConditions_( flightConditions ),
                             currentMass_( currentMass )
    {
        coefficientInterface_ = flightConditions_->getAerodynamicCoefficientInterface( );
        aerodynamicCoefficientFrame_ = coefficientInterface_->getForceCoefficientsFrame( );
        aerodynamicCompleteCoefficientFrame_ = getCompleteFrameForCoefficients( aerodynamicCoefficientFrame_ );
        coefficientMultiplier_ = areCoefficientsInNegativeDirection( aerodynamicCoefficientFrame_) == true ? -1.0 : 1.0;
    }

    //! Destructor
    virtual ~AerodynamicAcceleration( ) { }

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
            currentForceCoefficients_ = coefficientMultiplier_ *  ( flightConditions_->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                aerodynamicCompleteCoefficientFrame_, reference_frames::inertial_frame ) * currentForceCoefficients_ );

            currentAcceleration_ = computeAerodynamicAcceleration( flightConditions_->getCurrentDynamicPressure( ),
                                                                   coefficientInterface_->getReferenceArea( ),
                                                                   currentForceCoefficients_,
                                                                   currentMass_( ) );
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
                reference_frames::inertial_frame, reference_frames::aerodynamic_frame ) * ( -currentForceCoefficients_ );
    }

    double getCurrentMass( ) const
    {
        return currentMass_( );
    }
private:

    std::shared_ptr< AtmosphericFlightConditions > flightConditions_;

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface_;

    std::function< double( ) > currentMass_;

    Eigen::Vector3d currentForceCoefficients_;

    double coefficientMultiplier_;

    AerodynamicCoefficientFrames aerodynamicCoefficientFrame_;

    reference_frames::AerodynamicsReferenceFrames aerodynamicCompleteCoefficientFrame_;

};

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_AERODYNAMIC_ACCELERATION_H
