/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110617    D. Dirkx          File created.
 *      120324    K. Kumar          Minor Doxygen comment corrections, added astrodynamics
 *                                  namespace layer; added missing Eigen include-statement.
 *      121020    D. Dirkx          Update to new acceleration model architecture.
 *      130120    K. Kumar          Added shared pointer to AerodynamicAcceleration object.
 *      140129    D. Dirkx          Changed Doxygen descriptions
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_AERODYNAMIC_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ACCELERATION_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

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
Eigen::Vector3d computeAerodynamicAcceleration(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface,
        const double vehicleMass );

//! Class for calculation of aerodynamic accelerations.
/*!
 * Class for calculation of aerodynamic accelerations.
 * \sa AccelerationModel.
 */
class AerodynamicAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double ( ) > DoubleReturningFunction;

    //! Typedef for coefficient-returning function.
    typedef boost::function< Eigen::Vector3d( ) > CoefficientReturningFunction;

public:

    //! Acceleration model constructor, taking constant values of mass and reference area.
    /*!
     * Acceleration model constructor, taking constant values of mass and reference area.
     * \param coefficientFunction Function which retrieves current values of aerodynamic
     *          coefficients.
     * \param densityFunction Function which retrieves current value of the density.
     * \param airSpeedFunction Function which retrieves current value of the airspeed.
     * \param constantMass Value of vehicle mass that is used for all calls of this class.
     * \param constantReferenceArea Value of aerodynamic coefficient reference area that is used
     *          for all calls of this class.
     * \param areCoefficientsInNegativeDirection Boolean that determines whether to invert
     *          direction of aerodynamic coefficients. This is typically done for lift, drag and
     *          side force coefficients that point in negative direction in the local frame
     *          (default true).
     */
    AerodynamicAcceleration( const CoefficientReturningFunction coefficientFunction,
                             const DoubleReturningFunction densityFunction,
                             const DoubleReturningFunction airSpeedFunction,
                             const double constantMass,
                             const double constantReferenceArea,
                             const bool areCoefficientsInNegativeDirection = true ):
        coefficientFunction_( coefficientFunction ),
        densityFunction_( densityFunction ),
        airSpeedFunction_( airSpeedFunction ),
        massFunction_( boost::lambda::constant( constantMass ) ),
        referenceAreaFunction_( boost::lambda::constant( constantReferenceArea ) )
    {
        coefficientMultiplier_ = areCoefficientsInNegativeDirection == true ? -1.0 : 1.0;
    }

    //! Acceleration model constructor.
    /*!
     * Acceleration model constructor, taking function pointers for all member variables.
     * \param coefficientFunction Function which retrieves current values of aerodynamic
     *          coefficients.
     * \param densityFunction Function which retrieves current value of the density.
     * \param airSpeedFunction Function which retrieves current value of the airspeed.
     * \param massFunction Function which retrieves current value of the vehicle mass.
     * \param referenceAreaFunction Function which retrieves current value of the aerodynamic
     *          coefficient reference area.
     * \param areCoefficientsInNegativeDirection Boolean that determines whether to invert
     *          direction of aerodynamic coefficients. This is typically done for lift, drag and
     *          side force coefficients that point in negative direction in the local frame
     *          (default true).
     */
    AerodynamicAcceleration( const CoefficientReturningFunction coefficientFunction,
                             const DoubleReturningFunction densityFunction,
                             const DoubleReturningFunction airSpeedFunction,
                             const DoubleReturningFunction massFunction,
                             const DoubleReturningFunction referenceAreaFunction,
                             const bool areCoefficientsInNegativeDirection = true ):
        coefficientFunction_( coefficientFunction ),
        densityFunction_( densityFunction ),
        airSpeedFunction_( airSpeedFunction ),
        massFunction_( massFunction ),
        referenceAreaFunction_( referenceAreaFunction )
    {
        coefficientMultiplier_ = areCoefficientsInNegativeDirection == true ? -1.0 : 1.0;
    }

    //! Get acceleration.
    /*!
     * Returns the aerodynamic acceleration. All data required for the computation is taken
     * from member variables, which are set to their latest values by the last call of the
     * updateMembers function.
     * The returned acceleration is in the same reference frame as the aerodynamic coefficients,
     * with the coefficients assumed to be in  positive direction in the frame.
     * \return Acceleration.
     * \sa updateMembers().
     */
    Eigen::Vector3d getAcceleration( )
    {
        return computeAerodynamicAcceleration(
                    0.5 * currentDensity_ * currentAirspeed_ * currentAirspeed_,
                    currentReferenceArea_, currentForceCoefficients_, currentMass_ );

    }

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
            currentForceCoefficients_ = coefficientMultiplier_ * this->coefficientFunction_( );
            currentDensity_ = this->densityFunction_( );
            currentMass_ = this->massFunction_( );
            currentAirspeed_ = this->airSpeedFunction_( );
            currentReferenceArea_ = this->referenceAreaFunction_( );
        }
    }

private:

    //! Function to retrieve the current aerodynamic force coefficients.
    const CoefficientReturningFunction coefficientFunction_;

    //! Function to retrieve the current density.
    const DoubleReturningFunction densityFunction_;

    //! Function to retrieve the current airspeed.
    const DoubleReturningFunction airSpeedFunction_;

    //! Function to retrieve the current mass.
    const DoubleReturningFunction massFunction_;

    //! Function to retrieve the current reference area.
    const DoubleReturningFunction referenceAreaFunction_;

    //! Current aerodynamic force coefficients.
    Eigen::Vector3d currentForceCoefficients_;

    //! Current density.
    double currentDensity_;

    //! Current airspeed.
    double currentAirspeed_;

    //! Current mass as set by massFunction_.
    double currentMass_;

    //! Current reference area, as set by referenceAreaFunction_.
    double currentReferenceArea_;

    //! Multiplier to reverse direction of coefficients.
    double coefficientMultiplier_;
};

//! Typedef for shared-pointer to AerodynamicAcceleration object.
typedef boost::shared_ptr< AerodynamicAcceleration > AerodynamicAccelerationPointer;

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ACCELERATION_H
