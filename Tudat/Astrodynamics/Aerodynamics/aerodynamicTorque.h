/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_AERODYNAMIC_TORQUE_H
#define TUDAT_AERODYNAMIC_TORQUE_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic moment in same reference frame as input coefficients, with same reference lengths for each axis
/*!
 * This function calculates the aerodynamic moment. It takes primitive types as arguments to
 * perform the calculations. Therefor, these quantities (dynamic pressure, reference area, reference length and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param referenceLength Reference length of the aerodynamic coefficients. Note that this
 *          reference length is used for all three independent directions.
 * \param momentCoefficients Aerodynamic moment coefficients in right-handed reference frame.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given, but with opposite sign. i.e., a positive drag
 *          coefficient will give a negative force in -x direction (in the aerodynamic frame).
 */
Eigen::Vector3d computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const double referenceLength,
                                          const Eigen::Vector3d& momentCoefficients );

//! Compute the aerodynamic moment in same reference frame as input coefficients, with different reference lengths for each axis.
/*!
 * This function calculates the aerodynamic moment. It takes primitive types as arguments to
 * perform the calculations. Therefor, these quantities (dynamic pressure, reference area, reference lengths and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param referenceLengths Reference lengths of the aerodynamic coefficients, used in x, y and z directions.
 * \param momentCoefficients Aerodynamic moment coefficients in right-handed reference frame.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given, but with opposite sign. i.e., a positive drag
 *          coefficient will give a negative force in -x direction (in the aerodynamic frame).
 */
Eigen::Vector3d computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const Eigen::Vector3d& referenceLengths,
                                          const Eigen::Vector3d& momentCoefficients );

//! Compute the aerodynamic moment in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic moment. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area and length are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and length, and the moment coefficients are retrieved.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::Vector3d computeAerodynamicMoment(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface );

//! Class for calculation of aerodynamic torques.
class AerodynamicTorque : public basic_astrodynamics::TorqueModel
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double ( ) > DoubleReturningFunction;

    //! Typedef for coefficient-returning function.
    typedef boost::function< Eigen::Vector3d( ) > CoefficientReturningFunction;

public:

    //! Acceleration torque constructor.
    /*!
     * Acceleration torque constructor, taking function pointers for all member variables.
     * \param coefficientFunction Function which retrieves current values of aerodynamic moment coefficients.
     * \param densityFunction Function which retrieves current value of the density.
     * \param airSpeedFunction Function which retrieves current value of the airspeed.
     * \param referenceAreaFunction Function which retrieves current value of the aerodynamic coefficient reference area.
     * \param referenceLengthsFunction Function which retrieves current values of the aerodynamic  coefficient reference lengths.
     * \param areCoefficientsInNegativeDirection Boolean that determines whether to invert
     *          direction of aerodynamic coefficients. This is typically done for lift, drag and
     *          side force coefficients that point in negative direction in the local frame
     *          (default true).
     */
    AerodynamicTorque( const CoefficientReturningFunction coefficientFunction,
                       const DoubleReturningFunction densityFunction,
                       const DoubleReturningFunction airSpeedFunction,
                       const DoubleReturningFunction referenceAreaFunction,
                       const CoefficientReturningFunction referenceLengthsFunction,
                       const bool areCoefficientsInNegativeDirection = true ):
        coefficientFunction_( coefficientFunction ),
        densityFunction_( densityFunction ),
        airSpeedFunction_( airSpeedFunction ),
        referenceAreaFunction_( referenceAreaFunction ),
        referenceLengthsFunction_( referenceLengthsFunction )
    {
        coefficientMultiplier_ = areCoefficientsInNegativeDirection == true ? -1.0 : 1.0;
    }

    //! Get aerodynamic torque.
    /*!
     * Returns the aerodynamic torque. All data required for the computation is taken
     * from member variables, which are set to their latest values by the last call of the
     * updateMembers function.
     * The returned torque is in the same reference frame as the aerodynamic coefficients,
     * with the coefficients assumed to be in  positive direction in the frame.
     * \return Aerodynamic torque.
     * \sa updateMembers().
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the aerodynamic torque model.
    /*!
     * Updates member variables used by the aerodynamic accfeleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            currentMomentCoefficients_ = coefficientMultiplier_ * this->coefficientFunction_( );
            currentDensity_ = this->densityFunction_( );
            currentAirspeed_ = this->airSpeedFunction_( );
            currentReferenceArea_ = this->referenceAreaFunction_( );
            currentReferenceLengths_ =  this->referenceLengthsFunction_( );
            currentTorque_ = computeAerodynamicMoment(
                        0.5 * currentDensity_ * currentAirspeed_ * currentAirspeed_,
                        currentReferenceArea_, currentReferenceLengths_, currentMomentCoefficients_ );
            currentTime_ = currentTime;
        }
    }


private:

    //! Function to retrieve the current aerodynamic moment coefficients.
    const CoefficientReturningFunction coefficientFunction_;

    //! Function to retrieve the current density.
    const DoubleReturningFunction densityFunction_;

    //! Function to retrieve the current airspeed.
    const DoubleReturningFunction airSpeedFunction_;

    //! Function to retrieve the current reference area.
    const DoubleReturningFunction referenceAreaFunction_;

    const CoefficientReturningFunction referenceLengthsFunction_;

    //! Current aerodynamic moment coefficients, as set by updateMembers function.
    Eigen::Vector3d currentMomentCoefficients_;

    //! Current aerodynamic reference lengths, as set by updateMembers function.
    Eigen::Vector3d currentReferenceLengths_;

    //! Current aerodynamic torque, as set by updateMembers function.
    Eigen::Vector3d currentTorque_;

    //! Current density, as set by updateMembers function.
    double currentDensity_;

    //! Current airspeed, as set by updateMembers function.
    double currentAirspeed_;


    //! Current reference area, as set by updateMembers function.
    double currentReferenceArea_;

    //! Multiplier to reverse direction of coefficients.
    double coefficientMultiplier_;

};


} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_TORQUE_H
