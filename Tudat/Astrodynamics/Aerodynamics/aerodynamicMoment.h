/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110617    F.M. Engelen      First creation of code.
 *      110822    D. Dirkx          Removed pointer to double member, removed location of center
 *                                  of mass, minor changes.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_MOMENT_H
#define TUDAT_AERODYNAMIC_MOMENT_H

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/pureMomentModel.h"

namespace tudat
{
namespace astrodynamics
{
namespace moment_models
{

//! Compute the aerodynamic moment in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic moment. It takes primitive types as arguments to
 * perform the calculations. Therefor, these quantities (dynamic pressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param referenceLength Reference length of the aerodynamic coefficients. Note that this
 *          reference length is used for all three independent directions.
 * \param aerodynamicCoefficients. Aerodynamic moment coefficients in right-handed reference frame.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given, but with opposite sign. i.e., a positive drag
 *          coefficient will give a negative force in -x direction (in the aerodynamic frame).
 */
Eigen::MatrixXd computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const double referenceLength,
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
Eigen::MatrixXd computeAerodynamicMoment(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface );

class AerodynamicMoment : public PureMomentModel
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicMoment( AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface ) :
        pointerToAerodynamicCoefficientInterface_( pointerToAerodynamicCoefficientInterface ),
        dynamicPressure_( TUDAT_NAN )
    { }

    //! Set aerodynamic coefficient interface.
    /*!
     * Returns the pointer to the AerodunamicCoefficientDatabase object.
     * \return Aerodynamic coefficient interface used to retrieve aerodynamic coefficients.
     */
    AerodynamicCoefficientInterface* getAerodynamicCoefficientInterface( )
    {
        return pointerToAerodynamicCoefficientInterface_;
    }

    //! Set dynamic pressure.
    /*!
     * Sets the dynamic pressure.
     * \param dynamicPressure Dynamic pressure.
     */
    void setDynamicPressure( const double dynamicPressure ) { dynamicPressure_ = dynamicPressure; }

    //! Get dynamic pressure.
    /*!
     * Returns the dynamic pressure.
     * \return Dynamic pressure.
     */
    double getDynamicPressure( ) { return dynamicPressure_; }

    //! Compute aerodynamic moment.
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    void computeMoment( State* pointerToState = NULL, double time = 0.0 )
    {
        moment_ = computeAerodynamicMoment( dynamicPressure_,
                pointerToAerodynamicCoefficientInterface_->getReferenceArea( ),
                pointerToAerodynamicCoefficientInterface_->getReferenceLength( ),
                pointerToAerodynamicCoefficientInterface_->getCurrentMomentCoefficients( ) );
    }

protected:

private:

    //! Pointer to aerodynamic coefficient interface.
    /*!
     * Pointer to an aerodynamic coefficient interface.
     */
    AerodynamicCoefficientInterface* pointerToAerodynamicCoefficientInterface_;

    //! The dynamic pressure.
    /*!
     * The dynamic pressure.
     */
    double dynamicPressure_;
};

} // namespace moment_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_MOMENT_H
