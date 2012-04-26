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
 *      110822    D. Dirkx          Removed pointer to double member, minor changes.
 *      120324    K. Kumar          Added missing include-statements; used TUDAT_NAN in class
 *                                  constructor initialization list instead of -0.0; minor
 *                                  corrections in Doxygen comments; updated input parameters for
 *                                  free functions; added astrodynamics namespace layer.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_FORCE_H
#define TUDAT_AERODYNAMIC_FORCE_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include <TudatCore/Basics/utilityMacros.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/forceModel.h"
#include "Tudat/Astrodynamics/States/state.h"

namespace tudat
{
namespace astrodynamics
{
namespace force_models
{

//! Compute the aerodynamic force in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic force. It takes primitive types as arguments to
 * perform the calculations. Therefore, these quantities (dynamicPressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param aerodynamicCoefficients. Aerodynamic coefficients in right-handed reference frame.
 * \return Resultant aerodynamic force, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::VectorXd computeAerodynamicForce( const double dynamicPressure,
                                         const double referenceArea,
                                         const Eigen::Vector3d& aerodynamicCoefficients );

//! Compute the aerodynamic force in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic force. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and coefficients are retrieved.
 * \return Resultant aerodynamic force, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::MatrixXd computeAerodynamicForce(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface );

//! Aerodynamic force model.
/*!
 * Calculates the aerodynamic force based on C_D C_S and C_L, a reference lenght or area, the
 * local velocity and the local density.
 */
class AerodynamicForce : public ForceModel
{
public:

    //! Typedef for shared pointer to state.
    /*!
     * Typedef for shared pointer to state.
     */
    typedef boost::shared_ptr< states::State > StatePointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicForce( boost::shared_ptr< AerodynamicCoefficientInterface >
                      aerodynamicCoefficientInterface ):
        aerodynamicCoefficientInterface_( aerodynamicCoefficientInterface ),
        dynamicPressure_( TUDAT_NAN )
    { }

    //! Get aerodynamic coefficient interface.
    /*!
     * Returns the pointer to the AerodynamicCoefficientInterface object.
     * \return Aerodynamic coefficient interface used to retrieve aerodynamic coefficients.
     */
    boost::shared_ptr< AerodynamicCoefficientInterface > getAerodynamicCoefficientInterface( )
    {
        return aerodynamicCoefficientInterface_;
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

    //! Compute aerodynamic force.
    /*!
     * Computes the force due to the gravity field in Newtons.
     * \param pointerToState Pointer to an object of the State class containing current state.
     * \param time Current time.
     */
    void computeForce( StatePointer state, const double time = 0.0 )
    {
        TUDAT_UNUSED_PARAMETER( state );
        TUDAT_UNUSED_PARAMETER( time );

        force_ = computeAerodynamicForce(
                    dynamicPressure_,
                    aerodynamicCoefficientInterface_->getReferenceArea( ),
                    aerodynamicCoefficientInterface_->getCurrentForceCoefficients( ) );
    }

protected:

private:

    //! Pointer to aerodynamic coefficient interface.
    /*!
     * Pointer to an aerodynamic coefficient interface.
     */
    boost::shared_ptr< AerodynamicCoefficientInterface > aerodynamicCoefficientInterface_;

    //! The dynamic pressure.
    /*!
     * The dynamic pressure.
     */
    double dynamicPressure_;
};

} // namespace force_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_FORCE_H
