/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110608    F.M. Engelen      Creation of code.
 *      110714    D. Dirkx          Class name change and other minor changes during code check.
 *      110715    F.M. Engelen      Added the virtual compute function.
 *      110810    J. Leloux         Corrected doxygen documentation (function variable name).
 *
 *    References
 *
 *    The computeCoefficients() function is not yet implemented in any derived classes and is
 *    therefore not pure virtual in this base class.
 *
 */

#ifndef TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <Eigen/Core>

namespace tudat
{
namespace aerodynamics
{

//! Base class to hold an aerodynamic coefficient interface.
/*!
 * This interface can, for instance, be a database of coefficients or an aerodynamic analysis code
 * which generates coefficients.
 */
class AerodynamicCoefficientInterface
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    AerodynamicCoefficientInterface( ) : currentForceCoefficients_( Eigen::Vector3d::Zero( ) ),
        currentMomentCoefficients_( Eigen::Vector3d::Zero( ) ), referenceLength_( -0.0 ),
        referenceArea_( -0.0 ), lateralReferenceLength_( -0.0 ),
        momentReferencePoint_( Eigen::Vector3d::Zero( ) ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~AerodynamicCoefficientInterface( ) { }

    //! Set reference area.
    /*!
     * Sets reference area used to non-dimensionalize aerodynamic forces and moments.
     * \param referenceArea Aerodynamic reference area.
     */
    void setReferenceArea( double referenceArea ) { referenceArea_ = referenceArea; }

    //! Get reference area.
    /*!
     * Returns reference area used to non-dimensionalize aerodynamic forces and moments.
     * \return Aerodynamic reference area.
     */
    double getReferenceArea( ) { return referenceArea_; }

    //! Set reference length.
    /*!
     * Sets reference length used to non-dimensionalize aerodynamic moments.
     * \param referenceLength Aerodynamic reference length.
     */
    void setReferenceLength( double referenceLength )
    { referenceLength_ = referenceLength; }

    //! Get reference length.
    /*!
     * Returns reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic reference length.
     */
    double getReferenceLength( ) { return referenceLength_; }

    //! Set lateral reference length.
    /*!
     * Sets lateral reference length used to non-dimensionalize aerodynamic moments.
     * \param lateralReferenceLength Aerodynamic reference length.
     */
    void setLateralReferenceLength( double lateralReferenceLength )
    { lateralReferenceLength_ = lateralReferenceLength; }

    //! Get lateral reference length.
    /*!
     * Returns lateral reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic lateral reference length.
     */
    double getLateralReferenceLength( ) { return lateralReferenceLength_; }

    //! Set moment reference point.
    /*!
     * Sets the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \param momentReferencePoint Aerodynamic moment reference point.
     */
    void setMomentReferencePoint( const Eigen::Vector3d& momentReferencePoint )
    { momentReferencePoint_ = momentReferencePoint; }

    //! Get moment reference point.
    /*!
     * Returns the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \return Aerodynamic reference point.
     */
    Eigen::VectorXd getMomentReferencePoint( ) { return momentReferencePoint_; }

    //! Get the current force coefficients.
    /*!
     * Returns the force coefficients that have been set as the 'current' force coefficients,
     * i.e. at the current flight condition.
     * \return current force coefficients.
     */
    Eigen::Vector3d getCurrentForceCoefficients( ) { return currentForceCoefficients_; }

    //! Get the moment coefficients
    /*!
     * Return the moment coefficients that have been set as the 'current' moment coefficients,
     * i.e. at the current flight condition.
     * \return current moment coefficients.
     */
    Eigen::Vector3d getCurrentMomentCoefficients( ) { return currentMomentCoefficients_; }

    //! Set the force coefficients.
    /*!
     * Sets the current force coefficients, i.e. at the current flight condition.
     * \param currentForceCoefficients the current force coefficients.
     */
    void setCurrentForceCoefficients( const Eigen::Vector3d& currentForceCoefficients )
    { currentForceCoefficients_ = currentForceCoefficients; }

    //! Set the moment coefficients.
    /*!
     * Sets the current moment coefficients, i.e. at the current flight condition.
     * \param currentMomentCoefficients the current force coefficients.
     */
    void setCurrentMomentCoefficients( const Eigen::Vector3d& currentMomentCoefficients )
    { currentMomentCoefficients_ = currentMomentCoefficients; }

    //! Compute the aerodynamic coefficients at current flight condition.
    /*!
     * Computes the current force and moment coefficients and is to be
     * implemented in derived classes. Such a calculation would be performed by, for instance,
     * including a pointer to a Vehicle object, from which the current flight condition can
     * be retrieved. The function here is not made pure virtual pending the inclusion of
     * such functionality in derived classes.
     */
    virtual void computeCurrentCoefficients( ) { }

protected:

    //! The current force coefficients.
    /*!
     * The force coefficients at the current flight condition.
     */
    Eigen::Vector3d currentForceCoefficients_;

    //! The current moment coefficients.
    /*!
     * The moment coefficients at the current flight condition.
     */
    Eigen::Vector3d currentMomentCoefficients_;

    //! Aerodynamic reference length.
    /*!
     * Reference length with which aerodynamic moments are non-dimensionalized.
     */
    double referenceLength_;

    //! Aerodynamic reference area.
    /*!
     * Reference area with which aerodynamic forces and moments are non-dimensionalized.
     */
    double referenceArea_;

    //! Lateral aerodynamic reference length.
    /*!
     * Lateral reference length with which aerodynamic moments are non-dimensionalized.
     */
    double lateralReferenceLength_;

    //! Aerodynamic moment reference point.
    /*!
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;

private:

};

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H
