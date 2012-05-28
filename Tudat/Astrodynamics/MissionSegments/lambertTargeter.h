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
 *      101111    E. Iorfida        Creation of code.
 *      101126    E. Iorfida        Virtual void "solve" added.
 *      101206    E. Iorfida        setInitialLambertState added, protected variables added.
 *      101207    E. Iorfida        Set single variables, change variables names in more
 *                                  understandable ones.
 *      110113    E. Iorfida        Added necessary elements to build pointer-to-member-function to
 *                                  RootFinderAlgorithms and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to use the last
 *                                  version of Newton-Raphson code.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfVelocityVector' became
 *                                  'speed'. Also corrected 'tangential' to 'transverse'.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names
 *                                  (from heliocentric, to inertial).
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      Gooding, R.H. A procedure for the solution of Lambert's orbital boundary-value problem,
 *          Celestial Mechanics and Dynamical Astronomy, 48:145-165, 1990.
 *
 *    The number of revolutions from departure to arrival body is zero by definition in this
 *    routine. This can be made user-defined later on. The resulting trajectories are in
 *    anti-clockwise direction.
 *
 */

#ifndef TUDAT_LAMBERT_TARGETER_H
#define TUDAT_LAMBERT_TARGETER_H

#include <iostream>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/MissionSegments/trajectoryDesignMethod.h"
#include "Tudat/Astrodynamics/States/cartesianPositionElements.h"
#include "Tudat/Astrodynamics/States/cartesianVelocityElements.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

//! Lambert targeting algorithm class.
/*!
 * Implementation of Lambert targeting algorithm in Tudat.
 */
class LambertTargeter : public TrajectoryDesignMethod
{
public:

    //! Typedef for shared pointer to celestial body.
    /*!
     * Typedef for shared pointer to celestial body.
     */
    typedef boost::shared_ptr< bodies::CelestialBody > CelestialBodyPointer;

    //! Typedef for shared pointer to Newton-Raphson method.
    /*!
     * Typedef for shared pointer to Newton-Raphson method.
     */
    typedef boost::shared_ptr< NewtonRaphson > NewtonRaphsonPointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LambertTargeter( )
        : numberOfRevolutions_( TUDAT_NAN ),
          timeOfFlight_( TUDAT_NAN ),
          lambertSemiMajorAxis_( TUDAT_NAN ),
          xParameter_( TUDAT_NAN ),
          normalizedTimeOfFlight_( TUDAT_NAN ),
          qParameter_( TUDAT_NAN ),
          radialSpeedAtDeparture_( TUDAT_NAN ),
          radialSpeedAtArrival_( TUDAT_NAN ),
          transverseSpeedAtDeparture_( TUDAT_NAN ),
          transverseSpeedAtArrival_( TUDAT_NAN )
    {
        cartesianPositionAtDeparture_.setZero( 3 );
        cartesianPositionAtArrival_.setZero( 3 );
        cartesianVelocityAtDeparture_.setZero( 3 );
        cartesianVelocityAtArrival_.setZero( 3 );
    }

    //! Set position at departure.
    /*!
     * Sets the position at departure as cartesian Vector3d.
     * \param cartesianPositionAtDeparture Position at departure.
     */
    void setPositionAtDeparture( Eigen::Vector3d cartesianPositionAtDeparture )
    {
        cartesianPositionAtDeparture_ = cartesianPositionAtDeparture;
    }

    //! Set position at arrival.
    /*!
     * Sets position at arrival as cartesian Vector3d
     * \param cartesianPositionAtArrival Position at arrival.
     */
    void setPositionAtArrival( Eigen::Vector3d cartesianPositionAtArrival )
    {
        cartesianPositionAtArrival_ = cartesianPositionAtArrival;
    }

    //! Set number of revolutions.
    /*!
     * Sets the number of revolutions.
     * \param numberOfRevolutions Number of Revolutions.
     */
    void setNumberOfRevolutions( const int numberOfRevolutions )
    {
        numberOfRevolutions_ = numberOfRevolutions;
    }

    //! Set time-of-flight.
    /*!
     * Sets the time-of-flight.
     * \param timeOfFlight Time-of-flight.
     */
    void setTimeOfFlight( const double timeOfFlight ) { timeOfFlight_ = timeOfFlight; }

    //! Set central body.
    /*!
     * Sets pointer to central body.
     * \param celestialBody Central body
     */
    void setCentralBody( CelestialBodyPointer celestialBody ) { celestialBody_ = celestialBody; }

    //! Set pointer to Newton-Raphson method for Lambert targeting algorithm.
    /*!
     * Sets a pointer to the Newton-Raphson method for the Lambert targeting algorithm.
     * \param newtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphsonMethod( NewtonRaphsonPointer newtonRaphson )
    {
        newtonRaphson_ = newtonRaphson;
    }

    //! Get Lambert semi-major axis.
    /*!
     * Returns the Lambert semi-major axis.
     * \return Lambert semi-major axis.
     */
    double getLambertSemiMajorAxis( ) { return lambertSemiMajorAxis_; }

    //! Get radial speed at departure.
    /*!
     * Returns the radial speed at departure.
     * \return Radial speed at departure.
     */
    double getRadialSpeedAtDeparture( ) { return radialSpeedAtDeparture_; }

    //! Get radial speed at arrival.
    /*!
     * Returns the radial speed at arrival.
     * \return Radial speed at arrival.
     */
    double getRadialSpeedAtArrival( ) { return radialSpeedAtArrival_; }

    //! Get transverse speed at departure.
    /*!
     * Returns the transverse speed at departure.
     * \return Transverse speed at departure.
     */
    double getTransverseSpeedAtDeparture( ) { return transverseSpeedAtDeparture_; }

    //! Get transverse speed at arrival.
    /*!
     * Returns the transverse speed at arrival.
     * \return Transverse speed at arrival.
     */
    double getTransverseSpeedAtArrival( ) { return transverseSpeedAtArrival_; }

    //! Get inertial velocity at departure.
    /*!
     * Returns the inertial velocity at departure ( heliocentric or
     * planetocentric ).
     * \return Inertial velocity at departure.
     */
    Eigen::Vector3d getInertialVelocityAtDeparture( )
    {
        return cartesianVelocityAtDeparture_;
    }

    //! Get inertial velocity at arrival.
    /*!
     * Returns the inertial velocity at arrival ( heliocentric or
     * planetocentric ).
     * \return Inertial velocity at arrival.
     */
    Eigen::Vector3d getInertialVelocityAtArrival( ) { return cartesianVelocityAtArrival_; }

    //! Execute Lambert targeting algorithm.
    /*!
     * Executes the Lambert targeting algorithm.
     */
    void execute( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param lambertTargeter Lambert targeter object.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, LambertTargeter& lambertTargeter );

protected:

private:

    //! Number of Revolutions.
    /*!
     * Number of Revolutions.
     */
    int numberOfRevolutions_;

    //! Time-of-flight.
    /*!
     * Time-of-flight.
     */
    double timeOfFlight_;

    //! Pointer to CelestialBody class.
    /*!
     * Pointer to CelestialBody class.
     */
    CelestialBodyPointer celestialBody_;

    //! Cartesian position of object at departure.
    /*!
     * Cartesian position of object at departure.
     */
    Eigen::Vector3d cartesianPositionAtDeparture_;

    //! Cartesian position of object at arrival.
    /*!
     * Cartesian position of object at arrival.
     */
    Eigen::Vector3d cartesianPositionAtArrival_;

    //! Lambert semi-major axis.
    /*!
     * Lambert semi-major axis.
     */
    double lambertSemiMajorAxis_;

    //! Lambert x-parameter.
    /*!
     * Lambert x-parameter.
     */
    double xParameter_;

    //! Normalized time of flight for Lambert implementation.
    /*!
     * Normalized time of flight for Lambert implementation.
     */
    double normalizedTimeOfFlight_;

    //! Lambert q-parameter.
    /*!
     * Lambert q-parameter.
     */
    double qParameter_;

    //! Radial speed at departure.
    /*!
     * Radial speed at departure.
     */
    double radialSpeedAtDeparture_;

    //! Radial speed at arrival.
    /*!
     * Radial speed at arrival.
     */
    double radialSpeedAtArrival_;

    //! Transverse speed at departure.
    /*!
     * Transverse speed at departure.
     */
    double transverseSpeedAtDeparture_;

    //! Transverse speed at arrival.
    /*!
     * Transverse speed at arrival.
     */
    double transverseSpeedAtArrival_;

    //! Cartesian velocity of object at departure.
    /*!
     * Cartesian velocity of object at departure
     * that represents the inertial velocity at arrival (referred to central
     * body).
     */
    Eigen::Vector3d cartesianVelocityAtDeparture_;

    //! Cartesian velocity of object at arrival.
    /*!
     * Cartesian velocity of object at arrival
     * that represents the inertial velocity at arrival (referred to central
     * body).
     */
    Eigen::Vector3d cartesianVelocityAtArrival_;

    //! Pointer to object of NewtonRaphson class.
    /*!
     * Pointer to object of NewtonRaphson class.
     */
    NewtonRaphsonPointer newtonRaphson_;

    //! Pointer to adaptor object of NewtonRaphsonAdaptor class.
    /*!
     * Pointer to adaptor object of NewtonRaphsonAdaptor class. The template
     * parameter passed is this class.
     */
    NewtonRaphsonAdaptor< LambertTargeter > newtonRaphsonAdaptorForLambertTargeter_;

    //! Define Lambert function for positive lambertEccentricAnomaly_.
    /*!
     * Defines Lambert function for positive lambertEccentricAnomaly_.
     */
    double lambertFunctionPositive( double& xParameter_ );

    //! Define Lambert function for negative lambertEccentricAnomaly_.
    /*!
     * Defines Lambert function for negative lambertEccentricAnomaly_.
     */
    double lambertFunctionNegative( double& xParameter_ );

    //! Define first derivative of Lambert function for positive
    //! lambertEccentricAnomaly_.
    /*!
     * Defines first derivative of Lambert function for positive
     * lambertEccentricAnomaly_.
     */
    double lambertFirstDerivativeFunctionPositive( double& xParameter_ );

    //! Define first derivative of Lambert function for negative
    //! lambertEccentricAnomaly_.
    /*!
     * Defines first derivative of Lambert function for negative
     * lambertEccentricAnomaly_.
     */
    double lambertFirstDerivativeFunctionNegative( double& xParameter_ );

    //! Define general Lambert function.
    /*!
     * Defines general Lambert function.
     */
    double lambertFunction( double& xParameter_ );

    //! Define first derivative of general Lambert function.
    /*!
     * Defines first derivative of general Lambert function.
     */
    double lambertFirstDerivativeFunction( double& xParameter_ );
};

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_LAMBERT_TARGETER_H
