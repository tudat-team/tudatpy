/*! \file lambertTargeter.h
 *    Source file of the Lambert targeting solver implemented in Tudat.
 *
 *    Path              : /Astrodynamics/MissionSegments/LambertTargeter/
 *    Version           : 9
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2010
 *    Last modified     : 8 February, 2011
 *
 *    References        :
 *      Gooding, R.H. A procedure for the solution of Lambert's orbital
 *          boundary-value problem, Celestial Mechanics and Dynamical
 *          Astronomy, 48:145-165, 1990.
 *
 *    Notes
 *      The number of revolutions from departure to arrival body is zero
 *      by definition in this routine. This can be made user-defined later on.
 *      The resulting trajectories are in anti-clockwise direction.
 *      At the moment the CartesianVelocityElements are defined with their
 *      pointers, that are referenced to the objects of the class.
 *      In the future it should be possibile to have only objects of
 *      CartesianVelocityElements class with their direct reference.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      101111    E. Iorfida        First creation of code.
 *      101126    E. Iorfida        Virtual void "solve" added.
 *      101206    E. Iorfida        setInitialLambertState added, protected
 *                                  variables added.
 *      101207    E. Iorfida        Set single variables, change variables
 *                                  names in more understandable ones.
 *      110113    E. Iorfida        Added necessary elements to build
 *                                  pointer-to-member-function to
 *                                  RootFinderAlgorithms and
 *                                  NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to
 *                                  use the last version of Newton-Raphson code.
 *      110130    J. Melman         Simplified variable names, e.g.,
 *                                  'normOfVelocityVector' became 'speed'.
 *                                  Also corrected 'tangential' to 'transverse'.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified
 *                                  variable names (from heliocentric, to
 *                                  inertial).
 *      110208    E. Iorfida        Added CartesianPositionElements objects as
 *                                  input and CartesianVelocityElements objects
 *                                  as output.
 */

#ifndef LAMBERTTARGETER_H
#define LAMBERTTARGETER_H

// Include statements.
#include "trajectoryDesignMethod.h"
#include "newtonRaphson.h"
#include "newtonRaphsonAdaptor.h"
#include "basicMathematicsFunctions.h"
#include "linearAlgebra.h"
#include "celestialBody.h"
#include "cartesianPositionElements.h"
#include "cartesianVelocityElements.h"

//! Lambert targeting algorithm class.
/*!
 * Implementation of Lambert targeting algorithm in Tudat.
 */
class LambertTargeter : public TrajectoryDesignMethod
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    LambertTargeter( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~LambertTargeter( );

    //! Set position at departure.
    /*!
     * Sets the position at departure as pointer to object of
     * CartesianPositionElements class.
     * \param pointerToCartesianPositionAtDeparture Position at departure.
     */
    void setPositionAtDeparture( CartesianPositionElements*
                                 pointerToCartesianPositionAtDeparture );

    //! Set position at arrival.
    /*!
     * Sets position at arrival as pointer to object of
     * CartesianPositionElements class.
     * \param pointerToCartesianPositionAtArrival Position at arrival.
     */
    void setPositionAtArrival( CartesianPositionElements*
                               pointerToCartesianPositionAtArrival );

    //! Set number of revolutions.
    /*!
     * Sets the number of revolutions.
     * \param numberOfRevolutions Number of Revolutions.
     */
    void setNumberOfRevolutions( const int& numberOfRevolutions );

    //! Set time-of-flight.
    /*!
     * Sets the time-of-flight.
     * \param timeOfFlight Time-of-flight.
     */
    void setTimeOfFlight( const double& timeOfFlight );

    //! Set central body.
    /*!
     * Sets pointer to central body.
     * \param pointerToCelestialBody Central body
     */
    void setCentralBody( CelestialBody* pointerToCelestialBody );

    //! Set pointer to Newton-Raphson method for Lambert targeting algorithm.
    /*!
     * Sets a pointer to the Newton-Raphson method for the Lambert targeting
     * algorithm.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson
     *          object.
     */
    void setNewtonRaphsonMethod( NewtonRaphson* pointerToNewtonRaphson );

    //! Get Lambert semi-major axis.
    /*!
     * Returns the Lambert semi-major axis.
     * \return Lambert semi-major axis.
     */
    double& getLambertSemiMajorAxis( );

    //! Get radial speed at departure.
    /*!
     * Returns the radial speed at departure.
     * \return Radial speed at departure.
     */
    double& getRadialSpeedAtDeparture( );

    //! Get radial speed at arrival.
    /*!
     * Returns the radial speed at arrival.
     * \return Radial speed at arrival.
     */
    double& getRadialSpeedAtArrival( );

    //! Get transverse speed at departure.
    /*!
     * Returns the transverse speed at departure.
     * \return Transverse speed at departure.
     */
    double& getTransverseSpeedAtDeparture( );

    //! Get transverse speed at arrival.
    /*!
     * Returns the transverse speed at arrival.
     * \return Transverse speed at arrival.
     */
    double& getTransverseSpeedAtArrival( );

    //! Get inertial velocity at departure.
    /*!
     * Returns the inertial velocity at departure ( heliocentric or
     * planetocentric ).
     * \return Inertial velocity at departure.
     */
    CartesianVelocityElements* getInertialVelocityAtDeparture( );

    //! Get inertial velocity at arrival.
    /*!
     * Returns the inertial velocity at arrival ( heliocentric or
     * planetocentric ).
     * \return Inertial velocity at arrival.
     */
    CartesianVelocityElements* getInertialVelocityAtArrival( );

    //! Execute Lambert targeting algorithm.
    /*!
     * Executes the Lambert targeting algorithm.
     */
    void execute( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param pointerToLambertTargeter Pointer to Lambert Targeter.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     LambertTargeter& lambertTargeter );

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
    CelestialBody* pointerToCelestialBody_;

    //! Pointer to CartesianPositionElements object at departure.
    /*!
     * Pointer to CartesianPositionElements object at departure.
     */
    CartesianPositionElements* pointerToCartesianPositionAtDeparture_;

    //! Pointer to CartesianPositionElements object at arrival.
    /*!
     * Pointer to CartesianPositionElements object at arrival.
     */
    CartesianPositionElements* pointerToCartesianPositionAtArrival_;

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

    //! Pointer to CartesianVelocityElements object at departure.
    /*!
     * Pointer to CartesianVelocityElements object at departure
     * that represents the inertial velocity at departure (referred to central
     * body).
     */
    CartesianVelocityElements* pointerToCartesianVelocityAtDeparture_;

    //! Pointer to CartesianVelocityElements object at arrival.
    /*!
     * Pointer to CartesianVelocityElements object at arrival
     * that represents the inertial velocity at arrival (referred to central
     * body).
     */
    CartesianVelocityElements* pointerToCartesianVelocityAtArrival_;

    //! Pointer to object of NewtonRaphson class.
    /*!
     * Pointer to object of NewtonRaphson class.
     */
    NewtonRaphson* pointerToNewtonRaphson_;

    //! Pointer to adaptor object of NewtonRaphsonAdaptor class.
    /*!
     * Pointer to adaptor object of NewtonRaphsonAdaptor class. The template
     * parameter passed is this class.
     */
    NewtonRaphsonAdaptor< LambertTargeter >
            newtonRaphsonAdaptorForLambertTargeter_;

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

#endif // LAMBERTTARGETER_H

// End of file.
