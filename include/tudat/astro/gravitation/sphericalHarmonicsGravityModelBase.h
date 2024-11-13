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

#ifndef TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H
#define TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H

#include <functional>

namespace tudat
{
namespace gravitation
{

//! Template base class for spherical harmonics gravitational acceleration models.
/*!
 * This //template class serves as the base class for the
 * SphericalHarmonicsGravitationalAccelerationModel, CentralGravitationalAccelerationModel,
 * CentralJ2GravitationalAccelerationModel, CentralJ2J3GravitationalAccelerationModel, and
 * CentralJ2J3J4GravitationalAccelerationModel classes.
 * \tparam StateMatrix Type used to store a state matrix.
 */
template< typename StateMatrix >
class SphericalHarmonicsGravitationalAccelerationModelBase
{
protected:

    //! Typedef for a position-returning function.
    typedef std::function< void( StateMatrix& ) > StateFunction;

public:

    //! Default constructor taking position of body subject to acceleration, variable
    //! gravitational parameter, and position of body exerting acceleration.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, a pointer to a function returning the gravitational parameter of
     * the body exerting the acceleration, and a pointer to a function returning the position of
     * the body exerting the  gravitational acceleration (typically the central body). The
     * constructor also updates all the internal members.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter Gravitational parameter of body exerting gravitational acceleration.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    SphericalHarmonicsGravitationalAccelerationModelBase(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const StateFunction positionOfBodyExertingAccelerationFunction,
            const bool isMutualAttractionUsed,
            const bool updateGravitationalPotential = false )
        : subjectPositionFunction( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction( [ = ]( ){ return aGravitationalParameter; } ),
          sourcePositionFunction( positionOfBodyExertingAccelerationFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          currentPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential )
    { }

    //! Default constructor taking position of body subject to acceleration, variable
    //! gravitational parameter, and position of body exerting acceleration.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, a pointer to a function returning the gravitational parameter of
     * the body exerting the acceleration, and a pointer to a function returning the position of
     * the body exerting the  gravitational acceleration (typically the central body). The
     * constructor also updates all the internal members.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameterFunction Pointer to function returning gravitational parameter
     *          of body exerting gravitational acceleration.
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    SphericalHarmonicsGravitationalAccelerationModelBase(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const std::function< double( ) > aGravitationalParameterFunction,
            const StateFunction positionOfBodyExertingAccelerationFunction,
            const bool isMutualAttractionUsed,
            const bool updateGravitationalPotential = false )
        : subjectPositionFunction( positionOfBodySubjectToAccelerationFunction ),
          gravitationalParameterFunction( aGravitationalParameterFunction ),
          sourcePositionFunction( positionOfBodyExertingAccelerationFunction ),
          isMutualAttractionUsed_( isMutualAttractionUsed ),
          currentPotential_( TUDAT_NAN ),
          updatePotential_( updateGravitationalPotential )
    { }

    //! Virtual destructor.
    /*!
     * Base class virtual destructor.
     */
    virtual ~SphericalHarmonicsGravitationalAccelerationModelBase( ){ }

    //! Update base class members.
    /*!
     * Updates members of base class. In this case, the current position of the body subject to and
     * exerting the gravitational acceleration, as well as the gravitational parameter are updated
     * by calling the functions provided through the constructor.
     * \return True; this should be modified to return a flag indicating if the update was
     *          successful.
     */
    void updateBaseMembers( )
    {
        this->gravitationalParameter = this->gravitationalParameterFunction( );
        this->subjectPositionFunction( this->positionOfBodySubjectToAcceleration );
        this->sourcePositionFunction( this->positionOfBodyExertingAcceleration );
    }

    //! Function to return the function returning the relevant gravitational parameter.
    /*!
     * Function to return the function returning the relevant gravitational parameter.
     * \return Function returning the gravitational parameter used in the computations.
     */
    std::function< double( ) > getGravitationalParameterFunction( )
    { return gravitationalParameterFunction; }

    //! Function to return the function returning position of body exerting acceleration.
    /*!
     * Function to return the function returning position of body exerting acceleration.
     * \return Function returning position of body exerting acceleration.
     */
    StateFunction getStateFunctionOfBodyExertingAcceleration( )
    { return sourcePositionFunction; }

    //! Function to return the function returning position of body subject to acceleration.
    /*!
     * Function to return the function returning position of body subject to acceleration.
     * \return Function returning position of body subject to acceleration.
     */
    StateFunction getStateFunctionOfBodyUndergoingAcceleration( )
    { return subjectPositionFunction; }

    //! Function to return whether the mutual attraction is used.
    /*!
     *  Function to return whether the mutual attraction is used.
     * \return Boolean defining whether the mutual attraction is used.
     */
    bool getIsMutualAttractionUsed( )
    {
        return isMutualAttractionUsed_;
    }

    //! Function to return the current gravitational paramter.
    /*!
     *  Function to return the current gravitational paramter. Must be called after updateBaseMembers( ).
     * \return The value of the current gravitational parameter.
     */
    double getCurrentGravitationalParameter( )
    {
        return gravitationalParameter;
    }

    //! Function to return current position vector of body exerting gravitational acceleration in inertial frame.
    /*!
     *  Function to return current position vector of body exerting gravitational acceleration in inertial frame.
     * \return Current position vector of body exerting gravitational acceleration in inertial frame.
     */
    StateMatrix getCurrentPositionOfBodySubjectToAcceleration( )
    {
        return positionOfBodySubjectToAcceleration;
    }

    //! Function to return current position vector of body undergoing gravitational acceleration in inertial frame.
    /*!
     *  Function to return current position vector of body undergoing gravitational acceleration in inertial frame.
     * \return Current position vector of body undergoing gravitational acceleration in inertial frame.
     */
    StateMatrix getCurrentPositionOfBodyExertingAcceleration( )
    {
        return positionOfBodyExertingAcceleration;
    }

    //! Function to return the value of the current gravitational potential.
    double getCurrentPotential ( )
    { return currentPotential_; }

    //! Function to return the update potential flag.
    bool getUpdatePotential ( )
    { return updatePotential_; }

    //! Function to reset the update potential flag.
    void resetUpdatePotential ( bool updatePotential )
    { updatePotential_ = updatePotential; }


protected:

    //! Position of body subject to acceleration.
    /*!
     * Current position vector of body subject to gravitational acceleration in inertial frame.
     */
    StateMatrix positionOfBodySubjectToAcceleration;

    //! Pointer to function returning position of body subject to acceleration.
    /*!
     * Pointer to function that returns the current position of the body subject to the
     * gravitational acceleration.
     */
    const StateFunction subjectPositionFunction;

    //! Function returning a gravitational parameter [m^3 s^-2].
    /*!
     * Function returning current gravitational parameter of body exerting acceleration [m^3 s^-2].
     */
    const std::function< double( ) > gravitationalParameterFunction;

    //! Gravitational parameter [m^3 s^-2].
    /*!
     * Current gravitational parameter of body exerting acceleration [m^3 s^-2].
     */
    double gravitationalParameter;

    //! Position of body exerting acceleration.
    /*!
     * Current position vector of body exerting gravitational acceleration in inertial frame.
     */
    StateMatrix positionOfBodyExertingAcceleration;

    //! Pointer to function returning position of body exerting acceleration.
    /*!
     * Pointer to function that returns the current position of the body exerting the
     * gravitational acceleration.
     */
    const StateFunction sourcePositionFunction;

    //! Variable denoting whether mutual acceleration between bodies is included.
    /*!
     * Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    bool isMutualAttractionUsed_;

    //! Current gravitational potential acting on the body undergoing acceleration, as computed by last call to
    //! updateMembers function
    double currentPotential_;

    //! Flag indicating whether to update the gravitational potential when calling the updateMembers function.
    bool updatePotential_;


private:
};

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_GRAVITATIONAL_ACCELERATION_MODEL_BASE_H
