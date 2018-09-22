/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THRUSTACCELERATIONMODEL_H
#define TUDAT_THRUSTACCELERATIONMODEL_H

#include <limits>

#include <functional>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"
#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"

namespace tudat
{

namespace propulsion
{

//! Function used as interface to merge two update functions
/*!
 * Function used as interface to merge two update functions. Calling this function will update the two update functions
 * (in order)
 * \param updateFunction1 First update function.
 * \param updateFunction2 Second update function.
 * \param time Time to which both functions are to be updated
 */
inline void mergeUpdateFunctions(
        const std::function< void( const double ) > updateFunction1,
        const std::function< void( const double ) > updateFunction2,
        const double time )
{
    updateFunction1( time );
    updateFunction2( time );
}

//! Class used for computing an acceleration due to a continuous thrust.
/*!
 *  Class used for computing an acceleration due to a continuous thrust. The thrust magnitude and direction (in the
 *  propagation frame) are retrieved from separate functions provided by tye user.
 */
class ThrustAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    /*!
     * Constructor.
     * \param thrustMagnitudeFunction Function returning the current magnitude of the thrust. Any dependencies of the
     * thrust on (in)dependent variables is to be handled by the thrustUpdateFunction.
     * \param thrustDirectionFunction Function returning the direction of the thrust (as a unit vector).
     * Any dependencies of the thrust on (in)dependent variables is to be handled by the thrustUpdateFunction.
     * \param bodyMassFunction Function returning the current mass of the body being propagated.
     * \param massRateFunction Function returning total propellant mass rate from the thrust system.
     * \param associatedThrustSource ID associated with the source of the thrust (i.e. engine name).
     * \param thrustUpdateFunction Function used to update the thrust magnitude and direction to current time (default empty)
     * \param timeResetFunction Function to reset the time in the classes to which the thrustUpdateFunction function directs,
     * default empty.
     * \param requiredModelUpdates List of environment models that are to be updated before computing the acceleration,
     * list is included here to account for versatility of dependencies of thrust model (guidance) algorithms. Default empty.
     */
    ThrustAcceleration(
            const std::function< double( ) > thrustMagnitudeFunction,
            const std::function< Eigen::Vector3d( ) > thrustDirectionFunction,
            const std::function< double( ) > bodyMassFunction,
            const std::function< double( ) > massRateFunction,
            const std::string associatedThrustSource = "",
            const std::function< void( const double ) > thrustUpdateFunction = std::function< void( const double ) >( ),
            const std::function< void( const double ) > timeResetFunction = std::function< void( const double ) >( ),
            const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& requiredModelUpdates =
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >( ) ):
        AccelerationModel< Eigen::Vector3d >( ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        thrustDirectionFunction_( thrustDirectionFunction ),
        bodyMassFunction_( bodyMassFunction ),
        massRateFunction_( massRateFunction ),
        associatedThrustSource_( associatedThrustSource ),
        thrustUpdateFunction_( thrustUpdateFunction ),
        timeResetFunction_( timeResetFunction ),
        requiredModelUpdates_( requiredModelUpdates ){ }

    //! Destructor
    ~ThrustAcceleration( ){ }

    //! Function to retrieve the current acceleration, as set by last call to updateMembers function.
    /*!
     * Function to retrieve the current acceleration, as set by last call to updateMembers function.
     * \return Current thrust acceleration.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;

        if( !( timeResetFunction_ == nullptr ) )
        {
            timeResetFunction_( currentTime_ );
        }
    }

    //! Update member variables used by the thrust acceleration model.
    /*!
     * Updates member variables used by the thrust acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls and combines their output to
     * compute the acceleration.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        // Check if update is needed
        if( !( currentTime_ == currentTime ) )
        {
            // Update thrust dependencies if needed
            if( !( thrustUpdateFunction_ == nullptr ) )
            {
                thrustUpdateFunction_( currentTime );
            }

            // Retrieve thrust direction.
            currentAccelerationDirection_ = thrustDirectionFunction_( );

            if( ( std::fabs( currentAccelerationDirection_.norm( ) ) - 1.0 ) >
                    10.0 * std::numeric_limits< double >::epsilon( ) )
            {
                throw std::runtime_error( "Error in thrust acceleration, direction is not a unit vector"  );
            }

            // Retrieve magnitude of thrust and mass rate.
            currentThrustMagnitude_ = thrustMagnitudeFunction_( );
            currentMassRate_ = -massRateFunction_( );

            // Compute acceleration due to thrust
            currentAcceleration_ = currentAccelerationDirection_ * currentThrustMagnitude_ / bodyMassFunction_( );

            // Reset current time.
            currentTime_ = currentTime;

        }

    }

    //! Function to retreieve the current propellant mass rate, as computed by last call to updateMembers function.
    /*!
     * Function to retreieve the current propellant mass rate, as computed by last call to updateMembers function.
     * \return ICurrent propellant mass rate, as computed by last call to updateMembers function.
     */
    double getCurrentMassRate( )
    {
        return currentMassRate_;
    }

    //! Function to retreieve the ID associated with the source of the thrust (i.e. engine name).
    /*!
     * Function to retreieve the ID associated with the source of the thrust (i.e. engine name).
     * \return ID associated with the source of the thrust (i.e. engine name).
     */
    std::string getAssociatedThrustSource( )
    {
        return associatedThrustSource_;
    }

    //! Function to retreieve the list of environment models that are to be updated before computing the acceleration.
    /*!
     * Function to retreieve the list of environment models that are to be updated before computing the acceleration.
     * \return List of environment models that are to be updated before computing the acceleration.
     */
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > getRequiredModelUpdates( )
    {
        return requiredModelUpdates_;
    }

    //! Function to set or add a thrust update function
    /*!
     * Function to set or add a thrust update function
     * \param thrustUpdateFunction Update function that is to be added to class
     */
    void setThrustUpdateFunction( const std::function< void( const double ) > thrustUpdateFunction )
    {
        if( ( thrustUpdateFunction_ == nullptr ) )
        {
            thrustUpdateFunction_ = thrustUpdateFunction;
        }
        else
        {
            thrustUpdateFunction_ = std::bind( &mergeUpdateFunctions, thrustUpdateFunction, thrustUpdateFunction_, std::placeholders::_1 );
        }
    }

protected:

    //! Function returning the current magnitude of the thrust.
    /*!
     * Function returning the current magnitude of the thrust. Any dependencies of the
     * thrust on (in)dependent variables is to be handled by the thrustUpdateFunction.
     */
    std::function< double( ) > thrustMagnitudeFunction_;

    //! Function returning the direction of the thrust (as a unit vector).
    /*!
     * Function returning the direction of the thrust (as a unit vector). Any dependencies of the
     * thrust on (in)dependent variables is to be handled by the thrustUpdateFunction.
     */
    std::function< Eigen::Vector3d( ) > thrustDirectionFunction_;

    //! Function returning the current mass of the body being propagated.
    std::function< double( ) > bodyMassFunction_;

    //! Function returning total propellant mass rate from the thrust system.
    std::function< double( ) > massRateFunction_;

    //! ID associated with the source of the thrust (i.e. engine name).
    std::string associatedThrustSource_;

    //! Function used to update the thrust magnitude and direction to current time.
    std::function< void( const double ) > thrustUpdateFunction_;

    //! Current acceleration, as computed by last call to updateMembers function.
    Eigen::Vector3d currentAcceleration_;

    //! Current acceleration direction, as computed by last call to updateMembers function.
    Eigen::Vector3d currentAccelerationDirection_;

    //! Current thrust magnitude, as computed by last call to updateMembers function.
    double currentThrustMagnitude_;

    //! Current propellant mass rate, as computed by last call to updateMembers function.
    double currentMassRate_;

    //! Function to reset the time in the classes to which the thrustUpdateFunction function directs
    const std::function< void( const double ) > timeResetFunction_;

    //! List of environment models that are to be updated before computing the acceleration,
    /*!
     * List of environment models that are to be updated before computing the acceleration,
     *  list is included here to account for versatility of dependencies of thrust model (guidance) algorithms.
     */
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > requiredModelUpdates_;

};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_THRUSTACCELERATIONMODEL_H
