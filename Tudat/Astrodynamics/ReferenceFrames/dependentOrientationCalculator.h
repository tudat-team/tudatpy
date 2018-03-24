/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_DEPENDENTORIENTATIONCALCULATOR_H
#define TUDAT_DEPENDENTORIENTATIONCALCULATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{

namespace reference_frames
{

//! Class that computes the orientation of a body as a function of the current state and time
/*!
 *  Class that computes the orientation of a body as a function of the current state and time. The functionality of this
 *  class is distinct from the RotationalEphemeris, in the sense that it is only defined during the numerical propagation,
 *  as it relies on the current state of the environment for the computation of the rotations.
 */
class DependentOrientationCalculator
{
public:

    //! Constructor.
    DependentOrientationCalculator( ): currentTime_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~DependentOrientationCalculator( ){ }

    //! Function to get the current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
    /*!
     * Function to get the current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
     * Function is pure virtual and must be defined in derived class.
     * \return Current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
     */
    virtual Eigen::Quaterniond getRotationToLocalFrame( ) = 0;


    //! Function to get the current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
    /*!
     * Function to get the current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
     * Function is pure virtual and must be defined in derived class.
     * \return Current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
     */
    virtual Eigen::Quaterniond getRotationToGlobalFrame( ) = 0;

    //! Pure virtual function to update the object to the current state.
    /*!
     *  Pure virtual function to update the object to the current state.
     *  \param currentTime Time to which angle calculator is to be updated.
     */
    virtual void updateCalculator( const double currentTime ) = 0;

    //! Function to update the object to teh current time and retrieve the rotation to the local frame
    /*!
     * Function to update the object to teh current time and retrieve the rotation to the local frame
     * Note that the environment must have been updated to the current time before calling this function.
     * \param currentTime Time to which angle calculator is to be updated.
     * \return Current rotation from the global (propagation/inertial) to the local (body-fixed) frame.
     */
    Eigen::Quaterniond getRotationToLocalFrame( const double currentTime )
    {
        updateCalculator( currentTime );
        return getRotationToLocalFrame( );
    }

    //! Function to update the object to teh current time and retrieve the rotation to the global frame
    /*!
     * Function to update the object to teh current time and retrieve the rotation to the global frame
     * Note that the environment must have been updated to the current time before calling this function.
     * \param currentTime Time to which angle calculator is to be updated.
     * \return Current rotation from the local (body-fixed) to the global (propagation/inertial) frame.
     */
    Eigen::Quaterniond getRotationToGlobalFrame( const double currentTime )
    {
        updateCalculator( currentTime );
        return getRotationToGlobalFrame( );
    }

    //! Function to reset the current time in the derived class.
    /*!
     * Function to reset the current time in the derived class, must be implemented in the derived class.
     * \sa resetCurrentTime
     * \param currentTime New value of currentTime_ variable.
     */
    virtual void resetDerivedClassTime( const double currentTime = TUDAT_NAN ){ }

    //! Function to reset the current time of the object
    /*!
     * Function to reset the current time of the object. Note that the orientation is not recomputed by calling this function
     * only the currentTime_ (and any associated variables) are reset. Function is typically used to reset time to NaN,
     * signalling the need to recompute all variables.
     * \param currentTime New value of currentTime_ variable.
     */
    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        resetDerivedClassTime( currentTime );
    }

protected:

    //! Current simulation time.
    double currentTime_;
};

} // namespace reference_frames

} // namespace tudat

#endif // TUDAT_DEPENDENTORIENTATIONCALCULATOR_H
