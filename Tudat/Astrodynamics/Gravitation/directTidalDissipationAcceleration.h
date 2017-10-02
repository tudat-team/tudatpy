/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
#define TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

namespace tudat
{

namespace gravitation
{

Eigen::Vector3d computeDirectTidalDissipationAcceleration(
        const Eigen::Vector6d satelliteRelativeState, const Eigen::Vector3d planetAngularVelocityVector,
        const double satelliteMass, const double k2LoveNumber, const double timeLag, const double planetReferenceRadius,
        const bool includeDirectRadialComponent = true );

class DirectTidalDissipationAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    DirectTidalDissipationAcceleration(
            const boost::function< Eigen::Vector6d( ) > satelliteStateFunction,
            const boost::function< Eigen::Vector6d( ) > planetStateFunction,
            const boost::function< double( ) > satelliteMassFunction,
            const boost::function< Eigen::Vector3d( ) > planetAngularVelocityVector,
            const double k2LoveNumber,
            const double timeLag,
            const double planetEquatorialRadius,
            const bool includeDirectRadialComponent ):
        satelliteStateFunction_( satelliteStateFunction ), planetStateFunction_( planetStateFunction ),
        satelliteMassFunction_( satelliteMassFunction ),
        planetAngularVelocityVector_( planetAngularVelocityVector ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), planetEquatorialRadius_( planetEquatorialRadius ),
        includeDirectRadialComponent_( includeDirectRadialComponent )
    { }


    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
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
            currentAcceleration_ = computeDirectTidalDissipationAcceleration(
                        satelliteStateFunction_( ) - planetStateFunction_( ), planetAngularVelocityVector_( ),
                        satelliteMassFunction_( ),
                        k2LoveNumber_, timeLag_, planetEquatorialRadius_, includeDirectRadialComponent_ );

            currentTime_ = currentTime;
        }
    }

private:

    boost::function< Eigen::Vector6d( ) > satelliteStateFunction_;

    boost::function< Eigen::Vector6d( ) > planetStateFunction_;

    boost::function< double( ) > satelliteMassFunction_;

    boost::function< Eigen::Vector3d( ) > planetAngularVelocityVector_;

    double k2LoveNumber_;

    double timeLag_;

    double planetEquatorialRadius_;

    bool includeDirectRadialComponent_;

    Eigen::Vector3d currentAcceleration_;
};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
