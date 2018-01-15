/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FROMTHRUSTMASSRATEMODEL_H
#define TUDAT_FROMTHRUSTMASSRATEMODEL_H

#include <map>
#include <vector>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"

namespace tudat
{
namespace propulsion
{

//! This class provides the mass rate from a (set of) thrust model(s), used in numerically propagating a body's total mass.
/*!
 *  This class provides the mass rate from a (set of) thrust model(s, used in numerically propagating a body's total mass.
 *  The mass rate is directly obtained from the thrust accelerations. This class can be used to sum the mass rates due to
 *  any number of thrust acceleration models.
 */
class FromThrustMassRateModel: public basic_astrodynamics::MassRateModel
{
public:

    //! Constructor taking single acceleration model.
    /*!
     *  Constructor taking single acceleration model.
     *  \param thrustAcceleration Thrust acceleration model from which the mass rate is to be retrieved.
     */
    FromThrustMassRateModel(
            const boost::shared_ptr< ThrustAcceleration > thrustAcceleration )
    {
        thrustAccelerations_.push_back( thrustAcceleration );
    }

    //! Constructor taking list of acceleration models.
    /*!
     *  Constructor taking list of acceleration models.
     *  \param thrustAccelerations List of thrust acceleration models from which the mass rates are to be retrieved and
     *  summed.
     */
    /*!
     * Constructor
     */
    FromThrustMassRateModel(
            const std::vector< boost::shared_ptr< ThrustAcceleration > > thrustAccelerations ):
        thrustAccelerations_( thrustAccelerations ){ }

    //! Destructor.
    ~FromThrustMassRateModel( ){ }

    //! Update member variables used by the mass rate model and compute the mass rate
    /*!
     * Update member variables used by the mass rate model and compute the mass rate.
     * Updates the thrusta acceleration models and sums the resulting mass rates, updating the currentMassRate_ variable.
     * \param currentTime Time at which mass rate model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        // Check if update is needed.
        if( currentTime != currentTime_ )
        {
            // Update thrust accelerations
            for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
            {
                thrustAccelerations_.at( i )->updateMembers( currentTime );
            }

            // Compute total mass rate
            currentMassRate_ = 0.0;
            for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
            {
                currentMassRate_ += thrustAccelerations_.at( i )->getCurrentMassRate( );
            }

            currentTime_ = currentTime;
        }
    }

private:

    //! List of thrust accelerations from which the total mass rate is to be computed.
    std::vector< boost::shared_ptr< ThrustAcceleration > > thrustAccelerations_;
};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_FROMTHRUSTMASSRATEMODEL_H
