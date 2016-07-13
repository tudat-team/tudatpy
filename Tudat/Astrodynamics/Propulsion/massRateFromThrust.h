/*    Copyright (c) 2010-2016, Delft University of Technology
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
#include <iostream>
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

class FromThrustMassRateModel: public basic_astrodynamics::MassRateModel
{
public:

    //! Constructor.
    /*!
     * Constructor
     */
    FromThrustMassRateModel(
            const boost::shared_ptr< ThrustAcceleration > thrustAcceleration )
    {
        thrustAccelerations_.push_back( thrustAcceleration );
    }

    //! Constructor.
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
     * Update member variables used by the mass rate model and compute the mass rate
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
        {
            thrustAccelerations_.at( i )->updateMembers( currentTime );
        }

        currentMassRate_ = 0.0;
        for( unsigned int i = 0; i < thrustAccelerations_.size( ); i++ )
        {
            currentMassRate_ += thrustAccelerations_.at( i )->getCurrentMassRate( );
        }

    }

private:

    //! Function returning mass rate as a function of time.
    std::vector< boost::shared_ptr< ThrustAcceleration > > thrustAccelerations_;
};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_FROMTHRUSTMASSRATEMODEL_H
