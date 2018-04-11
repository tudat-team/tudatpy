/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/secondDegreeGravitationalTorquePartial.h"


namespace tudat
{

namespace acceleration_partials
{

//! Constructor
SecondDegreeGravitationalTorquePartial::SecondDegreeGravitationalTorquePartial(
        const boost::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > torqueModel,
        const std::string acceleratedBody,
        const std::string acceleratingBody ):
    TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::second_order_gravitational_torque ),
    torqueModel_( torqueModel ){ }

}

}
