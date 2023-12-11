/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_REFERENCEPOINTPOSITION_H
#define TUDAT_REFERENCEPOINTPOSITION_H

#include <Eigen/Core>

#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"


namespace tudat
{

namespace estimatable_parameters
{

class ReferencePointPosition: public EstimatableParameter< Eigen::VectorXd >
{

public:

    ReferencePointPosition( const std::shared_ptr< system_models::VehicleSystems > systemModels,
                            const std::string& associatedBody,
                            const std::string& associatedReferencePoint ):
        EstimatableParameter< Eigen::VectorXd  >( reference_point_position, associatedBody, associatedReferencePoint ),
        systemModels_( systemModels )
    {
        if( systemModels_->doesReferencePointExist( associatedReferencePoint ) == false )
        {
            throw std::runtime_error( "Error when making reference point position parameter for " +
                associatedBody + ", " + associatedReferencePoint + ", reference point not found" );
        }
    }

    //! Destructor    
    ~ReferencePointPosition( ) { }

    Eigen::VectorXd  getParameterValue( )
    {
        return systemModels_->getReferencePointPosition( parameterName_.second.second );
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        systemModels_->setReferencePointPosition( parameterName_.second.second, parameterValue );
    }

    int getParameterSize( )
    {
        return 3;
    }

protected:

private:

    std::shared_ptr< system_models::VehicleSystems > systemModels_;
};

}

}

#endif // TUDAT_REFERENCEPOINTPOSITION_H
