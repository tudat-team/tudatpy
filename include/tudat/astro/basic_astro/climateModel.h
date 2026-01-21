/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CLIMATEMODEL_H
#define TUDAT_CLIMATEMODEL_H

#include <memory>
#include <vector>

namespace tudat
{

namespace simulation_setup
{

class Body;

}

}


namespace tudat
{

namespace environment
{

class ClimateModel {

    public:

    ClimateModel( std::shared_ptr< simulation_setup::Body > bodyWithClimateModel )
    { 
        bodyWithClimateModel_ = std::weak_ptr< simulation_setup::Body >( bodyWithClimateModel );
    };

    virtual ~ClimateModel() = default;

    virtual void update( double ) = 0;

    void addBodyRequiringClimateModel( std::shared_ptr< simulation_setup::Body > bodyRequiringClimateModel )
    {
        listBodiesRequiringClimateModel_.push_back( std::weak_ptr< simulation_setup::Body >( bodyRequiringClimateModel ) ); 
    };

    protected:

    std::weak_ptr< simulation_setup::Body > bodyWithClimateModel_;

    std::vector< std::weak_ptr< simulation_setup::Body > > listBodiesRequiringClimateModel_;
        
};

}

}

#endif