/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIRECTTIDALTIMELAG_H
#define TUDAT_DIRECTTIDALTIMELAG_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace estimatable_parameters
{


class DirectTidalTimeLag: public EstimatableParameter< double >
{
public:

    DirectTidalTimeLag(
            const std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels,
            const std::string& deformedBody,
            const std::vector< std::string > bodiesCausingDeformation = std::vector< std::string >( ) ):
        EstimatableParameter< double >(
            direct_dissipation_tidal_time_lag, deformedBody ),
        tidalAccelerationModels_( tidalAccelerationModels ),
        bodiesCausingDeformation_( bodiesCausingDeformation )
    {
        for( unsigned int i = 1; i < tidalAccelerationModels_.size( ); i++ )
        {
            if( tidalAccelerationModels_.at( i )->getTimeLag( ) != tidalAccelerationModels_.at( 0 )->getTimeLag( ) )
            {
                std::cerr<<"Warning when making direct tidal time lag parameter. Time lags are different in model upon creation, but will be estimated to the same value"<<std::endl;
            }
        }
    }

    double getParameterValue( )
    {
        return tidalAccelerationModels_.at( 0 )->getTimeLag( );
    }

    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < tidalAccelerationModels_.size( ); i++ )
        {
            tidalAccelerationModels_.at( i )->resetTimeLag( parameterValue );
        }
    }


    int getParameterSize( )
    {
        return 1;
    }

    std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > getTidalAccelerationModels( )
    {
        return tidalAccelerationModels_;
    }

    std::vector< std::string > getBodiesCausingDeformation( )
    {
        return bodiesCausingDeformation_;
    }


private:

    std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels_;

    std::vector< std::string > bodiesCausingDeformation_;

};

}

}

#endif // TUDAT_DIRECTTIDALTIMELAG_H
