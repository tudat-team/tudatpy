/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity/relativisticLightTimeCorrection.h"

namespace tudat
{

namespace relativity
{

//! Function to calculate first order relativistic light time correction due to a gravitating point mass.
double calculateFirstOrderLightTimeCorrectionFromCentralBody( const double bodyGravitationalParameter,
                                                              const Eigen::Vector3d& transmitterPosition,
                                                              const Eigen::Vector3d& receiverPosition,
                                                              const Eigen::Vector3d& centralBodyPosition,
                                                              const double ppnParameterGamma )
{
    // Calculate Euclidean geometric distances between transmitter, receiver and gravitating body.
    double distanceToReceiver = ( receiverPosition - centralBodyPosition ).norm( );
    double distanceToTransmitter = ( transmitterPosition - centralBodyPosition ).norm( );    
    double linkEuclideanDistance = ( transmitterPosition - receiverPosition ).norm( );

    // Calculate and return light time correction.
    return ( 1.0 + ppnParameterGamma ) * bodyGravitationalParameter * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * std::log(
                ( distanceToReceiver + distanceToTransmitter + linkEuclideanDistance ) /
                ( distanceToReceiver + distanceToTransmitter - linkEuclideanDistance ) );

}

//! Function to calculate gradient of first order relativistic light time correction due to a gravitating point mass.
Eigen::Matrix< double, 1, 3 > calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
        const double bodyGravitationalParameter,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const Eigen::Vector3d& centralBodyPosition,
        const bool evaluateGradientAtReceiver,
        const double ppnParameterGamma )
{
    Eigen::Vector3d relativePositionVector = ( receiverPosition - transmitterPosition );
    double receiverDistance = ( receiverPosition - centralBodyPosition ).norm( );
    double transmitterDistance = ( transmitterPosition - centralBodyPosition ).norm( );
    double linkEndDistance = relativePositionVector.norm( );

    Eigen::Matrix< double, 1, 3 > gradient = ( receiverDistance + transmitterDistance ) *
            ( relativePositionVector.normalized( ) ).transpose( );
    if( evaluateGradientAtReceiver )
    {
       gradient -= relativePositionVector.norm( ) * ( ( receiverPosition - centralBodyPosition ).normalized( ) ).transpose( );
    }
    else
    {
        gradient += relativePositionVector.norm( ) * ( ( transmitterPosition - centralBodyPosition ).normalized( ) ).transpose( );
        gradient *= -1.0;
    }

    return 2.0 * ( 1.0 + ppnParameterGamma ) * bodyGravitationalParameter * physical_constants::INVERSE_CUBIC_SPEED_OF_LIGHT * gradient /
            ( ( receiverDistance + transmitterDistance ) * ( receiverDistance + transmitterDistance ) -
                  linkEndDistance * linkEndDistance );


}

}

}

