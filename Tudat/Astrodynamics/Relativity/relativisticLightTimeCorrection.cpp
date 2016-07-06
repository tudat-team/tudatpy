#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/relativisticLightTimeCorrection.h"

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

}

}

